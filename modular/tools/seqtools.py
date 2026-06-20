"""Modular enhancer-hijacking pipeline — pick a model, run the whole analysis.

The ONLY model-specific object is the `Predictor` (AlphaGenome / Borzoi / Enformer). Everything
else — building WT/der sequences, the dinucleotide-shuffle control, exon fetch, plotting, the full
TSS-centred analysis and the distance sweep — is model-agnostic and keys off the predictor's
`window` (input bp) and `bin_bp` (output resolution).

    import seqtools as st
    from candidates import CANDIDATES
    cand = CANDIDATES['GATA2-MECOM']; cand['exons'] = st.fetch_exons(cand['tx'])
    st.analyze(st.Borzoi(), cand)          # full TSS analysis
    st.distance_sweep(st.Enformer(), cand) # distance decay

Run with the repo root as the working directory (the AlphaGenome/Borzoi paths in `fusionseq` are
relative: `.env`, `raw/targets_human.txt`). Sequence primitives are reused from `fusionseq`;
Enformer is self-contained in `enformer.py`.
"""
import os, sys, json, urllib.request
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

for _p in ['.', '..', '../..', '../../..']:
    if os.path.isdir(os.path.join(_p, 'fusionseq')):
        sys.path.insert(0, os.path.abspath(_p)); break
from fusionseq import hijack as hj, pipeline as fp, motifs as mo   # sequence primitives


# ===================================================================== model adapters
class Predictor:
    """sequence -> 1-D expression coverage. Subclasses set name/window/bin_bp/expr."""
    name = window = bin_bp = expr = None
    def cell(self, cand):                      raise NotImplementedError   # model's cell spec
    def predict_rna(self, seq, cell, strand):  raise NotImplementedError


class AlphaGenome(Predictor):
    name, window, bin_bp, expr = 'AlphaGenome', 1_048_576, 1, 'RNA-seq'
    def cell(self, cand): return cand['cell']['AlphaGenome']                # ontology CURIE
    def predict_rna(self, seq, cell, strand):
        out, _ = fp.predict_alphagenome({'x': seq}, output_type='RNA_SEQ',
                                        ontology_terms=[cell], strand=strand)
        return out['x']


class Borzoi(Predictor):
    name, window, bin_bp, expr = 'Borzoi', 524_288, 32, 'RNA-seq'           # output: central ~196 kb
    def cell(self, cand): return cand['cell']['Borzoi']                     # targets_human.txt keyword
    def predict_rna(self, seq, cell, strand):
        idx = fp.borzoi_tracks('RNA', cell)
        assert idx, f"no Borzoi RNA track matches {cell!r} in raw/targets_human.txt"
        return fp.predict_borzoi({'x': seq}, idx)['x']


from enformer import Enformer       # self-contained; expr is CAGE (Enformer has no RNA-seq head)
Enformer.expr = 'CAGE'

MODELS = {'AlphaGenome': AlphaGenome, 'Borzoi': Borzoi, 'Enformer': Enformer}


# ===================================================================== sequence helpers
def build_wt(pred, cand):
    return hj.build_wt(cand['onco_chrom'], cand['onco_tss'], pred.window)

def build_der(pred, cand, d_kb):
    der, _ = hj.build_der(cand['onco_chrom'], cand['onco_tss'], cand['enh_chrom'],
                          cand['enh_start'], cand['enh_end'], d_kb=d_kb,
                          window=pred.window, side='high', anchor=cand.get('anchor', 'center'))
    return der

def scramble(wt_seq, der_seq, seed=0):
    """der with its inserted enhancer dinucleotide-shuffled (composition kept, motifs destroyed)."""
    d = np.where(np.frombuffer(wt_seq.encode(), 'S1') != np.frombuffer(der_seq.encode(), 'S1'))[0]
    if not len(d):
        return der_seq
    o, e = int(d[0]), int(d[-1]) + 1
    return der_seq[:o] + mo.dinuc_shuffle(der_seq[o:e], np.random.default_rng(seed)) + der_seq[e:]

def enh_genomic(wt_seq, der_seq, tss, window):
    """genomic span of the inserted element (= where der differs from WT)."""
    d = np.where(np.frombuffer(wt_seq.encode(), 'S1') != np.frombuffer(der_seq.encode(), 'S1'))[0]
    return (tss - window // 2 + int(d[0]), tss - window // 2 + int(d[-1]) + 1) if len(d) else None

def fetch_exons(tx):
    url = f'https://rest.ensembl.org/lookup/id/{tx}?expand=1;content-type=application/json'
    return [(e['start'], e['end']) for e in json.load(urllib.request.urlopen(url, timeout=60))['Exon']]


# ===================================================================== plotting (bin_bp-aware)
def coords(pred, track, tss, strand):
    """x (kb from TSS, transcription-oriented) and a genomic->kb mapper, at the model's resolution."""
    flip = -1 if strand == '-' else 1
    x = flip * (np.arange(len(track)) - len(track) // 2) * pred.bin_bp / 1000
    span = lambda s, e: sorted((flip * (s - tss) / 1000, flip * (e - tss) / 1000))
    return x, span

def _gene_half(pred, track, tss, strand, exon_list, pad=10):
    _, span = coords(pred, track, tss, strand)
    ga, gb = span(min(s for s, e in exon_list), max(e for s, e in exon_list))
    return max(abs(ga), abs(gb)) + pad

def plot_tracks(pred, gene, tss, strand, exon_list, series, title, desc,
                half=None, enh=None, ymax_fixed=None):
    """series = list of (1-D track, label, color, linestyle). One figure, any number of tracks."""
    t0 = series[0][0]
    x, span = coords(pred, t0, tss, strand)
    if half is None:
        half = len(t0) // 2 * pred.bin_bp / 1000          # full model output extent
    ga, gb = span(min(s for s, e in exon_list), max(e for s, e in exon_list))
    vis = np.abs(x) <= half
    ymax = ymax_fixed if ymax_fixed is not None else max(float(max(tr[vis].max() for tr, *_ in series)) * 1.15, 1e-4)
    eh = ymax * 0.05

    fig, ax = plt.subplots(figsize=(9, 3.3))
    legend = []
    for tr, label, color, ls in series:
        ax.plot(x, tr, lw=0.9, color=color, ls=ls, zorder=4)
        legend.append(Line2D([0], [0], color=color, lw=1.5, ls=ls, label=label))
    ax.axhline(0, color='0.7', lw=0.5); ax.axvline(0, color='0.7', ls=':', lw=0.7)
    for s, e in exon_list:                                  # exons (red) + enhancer (orange) on one row
        a, b = span(s, e)
        ax.add_patch(plt.Rectangle((min(a, b), -1.5 * eh), max(abs(b - a), 0.05), eh,
                                   color='firebrick', lw=0, zorder=5, clip_on=True))
    if enh is not None:
        ex0, ex1 = span(*enh)
        if ex1 >= -half and ex0 <= half:
            ax.add_patch(plt.Rectangle((max(ex0, -half), -1.5 * eh), min(ex1, half) - max(ex0, -half), eh,
                                       color='darkorange', lw=0, zorder=6, clip_on=True))
    ax.plot([max(ga, -half), min(gb, half)], [-2.7 * eh, -2.7 * eh],         # gene (green) below
            lw=4, color='seagreen', solid_capstyle='butt', zorder=5, clip_on=True)
    ax.text(np.clip((ga + gb) / 2, -half, half), -2.7 * eh, gene.split()[0],
            ha='center', va='center', fontsize=8, fontweight='bold', color='white', zorder=7,
            bbox=dict(boxstyle='round,pad=0.2', facecolor='seagreen', edgecolor='none'))
    ax.set_xlim(-half, half); ax.set_ylim(-3.7 * eh, ymax)
    ax.set_ylabel(f'pred {pred.expr}'); ax.set_xlabel('kb from TSS (gene body = positive)')
    legend += [Patch(facecolor='seagreen', label='gene'), Patch(facecolor='firebrick', label='exon')]
    if enh is not None:
        legend.append(Patch(facecolor='darkorange', label='enhancer'))
    ax.legend(handles=legend, loc='upper right', fontsize=7, frameon=True, framealpha=0.9)
    ax.set_title(f"{title}\n{desc}", fontsize=10, loc='left')
    fig.tight_layout()
    return fig


# ===================================================================== high-level pipeline
def analyze(pred, cand, base_kb=10, seed=0):
    """The full TSS-centred analysis for one candidate + model: WT baseline, WT-vs-der (wide +
    gene close-up), and the dinucleotide-scramble specificity control (vs WT and vs der)."""
    cell, strand, tss, ex, g = pred.cell(cand), cand['strand'], cand['onco_tss'], cand['exons'], cand['gene']
    wt_seq, der_seq = build_wt(pred, cand), build_der(pred, cand, base_kb)
    wt  = pred.predict_rna(wt_seq,  cell, strand)
    der = pred.predict_rna(der_seq, cell, strand)
    shf = pred.predict_rna(scramble(wt_seq, der_seq, seed), cell, strand)
    enh = enh_genomic(wt_seq, der_seq, tss, pred.window)
    gh  = _gene_half(pred, wt, tss, strand, ex)
    tag = f"{cell} · {pred.name} ({pred.window // 1000} kb in, {pred.bin_bp} bp bins)"

    plot_tracks(pred, g, tss, strand, ex, [(wt, 'WT', 'steelblue', '-')],
                f"{g} — WT baseline · {pred.name}", f"predicted {pred.expr} · {tag}")
    plot_tracks(pred, g, tss, strand, ex, [(wt, 'WT', '0.55', '--'), (der, 'der', 'steelblue', '-')],
                f"{g} — WT vs der (enhancer @{base_kb} kb) · {pred.name}", f"predicted {pred.expr} · {tag}", enh=enh)
    plot_tracks(pred, g, tss, strand, ex, [(wt, 'WT', '0.55', '--'), (der, 'der', 'steelblue', '-')],
                f"{g} — WT vs der, gene close-up · {pred.name}", f"predicted {pred.expr} · {tag}", half=gh, enh=enh)
    plot_tracks(pred, g, tss, strand, ex, [(wt, 'WT', '0.55', '--'), (shf, 'scrambled', 'darkviolet', '-')],
                f"{g} — scrambled enhancer vs WT · {pred.name}", f"specificity control · {tag}", half=gh)
    plot_tracks(pred, g, tss, strand, ex, [(der, 'der', 'steelblue', '-'), (shf, 'scrambled', 'darkviolet', '--')],
                f"{g} — scrambled enhancer vs der · {pred.name}", f"specificity control · {tag}", half=gh)
    return dict(wt=wt, der=der, scrambled=shf, enh=enh)


def distance_sweep(pred, cand, distances=(10, 30, 50, 100)):
    """Transplant the enhancer at each distance that fits the model's window; plot der vs WT on a
    fixed y so the decay with distance is visible. Distances too far for the window are skipped."""
    cell, strand, tss, ex, g = pred.cell(cand), cand['strand'], cand['onco_tss'], cand['exons'], cand['gene']
    w = cand['enh_end'] - cand['enh_start']
    fit = [d for d in distances if d * 1000 + w < pred.window // 2 - 1000]
    wt = pred.predict_rna(build_wt(pred, cand), cell, strand)
    der_at = {d: pred.predict_rna(build_der(pred, cand, d), cell, strand) for d in fit}
    gh = _gene_half(pred, wt, tss, strand, ex)
    x, _ = coords(pred, wt, tss, strand); gvis = np.abs(x) <= gh
    ymax = max(float(der_at[d][gvis].max()) for d in fit) * 1.15
    for d in fit:
        plot_tracks(pred, g, tss, strand, ex,
                    [(wt, 'WT', '0.55', '--'), (der_at[d], f'der @{d} kb', 'steelblue', '-')],
                    f"{g} — enhancer @ {d} kb · {pred.name}",
                    f"distance sweep (fixed y) · {cell} · {pred.name}", half=gh, ymax_fixed=ymax)
    skipped = [d for d in distances if d not in fit]
    if skipped:
        print(f"{pred.name}: skipped {skipped} kb (beyond the {pred.window // 1000} kb window)")
    return der_at
