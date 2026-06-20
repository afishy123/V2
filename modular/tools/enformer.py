"""Enformer predictor for the modular enhancer-hijacking pipeline.

Drops into the same `Predictor` interface used by AlphaGenome / Borzoi in `seqtools.py`:
a sequence -> 1-D predicted-expression (CAGE) coverage array on the gene's strand.

Enformer (DeepMind, ported to PyTorch by lucidrains as `enformer-pytorch`):
  - input  : 196,608 bp of DNA (one-hot (L,4) or token indices)
  - output : human head of shape (batch, 896, 5313); 896 bins x 128 bp covers the
             CENTRAL 896*128 = 114,688 bp of the 196,608-bp input.
  - expression proxy: CAGE tracks (there is NO RNA-seq head). The track set is the
    Basenji2/Enformer "targets_human.txt" (5313 rows, with a `description` column).

Usage:
    from enformer import Enformer, enformer_tracks
    pred = Enformer()
    y = pred.predict_rna(seq_196608bp, cell='K562', strand='+')   # -> np.ndarray (896,)

Requires: `pip install enformer-pytorch` (lucidrains), torch with CUDA.
The targets file `enformer_targets_human.txt` must sit next to this module; if it is
missing it is auto-downloaded on first use (see `_ensure_targets`).
"""
import os
import urllib.request
import functools
import numpy as np

# matches seqtools.Predictor without importing it (avoids pulling in fusionseq/.env deps
# when this module is used standalone). seqtools.Enformer can subclass either.
try:
    from seqtools import Predictor as _Predictor
except Exception:                                   # standalone use
    class _Predictor:
        name = window = bin_bp = None
        def cell(self, cand):                       raise NotImplementedError
        def predict_rna(self, seq, cell, strand):   raise NotImplementedError


_HERE = os.path.dirname(os.path.abspath(__file__))
TARGETS_PATH = os.path.join(_HERE, 'enformer_targets_human.txt')

# Canonical Basenji2/Enformer human targets (5313 rows). This is the exact ordering of
# the model's human head, NOT Borzoi's targets_human.txt (different model, different set).
_TARGETS_URL = ('https://raw.githubusercontent.com/calico/basenji/master/'
                'manuscripts/cross2020/targets_human.txt')

MODEL_ID = 'EleutherAI/enformer-official-rough'


# --------------------------------------------------------------------- targets table
def _ensure_targets():
    if not os.path.exists(TARGETS_PATH):
        urllib.request.urlretrieve(_TARGETS_URL, TARGETS_PATH)
    return TARGETS_PATH


@functools.lru_cache(maxsize=1)
def _targets():
    """Return list of (row_index, description) for all 5313 human tracks, in head order.

    targets_human.txt is tab-separated with a header row; columns include `index`,
    `identifier`, `description` (and others). The model head index == file data-row order
    (0..5312), which also equals the `index` column. We key off the `description` column.
    """
    _ensure_targets()
    rows = []
    with open(TARGETS_PATH) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        # locate the description column (fallback to last column)
        try:
            desc_col = header.index('description')
        except ValueError:
            desc_col = len(header) - 1
        for i, line in enumerate(fh):
            parts = line.rstrip('\n').split('\t')
            desc = parts[desc_col] if desc_col < len(parts) else parts[-1]
            rows.append((i, desc))
    assert len(rows) == 5313, f'expected 5313 human targets, got {len(rows)}'
    return rows


def enformer_tracks(keyword, assay='CAGE'):
    """Indices (into the 5313-wide human head) of `assay` tracks whose description
    contains `keyword` (case-insensitive). Default assay='CAGE' = Enformer's expression
    proxy. Descriptions look like 'CAGE:chronic myelogenous leukemia cell line:K562 ...'.
    Pass assay=None to match the keyword across all assays.
    """
    kw = keyword.lower()
    av = assay.lower() if assay else None
    out = []
    for i, desc in _targets():
        d = desc.lower()
        if av is not None and not d.startswith(av + ':') and (av + ':') not in d[:8]:
            # CAGE rows start with "CAGE:"; be lenient about leading whitespace/case
            if not d.lstrip().startswith(av + ':'):
                continue
        if kw in d:
            out.append(i)
    return out


def enformer_track_descs(keyword, assay='CAGE'):
    """(index, description) for tracks matching `enformer_tracks` — handy for inspection."""
    idx = set(enformer_tracks(keyword, assay))
    return [(i, d) for i, d in _targets() if i in idx]


# --------------------------------------------------------------------- model cache
_MODEL = None
_ONEHOT = None


def _load_model():
    """Load Enformer once, on GPU, in eval mode. Cached at module level."""
    global _MODEL, _ONEHOT
    if _MODEL is None:
        import torch
        from enformer_pytorch import from_pretrained, str_to_one_hot
        _ONEHOT = str_to_one_hot
        model = from_pretrained(MODEL_ID, use_checkpointing=True)
        model = model.to('cuda' if torch.cuda.is_available() else 'cpu').eval()
        _MODEL = model
    return _MODEL


# --------------------------------------------------------------------- predictor
class Enformer(_Predictor):
    name   = 'Enformer'
    window = 196_608          # Enformer input length (bp)
    bin_bp = 128              # output bin size; output covers CENTRAL 896*128 = 114,688 bp
    n_bins = 896

    def cell(self, cand):
        """Model-specific cell spec from a candidate dict (a CAGE-description keyword)."""
        return cand['cell']['Enformer']

    def predict_rna(self, seq, cell, strand):
        """seq: a `window`-length ACGT string; cell: a CAGE keyword; strand: '+'/'-'.
        Returns 1-D np.ndarray length 896 = mean predicted CAGE coverage over the cell's
        CAGE track(s).

        Strand note: Enformer's CAGE tracks are NOT cleanly +/- stranded (FANTOM5/ENCODE
        CAGE here are unstranded sample-level tracks), so `strand` is accepted for interface
        compatibility but does not select a strand-specific track; all CAGE tracks matching
        `cell` are averaged. Coverage is returned in genomic (input) orientation; the caller
        (seqtools.coords) flips the x-axis for '-' genes, matching the other predictors.
        """
        import torch

        assert len(seq) == self.window, (
            f'Enformer needs a {self.window}-bp sequence, got {len(seq)}')

        idx = enformer_tracks(cell, assay='CAGE')
        if not idx:
            raise ValueError(
                f"no CAGE track description contains {cell!r}; "
                f"try enformer_track_descs({cell!r}) to inspect candidates")

        model = _load_model()
        device = next(model.parameters()).device

        one_hot = _ONEHOT(seq).to(device)               # (L, 4) float
        if one_hot.ndim == 2:
            one_hot = one_hot.unsqueeze(0)              # (1, L, 4)

        with torch.no_grad():
            out = model(one_hot)                        # dict: 'human' (1,896,5313)
            human = out['human'] if isinstance(out, dict) else out
            human = human[0]                            # (896, 5313)
            track = human[:, idx].mean(dim=-1)         # mean over matching CAGE tracks
            return track.float().cpu().numpy()         # (896,)


__all__ = ['Enformer', 'enformer_tracks', 'enformer_track_descs',
           'TARGETS_PATH', 'MODEL_ID']
