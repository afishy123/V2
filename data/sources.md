# Sources & provenance — `fusion_controls.csv`

Where every value in `fusion_controls.csv` came from. Cases fall into: **positive controls** (enhancer
hijacking, wet-lab confirmed), a **negative control** (kinase fusion), **secondary**/**rejected** candidates,
and **clinical fusions** (Sergio's list). Genomic coordinates are **hg38** for the enhancer-hijacking driven-gene
TSS + enhancer; **breakpoints for fusion genes are hg19** (from the fusion catalog) and need liftover before use.

## Column legend
| column | meaning |
|---|---|
| `category` | positive_control / negative_control / secondary_candidate / rejected_candidate / clinical_fusion / worked_case |
| `mechanism` | enhancer_hijacking / promoter_donation / kinase_fusion / boundary_disruption / regulatory / uncertain |
| `driven_gene` + `driven_gene_TSS_hg38` + `strand` | the oncogene whose expression changes, and its TSS (hg38) |
| `regulatory_element` + `element_coords_hg38` | the hijacked enhancer (for enhancer-hijacking cases), hg38 |
| `fusion_5p_gene` / `fusion_3p_gene` + `breakpoint_5p/3p` + `breakpoint_assembly` | for fusion genes (hg19) |
| `wetlab_status` | the strength of causal proof (CRISPR_confirmed / associative / none) |
| `model_result` | what our in-silico pipeline found (captured / weak / null / untested) |
| `coord_source` | where the coordinates were obtained |

---

## Wet-lab / mechanism references (per case)
- **GATA2-MECOM** — Gröschel et al., *Cell* 2014, **PMID 24703711**, doi:10.1016/j.cell.2014.02.019 (CRISPR enhancer excision re-silences EVI1). Co-discovery: Yamazaki et al., *Cancer Cell* 2014, **PMID 24703906**.
- **TCR-TAL1** — t(1;14)(p32;q11): Begley et al., *PNAS* 1989, **PMID 2467296** (named the gene **SCL**; DU.528 line) + Chen et al., *EMBO J* 1990, **PMID 2303035** (breakpoint resolution). Hijacked enhancer = TCRδ **Eδ**: Redondo et al., *Science* 1990, **PMID 2156339** (GenBank M33967). ⚠️ Mansour et al., *Science* 2014, **PMID 25394790** (MuTE neo-enhancer) is a SEPARATE *cis* mechanism, NOT this control. Gene TAL1 a.k.a. SCL/TCL5.
- **ETV6-MNX1** — Naranjo/Heim et al., *Blood Advances* 2024, **PMID 39121370**, doi:10.1182/bloodadvances.2024013116 (CRISPR ΔEn of the ETV6 enhancer block abrogates MNX1 in a t(7;12) iPSC→HSPC model). bioRxiv 2023.09.13.557546.
- **MIPOL1-ETV1** (replaces LMO2) — Wang et al., *Nat Methods* 2021, **PMID 34092790** (t(7;14) enhancer hijacking in LNCaP; CRISPR deletion of two chr14q13/MIPOL1-intronic enhancers lowers ETV1 by 66–80%; ETV1 is NOT fused → genuine enhancer adoption). ⚠️ LMO2 removed: its cited paper (Rahman/Abraham 2017, **PMID 28270453**) is a neomorphic *promoter*, not enhancer hijacking; the t(11;14) translocation is loss-of-repression. See project memory [[tcr-controls-wrong-citation]].
- **BCL11B** — Montefiori et al., *Cancer Discovery* 2021, **PMID 34103329**, doi:10.1158/2159-8290.CD-21-0145 (HiChIP + allele-specific; *associative*). CRISPRi causal proof is for the reverse target TLX3: Botten et al., *Blood* 2023, **PMID 36947815**. SV data: EGA EGAS00001004810.
- **IGF2 / IRS4** — Weischenfeldt et al., *Nat Genet* 2017, **PMID 28604730**, doi:10.1038/ng.3722 (4C + luciferase; no enhancer deletion). IGF2 is imprinted.
- **TERT** — Peifer et al., *Nature* 2015, **PMID 26466568**, doi:10.1038/nature14980; Valentijn et al., *Nat Genet* 2015, **PMID 26523776** (correlative). The CRISPR-SV functional paper (Xu/Ren, *Nature* 2022, PMID 36477537) validated MYC, not TERT → rejected.
- **GFI1 / GFI1B** — Northcott et al., *Nature* 2014, **PMID 25043047**, doi:10.1038/nature13379 (luciferase + mouse co-transformation only; no enhancer deletion) → rejected.
- **EGFR-SEPT14** — Frattini et al., *Nat Genet* 2013, **PMID 23917401** (chimeric EGFR kinase → STAT3, in *trans*). The negative control.
- **TMPRSS2-ERG** — Tomlins et al., *Science* 2005, **PMID 16254181** (discovery); Kron et al., *Nat Genet* 2017, **PMID 28783165** (super-enhancer / promoter-swap mechanism).
- **PTPRZ1-MET** — Bao et al. 2016 (Chinese Glioma Genome Atlas; PTPRZ1 promoter drives MET overexpression in secondary GBM).
- **FGFR3-TACC3** — Singh et al., *Science* 2012, **PMID 22837387**.
- **BCR-ABL1**, **ABL1/BCR-NTRK2**, **EGFR-VOPP1** — recurrent fusions; mechanism from review literature (kinase fusions / overexpression).
- **TBC1D13-CDK4, QSER1-DCDC1, DCDC1-ELP4, IFITM3-IGF2R** — from Sergio's list; **mechanism + breakpoints not yet sourced** (flagged in the CSV).

## Coordinate provenance (genome assemblies & accessions)
- **Driven-gene TSS (hg38)** — Ensembl REST (`/lookup/id`) canonical/specified transcripts:
  EVI1 `ENST00000264674` (chr3:169,146,305, −); MNX1 `ENST00000252971` (chr7:157,010,663, −);
  BCL11B `ENST00000357195` (chr14:99,272,197, −); TMPRSS2 `ENST00000332149` (chr21:41,508,158, −).
  NCBI Gene / UCSC RefSeq: TAL1 Gene 6886; ETV1 Gene 2115 (ENST00000430479 / NM_004956.5, chr7:13,989,666, −); GFI1B Gene 8328 / GFI1 Gene 2672;
  TERT `NM_198253.3` (chr5:1,295,183, −); IGF2 `NM_000612` (chr11:2,139,389, −).
- **Enhancer coordinates (hg38)** — G2DHE: Gröschel 2014 + the MonoMAC −110 enhancer SNV (PMC10587712).
  TCRδ (Eδ) — the TAL1 t(1;14) hijacked enhancer: **chr14:22,460,450–22,460,810 (hg38, 360 bp)** = GenBank M33967
  (Redondo *Science* 1990, PMID 2156339); Jδ3–Cδ intron, early/DN enhancer; ENCODE pELS E1702206. This is NOT Eα.
  TCR-α (Eα): ~3′ of TRAC, chr14 (≈ chr14:22,555,000–22,557,000, approximate) — was wrongly used for TAL1 (~93 kb
  off, and Eα is the later DP-stage enhancer); now UNUSED (LMO2 was replaced by the ETV1 control).
  ETV1 (t(7;14) hijack, LNCaP): chr14q13 / MIPOL1-intronic DNase enhancers — M2 chr14:37,246,515–37,246,652 (primary,
  80.5% CRISPR effect) + M1 chr14:37,278,706–37,279,509 (66.3%); hg38 (Wang *Nat Methods* 2021, PMID 34092790).
  MNX1 ETV6 block: Naranjo 2024 (hg19 chr12:11,951,022–12,164,578 → hg38 chr12:11,798,088–12,011,644).
  This block IS the CRISPR-validated ΔEn region (Naranjo 2024 deletes all 213.5 kb at once); it contains
  **4 candidate enhancers (2 p300+), HSPC-specific** and the active sub-element(s) are **unresolved** — the
  paper does not localize them, so this is a faithful caveat, not a coordinate gap.
  **Both builds keep the FULL 213.5 kb block, edge-anchored** (`anchor='edge'`), so it never overwrites the
  promoter: the **TSS-centered** build edge-anchors it and **sweeps the near-edge distance**; the
  **event_site-centered** build keeps the block on the partner arm with the promoter intact on the oncogene
  arm. (The earlier "2 kb block-start proxy" approach was removed — it was unjustified, since the 4 enhancers
  are distributed across the block, not at the start.)
  BCL11B BETA, IGF2 SE, TERT/GFI1B enhancers: **coordinates unresolved** — must be pulled from each paper's
  supplementary tables (flagged "UNRESOLVED" in the CSV).
- **Fusion-gene breakpoints (hg19)** — the FusionGDB / **ChimerDB4 (TCGA)** breakpoint table
  (`raw/fusiongdb_breakpoints_hg19.txt`): EGFR-SEPT14 (TCGA-06-0750-01A), TMPRSS2-ERG, PTPRZ1-MET (LGG),
  EGFR-VOPP1 (ESCA), BCR-ABL1, ABL1-NTRK2 (ESCA). **hg19 → liftover required.**

## Datasets & tools used
- **FusionGDB / ChimerDB4** — fusion-transcript breakpoint catalog (TCGA). Source: compbio.uth.edu/FusionGDB ; kobic.re.kr/chimerdb. (Local: `raw/fusiongdb_breakpoints_hg19.txt`, hg19, **not in git — re-download**.)
- **JASPAR 2024** CORE vertebrates — TF motifs (`raw/jaspar2024_core_vert.meme`). jaspar.elixir.no.
- **AlphaGenome** — model + cell-type ontologies (K562 `EFO:0002067`; CD34+ CMP `CL:0001059`; hematopoietic MPP `CL:0000837`; CD4/CD8 T `CL:0000624`/`CL:0000625`; thymus CD4/CD8 `CL:0000895`/`CL:0000900`; colon `UBERON:0001155`; brain `UBERON:0000955`; prostate `UBERON:0002367`; neuroblastoma `EFO:0000621`). alphagenomedocs.com.
- **Borzoi** — open-weights model; human track set in `raw/targets_human.txt` (in git).
- **Ensembl REST** (rest.ensembl.org), **UCSC** (genome.ucsc.edu), **NCBI Gene** — coordinate lookups.

## Method precedent (why the project tests this)
- **MethylSeqNet** — Dixon-Luinenburg, Bajwa, …, Ioannidis, *bioRxiv* 2026, doi:10.64898/2026.06.02.729723. Fig. 4 (PDK3 enhancer adoption) showed plain Borzoi weakly detects a fusion-created distal enhancer — the seed of this project. (Local: `raw/methylseqnet-2026.md`.)

---
*Notes:* values marked `UNRESOLVED`, `NOT AVAILABLE`, or `(approx, VERIFY)` in the CSV still need exact
coordinates before use (paper supplements / ENCODE cCREs). Everything else is pipeline-ready.
