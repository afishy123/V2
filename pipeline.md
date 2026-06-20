# v2 Simple pipeline

```
fusion_controls.csv  ->  load_candidates.ipynb  ->  data CSVs  ->  readout notebooks
```

**3 framings** (each reads its own data CSV):

```
tss_centered      ->  rna_seq, h3k27ac,distance_sweep
event_centered    ->  rna_seq, h3k27ac
n_padded          ->  rna_seq
```
# rna_seq includes WT alone, der alone, and wt/der layered

**multi-model:**

```
modular/analysis.ipynb  ->  AlphaGenome / Borzoi / Enformer
```
