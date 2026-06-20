"""Enhancer-hijacking controls, with the cell spec each model needs.

`cell` maps model name -> its cell handle: AlphaGenome uses an ontology CURIE, Borzoi and Enformer
use a keyword matched against their track descriptions. GATA2-MECOM (EVI1/K562) is the cleanest
cross-model case (K562 exists in all three). Borzoi/Enformer keyword availability for CD4 / CD34
varies — inspect with `fusionseq.pipeline.borzoi_tracks('RNA', kw)` / `enformer.enformer_tracks(kw)`.
"""
CANDIDATES = {
    'GATA2-MECOM': dict(
        name='GATA2-MECOM', gene='EVI1 (MECOM)', onco_chrom='chr3', onco_tss=169_146_305, strand='-',
        enh_chrom='chr3', enh_start=128_603_548, enh_end=128_604_548, anchor='center',
        tx='ENST00000264674',
        cell={'AlphaGenome': 'EFO:0002067', 'Borzoi': 'K562', 'Enformer': 'K562'}),

    'TCR-TAL1': dict(  # TAL1 a.k.a. SCL/TCL5; t(1;14) hijacks the TCRδ enhancer Eδ (Redondo 1990, GenBank M33967)
        name='TCR-TAL1', gene='TAL1', onco_chrom='chr1', onco_tss=47_232_335, strand='-',
        enh_chrom='chr14', enh_start=22_460_450, enh_end=22_460_810, anchor='center',   # Eδ (NOT Eα, ~93kb away)
        tx='ENST00000691006',
        cell={'AlphaGenome': 'CL:0000624', 'Borzoi': 'CD4', 'Enformer': 'CD4+ T'}),

    'MIPOL1-ETV1': dict(  # t(7;14) enhancer hijacking (Wang Nat Methods 2021, PMID 34092790); replaces discredited LMO2
        name='MIPOL1-ETV1', gene='ETV1', onco_chrom='chr7', onco_tss=13_989_666, strand='-',
        enh_chrom='chr14', enh_start=37_246_515, enh_end=37_246_652, anchor='center',  # M2 enhancer (80.5% CRISPR); M1=37_278_706-37_279_509 (66.3%)
        tx='ENST00000430479',
        cell={'AlphaGenome': 'UBERON:0002367', 'Borzoi': 'LNCaP', 'Enformer': 'prostate'}),  # AlphaGenome has NO LNCaP RNA track -> use prostate (UBERON:0002367); LNCaP kept for Borzoi DNase/ChIP

    'ETV6-MNX1': dict(
        name='ETV6-MNX1', gene='MNX1', onco_chrom='chr7', onco_tss=157_010_663, strand='-',
        enh_chrom='chr12', enh_start=11_798_088, enh_end=12_011_644, anchor='edge',   # 213 kb block
        tx='ENST00000252971',
        cell={'AlphaGenome': 'CL:0001059', 'Borzoi': 'CD34', 'Enformer': 'CD34'}),
}
