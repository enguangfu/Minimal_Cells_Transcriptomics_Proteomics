#!/usr/bin/env python3
"""Assemble per-gene sense-TPM tables across RNA-seq platforms for Fig S1.

Background
----------
Fig S1 compares transcript quantification across the sequencing platforms used
for each organism: Syn1 has Illumina, PacBio Iso-Seq, and two ONT direct-RNA
runs; Syn3A has Illumina and one ONT run. This script builds the per-gene
sense-TPM tables the figure reads (compute-then-plot; the figure never recomputes).

Algorithm
---------
- Syn1 Illumina (avg over 3 libraries) and PacBio sense TPM are taken directly
  from Syn1_Transcriptomics/Gene_TPM/syn1_Illumina_PacBio_TPM_profiles.csv.
- Syn1 ONT run1 / run2 / merged sense TPM are computed here from the strand-split
  depth bedGraphs (sense = plus track for + genes, minus for - genes): per-gene
  mean per-base sense depth over the gene body, then length-implicit TPM
  normalization (depth / sum(depth) * 1e6), matching the depth-based TPM used
  for the Illumina/PacBio tracks.
- Syn3A Illumina and ONT sense TPM are taken from
  Syn3A_Transcriptomics/Gene_TPM/syn3A_TPM_Illumina_ONT.tsv.

Outputs
-------
- platform_TPM_syn1.tsv  : locus_tag, gene_name, gene_len, strand,
                           Illumina_TPM, PacBio_TPM, ONT1_TPM, ONT2_TPM, ONTmerged_TPM
- platform_TPM_syn3A.tsv : locus_tag, gene_name, gene_len, Illumina_TPM, ONT_TPM
- compute_platform_TPM.txt : pairwise Pearson r (log10, genes expressed in both, TPM>0.5)

Run in the RNAseq conda env (numpy/pandas/scipy).
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SYN1_LEN = 1_078_809
DBG = os.path.join(ROOT, "Syn1_Transcriptomics/ONT/ONT_Processing/depth_bedgraph")
THR = 0.5   # TPM threshold for "expressed in both" when correlating


def load_track(fname):
    a = np.zeros(SYN1_LEN, dtype=np.float64)
    with open(os.path.join(DBG, fname)) as fh:
        for line in fh:
            _, s, _e, v = line.split()
            a[int(s)] = float(v)
    return a


def ont_sense_tpm(genes, plus, minus):
    d = np.empty(len(genes))
    for i, (_, r) in enumerate(genes.iterrows()):
        s, e = int(r.start0), int(r.end0)
        arr = plus if r.strand == "+" else minus
        d[i] = arr[s:e].mean() if e > s else 0.0
    return d / d.sum() * 1e6


# ---------- Syn1 ----------
syn1 = pd.read_csv(os.path.join(ROOT, "Syn1_Transcriptomics/Gene_TPM/syn1_Illumina_PacBio_TPM_profiles.csv"))
syn1 = syn1[syn1.feature.isin(["gene", "pseudogene"])].reset_index(drop=True)

tracks = {r: (load_track(f"syn1.ONT.{r}.plus.bedGraph"),
              load_track(f"syn1.ONT.{r}.minus.bedGraph")) for r in ("rep1", "rep2", "merged")}

out1 = pd.DataFrame({
    "locus_tag": syn1.locus_tag,
    "gene_name": syn1.gene_name,
    "gene_len": syn1.gene_len,
    "strand": syn1.strand,
    "Illumina_TPM": syn1.avg_sense_TPM,
    "PacBio_TPM": syn1.PacBio_sense_TPM,
    "ONT1_TPM": ont_sense_tpm(syn1, *tracks["rep1"]),
    "ONT2_TPM": ont_sense_tpm(syn1, *tracks["rep2"]),
    "ONTmerged_TPM": ont_sense_tpm(syn1, *tracks["merged"]),
})
out1.to_csv(os.path.join(HERE, "platform_TPM_syn1.tsv"), sep="\t", index=False, float_format="%.4f")

# ---------- Syn3A ----------
syn3 = pd.read_csv(os.path.join(ROOT, "Syn3A_Transcriptomics/Gene_TPM/syn3A_TPM_Illumina_ONT.tsv"), sep="\t")
out3 = pd.DataFrame({
    "locus_tag": syn3.locus_tag,
    "gene_name": syn3.gene_name,
    "gene_len": syn3.gene_len,
    "Illumina_TPM": syn3.Illumina_sense_TPM,
    "ONT_TPM": syn3.ONT_sense_TPM,
})
out3.to_csv(os.path.join(HERE, "platform_TPM_syn3A.tsv"), sep="\t", index=False, float_format="%.4f")


# ---------- pairwise correlations ----------
def rlog(a, b):
    m = (a > THR) & (b > THR)
    if m.sum() < 3:
        return float("nan"), int(m.sum())
    return pearsonr(np.log10(a[m]), np.log10(b[m]))[0], int(m.sum())


lines = ["Pairwise Pearson r on log10 sense TPM (genes with both TPM > %.1f)\n" % THR,
         "== Syn1 (n genes = %d) ==" % len(out1)]
cols1 = [("Illumina", "Illumina_TPM"), ("PacBio", "PacBio_TPM"),
         ("ONT1", "ONT1_TPM"), ("ONT2", "ONT2_TPM")]
for i in range(len(cols1)):
    for j in range(i + 1, len(cols1)):
        (na, ca), (nb, cb) = cols1[i], cols1[j]
        r, n = rlog(out1[ca].values, out1[cb].values)
        lines.append(f"  {na:9s} vs {nb:9s}: r = {r:.3f}  (n = {n})")
lines.append("  [ONTmerged vs Illumina]: r = %.3f  (n = %d)" % rlog(out1.ONTmerged_TPM.values, out1.Illumina_TPM.values))
lines.append("")
lines.append("== Syn3A (n genes = %d) ==" % len(out3))
r, n = rlog(out3.ONT_TPM.values, out3.Illumina_TPM.values)
lines.append(f"  ONT       vs Illumina : r = {r:.3f}  (n = {n})")

txt = "\n".join(lines) + "\n"
with open(os.path.join(HERE, "compute_platform_TPM.txt"), "w") as fh:
    fh.write(txt)
print(txt)
print("wrote platform_TPM_syn1.tsv (%d genes), platform_TPM_syn3A.tsv (%d genes)" % (len(out1), len(out3)))
