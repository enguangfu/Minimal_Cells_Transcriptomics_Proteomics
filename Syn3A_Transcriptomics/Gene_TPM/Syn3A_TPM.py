#!/usr/bin/env python
# coding: utf-8

# ## Gene-level Transcriptomics
# 
# Illumina of syn3A

# In[1]:


import pandas as pd
import numpy as np
import os
import pysam
from typing import Dict, List, Tuple, Optional


# In[2]:


GENES_GFF    = "../../Genomes_Input/syn3a_genome.gff3"

# Gene annotation parsing — syn3A annotates 3 pseudogenes (0051, 0546, 0602)
# only under the `pseudogene` feature type, so include both.
PRIMARY_FEATURES = {"gene", "pseudogene"}
INCLUDE_CDS_FALLBACK = False

def parse_gff_attributes(attr: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in attr.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out

def read_genes_gff(path: str) -> pd.DataFrame:
    rows = []
    with open(path, "r") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start1, end1, score, strand, phase, attrs = parts
            if feature == "CDS" and not INCLUDE_CDS_FALLBACK:
                continue
            if feature not in PRIMARY_FEATURES and not (INCLUDE_CDS_FALLBACK and feature == "CDS"):
                continue
            a = parse_gff_attributes(attrs)
            locus_tag = a.get("locus_tag", "")
            gene_name = a.get("gene", "")
            rna_type = a.get("rna_type", "")
            gene_product = a.get("product", "")
            s1 = int(start1); e1 = int(end1)
            start0 = s1 - 1
            end0 = e1
            rows.append((chrom, feature, start0, end0, strand, locus_tag, gene_name, rna_type, gene_product))
    df = pd.DataFrame(rows, columns=["chrom", "feature", "start0", "end0", "strand", "locus_tag", "gene_name", "rna_type", "gene_product"])
    df = df.sort_values(["chrom", "start0", "end0"]).reset_index(drop=True)
    return df

GENES = read_genes_gff(GENES_GFF)

GENES.head()


# ### Helpers + syn3A gene table

import pandas as pd
import numpy as np
from pathlib import Path

HOME_DIR = "../.."  # project root; this script lives at Syn3A_Transcriptomics/Gene_TPM/

# Use TPM = (avg_depth_over_gene_body) / (sum_of_avg_depths_over_all_genes) * 1e6.
# Length normalisation is implicit in the avg-depth measurement (depth = read-bases / length).


def load_depth_array_from_bedgraph(path):
    """
    Load a bedGraph and return (chrom -> per-base depth array, chrom -> prefix sum).
    bedGraph format: chrom, start0, end0, depth — each row covers start0..end0-1.
    """
    df = pd.read_csv(path, sep="\t", header=None, names=["chrom", "start0", "end0", "depth"])
    df["start0"] = df["start0"].astype(int)
    df["end0"]   = df["end0"].astype(int)
    df["depth"]  = df["depth"].astype(float)

    chrom_arrays = {}
    chrom_prefix = {}
    for chrom, sub in df.groupby("chrom", sort=False):
        max_end = int(sub["end0"].max())
        arr = np.zeros(max_end, dtype=np.float32)
        for s, e, d in zip(sub["start0"].to_numpy(dtype=int),
                           sub["end0"].to_numpy(dtype=int),
                           sub["depth"].to_numpy(dtype=np.float32)):
            arr[s:e] = d
        chrom_arrays[chrom] = arr
        chrom_prefix[chrom] = np.concatenate([[0.0], np.cumsum(arr, dtype=np.float64)])
    return chrom_arrays, chrom_prefix


def compute_strand_tpm(genes_df: pd.DataFrame,
                       plus_bedgraph: str, minus_bedgraph: str,
                       prefix: str) -> pd.DataFrame:
    """Compute sense / antisense average depth + TPM per gene from two strand
    bedGraphs. Adds `{prefix}_sense_avg_depth`, `..._antisense_avg_depth`,
    `..._sense_TPM`, `..._antisense_TPM` columns to a copy of `genes_df`
    and returns the new dataframe."""
    out = genes_df.copy().reset_index(drop=True)
    _, plus_prefix  = load_depth_array_from_bedgraph(plus_bedgraph)
    _, minus_prefix = load_depth_array_from_bedgraph(minus_bedgraph)

    sense_vals     = np.zeros(len(out), dtype=np.float64)
    antisense_vals = np.zeros(len(out), dtype=np.float64)
    for i, row in out.iterrows():
        c  = row["chrom"]
        s0 = int(row["start0"])
        e0 = int(row["end0"])
        glen = max(1, e0 - s0)
        pm = (plus_prefix[c][e0]  - plus_prefix[c][s0])  / glen
        mm = (minus_prefix[c][e0] - minus_prefix[c][s0]) / glen
        if row["strand"] == "+":
            sense_vals[i], antisense_vals[i] = pm, mm
        else:
            sense_vals[i], antisense_vals[i] = mm, pm
    out[f"{prefix}_sense_avg_depth"]     = sense_vals
    out[f"{prefix}_antisense_avg_depth"] = antisense_vals
    denom = sense_vals.sum() + antisense_vals.sum()
    if denom > 0:
        out[f"{prefix}_sense_TPM"]     = sense_vals     / denom * 1e6
        out[f"{prefix}_antisense_TPM"] = antisense_vals / denom * 1e6
    else:
        out[f"{prefix}_sense_TPM"]     = 0.0
        out[f"{prefix}_antisense_TPM"] = 0.0
    return out


# ### syn3A Illumina TPM
#
# Inputs: ./depth_bedgraph/syn3A_rep1.{plus,minus}.bedGraph (from 02 alignment).
# Output: per-gene sense and antisense TPM, prefix "Illumina".

ILLUMINA_PREFIX = "Illumina"
ILLUMINA_PLUS  = "../Illumina/Illumina_Processing/depth_bedgraph/syn3A_rep1.plus.bedGraph"
ILLUMINA_MINUS = "../Illumina/Illumina_Processing/depth_bedgraph/syn3A_rep1.minus.bedGraph"

syn3a_genes_base = GENES.copy().reset_index(drop=True)
syn3a_genes_base["start0"]   = syn3a_genes_base["start0"].astype(int)
syn3a_genes_base["end0"]     = syn3a_genes_base["end0"].astype(int)
syn3a_genes_base["gene_len"] = syn3a_genes_base["end0"] - syn3a_genes_base["start0"]

illumina_df = compute_strand_tpm(syn3a_genes_base, ILLUMINA_PLUS, ILLUMINA_MINUS,
                                 prefix=ILLUMINA_PREFIX)
print(f"syn3A Illumina TPM computed for {len(illumina_df)} loci")


# ### syn3A ONT TPM
#
# Inputs: ../ONT_Processing/depth_bedgraph/syn3A.ONT.rep1.{plus,minus}.bedGraph.
# Output: per-gene sense and antisense TPM, prefix "ONT".

ONT_PLUS  = HOME_DIR + "/Syn3A_Transcriptomics/ONT/ONT_Processing/depth_bedgraph/syn3A.ONT.rep1.plus.bedGraph"
ONT_MINUS = HOME_DIR + "/Syn3A_Transcriptomics/ONT/ONT_Processing/depth_bedgraph/syn3A.ONT.rep1.minus.bedGraph"
ont_df = compute_strand_tpm(syn3a_genes_base, ONT_PLUS, ONT_MINUS, prefix="ONT")
print(f"syn3A ONT      TPM computed for {len(ont_df)} loci")

# In syn3a.gff3 the `product=` attribute lives on CDS rows, not gene rows.
# Build a locus_tag -> product map from CDS lines and merge it in.
SYN3A_GFF = HOME_DIR + "/Genomes_Input/syn3a_genome.gff3"
_cds_product = {}
with open(SYN3A_GFF, "r") as _fh:
    for _line in _fh:
        if not _line or _line.startswith("#"):
            continue
        _parts = _line.rstrip("\n").split("\t")
        if len(_parts) < 9 or _parts[2] != "CDS":
            continue
        _attrs = parse_gff_attributes(_parts[8])
        _lt = _attrs.get("locus_tag", "")
        if _lt and _lt not in _cds_product:
            _cds_product[_lt] = _attrs.get("product", "")


# ### Merge Illumina + ONT into one syn3A TPM table

ill_keep = ["locus_tag", "Illumina_sense_avg_depth", "Illumina_antisense_avg_depth",
            "Illumina_sense_TPM", "Illumina_antisense_TPM"]
ont_keep = ["locus_tag", "ONT_sense_avg_depth", "ONT_antisense_avg_depth",
            "ONT_sense_TPM", "ONT_antisense_TPM"]

merged = (syn3a_genes_base[["locus_tag", "gene_name", "gene_product", "chrom",
                            "start0", "end0", "strand", "gene_len"]]
          .merge(illumina_df[ill_keep], on="locus_tag", how="left")
          .merge(ont_df[ont_keep],     on="locus_tag", how="left"))
merged["gene_product"] = merged["locus_tag"].map(_cds_product).fillna(merged["gene_product"])

OUT_TSV = "./syn3A_TPM_Illumina_ONT.tsv"
merged.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.4f")
print(f"Exported to {OUT_TSV}  ({len(merged)} loci)")


# ### Correlation: Illumina vs ONT sense TPM (gene level)
#
# Use sense TPM on log10 with a small pseudocount to handle zeros. Report
# Pearson + Spearman; save a scatter plot. Restrict to genes with non-zero
# sense TPM in at least one assay (so pseudo-count-dominated comparisons
# don't drag correlations).

from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt

EPS = 1e-2  # pseudo-count for log10
def _corr_pair(x, y):
    keep = (x > 0) | (y > 0)
    xk, yk = x[keep], y[keep]
    lx = np.log10(xk + EPS)
    ly = np.log10(yk + EPS)
    pr, pp = pearsonr(lx, ly)
    sr, sp = spearmanr(xk, yk)
    return len(xk), pr, sr

n_io, pr_io, sr_io = _corr_pair(
    merged["Illumina_sense_TPM"].to_numpy(),
    merged["ONT_sense_TPM"].to_numpy(),
)
print(f"\nIllumina vs ONT (syn3A, sense TPM, n={n_io} non-zero loci):")
print(f"  Pearson  (log10) : r = {pr_io:.3f}")
print(f"  Spearman (rank)  : rho = {sr_io:.3f}")


# ### Correlation: our Illumina TPM vs Palsson's reported Illumina TPM
#
# Palsson GSM6204176_3A.csv has columns Geneid, syn3A_3A, Illumina_TPM.
# We compare against the reported Illumina_TPM (per gene).

PALSSON_CSV = "./Processed_TPM_Palsson/GSM6204176_3A.csv"
palsson = pd.read_csv(PALSSON_CSV)
palsson = palsson.rename(columns={"Geneid": "locus_tag", "Illumina_TPM": "Palsson_Illumina_TPM"})
print(f"\nPalsson Illumina TPM rows : {len(palsson)}")

cmp_p = merged.merge(palsson[["locus_tag", "Palsson_Illumina_TPM"]],
                     on="locus_tag", how="inner")
print(f"Overlap with our Illumina TPM : {len(cmp_p)} loci")
n_pal, pr_pal, sr_pal = _corr_pair(
    cmp_p["Illumina_sense_TPM"].to_numpy(),
    cmp_p["Palsson_Illumina_TPM"].to_numpy(),
)
print(f"Our Illumina vs Palsson Illumina (sense TPM, n={n_pal} non-zero):")
print(f"  Pearson  (log10) : r = {pr_pal:.3f}")
print(f"  Spearman (rank)  : rho = {sr_pal:.3f}")


# ### Scatter PDFs (one figure per correlation)

def _scatter_one(x, y, xlabel, ylabel, title, n, pr, sr, out_pdf: str):
    fig, ax = plt.subplots(figsize=(5.5, 5))
    keep = (x > 0) | (y > 0)
    xk = x[keep] + EPS
    yk = y[keep] + EPS
    ax.scatter(xk, yk, s=10, alpha=0.4, edgecolor="none")
    lo = min(xk.min(), yk.min()) * 0.5
    hi = max(xk.max(), yk.max()) * 2
    ax.plot([lo, hi], [lo, hi], ls=":", color="grey", lw=0.8)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_title(f"{title}\nn={n}, Pearson(log10) r={pr:.2f}, Spearman ρ={sr:.2f}",
                 fontsize=11)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_pdf, dpi=200)
    plt.close(fig)


OUT_PDF_IO  = "./syn3A_TPM_correlation_Illumina_vs_ONT.pdf"
OUT_PDF_PAL = "./syn3A_TPM_correlation_Illumina_vs_Palsson.pdf"

_scatter_one(
    merged["Illumina_sense_TPM"].to_numpy(),
    merged["ONT_sense_TPM"].to_numpy(),
    "syn3A Illumina sense TPM", "syn3A ONT sense TPM",
    "Illumina vs ONT (syn3A)", n_io, pr_io, sr_io,
    OUT_PDF_IO,
)
_scatter_one(
    cmp_p["Illumina_sense_TPM"].to_numpy(),
    cmp_p["Palsson_Illumina_TPM"].to_numpy(),
    "Our Illumina sense TPM", "Palsson reported Illumina TPM",
    "Our Illumina vs Palsson Illumina", n_pal, pr_pal, sr_pal,
    OUT_PDF_PAL,
)

print(f"\nWrote: {OUT_TSV}")
print(f"Wrote: {OUT_PDF_IO}")
print(f"Wrote: {OUT_PDF_PAL}")


