#!/usr/bin/env python
# coding: utf-8

# ## Gene-level Transcriptomics and Proteomics
# 
# PacBio and Illumina of Syn1

# In[1]:


import pandas as pd
import numpy as np
import os
import pysam
from typing import Dict, List, Tuple, Optional


# ## TPM and correlation of Illumina
# 
# Check the consistency of three Illumina reps Reads
#  
# SRR35996296 and SRR35996297 are the same RNA sample; technical replicates to each other.  
# **Technical replicates with high Pearson r of 0.98.**
# 
# 
# SRR35996298 is a different RNA biological sample.  
# **Person r of 0.92 and 0.94 with the other two reps**
# 
# We don't know if PacBio RNA sample the same or not with Illumia.
# 
# Thus, merge three Illumina to get the most representative.
# 
# | Sample Name | Read 1 OR Read 2 | # of Reads | Strand rel to RNA | Time | File Name |
# |---|---|---|---|---|---|
# | Syn1_enr | R1 | 1,085,803 | Reverse-Complementary | 2023.07 | SRR35996298_1 |
# | Syn1_enr | R2 | 1,085,803 | Same | 2023.07 | SRR35996298_2 |
# | 95A | R1 | 510,910 | Reverse-Complementary | 2023.09 | SRR35996297_1 |
# | 95A | R2 | 510,910 | Same | 2023.09 | SRR35996297_2 |
# | 95B | R1 | 504,891 | Reverse-Complementary | 2023.09 | SRR35996296_1 |
# | 95B | R2 | 504,891 | Same | 2023.09 | SRR35996296_2 |

# In[2]:


GENES_GFF    = "../Genomes_Input/syn1.genes.gff3"

# Gene annotation parsing
PRIMARY_FEATURES = {"gene"}
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


# ### Correlation between Illuminas

# In[3]:


import pandas as pd
import numpy as np
from pathlib import Path


HOME_DIR = ".."
SEQDEPTH_FOLDER = HOME_DIR + "/Illumina_Processing/depth_bedgraph"

SAMPLES = ["SRR35996296", "SRR35996297", "SRR35996298"]

genes = GENES.copy().reset_index(drop=True)
genes["start0"] = genes["start0"].astype(int)
genes["end0"] = genes["end0"].astype(int)
genes["gene_len"] = genes["end0"] - genes["start0"]


def load_depth_array_from_bedgraph(path):
    """
    Load a bedGraph file and build per-chromosome depth arrays + prefix sums.
    bedGraph format: chrom, start0, end0, depth (each row spans start0..end0-1).
    """
    df = pd.read_csv(
        path, sep="\t", header=None, names=["chrom", "start0", "end0", "depth"]
    )
    df["start0"] = df["start0"].astype(int)
    df["end0"] = df["end0"].astype(int)
    df["depth"] = df["depth"].astype(float)

    chrom_arrays = {}
    chrom_prefix = {}

    for chrom, sub in df.groupby("chrom", sort=False):
        max_end = int(sub["end0"].max())
        arr = np.zeros(max_end, dtype=np.float32)

        # --- FIX: fill the entire span, not just the start position ---
        starts = sub["start0"].to_numpy(dtype=int)
        ends = sub["end0"].to_numpy(dtype=int)
        depths = sub["depth"].to_numpy(dtype=np.float32)
        for s, e, d in zip(starts, ends, depths):
            arr[s:e] = d

        chrom_arrays[chrom] = arr
        chrom_prefix[chrom] = np.concatenate(
            [[0.0], np.cumsum(arr, dtype=np.float64)]
        )

    return chrom_arrays, chrom_prefix


for sample in SAMPLES:
    plus_file = SEQDEPTH_FOLDER + f"/{sample}.plus.bedGraph"
    minus_file = SEQDEPTH_FOLDER + f"/{sample}.minus.bedGraph"

    _, plus_prefix = load_depth_array_from_bedgraph(plus_file)
    _, minus_prefix = load_depth_array_from_bedgraph(minus_file)

    sense_vals = np.zeros(len(genes), dtype=np.float64)
    antisense_vals = np.zeros(len(genes), dtype=np.float64)

    for i, row in genes.iterrows():
        chrom = row["chrom"]
        start0 = int(row["start0"])
        end0 = int(row["end0"])
        glen = int(row["gene_len"])
        strand = row["strand"]

        plus_mean = (
            plus_prefix[chrom][end0] - plus_prefix[chrom][start0]
        ) / glen
        minus_mean = (
            minus_prefix[chrom][end0] - minus_prefix[chrom][start0]
        ) / glen

        if strand == "+":
            sense_vals[i] = plus_mean
            antisense_vals[i] = minus_mean
        else:
            sense_vals[i] = minus_mean
            antisense_vals[i] = plus_mean

    genes[f"{sample}_sense_avg_depth"] = sense_vals
    genes[f"{sample}_antisense_avg_depth"] = antisense_vals

    # TPM: shared denominator across both strands
    denom = sense_vals.sum() + antisense_vals.sum()
    if denom == 0:
        genes[f"{sample}_sense_TPM"] = 0.0
        genes[f"{sample}_antisense_TPM"] = 0.0
    else:
        genes[f"{sample}_sense_TPM"] = sense_vals / denom * 1e6
        genes[f"{sample}_antisense_TPM"] = antisense_vals / denom * 1e6

cols = (
    ["locus_tag", "strand"]
    + [f"{s}_sense_TPM" for s in SAMPLES]
    + [f"{s}_antisense_TPM" for s in SAMPLES]
)
genes[cols].head()


# In[4]:


import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations

tpm_cols = [f"{sample}_sense_TPM" for sample in SAMPLES]

for a, b in combinations(tpm_cols, 2):
    x = np.log10(genes[a].astype(float) + 1.0)
    y = np.log10(genes[b].astype(float) + 1.0)

    r = np.corrcoef(x, y)[0, 1]

    plt.figure(figsize=(6, 6))
    plt.scatter(x, y, s=10, alpha=0.5)
    lims = [
        min(x.min(), y.min()),
        max(x.max(), y.max())
    ]
    plt.plot(lims, lims, linestyle="--")
    plt.xlim(lims)
    plt.ylim(lims)
    plt.xlabel(f"log10(TPM+1)\n{a}")
    plt.ylabel(f"log10(TPM+1)\n{b}")
    plt.title(f"{a} vs {b}\nPearson r = {r:.4f}")
    plt.tight_layout()
    # plt.show()


# ### Merge into one Illumina Profile

# In[5]:


SAMPLE_TO_BIO = {
    "SRR35996296": "sample_95",
    "SRR35996297": "sample_95",
    "SRR35996298": "sample_enr",
}

OUTPUT_TSV = HOME_DIR + "/Transcriptomics_Quantification/syn1_illumina_TPM_profiles.tsv"

 #---------------------------------------------------------------------------
# Two-step averaging: tech reps within bio sample first, then bio samples
# ---------------------------------------------------------------------------
# Step 1: average technical replicates within each biological sample
bio_groups = {}
for sample in SAMPLES:
    bio = SAMPLE_TO_BIO[sample]
    bio_groups.setdefault(bio, []).append(sample)
 
for bio, reps in bio_groups.items():
    for strand_type in ["sense", "antisense"]:
        cols = [f"{s}_{strand_type}_TPM" for s in reps]
        genes[f"{bio}_{strand_type}_TPM"] = genes[cols].mean(axis=1)
 
# Step 2: average across biological samples (equal weight to sample1 and sample3)
bio_samples = sorted(bio_groups.keys())
for strand_type in ["sense", "antisense"]:
    bio_cols = [f"{bio}_{strand_type}_TPM" for bio in bio_samples]
    genes[f"avg_{strand_type}_TPM"] = genes[bio_cols].mean(axis=1)
 
# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
export_cols = (
    ["locus_tag", "chrom", "start0", "end0", "strand", "gene_len"]
    # Per-sample TPM
    + [f"{s}_sense_TPM" for s in SAMPLES]
    + [f"{s}_antisense_TPM" for s in SAMPLES]
    # Per-biological-sample averaged TPM
    + [f"{bio}_sense_TPM" for bio in bio_samples]
    + [f"{bio}_antisense_TPM" for bio in bio_samples]
    # Final averaged TPM
    + ["avg_sense_TPM", "avg_antisense_TPM"]
)

genes[export_cols].to_csv(OUTPUT_TSV, sep="\t", index=False, float_format="%.4f")
print(f"Exported to {OUTPUT_TSV}")
print(f"  Genes: {len(genes)}")
print(f"  Columns: {len(export_cols)}")
print(genes[["locus_tag", "avg_sense_TPM", "avg_antisense_TPM"]].head(10))


# ## Calculate TPM of PacBio
# 
# 

# In[6]:


from pathlib import Path
import numpy as np
import pandas as pd

PACBIO_DEPTH_FOLDER = HOME_DIR + "/PacBio_Processing/depth_bedgraph"
PACBIO_PREFIX = "syn1.PacBio.FLNC.HQ"

pacbio_genes = GENES.copy().reset_index(drop=True)
pacbio_genes["start0"] = pacbio_genes["start0"].astype(int)
pacbio_genes["end0"] = pacbio_genes["end0"].astype(int)
pacbio_genes["gene_len"] = pacbio_genes["end0"] - pacbio_genes["start0"]

plus_file  = PACBIO_DEPTH_FOLDER + f"/{PACBIO_PREFIX}.plus.bedGraph"
minus_file = PACBIO_DEPTH_FOLDER + f"/{PACBIO_PREFIX}.minus.bedGraph"

_, plus_prefix  = load_depth_array_from_bedgraph(plus_file)
_, minus_prefix = load_depth_array_from_bedgraph(minus_file)

sense_vals     = np.zeros(len(pacbio_genes), dtype=np.float64)
antisense_vals = np.zeros(len(pacbio_genes), dtype=np.float64)

for i, row in pacbio_genes.iterrows():
    chrom  = row["chrom"]
    start0 = int(row["start0"])
    end0   = int(row["end0"])
    glen   = int(row["gene_len"])
    strand = row["strand"]

    plus_mean  = (plus_prefix[chrom][end0]  - plus_prefix[chrom][start0])  / glen
    minus_mean = (minus_prefix[chrom][end0] - minus_prefix[chrom][start0]) / glen

    if strand == "+":
        sense_vals[i]     = plus_mean
        antisense_vals[i] = minus_mean
    else:
        sense_vals[i]     = minus_mean
        antisense_vals[i] = plus_mean

pacbio_genes["PacBio_sense_avg_depth"]     = sense_vals
pacbio_genes["PacBio_antisense_avg_depth"] = antisense_vals

# TPM: normalise by gene length (already done via avg depth) then scale to 1e6
denom = sense_vals.sum() + antisense_vals.sum()
if denom == 0:
    pacbio_genes["PacBio_sense_TPM"]     = 0.0
    pacbio_genes["PacBio_antisense_TPM"] = 0.0
else:
    pacbio_genes["PacBio_sense_TPM"]     = sense_vals     / denom * 1e6
    pacbio_genes["PacBio_antisense_TPM"] = antisense_vals / denom * 1e6

OUTPUT_PACBIO_TSV = HOME_DIR + "/Transcriptomics_Quantification/syn1_pacbio_TPM_profiles.tsv"

export_cols = [
    "locus_tag", "chrom", "start0", "end0", "strand", "gene_len",
    "PacBio_sense_avg_depth", "PacBio_antisense_avg_depth",
    "PacBio_sense_TPM", "PacBio_antisense_TPM",
]
pacbio_genes[export_cols].to_csv(OUTPUT_PACBIO_TSV, sep="\t", index=False, float_format="%.4f")
print(f"Exported to {OUTPUT_PACBIO_TSV}")
print(f"  Genes: {len(pacbio_genes)}")
pacbio_genes[["locus_tag", "strand", "PacBio_sense_TPM", "PacBio_antisense_TPM"]].head(10)


# In[7]:


genes['PacBio_sense_TPM'] = pacbio_genes['PacBio_sense_TPM']
genes["PacBio_antisense_TPM"] = pacbio_genes["PacBio_antisense_TPM"]


# In[8]:


genes.to_csv("./syn1_Illumina_PacBio_TPM_profiles.csv",index=False, float_format="%.4f")


# ## Correlate PacBio w. Illumina
# 
# - Relatively low correlation found:  
# Pearson r = 0.6 for log10 between PacBio and Illumina.   
# No significant read length biase observed.  
# No significant read count biase found.
# 
# 

# In[9]:


import numpy as np
import matplotlib.pyplot as plt

low_threshold = .5

FONT_SIZE = 14
# ── Scatter plot (log10 TPM) ────────────────────────────────────────────────
x_all = np.log10(genes["avg_sense_TPM"].astype(float))
y_all = np.log10(pacbio_genes["PacBio_sense_TPM"].astype(float))

mask = (x_all >= np.log10(low_threshold)) & (y_all >= np.log10(low_threshold))
x = x_all[mask]
y = y_all[mask]

r = np.corrcoef(x, y)[0, 1]

plt.figure(figsize=(6, 6))
plt.scatter(x, y, s=10, alpha=0.5)
lims = [min(x.min(), y.min()), max(x.max(), y.max())]
plt.plot(lims, lims, linestyle="--")
plt.xlim(lims)
plt.ylim(lims)
plt.xlabel("log10(TPM)\nIllumina", fontsize=FONT_SIZE)
plt.ylabel("log10(TPM)\nPacBio", fontsize=FONT_SIZE)
plt.title(f"Illumina vs PacBio sense TPM", fontsize=FONT_SIZE)
plt.text(0.05, 0.95, f"n = {mask.sum()} genes w. TPM ≥ {low_threshold}", transform=plt.gca().transAxes, fontsize=FONT_SIZE-2, verticalalignment="top")
plt.text(0.05, 0.90, f"Pearson r = {r:.4f}", transform=plt.gca().transAxes, fontsize=FONT_SIZE-2, verticalalignment="top")
plt.tight_layout()
plt.savefig("./Illumina_vs_PacBio_sense_TPM_scatter.pdf", dpi=300)
# plt.show()

# ── Pearson correlation tables ────────────────────────────────────────────────
corr_df = pd.DataFrame({
    "Illumina_avg_sense_TPM": genes["avg_sense_TPM"].values[mask],
    "PacBio_sense_TPM": pacbio_genes["PacBio_sense_TPM"].values[mask],
})
print(f"Genes passing threshold (Illumina TPM ≥ {low_threshold}): {mask.sum()}")
print("Pearson correlation of sense TPM:")
print(corr_df.corr(method="pearson"))

log_corr_df = np.log10(corr_df)
print("Pearson correlation of log10(TPM):")
print(log_corr_df.corr(method="pearson"))

# ── MA plot ───────────────────────────────────────────────────────────────────
illumina = genes["avg_sense_TPM"].astype(float).values[mask]
pacbio   = pacbio_genes["PacBio_sense_TPM"].astype(float).values[mask]

M = np.log10((pacbio) / (illumina))
A = 0.5 * np.log10((pacbio) * (illumina))

plt.figure(figsize=(6, 6))
plt.scatter(A, M, s=10, alpha=0.5)
plt.axhline(0, color="red", linewidth=1, linestyle="--")
plt.xlabel("A  =  0.5 × log₁₀(PacBio × Illumina)", fontsize=FONT_SIZE)
plt.ylabel("M  =  log₁₀(PacBio / Illumina)", fontsize=FONT_SIZE)
plt.text(0.05, 0.05, f"n = {mask.sum()} genes w. TPM ≥ {low_threshold}", transform=plt.gca().transAxes, fontsize=FONT_SIZE-2, verticalalignment="top")
plt.title(f"MA Plot: PacBio vs Illumina sense TPM", fontsize=FONT_SIZE)
plt.tight_layout()
plt.savefig("./Illumina_vs_PacBio_sense_TPM_MA_plot.pdf", dpi=300)
# plt.show()


# In[10]:


import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Apply same threshold mask
illumina_tpm = genes["avg_sense_TPM"].astype(float).values
pacbio_tpm   = pacbio_genes["PacBio_sense_TPM"].astype(float).values
gene_len     = genes["gene_len"].astype(float).values

mask = (illumina_tpm >= low_threshold) & (pacbio_tpm >= low_threshold)

illumina_m = illumina_tpm[mask]
pacbio_m   = pacbio_tpm[mask]
gene_len_m = gene_len[mask]

log10_ratio = np.log10((pacbio_m) / (illumina_m))

# Linear regression
slope, intercept, r, p, _ = stats.linregress(gene_len_m, log10_ratio)

x_line = np.array([gene_len_m.min(), gene_len_m.max()])
y_line = slope * x_line + intercept

plt.figure(figsize=(6, 6))
plt.scatter(gene_len_m, log10_ratio, s=10, alpha=0.5)
plt.plot(x_line, y_line, color="red", linewidth=1.5,
         label=f"slope={slope:.4f}, r={r:.3f}, p={p:.2e}")
plt.axhline(0, color="gray", linewidth=0.8, linestyle="--")
plt.xlabel("Gene length (bp)", fontsize=FONT_SIZE)
plt.ylabel("log₁₀(PacBio TPM / Illumina TPM)", fontsize=FONT_SIZE)
plt.title(f"Length bias: PacBio vs Illumina ", fontsize=FONT_SIZE)
plt.text(0.4, 0.05, f"n = {mask.sum()} genes with TPM ≥ {low_threshold}", transform=plt.gca().transAxes, fontsize=FONT_SIZE-2, verticalalignment="top")
plt.legend(fontsize=9)
plt.tight_layout()
plt.savefig("./Illumina_vs_PacBio_length_bias.pdf", dpi=300)
# plt.show()

print(f"Linear regression: slope={slope:.5f}, Pearson r={r:.4f}, p={p:.2e}")
print("Negative slope → shorter genes enriched in PacBio (unexpected for full-length Iso-Seq).")
print("Positive slope → longer genes enriched in PacBio (consistent with Iso-Seq length preference).")


# <!-- ## Draw Sequencing Depth -->

# ## Syn3A
#
# For now, only one set of ONT for Syn3A.

# In[11]:


SYN3A_GFF              = HOME_DIR + "/Genomes_Input/syn3a_genome.gff3"
SYN3A_DEPTH_FOLDER     = HOME_DIR + "/ONT_Processing/depth_bedgraph"
SYN3A_PREFIX           = "syn3A.ONT.rep1"

syn3a_genes = read_genes_gff(SYN3A_GFF).reset_index(drop=True)
syn3a_genes["start0"]   = syn3a_genes["start0"].astype(int)
syn3a_genes["end0"]     = syn3a_genes["end0"].astype(int)
syn3a_genes["gene_len"] = syn3a_genes["end0"] - syn3a_genes["start0"]

plus_file  = SYN3A_DEPTH_FOLDER + f"/{SYN3A_PREFIX}.plus.bedGraph"
minus_file = SYN3A_DEPTH_FOLDER + f"/{SYN3A_PREFIX}.minus.bedGraph"

_, plus_prefix  = load_depth_array_from_bedgraph(plus_file)
_, minus_prefix = load_depth_array_from_bedgraph(minus_file)

sense_vals     = np.zeros(len(syn3a_genes), dtype=np.float64)
antisense_vals = np.zeros(len(syn3a_genes), dtype=np.float64)

for i, row in syn3a_genes.iterrows():
    chrom  = row["chrom"]
    start0 = int(row["start0"])
    end0   = int(row["end0"])
    glen   = int(row["gene_len"])
    strand = row["strand"]

    plus_mean  = (plus_prefix[chrom][end0]  - plus_prefix[chrom][start0])  / glen
    minus_mean = (minus_prefix[chrom][end0] - minus_prefix[chrom][start0]) / glen

    if strand == "+":
        sense_vals[i]     = plus_mean
        antisense_vals[i] = minus_mean
    else:
        sense_vals[i]     = minus_mean
        antisense_vals[i] = plus_mean

syn3a_genes["ONT_sense_avg_depth"]     = sense_vals
syn3a_genes["ONT_antisense_avg_depth"] = antisense_vals

denom = sense_vals.sum() + antisense_vals.sum()
if denom == 0:
    syn3a_genes["ONT_sense_TPM"]     = 0.0
    syn3a_genes["ONT_antisense_TPM"] = 0.0
else:
    syn3a_genes["ONT_sense_TPM"]     = sense_vals     / denom * 1e6
    syn3a_genes["ONT_antisense_TPM"] = antisense_vals / denom * 1e6

OUTPUT_SYN3A_TSV = HOME_DIR + "/Transcriptomics_Quantification/syn3a_ONT_TPM_profiles.tsv"
export_cols = [
    "locus_tag", "chrom", "start0", "end0", "strand", "gene_len",
    "ONT_sense_avg_depth", "ONT_antisense_avg_depth",
    "ONT_sense_TPM", "ONT_antisense_TPM",
]
syn3a_genes[export_cols].to_csv(OUTPUT_SYN3A_TSV, sep="\t", index=False, float_format="%.4f")
print(f"Exported to {OUTPUT_SYN3A_TSV}")
print(f"  Genes: {len(syn3a_genes)}")
syn3a_genes[["locus_tag", "strand", "ONT_sense_TPM", "ONT_antisense_TPM"]].head(10)


