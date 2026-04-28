"""
Combined syn1 vs syn3A comparison at RNA (TPM) and protein (iPM) levels.

Produces:
  - distribution + correlation PDFs for both RNA and protein layers
  - operon summary txt files + copied operon plots for extreme-FC groups
  - one combined CSV: syn1_vs_syn3a_RNA_protein_combined.csv
"""

import os
import shutil
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

# ── Inputs ───────────────────────────────────────────────────────────────────
SYN1_TPM_CSV          = "../Transcriptomics_Quantification/syn1_Illumina_PacBio_TPM_profiles.csv"
SYN3A_TPM_TSV         = "../Transcriptomics_Quantification/syn3a_ONT_TPM_profiles.tsv"
SYN3A_ILLUMINA_TPM_CSV = "../Syn3A_Illumina/Processed_TPM_Palsson/GSM6204176_3A.csv"
SYN1_PTN_CSV   = "../Proteomics_Quantification/syn1_proteomics_localization_2026.csv"
SYN3A_PTN_CSV  = "../Proteomics_Quantification/syn3a_proteomics_summary_2026.csv"
OPERON_COV     = "../Operon_Annotation_Visualization/gene_operon_coverage.tsv"
OPERON_PLOTS   = "../Operon_Annotation_Visualization/operon_plots"

SYN1_TPM_COL   = "avg_sense_TPM"   # Illumina merged across bio samples
SYN3A_TPM_COL  = "ONT_sense_TPM"
SYN3A_ILLUMINA_TPM_COL = "Illumina_TPM"

OUT_CODING     = "syn1_vs_syn3a_RNA_protein.csv"
OUT_NONCODING  = "syn1_vs_syn3A_noncoding_RNA.csv"

CODING_RNA_TYPES    = {"mRNA", "pseudo"}
NONCODING_RNA_TYPES = {"rRNA", "tRNA", "ncRNA", "tmRNA"}

# ── Operon coverage (used by both layers) ────────────────────────────────────
operon_cov = pd.read_csv(OPERON_COV, sep="\t")
operon_cov["locus_num"] = operon_cov["locus_tag"].str.extract(r"(\d+)$").astype(int)


def _extract_locus_num(s: pd.Series) -> pd.Series:
    return s.str.extract(r"(\d+)$").astype(int)


def _fc_plots(result: pd.DataFrame, value_col: str, fc_col: str, log10fc_col: str,
              tag: str, ylabel_prefix: str) -> None:
    """Make FC + log10FC distribution histograms and the syn1-vs-syn3a scatter."""
    fc = result[fc_col]
    fc_stats = fc.describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])
    print(f"\n[{tag}] fold change statistics (n={int(fc_stats['count'])}):")
    print(f"  mean={fc_stats['mean']:.3f}  median={fc_stats['50%']:.3f}  std={fc_stats['std']:.3f}")
    print(f"  5th={fc_stats['5%']:.3f}  25th={fc_stats['25%']:.3f}  "
          f"75th={fc_stats['75%']:.3f}  95th={fc_stats['95%']:.3f}")

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.hist(fc, bins=50, color="steelblue", edgecolor="white", linewidth=0.4)
    ax.axvline(10,  color="red",    linestyle="--", linewidth=1.2, label="FC = 10")
    ax.axvline(0.1, color="orange", linestyle="--", linewidth=1.2, label="FC = 0.1")
    ax.set_xlabel(f"{tag} fold change  (syn3A / syn1)", fontsize=11)
    ax.set_ylabel(f"Number of {ylabel_prefix}", fontsize=11)
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(f"{tag}_fold_change_distribution.pdf")
    plt.close(fig)

    log10fc = result[log10fc_col]
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.hist(log10fc, bins=50, color="steelblue", edgecolor="white", linewidth=0.4)
    ax.axvline(np.log10(10),  color="red",    linestyle="--", linewidth=1.2, label="log10(10) = 1")
    ax.axvline(np.log10(0.1), color="orange", linestyle="--", linewidth=1.2, label="log10(0.1) = −1")
    ax.set_xlabel(f"log10({tag} fold change)  (syn3A / syn1)", fontsize=12)
    ax.set_ylabel(f"Number of {ylabel_prefix}", fontsize=11)
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(f"{tag}_log10FC_distribution.pdf")
    plt.close(fig)

    # Restrict to rows where both syn1 and syn3a have positive, non-NaN values
    # (outer merges may produce NaN for genes detected in only one condition).
    plot_df = result[
        result[f"{value_col}_syn1"].notna() & (result[f"{value_col}_syn1"] > 0) &
        result[f"{value_col}_syn3a"].notna() & (result[f"{value_col}_syn3a"] > 0)
    ]
    x = plot_df[f"{value_col}_syn1"]
    y = plot_df[f"{value_col}_syn3a"]
    pearson_r, pearson_p   = stats.pearsonr(np.log10(x), np.log10(y))
    spearman_r, spearman_p = stats.spearmanr(x, y)

    abundant  = plot_df[f"{value_col}_syn1"] > 100
    mask_high = (plot_df[fc_col] > 10)  & abundant
    mask_low  = (plot_df[fc_col] < 0.1) & abundant
    mask_mid  = ~mask_high & ~mask_low

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(x[mask_mid],  y[mask_mid],  s=12, alpha=0.55, color="steelblue", edgecolors="none", label="other")
    ax.scatter(x[mask_high], y[mask_high], s=25, alpha=0.90, color="red",       edgecolors="none",
               label=f"FC > 10  & syn1 {tag} > 100")
    ax.scatter(x[mask_low],  y[mask_low],  s=25, alpha=0.90, color="orange",    edgecolors="none",
               label=f"FC < 0.1 & syn1 {tag} > 100")

    lims = [min(x.min(), y.min()) * 0.5, max(x.max(), y.max()) * 2]
    ax.plot(lims, lims,                    color="black",  linewidth=0.8, linestyle="--")
    ax.plot(lims, [v * 10  for v in lims], color="red",    linewidth=0.6, linestyle=":")
    ax.plot(lims, [v * 0.1 for v in lims], color="orange", linewidth=0.6, linestyle=":")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel(f"syn1  {tag} (log)",  fontsize=12)
    ax.set_ylabel(f"syn3A {tag} (log)",  fontsize=12)
    ax.text(0.04, 0.96,
            f"Pearson r = {pearson_r:.3f}  (log10)\nSpearman r = {spearman_r:.3f}\nn = {len(plot_df)}",
            transform=ax.transAxes, fontsize=12, verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))
    ax.legend(fontsize=10, loc="lower right")
    plt.tight_layout()
    fig.savefig(f"{tag}_correlation_syn1_vs_syn3a.pdf")
    plt.close(fig)
    print(f"[{tag}] Pearson r (log10) = {pearson_r:.3f} (p={pearson_p:.2e}); "
          f"Spearman r = {spearman_r:.3f} (p={spearman_p:.2e})")


def _operon_group_report(df_extreme: pd.DataFrame, folder: str, label: str,
                         tag: str, value_col: str, fc_col: str) -> None:
    """Write per-operon summary txt and copy operon plots (no per-group CSVs)."""
    os.makedirs(folder, exist_ok=True)

    all_operons = set()
    for ops_str in df_extreme["sense_covering_ops"].dropna():
        for op in str(ops_str).split(","):
            all_operons.add(op.strip())

    lines = [f"# {tag} {label.upper()} genes  (n={len(df_extreme)})\n"]
    for opid in sorted(all_operons):
        op_genes = operon_cov[operon_cov["sense_covering_ops"].str.contains(opid, na=False)]
        extreme_in_op = df_extreme[df_extreme["sense_covering_ops"].str.contains(opid, na=False)]
        lines.append(f"\n## Operon {opid}  —  {len(op_genes)} total genes, "
                     f"{len(extreme_in_op)} with extreme FC\n")
        lines.append(f"{'locus_tag':<18} {'gene_name':<14} {f'{tag}_syn1':>12} "
                     f"{f'{tag}_syn3a':>12} {'FC':>8}\n")
        lines.append("-" * 70 + "\n")
        for _, row in op_genes.iterrows():
            r = df_extreme[df_extreme["locus_num"] == row["locus_num"]]
            if r.empty:
                # fall back to all-genes table to fill values
                r = result_lookup.get(tag, pd.DataFrame())
                r = r[r["locus_num"] == row["locus_num"]] if not r.empty else r
            if r.empty:
                v1 = v3 = fc_str = "N/A"
            else:
                rr = r.iloc[0]
                v1 = f"{rr[f'{value_col}_syn1']:.1f}"
                v3 = f"{rr[f'{value_col}_syn3a']:.1f}"
                fc_str = f"{rr[fc_col]:.3f}"
            marker = " <<<" if row["locus_num"] in extreme_in_op["locus_num"].values else ""
            gname  = str(row["gene_name"]) if pd.notna(row.get("gene_name")) else ""
            lines.append(f"{row['locus_tag']:<18} {gname:<14} {v1:>12} {v3:>12} {fc_str:>8}{marker}\n")

    with open(f"{folder}/operon_summary_{label}.txt", "w") as fh:
        fh.writelines(lines)

    copied, missing = [], []
    for opid in sorted(all_operons):
        for suffix in ("", "_wdepth"):
            fname = f"{opid}{suffix}.pdf"
            src = os.path.join(OPERON_PLOTS, fname)
            if os.path.exists(src):
                shutil.copy2(src, os.path.join(folder, fname))
                copied.append(fname)
            else:
                missing.append(fname)
    print(f"[{tag} {label}] operon_summary written, {len(copied)} plot(s) copied to {folder}/")
    if missing:
        print(f"  missing: {missing}")


# ─────────────────────────────────────────────────────────────────────────────
# Load RNA tables once. rna_type comes from the syn1 TPM CSV; we use syn1's
# classification to split syn3A genes (paired by locus_num) into coding vs
# non-coding subsets. TPMs are then renormalized within each subset so the
# fold change reflects abundance among that gene class only.
# ─────────────────────────────────────────────────────────────────────────────
syn1_tpm_full  = pd.read_csv(SYN1_TPM_CSV)
syn3a_tpm_full = pd.read_csv(SYN3A_TPM_TSV, sep="\t")
syn3a_ill_full = pd.read_csv(SYN3A_ILLUMINA_TPM_CSV).rename(
    columns={"Geneid": "locus_tag", SYN3A_ILLUMINA_TPM_COL: "TPM_mean"})

# Mask zero TPM as NaN so missing measurements stay NaN (not inf) through the
# division. Rows with NaN TPM are kept so we can report the gene with empty
# values rather than dropping it.
syn1_tpm_full[SYN1_TPM_COL]   = syn1_tpm_full[SYN1_TPM_COL].mask(syn1_tpm_full[SYN1_TPM_COL]   <= 0)
syn3a_tpm_full[SYN3A_TPM_COL] = syn3a_tpm_full[SYN3A_TPM_COL].mask(syn3a_tpm_full[SYN3A_TPM_COL] <= 0)
syn3a_ill_full["TPM_mean"]    = syn3a_ill_full["TPM_mean"].mask(syn3a_ill_full["TPM_mean"]    <= 0)

# Add locus_num to each
for _df in (syn1_tpm_full, syn3a_tpm_full, syn3a_ill_full):
    _df["locus_num"] = _extract_locus_num(_df["locus_tag"])

# rename TPM columns to the common name
syn1_tpm_full  = syn1_tpm_full.rename(columns={SYN1_TPM_COL:  "TPM_mean"})
syn3a_tpm_full = syn3a_tpm_full.rename(columns={SYN3A_TPM_COL: "TPM_mean"})

# Reference metadata: locus_num -> rna_type / gene_name / gene_product, taken
# from syn1 (authoritative since the locus_num scheme originates there).
syn1_meta = syn1_tpm_full[["locus_num", "rna_type", "gene_name", "gene_product"]].copy()
locus_to_rnatype = dict(zip(syn1_meta["locus_num"], syn1_meta["rna_type"]))


def _build_rna_combined(rna_types: set, label: str):
    """Filter all three TPM tables by rna_type, renormalize TPM (sum=1e6) within
    the subset on detected genes, rank within subset, and outer-merge into one
    wide table. Genes detected in only one source still appear with NaN.
    Returns (rna_ont_df, rna_illumina_df)."""
    keep_loci = {ln for ln, rt in locus_to_rnatype.items() if rt in rna_types}

    def _subset_renorm(df: pd.DataFrame) -> pd.DataFrame:
        sub = df[df["locus_num"].isin(keep_loci)].copy()
        total = sub["TPM_mean"].sum(skipna=True)
        if total > 0:
            sub["TPM_mean"] = sub["TPM_mean"] / total * 1e6
        return sub

    s1   = _subset_renorm(syn1_tpm_full)
    s3o  = _subset_renorm(syn3a_tpm_full)
    s3i  = _subset_renorm(syn3a_ill_full)

    # Rank only on detected (non-NaN) values; NaN gets NaN rank automatically.
    s1["rank_syn1_TPM"]            = s1["TPM_mean"].rank(ascending=False, method="min")
    s3o["rank_syn3a_TPM"]          = s3o["TPM_mean"].rank(ascending=False, method="min")
    s3i["rank_syn3a_TPM_Illumina"] = s3i["TPM_mean"].rank(ascending=False, method="min")

    # ONT side — outer merge so syn1-only or syn3A-only loci still appear
    rna_ont = pd.merge(
        s1[["locus_num", "locus_tag", "TPM_mean", "rank_syn1_TPM"]],
        s3o[["locus_num", "locus_tag", "gene_name", "gene_product",
             "TPM_mean", "rank_syn3a_TPM"]],
        on="locus_num", how="outer", suffixes=("_syn1", "_syn3a"),
    )
    # Coalesce annotation: prefer syn3a, fall back to syn1 metadata (always present)
    rna_ont = rna_ont.merge(syn1_meta, on="locus_num", how="left", suffixes=("", "_syn1meta"))
    rna_ont["gene_name"]    = rna_ont["gene_name"].fillna(rna_ont["gene_name_syn1meta"])
    rna_ont["gene_product"] = rna_ont["gene_product"].fillna(rna_ont["gene_product_syn1meta"])
    rna_ont = rna_ont.drop(columns=["gene_name_syn1meta", "gene_product_syn1meta"])

    rna_ont["TPM_fold_change"] = rna_ont["TPM_mean_syn3a"] / rna_ont["TPM_mean_syn1"]
    rna_ont["TPM_log10FC"]     = np.log10(rna_ont["TPM_fold_change"].replace(0, np.nan))
    rna_ont = rna_ont.merge(
        operon_cov[["locus_num", "sense_covering_ops", "antisense_covering_ops", "coverage_type"]],
        on="locus_num", how="left",
    )

    # Illumina side — outer merge as well
    rna_ill = pd.merge(
        s1[["locus_num", "TPM_mean", "rank_syn1_TPM"]],
        s3i[["locus_num", "TPM_mean", "rank_syn3a_TPM_Illumina"]],
        on="locus_num", how="outer", suffixes=("_syn1", "_syn3a_Illumina"),
    )
    rna_ill["TPM_fold_change_Illumina"] = rna_ill["TPM_mean_syn3a_Illumina"] / rna_ill["TPM_mean_syn1"]
    rna_ill["TPM_log10FC_Illumina"]     = np.log10(rna_ill["TPM_fold_change_Illumina"].replace(0, np.nan))
    rna_ill = rna_ill.merge(
        operon_cov[["locus_num", "sense_covering_ops"]],
        on="locus_num", how="left",
    )

    print(f"[{label}] kept loci={len(keep_loci)}; "
          f"syn1 detected={s1['TPM_mean'].notna().sum()}, "
          f"syn3a_ONT detected={s3o['TPM_mean'].notna().sum()}, "
          f"syn3a_Illumina detected={s3i['TPM_mean'].notna().sum()}; "
          f"ONT rows={len(rna_ont)}, Illumina rows={len(rna_ill)}")
    return rna_ont, rna_ill


rna,          rna_illumina      = _build_rna_combined(CODING_RNA_TYPES,    "coding")
rna_nc,       rna_illumina_nc   = _build_rna_combined(NONCODING_RNA_TYPES, "non-coding")

# Distribution + scatter PDFs only for the coding subset (FC on rRNA/tRNA is unreliable)
_fc_plots(rna, value_col="TPM_mean", fc_col="TPM_fold_change",
          log10fc_col="TPM_log10FC", tag="TPM", ylabel_prefix="genes")
_rna_ill_plot = rna_illumina.rename(columns={
    "TPM_mean_syn3a_Illumina":  "TPM_mean_syn3a",
    "TPM_fold_change_Illumina": "TPM_fold_change",
    "TPM_log10FC_Illumina":     "TPM_log10FC",
})
_fc_plots(_rna_ill_plot, value_col="TPM_mean", fc_col="TPM_fold_change",
          log10fc_col="TPM_log10FC", tag="TPM_Illumina", ylabel_prefix="genes")


# ─────────────────────────────────────────────────────────────────────────────
# Protein layer (iPM)
# ─────────────────────────────────────────────────────────────────────────────
syn1_ptn  = pd.read_csv(SYN1_PTN_CSV)
syn3a_ptn = pd.read_csv(SYN3A_PTN_CSV)

syn1_ptn  = syn1_ptn[syn1_ptn["iPM_mean"].notna()   & (syn1_ptn["iPM_mean"]   > 0)].copy()
syn3a_ptn = syn3a_ptn[syn3a_ptn["iPM_mean"].notna() & (syn3a_ptn["iPM_mean"] > 0)].copy()

n_syn1_ptn, n_syn3a_ptn = len(syn1_ptn), len(syn3a_ptn)
syn1_ptn["rank_syn1_iPM"]   = syn1_ptn["iPM_mean"].rank(ascending=False, method="min")
syn3a_ptn["rank_syn3a_iPM"] = syn3a_ptn["iPM_mean"].rank(ascending=False, method="min")

syn1_ptn["locus_num"]  = _extract_locus_num(syn1_ptn["locus_tag"])
syn3a_ptn["locus_num"] = _extract_locus_num(syn3a_ptn["locus_tag"])

# Outer merge so syn1-only proteins (genes deleted in syn3A) still appear with
# their iPM_mean_syn1 / rank_syn1_iPM populated and syn3A side as NaN.
ptn = pd.merge(
    syn1_ptn[["locus_num", "locus_tag", "iPM_mean", "rank_syn1_iPM"]],
    syn3a_ptn[["locus_num", "locus_tag", "iPM_mean", "rank_syn3a_iPM"]],
    on="locus_num", how="outer", suffixes=("_syn1", "_syn3a"),
)
ptn["iPM_fold_change"] = ptn["iPM_mean_syn3a"] / ptn["iPM_mean_syn1"]
ptn["iPM_log10FC"]     = np.log10(ptn["iPM_fold_change"])

ptn = ptn.merge(
    operon_cov[["locus_num", "sense_covering_ops", "antisense_covering_ops", "coverage_type"]],
    on="locus_num", how="left",
)
print(f"[iPM] syn1 detected={n_syn1_ptn}, syn3a detected={n_syn3a_ptn}, merged={len(ptn)}")

_fc_plots(ptn, value_col="iPM_mean", fc_col="iPM_fold_change",
          log10fc_col="iPM_log10FC", tag="iPM", ylabel_prefix="proteins")


# ─────────────────────────────────────────────────────────────────────────────
# Operon summaries for extreme-FC groups (no per-group CSVs written)
# ─────────────────────────────────────────────────────────────────────────────
result_lookup = {"TPM": rna, "iPM": ptn}

_operon_group_report(rna[rna["TPM_fold_change"] > 10].sort_values("TPM_fold_change", ascending=False),
                     "RNA_upgrade",   "FC_gt10",  "TPM", "TPM_mean", "TPM_fold_change")
_operon_group_report(rna[rna["TPM_fold_change"] < 0.1].sort_values("TPM_fold_change"),
                     "RNA_downgrade", "FC_lt0.1", "TPM", "TPM_mean", "TPM_fold_change")
_operon_group_report(ptn[ptn["iPM_fold_change"] > 10].sort_values("iPM_fold_change", ascending=False),
                     "protein_upgrade",   "FC_gt10",  "iPM", "iPM_mean", "iPM_fold_change")
_operon_group_report(ptn[ptn["iPM_fold_change"] < 0.1].sort_values("iPM_fold_change"),
                     "protein_downgrade", "FC_lt0.1", "iPM", "iPM_mean", "iPM_fold_change")


# ─────────────────────────────────────────────────────────────────────────────
# Combined CSV outputs
#   - syn1_vs_syn3a_RNA_protein.csv      : coding (mRNA + pseudo) + protein layer
#   - syn1_vs_syn3A_noncoding_RNA.csv    : rRNA / tRNA / ncRNA / tmRNA (RNA only)
# ─────────────────────────────────────────────────────────────────────────────
def _build_combined(rna_df, rna_ill_df, ptn_df=None):
    rna_df = rna_df.copy()
    rna_df["rna_type"] = rna_df["locus_num"].map(locus_to_rnatype)
    rna_keep = rna_df.rename(columns={
        "locus_tag_syn1":   "locus_syn1",
        "locus_tag_syn3a":  "locus_syn3a",
        "TPM_mean_syn3a":   "TPM_mean_syn3A_ONT",
        "rank_syn3a_TPM":   "rank_syn3a_TPM_ONT",
        "TPM_fold_change":  "TPM_fold_change_ONT",
    })[[
        "locus_num", "locus_syn1", "locus_syn3a", "gene_name", "gene_product", "rna_type",
        "TPM_mean_syn1", "rank_syn1_TPM",
        "TPM_mean_syn3A_ONT", "rank_syn3a_TPM_ONT",
        "TPM_fold_change_ONT",
        "sense_covering_ops",
    ]]
    rna_ill_keep = rna_ill_df[[
        "locus_num",
        "TPM_mean_syn3a_Illumina", "rank_syn3a_TPM_Illumina",
        "TPM_fold_change_Illumina",
    ]]
    out = rna_keep.merge(rna_ill_keep, on="locus_num", how="outer")
    if ptn_df is not None:
        ptn_keep = ptn_df[[
            "locus_num",
            "iPM_mean_syn1", "rank_syn1_iPM",
            "iPM_mean_syn3a", "rank_syn3a_iPM",
            "iPM_fold_change",
        ]]
        out = out.merge(ptn_keep, on="locus_num", how="outer")
    return out.sort_values("locus_num")


coding = _build_combined(rna, rna_illumina, ptn)
coding_cols = [
    "locus_syn1", "locus_syn3a", "gene_name", "gene_product", "rna_type",
    "TPM_mean_syn1", "rank_syn1_TPM",
    "TPM_mean_syn3A_ONT", "rank_syn3a_TPM_ONT",
    "TPM_fold_change_ONT",
    "TPM_mean_syn3a_Illumina", "rank_syn3a_TPM_Illumina",
    "TPM_fold_change_Illumina",
    "iPM_mean_syn1", "rank_syn1_iPM",
    "iPM_mean_syn3a", "rank_syn3a_iPM",
    "iPM_fold_change",
    "sense_covering_ops",
]
coding[coding_cols].to_csv(OUT_CODING, index=False)
print(f"\nSaved: {OUT_CODING}  ({len(coding)} coding genes)")

noncoding = _build_combined(rna_nc, rna_illumina_nc, ptn_df=None)
noncoding_cols = [
    "locus_syn1", "locus_syn3a", "gene_name", "gene_product", "rna_type",
    "TPM_mean_syn1", "rank_syn1_TPM",
    "TPM_mean_syn3A_ONT", "rank_syn3a_TPM_ONT",
    "TPM_fold_change_ONT",
    "TPM_mean_syn3a_Illumina", "rank_syn3a_TPM_Illumina",
    "TPM_fold_change_Illumina",
    "sense_covering_ops",
]
noncoding[noncoding_cols].to_csv(OUT_NONCODING, index=False)
print(f"Saved: {OUT_NONCODING}  ({len(noncoding)} non-coding genes)")


# ─────────────────────────────────────────────────────────────────────────────
# Extreme fold-change subsets (protein-coding only) — printed, not saved
#   - TPM extremes: require BOTH ONT and Illumina FC to agree (>10 or <0.1)
#   - iPM extremes: FC > 10 or < 0.1
# ─────────────────────────────────────────────────────────────────────────────
TPM_PRINT_COLS = ["locus_syn1", "locus_syn3a", "gene_name", "gene_product",
                  "TPM_mean_syn1", "TPM_mean_syn3A_ONT", "TPM_fold_change_ONT",
                  "TPM_mean_syn3a_Illumina", "TPM_fold_change_Illumina"]
IPM_PRINT_COLS = ["locus_syn1", "locus_syn3a", "gene_name", "gene_product",
                  "iPM_mean_syn1", "iPM_mean_syn3a", "iPM_fold_change"]

fc_ont = coding["TPM_fold_change_ONT"]
fc_ill = coding["TPM_fold_change_Illumina"]
tpm_up   = coding[(fc_ont > 10)  & (fc_ill > 10)].sort_values("TPM_fold_change_ONT", ascending=False)
tpm_down = coding[(fc_ont < 0.1) & (fc_ill < 0.1)].sort_values("TPM_fold_change_ONT")

print(f"\n=== TPM extreme FC (both ONT & Illumina agree) — UP (n={len(tpm_up)}) ===")
print(tpm_up[TPM_PRINT_COLS].to_string(index=False))
print(f"\n=== TPM extreme FC (both ONT & Illumina agree) — DOWN (n={len(tpm_down)}) ===")
print(tpm_down[TPM_PRINT_COLS].to_string(index=False))

fc_ipm = coding["iPM_fold_change"]
ipm_up   = coding[fc_ipm > 10].sort_values("iPM_fold_change", ascending=False)
ipm_down = coding[fc_ipm < 0.1].sort_values("iPM_fold_change")

print(f"\n=== iPM extreme FC — UP (n={len(ipm_up)}) ===")
print(ipm_up[IPM_PRINT_COLS].to_string(index=False))
print(f"\n=== iPM extreme FC — DOWN (n={len(ipm_down)}) ===")
print(ipm_down[IPM_PRINT_COLS].to_string(index=False))
