import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Load ──────────────────────────────────────────────────────────────────────
SYN1_CSV  = "../Proteomics/syn1_proteomics_localization_2026.csv"
SYN3A_CSV = "../Proteomics/syn3a_proteomics_summary_2026.csv"

syn1  = pd.read_csv(SYN1_CSV)
syn3a = pd.read_csv(SYN3A_CSV)

# Keep only rows with measured iPM (non-zero, non-NaN)
syn1  = syn1[syn1["iPM_mean"].notna()  & (syn1["iPM_mean"]  > 0)].copy()
syn3a = syn3a[syn3a["iPM_mean"].notna() & (syn3a["iPM_mean"] > 0)].copy()

# Rank within each organism's own detected proteome
n_syn1  = len(syn1)
n_syn3a = len(syn3a)
syn1["rank_syn1"]   = syn1["iPM_mean"].rank(ascending=False, method="min").astype(int)
syn3a["rank_syn3a"] = syn3a["iPM_mean"].rank(ascending=False, method="min").astype(int)

# ── Match genes by locus-tag number ──────────────────────────────────────────
# syn1  locus_tag: MMSYN1_0001    → key 1
# syn3a locus_tag: JCVISYN3A_0001 → key 1
syn1["locus_num"]  = syn1["locus_tag"].str.extract(r"(\d+)$").astype(int)
syn3a["locus_num"] = syn3a["locus_tag"].str.extract(r"(\d+)$").astype(int)

# ── Merge ─────────────────────────────────────────────────────────────────────
merged = pd.merge(
    syn1[["locus_num", "locus_tag", "gene_name", "iPM_mean",
          "rank_syn1", "ptn_copy_number", "ptn_localization"]],
    syn3a[["locus_num", "locus_tag", "gene_name", "iPM_mean",
           "rank_syn3a", "copy_number_2026", "localization"]],
    on="locus_num",
    suffixes=("_syn1", "_syn3a"),
)

# ── Fold change ───────────────────────────────────────────────────────────────
merged["iPM_fold_change"] = merged["iPM_mean_syn3a"] / merged["iPM_mean_syn1"]
merged["iPM_log10FC"]     = np.log10(merged["iPM_fold_change"])

# ── Tidy output columns ───────────────────────────────────────────────────────
result = merged[[
    "locus_num",
    "locus_tag_syn1", "locus_tag_syn3a",
    "gene_name_syn1", "gene_name_syn3a",
    "iPM_mean_syn1", "rank_syn1", "iPM_mean_syn3a", "rank_syn3a",
    "iPM_fold_change", "iPM_log10FC",
    "ptn_copy_number", "copy_number_2026",
    "ptn_localization", "localization",
]].sort_values("locus_num").reset_index(drop=True)

# ── Attach syn1 operon membership early (used for reporting and later for plots)
import os, shutil
OPERON_COV   = "../Operon_Annotation_Visualization/gene_operon_coverage.tsv"
OPERON_PLOTS = "../Operon_Annotation_Visualization/operon_plots"
operon_cov = pd.read_csv(OPERON_COV, sep="\t")
operon_cov["locus_num"] = operon_cov["locus_tag"].str.extract(r"(\d+)$").astype(int)
result = result.merge(
    operon_cov[["locus_num", "sense_covering_ops", "antisense_covering_ops", "coverage_type"]],
    on="locus_num", how="left",
)

# ── Save CSV ──────────────────────────────────────────────────────────────────
OUT_CSV = "syn1_vs_syn3a_protein_fold_change.csv"
result.to_csv(OUT_CSV, index=False)
print(f"Saved: {OUT_CSV}  ({len(result)} genes)")

# ── Report extreme fold changes ───────────────────────────────────────────────
result["rank_syn1_str"]  = result["rank_syn1"].astype(str)  + f"/{n_syn1}"
result["rank_syn3a_str"] = result["rank_syn3a"].astype(str) + f"/{n_syn3a}"

high = result[result["iPM_fold_change"] > 10].sort_values("iPM_fold_change", ascending=False)
low  = result[result["iPM_fold_change"] < 0.1].sort_values("iPM_fold_change")

print(f"\nFC > 10  ({len(high)} genes):")
print(high[["locus_num", "gene_name_syn1", "gene_name_syn3a",
            "iPM_mean_syn1", "rank_syn1_str",
            "iPM_mean_syn3a", "rank_syn3a_str",
            "iPM_fold_change", "sense_covering_ops"]].to_string(index=False))

print(f"\nFC < 0.1  ({len(low)} genes):")
print(low[["locus_num", "gene_name_syn1", "gene_name_syn3a",
           "iPM_mean_syn1", "rank_syn1_str",
           "iPM_mean_syn3a", "rank_syn3a_str",
           "iPM_fold_change", "sense_covering_ops"]].to_string(index=False))

# ── Plot: distribution of iPM_fold_change ────────────────────────────────────
fc = result["iPM_fold_change"]
fc_stats = fc.describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])
print(f"\niPM fold change statistics (n={int(fc_stats['count'])}):")
print(f"  mean={fc_stats['mean']:.3f}  median={fc_stats['50%']:.3f}  std={fc_stats['std']:.3f}")
print(f"  5th={fc_stats['5%']:.3f}  25th={fc_stats['25%']:.3f}  75th={fc_stats['75%']:.3f}  95th={fc_stats['95%']:.3f}")
print(f"  min={fc_stats['min']:.3f}  max={fc_stats['max']:.3f}")

fig, ax = plt.subplots(figsize=(4, 4))
ax.hist(fc, bins=50, color="steelblue", edgecolor="white", linewidth=0.4)
ax.axvline(10,  color="red",    linestyle="--", linewidth=1.2, label="FC = 10")
ax.axvline(0.1, color="orange", linestyle="--", linewidth=1.2, label="FC = 0.1")
ax.set_xlabel("iPM fold change  (syn3A / syn1)", fontsize=11)
ax.set_ylabel("Number of proteins", fontsize=11)
ax.legend(fontsize=9)
plt.tight_layout()
fig.savefig("iPM_fold_change_distribution.pdf")
plt.close(fig)
print("Saved: iPM_fold_change_distribution.pdf")

# ── Plot: distribution of iPM_log10FC ────────────────────────────────────────
log10fc = result["iPM_log10FC"]
log10fc_stats = log10fc.describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])
print(f"\niPM log10FC statistics (n={int(log10fc_stats['count'])}):")
print(f"  mean={log10fc_stats['mean']:.3f}  median={log10fc_stats['50%']:.3f}  std={log10fc_stats['std']:.3f}")
print(f"  5th={log10fc_stats['5%']:.3f}  25th={log10fc_stats['25%']:.3f}  75th={log10fc_stats['75%']:.3f}  95th={log10fc_stats['95%']:.3f}")
print(f"  min={log10fc_stats['min']:.3f}  max={log10fc_stats['max']:.3f}")

fig, ax = plt.subplots(figsize=(4, 4))
ax.hist(log10fc, bins=50, color="steelblue", edgecolor="white", linewidth=0.4)
ax.axvline( np.log10(10),  color="red",    linestyle="--", linewidth=1.2, label="log10(10) = 1")
ax.axvline( np.log10(0.1), color="orange", linestyle="--", linewidth=1.2, label="log10(0.1) = −1")
ax.set_xlabel("log10(iPM fold change)  (syn3A / syn1)", fontsize=12)
ax.set_ylabel("Number of proteins", fontsize=11)
ax.legend(fontsize=9)
plt.tight_layout()
fig.savefig("iPM_log10FC_distribution.pdf")
plt.close(fig)
print("Saved: iPM_log10FC_distribution.pdf")

# ── Plot: iPM correlation (syn1 vs syn3A) ─────────────────────────────────────
from scipy import stats

x = result["iPM_mean_syn1"]
y = result["iPM_mean_syn3a"]

pearson_r, pearson_p = stats.pearsonr(np.log10(x), np.log10(y))
spearman_r, spearman_p = stats.spearmanr(x, y)

# Colour groups: highlight only if FC crosses boundary AND syn1 iPM > 100
abundant = result["iPM_mean_syn1"] > 100
mask_high = (result["iPM_fold_change"] > 10)  & abundant
mask_low  = (result["iPM_fold_change"] < 0.1) & abundant
mask_mid  = ~mask_high & ~mask_low   # everything else stays blue

fig, ax = plt.subplots(figsize=(8, 5))
ax.scatter(x[mask_mid],  y[mask_mid],  s=12, alpha=0.55, color="steelblue", edgecolors="none", label="other")
ax.scatter(x[mask_high], y[mask_high], s=25, alpha=0.90, color="red",       edgecolors="none", label="FC > 10  & syn1 iPM > 100")
ax.scatter(x[mask_low],  y[mask_low],  s=25, alpha=0.90, color="orange",    edgecolors="none", label="FC < 0.1 & syn1 iPM > 100")

lims = [min(x.min(), y.min()) * 0.5, max(x.max(), y.max()) * 2]
ax.plot(lims, lims,                    color="black",  linewidth=0.8, linestyle="--")
ax.plot(lims, [v * 10  for v in lims], color="red",    linewidth=0.6, linestyle=":")
ax.plot(lims, [v * 0.1 for v in lims], color="orange", linewidth=0.6, linestyle=":")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_xlabel("syn1  iPM mean (log scale)", fontsize=12)
ax.set_ylabel("syn3A  iPM mean (log scale)", fontsize=12)
ax.text(0.04, 0.96,
        f"Pearson r = {pearson_r:.3f}  (log10)\nSpearman r = {spearman_r:.3f}\nn = {len(result)}",
        transform=ax.transAxes, fontsize=12, verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))
ax.legend(fontsize=10, loc="lower right")
plt.tight_layout()
fig.savefig("iPM_correlation_syn1_vs_syn3a.pdf")
plt.close(fig)
print("Saved: iPM_correlation_syn1_vs_syn3a.pdf")
print(f"\nCorrelation (log10 iPM):  Pearson r = {pearson_r:.3f}  (p = {pearson_p:.2e})")
print(f"Correlation (raw iPM):    Spearman r = {spearman_r:.3f}  (p = {spearman_p:.2e})")

# ── Operon plots and summaries for extreme FC genes ──────────────────────────

def _write_group(df_extreme: pd.DataFrame, folder: str, label: str):
    os.makedirs(folder, exist_ok=True)

    # Columns to save
    cols = ["locus_num", "locus_tag_syn1", "locus_tag_syn3a",
            "gene_name_syn1", "gene_name_syn3a",
            "iPM_mean_syn1", "iPM_mean_syn3a", "iPM_fold_change", "iPM_log10FC",
            "ptn_localization", "sense_covering_ops", "antisense_covering_ops", "coverage_type"]
    df_extreme[cols].to_csv(f"{folder}/genes_{label}.csv", index=False)

    # Collect all unique sense operons and write a per-operon summary
    all_operons = set()
    for ops_str in df_extreme["sense_covering_ops"].dropna():
        for op in str(ops_str).split(","):
            all_operons.add(op.strip())

    summary_lines = [f"# {label.upper()} proteins  (n={len(df_extreme)})\n"]
    for opid in sorted(all_operons):
        # genes belonging to this operon in the full operon coverage table
        op_genes = operon_cov[operon_cov["sense_covering_ops"].str.contains(opid, na=False)]
        # which of those are in the extreme group
        extreme_in_op = df_extreme[df_extreme["sense_covering_ops"].str.contains(opid, na=False)]
        summary_lines.append(f"\n## Operon {opid}  —  {len(op_genes)} total genes, "
                             f"{len(extreme_in_op)} with extreme FC\n")
        summary_lines.append(f"{'locus_tag':<18} {'gene_name':<14} {'iPM_syn1':>10} "
                             f"{'iPM_syn3a':>10} {'FC':>8}\n")
        summary_lines.append("-" * 66 + "\n")
        for _, row in op_genes.iterrows():
            # find FC info if available
            fc_row = result[result["locus_num"] == row["locus_num"]]
            if fc_row.empty:
                fc_str = ipm1_str = ipm3_str = "N/A"
            else:
                r = fc_row.iloc[0]
                ipm1_str  = f"{r['iPM_mean_syn1']:.1f}"
                ipm3_str  = f"{r['iPM_mean_syn3a']:.1f}"
                fc_str    = f"{r['iPM_fold_change']:.3f}"
            marker = " <<<" if row["locus_num"] in extreme_in_op["locus_num"].values else ""
            gname  = str(row["gene_name"]) if pd.notna(row["gene_name"]) else ""
            summary_lines.append(f"{row['locus_tag']:<18} {gname:<14} {ipm1_str:>10} "
                                 f"{ipm3_str:>10} {fc_str:>8}{marker}\n")

    report_path = f"{folder}/operon_summary_{label}.txt"
    with open(report_path, "w") as fh:
        fh.writelines(summary_lines)
    print(f"Saved: {folder}/genes_{label}.csv  +  operon_summary_{label}.txt")

    # Copy operon plots
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
    print(f"  Copied {len(copied)} plot(s) to {folder}/")
    if missing:
        print(f"  Not found: {missing}")

# Attach operon cols to high/low subsets
high_ops = result[result["iPM_fold_change"] > 10].sort_values("iPM_fold_change", ascending=False)
low_ops  = result[result["iPM_fold_change"] < 0.1].sort_values("iPM_fold_change")

_write_group(high_ops, "upgrade",   "FC_gt10")
_write_group(low_ops,  "downgrade", "FC_lt0.1")
