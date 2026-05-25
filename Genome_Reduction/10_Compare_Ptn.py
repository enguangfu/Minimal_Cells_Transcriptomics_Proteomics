"""
Protein (iPM) syn1 vs syn3A comparison — the proteome counterpart of 09.

Split out of 09_Compare_RNA_Protein.py: 09 owns the RNA/TPM story (and still
writes the joined TPM+iPM CSV); this script owns the protein-specific figures and
outlier report. It READS 09's combined table
(Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv), which already carries the
mean-normalized protein columns relIPM_syn1, relIPM_syn3a, iPM_fold_change,
iPM_abs_change. Outputs go to the same Compare_RNA_Protein/ folder.

Normalization / metrics are identical to 09 (see its docstring): rel = value /
per-gene mean over detected proteins; fold change = ratio, absolute change =
rel_syn3A - rel_syn1.

Proteome reallocation story plots (iPM) mirror 09's mRNA pair, using the curated
Secondary/Tertiary function annotation and the deletion-corrected (retained-pool)
normalization. NOTE: ribosomal-protein iPM is digestion-biased — read Ribosome /
Translation at the protein level with caution.

Outputs (Compare_RNA_Protein/)
  - iPM_fold_change_distribution.pdf, iPM_log10FC_distribution.pdf, iPM_correlation_syn1_vs_syn3a.pdf
  - iPM_FC_vs_absChange.pdf, iPM_FC_vs_absChange_rprotein.pdf
  - PTR_TPMfc_vs_iPMfc.pdf, PTR_by_category_boxplot.pdf
  - iPM_pool_composition_by_secondary.pdf : tall (x:2x) two full-pool stacked bars
        (syn1 w/ hatched deleted block | syn3A), blocks = secondary (Unclear & Cellular
        collapsed to Primary), Roman-indexed, colored by Primary.
  - iPM_tertiary_share_change_dumbbell.pdf : tall (x:2x) retained-pool deletion-corrected
        tertiary dumbbell (|delta|>0.5pp), broken x-axis; labels/dots colored like the bars.
  - protein_upgrade/, protein_downgrade/ (operon summaries + copied operon plots)
  - Compare_Ptn.txt (iPM FC / abs-change outliers + ribosomal-protein iPM table)
"""

import os
import shutil
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from adjustText import adjust_text
from scipy import stats

# ── Inputs / outputs ─────────────────────────────────────────────────────────
IN_CSV       = "Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv"   # written by 09
SYN3A_PTN_CSV = "../Syn1_Syn3A_Proteomics/syn3a_proteomics_summary_2026.csv"
OPERON_COV   = "../Syn1_Operon/gene_operon_coverage.tsv"
OPERON_PLOTS = "../Syn1_Operon/operon_plots"
OUTDIR       = "Compare_RNA_Protein"
OUT_REPORT   = f"{OUTDIR}/Compare_Ptn.txt"
os.makedirs(OUTDIR, exist_ok=True)

ABUNDANT_THRESH = 1.0

operon_cov = pd.read_csv(OPERON_COV, sep="\t")
operon_cov["locus_num"] = operon_cov["locus_tag"].str.extract(r"(\d+)$").astype(int)


def _extract_locus_num(s: pd.Series) -> pd.Series:
    return s.str.extract(r"(\d+)$").astype(int)


def _fc_plots(result, value_col, fc_col, log10fc_col, tag, ylabel_prefix,
              abundant_thresh=ABUNDANT_THRESH):
    """FC + log10FC histograms and the syn1-vs-syn3a scatter (mean-normalized)."""
    fc = result[fc_col]
    fc_stats = fc.describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])
    print(f"\n[{tag}] fold change statistics (n={int(fc_stats['count'])}):")
    print(f"  mean={fc_stats['mean']:.3f}  median={fc_stats['50%']:.3f}  std={fc_stats['std']:.3f}")

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.hist(fc, bins=50, color="steelblue", edgecolor="white", linewidth=0.4)
    ax.axvline(1,   color="black",  linestyle="-",  linewidth=1.0, label="FC = 1 (baseline)")
    ax.axvline(10,  color="red",    linestyle="--", linewidth=1.2, label="FC = 10")
    ax.axvline(0.1, color="orange", linestyle="--", linewidth=1.2, label="FC = 0.1")
    ax.set_xlabel(f"{tag} fold change  (syn3A / syn1)", fontsize=11)
    ax.set_ylabel(f"Number of {ylabel_prefix}", fontsize=11)
    ax.legend(fontsize=8)
    plt.tight_layout()
    fig.savefig(f"{OUTDIR}/{tag}_fold_change_distribution.pdf")
    plt.close(fig)

    log10fc = result[log10fc_col]
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.hist(log10fc, bins=50, color="steelblue", edgecolor="white", linewidth=0.4)
    ax.axvline(0,             color="black",  linestyle="-",  linewidth=1.0, label="log10FC = 0")
    ax.axvline(np.log10(10),  color="red",    linestyle="--", linewidth=1.2, label="log10(10) = 1")
    ax.axvline(np.log10(0.1), color="orange", linestyle="--", linewidth=1.2, label="log10(0.1) = -1")
    ax.set_xlabel(f"log10({tag} fold change)  (syn3A / syn1)", fontsize=12)
    ax.set_ylabel(f"Number of {ylabel_prefix}", fontsize=11)
    ax.legend(fontsize=8)
    plt.tight_layout()
    fig.savefig(f"{OUTDIR}/{tag}_log10FC_distribution.pdf")
    plt.close(fig)

    plot_df = result[
        result[f"{value_col}_syn1"].notna() & (result[f"{value_col}_syn1"] > 0) &
        result[f"{value_col}_syn3a"].notna() & (result[f"{value_col}_syn3a"] > 0)
    ]
    x = plot_df[f"{value_col}_syn1"]
    y = plot_df[f"{value_col}_syn3a"]
    pearson_r, pearson_p   = stats.pearsonr(np.log10(x), np.log10(y))
    spearman_r, spearman_p = stats.spearmanr(x, y)

    abundant  = plot_df[f"{value_col}_syn1"] > abundant_thresh
    mask_high = (plot_df[fc_col] > 10)  & abundant
    mask_low  = (plot_df[fc_col] < 0.1) & abundant
    mask_mid  = ~mask_high & ~mask_low

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(x[mask_mid],  y[mask_mid],  s=12, alpha=0.55, color="steelblue", edgecolors="none", label="other")
    ax.scatter(x[mask_high], y[mask_high], s=25, alpha=0.90, color="red",       edgecolors="none",
               label=f"FC > 10  & syn1 rel-{tag} > {abundant_thresh:g}")
    ax.scatter(x[mask_low],  y[mask_low],  s=25, alpha=0.90, color="orange",    edgecolors="none",
               label=f"FC < 0.1 & syn1 rel-{tag} > {abundant_thresh:g}")
    lims = [min(x.min(), y.min()) * 0.5, max(x.max(), y.max()) * 2]
    ax.plot(lims, lims, color="black", linewidth=0.8, linestyle="--")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel(f"syn1  rel-{tag} (log)", fontsize=12)
    ax.set_ylabel(f"syn3A rel-{tag} (log)", fontsize=12)
    ax.text(0.04, 0.96,
            f"Pearson r = {pearson_r:.3f}  (log10)\nSpearman r = {spearman_r:.3f}\nn = {len(plot_df)}",
            transform=ax.transAxes, fontsize=12, verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))
    ax.legend(fontsize=10, loc="lower right")
    plt.tight_layout()
    fig.savefig(f"{OUTDIR}/{tag}_correlation_syn1_vs_syn3a.pdf")
    plt.close(fig)
    print(f"[{tag}] Pearson r (log10) = {pearson_r:.3f} (p={pearson_p:.2e}); "
          f"Spearman r = {spearman_r:.3f} (p={spearman_p:.2e})")


def _fc_vs_abs_plot(df, fc_col, abs_col, base_col, val3_col, tag, outname,
                    highlight_col=None, highlight_label="highlighted"):
    """FC (x, log) vs absolute change (y, symlog); baseline by alpha; optional
    green highlight. Mirrors 09's _fc_vs_abs_plot."""
    sub = df[df[fc_col].notna() & (df[fc_col] > 0) &
             df[abs_col].notna() &
             df[base_col].notna() & (df[base_col] > 0)].copy()
    x, y = sub[fc_col], sub[abs_col]

    nz = y[y != 0].abs()
    linthresh = max(1e-3, float(np.nanmedian(nz))) if len(nz) else 1e-3

    DOT_SIZE = 28
    logb = np.log10(sub[base_col])
    lo, hi = float(logb.min()), float(logb.max())
    norm = (logb - lo) / (hi - lo) if hi > lo else pd.Series(0.5, index=logb.index)
    hue = (165 / 255, 15 / 255, 21 / 255)
    green = (16 / 255, 130 / 255, 60 / 255)
    rgba = np.zeros((len(sub), 4))
    rgba[:, 0], rgba[:, 1], rgba[:, 2] = hue
    rgba[:, 3] = 0.12 + 0.78 * norm.to_numpy()

    hlmask = None
    if highlight_col is not None and highlight_col in sub:
        hlmask = sub[highlight_col].fillna(False).to_numpy().astype(bool)
        rgba[hlmask, 0], rgba[hlmask, 1], rgba[hlmask, 2] = green

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.scatter(x, y, s=DOT_SIZE, facecolors=rgba, edgecolors="none")
    if hlmask is not None:
        ax.legend(handles=[Line2D([0], [0], marker="o", linestyle="", color=green,
                                  markersize=7, label=f"{highlight_label} (n={int(hlmask.sum())})")],
                  loc="lower left", fontsize=9)

    ax.axvline(1, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(0, color="black", linestyle="--", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=linthresh)
    ax.set_xlabel(rf"{tag} fold change  ($\mathrm{{rel}}_{{\mathrm{{syn3A}}}}/\mathrm{{rel}}_{{\mathrm{{syn1}}}}$)", fontsize=12)
    ax.set_ylabel(rf"{tag} absolute change  ($\mathrm{{rel}}_{{\mathrm{{syn3A}}}}-\mathrm{{rel}}_{{\mathrm{{syn1}}}}$)", fontsize=12)

    mlog = float(np.log10(x).abs().max()) * 1.05
    ax.set_xlim(10 ** (-mlog), 10 ** mlog)
    ymax = float(y.abs().max()) * 1.15
    ax.set_ylim(-ymax, ymax)

    both = sub[(sub[base_col] > 0) & (sub[val3_col] > 0)]
    pr, _ = stats.pearsonr(np.log10(both[base_col]), np.log10(both[val3_col]))
    ax.text(0.03, 0.97, f"n = {len(sub)}\nPearson r ({tag} syn1 vs syn3A) = {pr:.3f}",
            transform=ax.transAxes, va="top", ha="left", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85))

    cmap = LinearSegmentedColormap.from_list("hue_alpha", [(*hue, 0.12), (*hue, 0.90)])
    sm = plt.cm.ScalarMappable(norm=Normalize(vmin=lo, vmax=hi), cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label(f"log10 syn1 baseline (rel-{tag})", fontsize=10)

    top = pd.concat([sub[sub[abs_col] > 0].nlargest(5, abs_col),
                     sub[sub[abs_col] < 0].nsmallest(5, abs_col)])
    texts = [ax.text(r[fc_col], r[abs_col], f"{int(r['locus_num']):04d}", fontsize=7)
             for _, r in top.iterrows()]
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))

    plt.tight_layout()
    fig.savefig(f"{OUTDIR}/{outname}.pdf")
    plt.close(fig)
    print(f"[{outname}] plotted n={len(sub)}; Pearson(syn1 vs syn3A)={pr:.3f}")


def _operon_group_report(df_extreme, folder, label, tag, value_col, fc_col):
    """Per-operon summary txt + copied operon plots for an extreme-FC group."""
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
                r = result_lookup.get(tag, pd.DataFrame())
                r = r[r["locus_num"] == row["locus_num"]] if not r.empty else r
            if r.empty:
                v1 = v3 = fc_str = "N/A"
            else:
                rr = r.iloc[0]
                v1 = f"{rr[f'{value_col}_syn1']:.3f}" if pd.notna(rr[f'{value_col}_syn1']) else "N/A"
                v3 = f"{rr[f'{value_col}_syn3a']:.3f}" if pd.notna(rr[f'{value_col}_syn3a']) else "N/A"
                fc_str = f"{rr[fc_col]:.3f}" if pd.notna(rr[fc_col]) else "N/A"
            marker = " <<<" if row["locus_num"] in extreme_in_op["locus_num"].values else ""
            gname = str(row["gene_name"]) if pd.notna(row.get("gene_name")) else ""
            lines.append(f"{row['locus_tag']:<18} {gname:<14} {v1:>12} {v3:>12} {fc_str:>8}{marker}\n")

    with open(f"{folder}/operon_summary_{label}.txt", "w") as fh:
        fh.writelines(lines)
    copied = 0
    for opid in sorted(all_operons):
        for suffix in ("", "_wdepth"):
            src = os.path.join(OPERON_PLOTS, f"{opid}{suffix}.pdf")
            if os.path.exists(src):
                shutil.copy2(src, os.path.join(folder, f"{opid}{suffix}.pdf"))
                copied += 1
    print(f"[{tag} {label}] operon_summary written, {copied} plot(s) copied to {folder}/")


def _tukey_outliers(df, col, log):
    s = df[col].dropna()
    if log:
        s = s[s > 0]
        vals = np.log10(s)
    else:
        vals = s
    q1, q3 = vals.quantile(0.25), vals.quantile(0.75)
    iqr = q3 - q1
    lo, hi = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    if log:
        lo, hi = 10 ** lo, 10 ** hi
    up = df[df[col] > hi].sort_values(col, ascending=False)
    dn = df[df[col] < lo].sort_values(col, ascending=True)
    return up, dn, lo, hi


def _outlier_section(df, title, col, log, show_cols):
    up, dn, lo, hi = _tukey_outliers(df, col, log)
    out = [f"\n{'=' * 78}\n## {title}\n",
           f"   Tukey fences ({'log10 ' if log else ''}1.5*IQR): low < {lo:.3g}, high > {hi:.3g}\n",
           f"\n-- HIGH outliers (n={len(up)}) --\n",
           (up[show_cols].to_string(index=False) + "\n") if len(up) else "   (none)\n",
           f"\n-- LOW outliers (n={len(dn)}) --\n",
           (dn[show_cols].to_string(index=False) + "\n") if len(dn) else "   (none)\n"]
    return "".join(out)


# ─────────────────────────────────────────────────────────────────────────────
# Load 09's combined table and build the protein frame.
# ─────────────────────────────────────────────────────────────────────────────
coding = pd.read_csv(IN_CSV, sep="\t")
_ln = coding["locus_syn1"].fillna(coding["locus_syn3a"]).astype(str).str.extract(r"(\d+)$")[0]
coding["locus_num"] = pd.to_numeric(_ln, errors="coerce").astype("Int64")
_gn = coding["gene_name"].fillna("").astype(str)
_gp = coding["gene_product"].fillna("").astype(str).str.lower()
coding["is_rprotein"] = _gn.str.match(r"rp[slm]") | _gp.str.contains("ribosomal protein")

# Essentiality (syn3A design) + truncated product for the report tables.
_ess = pd.read_csv(SYN3A_PTN_CSV, usecols=["locus_tag", "essentiality"])
_ess["locus_num"] = _extract_locus_num(_ess["locus_tag"]).astype("Int64")
coding = coding.merge(_ess[["locus_num", "essentiality"]], on="locus_num", how="left")
coding["gene_product_disp"] = coding["gene_product"].fillna("").astype(str).str.slice(0, 45)

# Functional category (keyword on product). Priority: ribosomal > glycolysis >
# chaperone > protease > other (apply lowest priority first, override upward).
GLYCOLYSIS_KW = ["phosphofructokinase", "fructose-bisphosphate aldolase",
                 "triosephosphate isomerase", "glyceraldehyde-3-phosphate dehydrogenase",
                 "phosphoglycerate kinase", "phosphoglycerate mutase", "enolase",
                 "pyruvate kinase", "glucose-6-phosphate isomerase",
                 "phosphoglucose isomerase", "lactate dehydrogenase"]
CHAPERONE_KW = ["chaperone", "chaperonin", "trigger factor", "heat shock"]
PROTEASE_KW = ["protease", "peptidase"]
_PROD = coding["gene_product"].fillna("").astype(str).str.lower()


def _has(kws):
    m = pd.Series(False, index=coding.index)
    for k in kws:
        m |= _PROD.str.contains(k, regex=False)
    return m


_cat = pd.Series("other", index=coding.index)
_cat = _cat.mask(_has(PROTEASE_KW), "protease")
_cat = _cat.mask(_has(CHAPERONE_KW), "chaperone")
_cat = _cat.mask(_has(GLYCOLYSIS_KW), "glycolysis")
_cat = _cat.mask(coding["is_rprotein"].fillna(False), "ribosomal")
coding["func_category"] = _cat

ptn = coding.rename(columns={"relIPM_syn1": "iPM_mean_syn1",
                             "relIPM_syn3a": "iPM_mean_syn3a"}).copy()
ptn["iPM_log10FC"] = np.log10(ptn["iPM_fold_change"].replace(0, np.nan))
n_det = int((ptn["iPM_mean_syn1"].notna() & ptn["iPM_mean_syn3a"].notna()).sum())
print(f"[iPM] proteins with syn1&syn3A iPM: {n_det}; rows: {len(ptn)}")

# ── Distribution / correlation + FC-vs-abs figures ───────────────────────────
_fc_plots(ptn, value_col="iPM_mean", fc_col="iPM_fold_change",
          log10fc_col="iPM_log10FC", tag="iPM", ylabel_prefix="proteins")
_fc_vs_abs_plot(ptn, "iPM_fold_change", "iPM_abs_change", "iPM_mean_syn1",
                "iPM_mean_syn3a", "iPM", "iPM_FC_vs_absChange")
_fc_vs_abs_plot(ptn, "iPM_fold_change", "iPM_abs_change", "iPM_mean_syn1",
                "iPM_mean_syn3a", "iPM", "iPM_FC_vs_absChange_rprotein",
                highlight_col="is_rprotein", highlight_label="ribosomal protein")

# ── Operon summaries for extreme-FC protein groups ───────────────────────────
result_lookup = {"iPM": ptn}
_operon_group_report(ptn[ptn["iPM_fold_change"] > 10].sort_values("iPM_fold_change", ascending=False),
                     f"{OUTDIR}/protein_upgrade",   "FC_gt10",  "iPM", "iPM_mean", "iPM_fold_change")
_operon_group_report(ptn[ptn["iPM_fold_change"] < 0.1].sort_values("iPM_fold_change"),
                     f"{OUTDIR}/protein_downgrade", "FC_lt0.1", "iPM", "iPM_mean", "iPM_fold_change")

# ── PTR (protein-to-transcript ratio) = iPM_FC / TPM_FC ──────────────────────
# PTR > baseline => protein rose more than its mRNA (post-transcriptional gain).
CAT_COLORS = {"glycolysis": "#1f77b4", "chaperone": "#d62728", "protease": "#ff7f0e",
              "ribosomal": "#2ca02c", "other": "lightgray"}
CAT_ORDER = ["glycolysis", "other", "chaperone", "protease", "ribosomal"]


def _ptr_scatter():
    d = coding[coding["TPM_fold_change"].notna() & (coding["TPM_fold_change"] > 0) &
               coding["iPM_fold_change"].notna() & (coding["iPM_fold_change"] > 0)].copy()
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    texts = []
    for c in ["other", "ribosomal", "glycolysis", "chaperone", "protease"]:
        sub = d[d["func_category"] == c]
        ax.scatter(sub["TPM_fold_change"], sub["iPM_fold_change"],
                   s=(14 if c == "other" else 34), alpha=(0.4 if c in ("other", "ribosomal") else 0.9),
                   color=CAT_COLORS[c], edgecolors="none", label=f"{c} (n={len(sub)})")
        if c in ("glycolysis", "chaperone", "protease"):
            for _, r in sub.iterrows():
                texts.append(ax.text(r["TPM_fold_change"], r["iPM_fold_change"],
                                     f"{int(r['locus_num']):04d}", fontsize=6))
    lims = [min(d["TPM_fold_change"].min(), d["iPM_fold_change"].min()) * 0.5,
            max(d["TPM_fold_change"].max(), d["iPM_fold_change"].max()) * 2]
    ax.plot(lims, lims, "k--", lw=0.9, label="y = x (PTR unchanged)")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel(r"TPM fold change  ($\mathrm{rel}_{\mathrm{syn3A}}/\mathrm{rel}_{\mathrm{syn1}}$)", fontsize=12)
    ax.set_ylabel(r"iPM fold change  ($\mathrm{rel}_{\mathrm{syn3A}}/\mathrm{rel}_{\mathrm{syn1}}$)", fontsize=12)
    ax.set_title("Transcript vs protein fold change (PTR; above y=x = protein-favored)", fontsize=11)
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.4))
    ax.legend(fontsize=8, loc="lower right")
    plt.tight_layout()
    fig.savefig(f"{OUTDIR}/PTR_TPMfc_vs_iPMfc.pdf")
    plt.close(fig)
    print(f"[PTR scatter] n={len(d)}")


def _ptr_boxplot():
    d = coding[coding["PTR_fold_change"].notna() & (coding["PTR_fold_change"] > 0)].copy()
    d["log10PTR"] = np.log10(d["PTR_fold_change"])
    cats = [c for c in CAT_ORDER if (d["func_category"] == c).sum() >= 1]
    data = [d.loc[d["func_category"] == c, "log10PTR"].values for c in cats]
    base = d.loc[d["func_category"] == "other", "log10PTR"].values
    other_med = float(np.median(base)) if len(base) else 0.0

    fig, ax = plt.subplots(figsize=(8, 5.5))
    bp = ax.boxplot(data, positions=range(len(cats)), widths=0.6, showfliers=False, patch_artist=True,
                    medianprops=dict(color="black"))
    for patch, c in zip(bp["boxes"], cats):
        patch.set_facecolor(CAT_COLORS[c]); patch.set_alpha(0.5)
    rng = np.random.default_rng(0)
    for i, dd in enumerate(data):
        ax.scatter(i + rng.uniform(-0.18, 0.18, size=len(dd)), dd, s=10, color="black", alpha=0.4, edgecolors="none")
    ax.axhline(other_med, color="red", linestyle="--", linewidth=1, label="'other' median (baseline)")

    labels = []
    print("\n[PTR by functional category]")
    for c, dd in zip(cats, data):
        if c == "other" or len(dd) < 2 or len(base) < 2:
            p, sig = np.nan, "baseline" if c == "other" else "n/a"
        else:
            p = stats.mannwhitneyu(dd, base, alternative="two-sided").pvalue
            sig = f"p={p:.1e} ({'sig' if p < 0.05 else 'n.s.'})"
        labels.append(f"{c}\nn={len(dd)}\n{sig}")
        print(f"  {c:<12} n={len(dd):3d}  median log10PTR={np.median(dd):+.3f}  p(vs other)={p:.2e}")
    ax.set_xticks(range(len(cats))); ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("log10 PTR fold change  (iPM_FC / TPM_FC)", fontsize=12)
    ax.set_title("Protein-to-transcript ratio change by functional category", fontsize=11)
    ax.legend(fontsize=8, loc="upper right")
    plt.tight_layout()
    fig.savefig(f"{OUTDIR}/PTR_by_category_boxplot.pdf")
    plt.close(fig)


_ptr_scatter()
_ptr_boxplot()


# ─────────────────────────────────────────────────────────────────────────────
# Proteome reallocation story plots (iPM) — the protein counterpart of 09's two
# mRNA plots. Same design: deletion-corrected (retained-pool) normalization, the
# curated Secondary/Tertiary function annotation, Primary-family colors.
#   1) iPM_pool_composition_by_secondary.pdf — two FULL-pool stacked bars (syn1
#      with hatched 'deleted' block | syn3A); blocks = secondary (Unclear &
#      Cellular Processes collapsed to Primary; retained-unannotated dropped),
#      Roman-indexed, colored by Primary; no spines/title/y-axis; legend below.
#   2) iPM_tertiary_share_change_dumbbell.pdf — retained-pool tertiary dumbbell
#      (|Δ|>0.5pp), broken x-axis, labels/dots colored like (1), share at dots.
# NOTE: ribosomal-protein iPM is digestion-biased; read Ribosome/Translation at
# the protein level with caution.
# ─────────────────────────────────────────────────────────────────────────────
import re as _re
SYN3A_GFF = "../Genomes_Input/syn3a_genome.gff3"
SYN3A_FUNC_XLSX = "../Syn1_Syn3A_Proteomics/syn3A_proteome_annotated.xlsx"
PRIM_COLORS = {
    "Genetic Information Processing":       "#3b6db3",
    "Metabolism":                          "#3f9e5a",
    "Unclear":                             "#9aa0a6",
    "Cellular Processes":                  "#8e6bb1",
    "Environmental Information Processing": "#2aa6a0",
    "Exogenous":                           "#c0654e",
}


def _shade(base, frac):
    return tuple(c + (1 - c) * frac for c in mcolors.to_rgb(base))


def _text_color(rgb):
    lum = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
    return "white" if lum < 0.55 else "black"


def _roman(num):
    out = ""
    for v, s in [(10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I")]:
        while num >= v:
            out += s; num -= v
    return out


def _load_gff_locus_nums(path):
    nums = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
                continue
            m = _re.search(r"locus_tag=([^;]+)", f[8])
            mm = _re.search(r"(\d+)$", m.group(1)) if m else None
            if mm:
                nums.add(int(mm.group(1)))
    return nums


def _load_function_for_plots():
    fa = pd.read_excel(SYN3A_FUNC_XLSX, sheet_name=0)
    fa["locus_num"] = _extract_locus_num(fa["Locus Tag"]).astype("Int64")
    cm = coding.merge(fa[["locus_num", "Primary Function", "Secondary Function",
                          "Tertiary Function"]], on="locus_num", how="left")
    cm["is_deleted"] = ~cm["locus_num"].isin(_load_gff_locus_nums(SYN3A_GFF))
    return cm


def _block_palette(cm):
    """Add cm['block'] (secondary; Unclear/Cellular -> Primary) and return the
    Primary-family-shaded color per block. Ordering/shading use the syn1 mRNA
    (TPM) share so the Roman index + colors match 09's composition bar exactly
    (left-bar consistency across 09/10); bar heights here still use iPM."""
    from collections import defaultdict
    cm["block"] = cm["Secondary Function"]
    coll = cm["Primary Function"].isin(["Unclear", "Cellular Processes"])
    cm.loc[coll, "block"] = cm.loc[coll, "Primary Function"]
    annm = cm["block"].notna()
    t1, t3 = cm["relTPM_syn1"] > 0, cm["relTPM_syn3a"] > 0
    g1b = cm.loc[t1 & annm].groupby("block")["relTPM_syn1"].sum()
    g3k = set(cm.loc[t3 & annm, "block"])
    block_prim = cm.dropna(subset=["block"]).groupby("block")["Primary Function"].first().to_dict()
    prim_order = list(PRIM_COLORS)
    blocks = sorted(set(g1b.index) | g3k,
                    key=lambda b: (prim_order.index(block_prim.get(b, "Unclear"))
                                   if block_prim.get(b) in prim_order else 99, -float(g1b.get(b, 0))))
    prim_blocks = defaultdict(list)
    for b in blocks:
        prim_blocks[block_prim.get(b, "Unclear")].append(b)
    block_color = {}
    for p, bl in prim_blocks.items():
        base, nb = PRIM_COLORS.get(p, "#9aa0a6"), len(bl)
        for i, b in enumerate(bl):
            block_color[b] = _shade(base, 0 if nb == 1 else (i / (nb - 1)) * 0.5)
    return blocks, block_color


def _ipm_composition_secondary_bar():
    """Plot 1: two full-pool proteome (iPM) stacked bars by secondary function."""
    cm = _load_function_for_plots()
    blocks, block_color = _block_palette(cm)
    det1, det3 = cm["relIPM_syn1"] > 0, cm["relIPM_syn3a"] > 0
    annm = cm["block"].notna()
    s_sec1 = cm.loc[det1 & annm].groupby("block")["relIPM_syn1"].sum()
    s_del1 = cm.loc[det1 & cm["is_deleted"], "relIPM_syn1"].sum()
    tot1 = s_sec1.sum() + s_del1
    g1, deleted1 = s_sec1 / tot1 * 100, s_del1 / tot1 * 100
    s_sec3 = cm.loc[det3 & annm].groupby("block")["relIPM_syn3a"].sum()
    g3 = s_sec3 / s_sec3.sum() * 100
    roman = {b: _roman(i + 1) for i, b in enumerate(blocks)}

    W, MIN_INSIDE = 0.84, 2.0
    fig, ax = plt.subplots(figsize=(4.5, 9))
    side = {0: [], 1: []}
    for xi, (g, deleted) in enumerate([(g1, deleted1), (g3, 0.0)]):
        bottom = 0.0
        for b in blocks:
            h = float(g.get(b, 0))
            if h <= 0:
                continue
            ax.bar(xi, h, W, bottom=bottom, color=block_color[b], edgecolor="white", linewidth=0.5)
            if h >= MIN_INSIDE:
                ax.text(xi, bottom + h / 2, f"{roman[b]} {h:.0f}%", ha="center", va="center",
                        fontsize=8, color=_text_color(block_color[b]))
            else:
                side[xi].append((bottom + h / 2, f"{roman[b]} {h:.1f}%"))
            bottom += h
        if deleted > 0:
            ax.bar(xi, deleted, W, bottom=bottom, color=(0.84, 0.19, 0.15),
                   hatch="///", edgecolor="white", linewidth=0)
            ax.text(xi, bottom + deleted / 2, f"deleted {deleted:.0f}%", ha="center",
                    va="center", color="white", fontweight="bold", fontsize=8)

    def _place(items, x_text, ha, x_edge):
        items = sorted(items)
        ys = [y for y, _ in items]
        for i in range(1, len(ys)):
            ys[i] = max(ys[i], ys[i - 1] + 2.8)
        for (y0, txt), y in zip(items, ys):
            ax.annotate(txt, xy=(x_edge, y0), xytext=(x_text, y), ha=ha, va="center",
                        fontsize=7, color="#333", arrowprops=dict(arrowstyle="-", color="gray", lw=0.4))
    _place(side[0], -0.72, "right", -W / 2)
    _place(side[1], 1.72, "left", 1 + W / 2)

    ax.set_xticks([0, 1]); ax.set_xticklabels(["syn1", "syn3A"], fontsize=13)
    ax.set_ylim(0, 100); ax.set_xlim(-1.35, 2.1)
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    handles = [Patch(facecolor=block_color[b], edgecolor="white", label=f"{roman[b]} — {b}")
               for b in blocks]
    handles.append(Patch(facecolor=(0.84, 0.19, 0.15), hatch="///", edgecolor="white",
                         label="deleted (syn1 only)"))
    ax.legend(handles=handles, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.02),
              ncol=2, frameon=False)
    fig.savefig(f"{OUTDIR}/iPM_pool_composition_by_secondary.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUTDIR}/iPM_pool_composition_by_secondary.pdf "
          f"(deleted={deleted1:.1f}%, Translation {g1.get('Translation',0):.0f}->{g3.get('Translation',0):.0f}%)")


def _ipm_tertiary_share_dumbbell(min_change=0.5):
    """Plot 2: retained-pool (deletion-corrected) tertiary iPM dumbbell."""
    import math
    cm = _load_function_for_plots()
    _, block_color = _block_palette(cm)
    det1, det3 = cm["relIPM_syn1"] > 0, cm["relIPM_syn3a"] > 0
    ret = ~cm["is_deleted"]
    tot1 = cm.loc[ret & det1, "relIPM_syn1"].sum()
    tot3 = cm.loc[det3, "relIPM_syn3a"].sum()
    ann = cm.dropna(subset=["Tertiary Function"])
    s1 = ann.loc[det1.loc[ann.index]].groupby("Tertiary Function")["relIPM_syn1"].sum() / tot1 * 100
    s3 = ann.loc[det3.loc[ann.index]].groupby("Tertiary Function")["relIPM_syn3a"].sum() / tot3 * 100
    ter_prim = ann.groupby("Tertiary Function")["Primary Function"].first().to_dict()
    ter_sec = ann.groupby("Tertiary Function")["Secondary Function"].first().to_dict()
    d = pd.DataFrame({"syn1": s1, "syn3a": s3}).fillna(0.0)
    d["change"] = d["syn3a"] - d["syn1"]
    d = d[d["change"].abs() > min_change].sort_values("change")
    n = len(d)

    def _color(name):
        p, s = ter_prim.get(name), ter_sec.get(name)
        return block_color.get(p if p in ("Unclear", "Cellular Processes") else s, "#555555")

    hi_row = d[["syn1", "syn3a"]].min(axis=1) >= 20
    assert hi_row.any() and (~hi_row).any(), "no clear high/low cluster for the broken axis"
    lo_max = float(d.loc[~hi_row, ["syn1", "syn3a"]].values.max())
    hi_min = float(d.loc[hi_row, ["syn1", "syn3a"]].values.min())
    hi_max = float(d.loc[hi_row, ["syn1", "syn3a"]].values.max())
    lo_xlim = (0, math.ceil(lo_max) + 1)
    hi_xlim = (math.floor(hi_min) - 2, math.ceil(hi_max) + 2)
    dlo_w, dhi_w = lo_xlim[1] - lo_xlim[0], hi_xlim[1] - hi_xlim[0]

    fig = plt.figure(figsize=(4.5, 9))
    gs = fig.add_gridspec(1, 2, width_ratios=[dlo_w, dhi_w], wspace=0.06)
    axbg = fig.add_subplot(gs[0, :])
    axbg.set_xlim(0, 1); axbg.set_ylim(-0.6, n - 0.4)
    for i in range(n):
        axbg.axhline(i, color="0.9", lw=0.8, zorder=0)
    axbg.set_xticks([]); axbg.set_yticks([])
    for sp in axbg.spines.values():
        sp.set_visible(False)
    axbg.patch.set_visible(False)
    axl = fig.add_subplot(gs[0]); axr = fig.add_subplot(gs[1], sharey=axl)
    axl.patch.set_visible(False); axr.patch.set_visible(False)
    axl.set_yticks(range(n))
    for i, (name, row) in enumerate(d.iterrows()):
        ax = axr if max(row["syn1"], row["syn3a"]) > lo_xlim[1] else axl
        col = _color(name)
        ax.plot([row["syn1"], row["syn3a"]], [i, i], color="lightgray", lw=1.8, zorder=2)
        ax.scatter(row["syn1"], i, s=52, facecolor="white", edgecolor="gray", linewidth=1.2, zorder=3)
        ax.scatter(row["syn3a"], i, s=66, color=col, edgecolor="white", linewidth=0.6, zorder=4)
        ax.text(row["syn3a"], i + 0.28, f"{row['syn3a']:.1f}", ha="center", va="center", fontsize=8, color="#333")
        ax.text(row["syn1"], i - 0.28, f"{row['syn1']:.1f}", ha="center", va="center", fontsize=8, color="#333")

    axl.set_xlim(*lo_xlim); axr.set_xlim(*hi_xlim)
    axl.set_ylim(-0.6, n - 0.4)
    axl.set_yticklabels([])
    tr = axl.get_yaxis_transform()
    for i, name in enumerate(d.index):
        c = _color(name)
        axl.text(-0.04, i + 0.22, ter_sec.get(name, ""), transform=tr, ha="right", va="center",
                 fontsize=8, color=c, style="italic")
        axl.text(-0.04, i - 0.18, name, transform=tr, ha="right", va="center",
                 fontsize=10, color=c)
    for sp in ("top", "left", "right"):
        axl.spines[sp].set_visible(False); axr.spines[sp].set_visible(False)
    axl.tick_params(left=False, labelsize=10)
    axr.tick_params(left=False, labelsize=10)
    plt.setp(axr.get_yticklabels(), visible=False)

    dd = 0.012
    kw = dict(transform=axl.transAxes, color="k", clip_on=False, lw=0.9)
    axl.plot((1 - dd, 1 + dd), (-dd, dd), **kw)
    kw.update(transform=axr.transAxes)
    axr.plot((-dd, dd), (-dd, dd), **kw)

    dumb_mid = (axl.get_position().x0 + axr.get_position().x1) / 2
    fig.text(dumb_mid, 0.05, "Proteome Pool Share Change (%)", ha="center", fontsize=12)
    mk = [Line2D([0], [0], marker="o", linestyle="", markersize=9, markerfacecolor="white",
                 markeredgecolor="gray", label="syn1"),
          Line2D([0], [0], marker="o", linestyle="", markersize=9, markerfacecolor="#555",
                 markeredgecolor="white", label="syn3A")]
    axr.legend(handles=mk, fontsize=10, loc="lower right", frameon=False)
    fig.savefig(f"{OUTDIR}/iPM_tertiary_share_change_dumbbell.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUTDIR}/iPM_tertiary_share_change_dumbbell.pdf "
          f"({n} tertiary, |Δ|>{min_change}pp; x-break {lo_xlim[1]:.0f}-{hi_xlim[0]:.0f}% hidden)")


_ipm_composition_secondary_bar()
_ipm_tertiary_share_dumbbell()


# ── Outlier report (Compare_Ptn.txt) ─────────────────────────────────────────
IPM_SHOW = ["locus_syn1", "gene_name", "essentiality", "gene_product_disp",
            "relIPM_syn1", "relIPM_syn3a", "iPM_fold_change", "iPM_abs_change"]
report = ["#" * 78, "\n# PROTEIN (iPM) FOLD-CHANGE AND ABSOLUTE-CHANGE OUTLIERS (coding genes)\n",
          "#" * 78, "\n",
          "# iPM = proteomics; rel* = mean-normalized. Source: " + IN_CSV + "\n"]
report.append(_outlier_section(coding, "iPM fold change", "iPM_fold_change", True,  IPM_SHOW))
report.append(_outlier_section(coding, "iPM absolute change (relIPM delta)", "iPM_abs_change", False, IPM_SHOW))

# PTR (protein-to-transcript ratio) outliers — r-proteins excluded (unreliable iPM).
PTR_SHOW = ["locus_syn1", "gene_name", "func_category", "essentiality", "gene_product_disp",
            "TPM_fold_change", "iPM_fold_change", "PTR_fold_change"]
report.append(_outlier_section(coding[~coding["is_rprotein"].fillna(False)],
                               "PTR fold change = iPM_FC / TPM_FC  (r-proteins excluded)",
                               "PTR_fold_change", True, PTR_SHOW))

RPROT_COLS = ["locus_syn1", "gene_name", "essentiality", "gene_product_disp",
              "relTPM_syn1", "relTPM_syn3a", "TPM_fold_change", "TPM_abs_change",
              "relIPM_syn1", "relIPM_syn3a", "iPM_fold_change", "iPM_abs_change"]
rprot_tbl = coding[coding["is_rprotein"]].sort_values("locus_num")
report += ["\n\n", "#" * 78,
           f"\n# RIBOSOMAL PROTEIN GENES (n={len(rprot_tbl)}) -- TPM + iPM, ascending locus\n",
           "#" * 78, "\n",
           rprot_tbl[RPROT_COLS].to_string(index=False) + "\n"]

with open(OUT_REPORT, "w") as fh:
    fh.write("".join(report))
print(f"\nSaved: {OUT_REPORT}")
for title, col, log in [("iPM fold change", "iPM_fold_change", True),
                        ("iPM abs change",  "iPM_abs_change",  False)]:
    up, dn, lo, hi = _tukey_outliers(coding, col, log)
    print(f"  [{title}] high={len(up)}, low={len(dn)}  (fences {lo:.3g} .. {hi:.3g})")
