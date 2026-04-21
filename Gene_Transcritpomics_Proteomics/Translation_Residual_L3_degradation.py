'''
Residual analysis of Transcriptomics and Proteomics for Syn1 - Level 3: protein degradation

Map Mpn protein half-lives (Burgos et al. 2020) onto Syn1 via the Mpn<->Syn1
reciprocal-best-hit homology table.
'''

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from Bio import SeqIO

MOTHER_DIR = ".."
RBH_TSV = MOTHER_DIR + "/Genomes_Input/homology_syn1_mpn/mpn_syn1_rbh_homology_table.tsv"
SYN1_MPN_TSV = MOTHER_DIR + "/Genomes_Input/homology_syn1_mpn/syn1_mpn.tsv"
MPN_GB = MOTHER_DIR + "/Genomes_Input/Mpn.gb"

HALFLIFE_XLSX = Path("Burgos_2020_SI_Ptn_Halflives.xlsx")
OMICS_CSV = Path("./syn1_genes_transcriptomics_proteomics.csv")

OUT_DIR = Path("./residual_analysis/")
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUT_DIR / "syn1_ptn_degradation_from_mpn.csv"

BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]

Lon_Syn1_Num = 128 # This study
FtsH_Syn1_Num = 292 # This study
Syn1_Volume = 0.033 # fL, 0.2 um radius

Lon_Mpn_Num = 122 # Quantification of mRNA and protein and integration with protein turnover in a bacterium, 2011
FtsH_Mpn_Num = 689 # Quantification of mRNA and protein and integration with protein turnover in a bacterium, 2011
Mpn_Volume = 0.05 # fL, Quantification of mRNA and protein and integration with protein turnover in a bacterium, 2011

Cyto_Halflive_factor = (Lon_Mpn_Num / Mpn_Volume) / (Lon_Syn1_Num / Syn1_Volume) # halflive in Syn1 = Catyo_halflive in Mpn * Cyto_Halflive_factor
Mem_Halflive_factor = (FtsH_Mpn_Num / Mpn_Volume) / (FtsH_Syn1_Num / Syn1_Volume) # Assume same surface area to volume ratio between Mpn and Syn1

def build_old_to_new_locus_map(gb_path: Path) -> dict[str, str]:
    """Map Burgos-style old locus tags (e.g. MPN001) to RefSeq locus tags (e.g. MPN_RS00005)."""
    mapping: dict[str, str] = {}
    for rec in SeqIO.parse(str(gb_path), "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            new_tag = feat.qualifiers.get("locus_tag", [None])[0]
            if not new_tag:
                continue
            for old in feat.qualifiers.get("old_locus_tag", []):
                mapping[old] = new_tag
    return mapping


def main() -> None:
    halflife = pd.read_excel(HALFLIFE_XLSX, sheet_name="Sheet1")
    halflife = halflife.rename(columns={
        "protein half-life_h": "halflife_h",
        "protein half-life_clipped": "halflife_clipped",
    })
    halflife = halflife[["gene", "halflife_h", "halflife_clipped"]].dropna(subset=["halflife_h"])

    old2new = build_old_to_new_locus_map(MPN_GB)
    halflife["mpn_locus_tag"] = halflife["gene"].map(old2new)

    n_burgos = len(halflife)
    n_mapped = halflife["mpn_locus_tag"].notna().sum()
    print(f"Burgos entries with half-life: {n_burgos}")
    print(f"Mapped to RefSeq locus_tag:    {n_mapped}")

    rbh = pd.read_csv(RBH_TSV, sep="\t")
    rbh_merged = rbh.merge(
        halflife.dropna(subset=["mpn_locus_tag"]),
        on="mpn_locus_tag",
        how="inner",
    )
    n_syn1_rbh_total = rbh["syn1_locus_tag"].nunique()
    n_syn1_rbh_hl = rbh_merged["syn1_locus_tag"].nunique()
    print(f"Syn1 proteins with a Mpn RBH ortholog:               {n_syn1_rbh_total}")
    print(f"Syn1 proteins with Mpn RBH-derived half-life mapped: {n_syn1_rbh_hl}")

    # Keep only RBH pairs for high-confidence half-life transfer.
    out = rbh_merged.rename(columns={"gene": "mpn_old_locus_tag"}).copy()
    keep = [
        "syn1_locus_tag", "syn1_gene", "mpn_locus_tag", "mpn_old_locus_tag", "mpn_gene",
        "pident", "evalue", "bitscore", "halflife_h", "halflife_clipped",
    ]
    out = out[keep].drop_duplicates("syn1_locus_tag").sort_values("bitscore", ascending=False)

    # Attach Syn1 localization from the omics table.
    omics_loc = pd.read_csv(OMICS_CSV)[["locus_tag", "ptn_localization"]]
    out = out.merge(omics_loc, left_on="syn1_locus_tag", right_on="locus_tag",
                    how="left").drop(columns=["locus_tag"])

    loc_counts = out["ptn_localization"].value_counts(dropna=False)
    print("\nLocalization of the 245 RBH-mapped Syn1 proteins:")
    print(loc_counts.to_string())
    # loc_counts.rename_axis("ptn_localization").to_csv(
    #     OUT_DIR / "syn1_halflife_localization_stats.csv", header=["n"],
    # )

    # Correct Mpn half-life -> Syn1 half-life using protease-abundance ratios.
    # Cyto (Lon) for cytoplasmic; Mem (FtsH) for membrane and lipoprotein.
    def _factor(loc):
        if loc == "cytoplasmic":
            return Cyto_Halflive_factor
        if loc in ("membrane", "lipoprotein"):
            return Mem_Halflive_factor
        return np.nan

    out["halflife_factor"] = out["ptn_localization"].map(_factor)
    out["halflife_h_syn1"] = out["halflife_h"] * out["halflife_factor"]

    print(f"\nCyto_Halflive_factor (Lon):  {Cyto_Halflive_factor:.3f}")
    print(f"Mem_Halflive_factor  (FtsH): {Mem_Halflive_factor:.3f}")

    def _stats(s, tag):
        s = s.dropna()
        print(f"  [{tag}] n={len(s)}  median={s.median():.2f}  mean={s.mean():.2f}  "
              f"min={s.min():.2f}  max={s.max():.2f}")

    print("\nHalf-life stats (hours):")
    _stats(out["halflife_h"], "Mpn, all 245")
    _stats(out.loc[out["ptn_localization"] == "cytoplasmic", "halflife_h"], "Mpn cytoplasmic")
    _stats(out.loc[out["ptn_localization"].isin(["membrane", "lipoprotein"]),
                   "halflife_h"], "Mpn membrane+lipo")
    _stats(out["halflife_h_syn1"], "Syn1-corrected (all w/ factor)")
    _stats(out.loc[out["ptn_localization"] == "cytoplasmic", "halflife_h_syn1"],
           "Syn1-corrected cytoplasmic")
    _stats(out.loc[out["ptn_localization"].isin(["membrane", "lipoprotein"]),
                   "halflife_h_syn1"], "Syn1-corrected membrane+lipo")

    out.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV} ({len(out)} rows)")

    # ---- Half-life distribution over the RBH-mapped Syn1 proteins ----
    hl = out.dropna(subset=["halflife_h_syn1"]).drop_duplicates("syn1_locus_tag").copy()
    min_row = hl.loc[hl["halflife_h_syn1"].idxmin()]
    # min_gene = min_row["syn1_locus_tag"]
    min_val = min_row["halflife_h_syn1"]

    fig, ax = plt.subplots(figsize=(3, 3))
    ax.hist(hl["halflife_h_syn1"], bins=40, color="#4C8BB5", edgecolor="white")
    ax.axvline(hl["halflife_h_syn1"].median(), color="k", ls="--", lw=1,
               label=f"median = {hl['halflife_h_syn1'].median():.1f} h")
    ax.axvline(min_val, color="#c0392b", ls=":", lw=1.2,
               label=f"shortest: {min_val:.2f} h")
    ax.set_xlabel(r" $t_{1/2}$ from Mpn" + f" \n corrected by Lon/FtsH Conc.", fontsize=10)
    ax.set_ylabel("Syn1 proteins (count)", fontsize=10)
    ax.set_title(f"Syn1 protein half-life (h)", fontsize=10)
    ax.text(0.95, 0.6, f"n = {len(hl)}", fontsize=10,
            ha="right", va="center", transform=ax.transAxes)
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "syn1_halflife_from_mpn_dist.pdf")
    plt.close(fig)

    # ---- Residual vs half-life ----
    omics = pd.read_csv(OMICS_CSV)
    omics = omics[["locus_tag", "avg_sense_TPM", "iPM_mean"]].copy()
    omics = omics[(omics["avg_sense_TPM"] > 0) & (omics["iPM_mean"] > 0)].dropna()
    omics["log10_TPM"] = np.log10(omics["avg_sense_TPM"])
    omics["log10_iPM"] = np.log10(omics["iPM_mean"])
    slope, intercept = np.polyfit(omics["log10_TPM"], omics["log10_iPM"], 1)
    omics["proteome_residual"] = omics["log10_iPM"] - (slope * omics["log10_TPM"] + intercept)
    r_base, p_base = pearsonr(omics["log10_TPM"], omics["log10_iPM"])
    print(f"\nBase correlation log10(TPM) vs log10(iPM): Pearson r = {r_base:.3f} (p={p_base:.2e}), "
          f"n={len(omics)}")
    print(f"Fit: log10(iPM) = {slope:.3f} * log10(TPM) + {intercept:.3f}")

    merged = hl.merge(
        omics[["locus_tag", "proteome_residual", "log10_TPM", "log10_iPM"]],
        left_on="syn1_locus_tag", right_on="locus_tag", how="inner",
    )
    merged = merged[merged["halflife_h_syn1"] > 0].copy()
    merged["log10_halflife"] = np.log10(merged["halflife_h_syn1"])

    r_p, p_p = pearsonr(merged["log10_halflife"], merged["proteome_residual"])
    r_s, p_s = spearmanr(merged["halflife_h_syn1"], merged["proteome_residual"])
    print(f"Residual vs log10(half-life): Pearson r = {r_p:.3f} (p={p_p:.2e}), "
          f"Spearman rho = {r_s:.3f} (p={p_s:.2e}), n={len(merged)}")

    fig, ax = plt.subplots(figsize=(3, 3))
    x = merged["log10_halflife"].values
    y = merged["proteome_residual"].values
    ax.scatter(x, y, s=12, alpha=0.55, color="#4C8BB5", edgecolors="none")
    m, b = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, m * xs + b, color="#c0392b", linewidth=1.2)
    ax.set_xlabel("log₁₀(Syn1-corrected t½, h)", fontsize=10)
    ax.set_ylabel("log₁₀(iPM) − fit(log₁₀TPM)", fontsize=10)
    ax.set_title(f"r = {r_p:.3f}  p = {p_p:.1e}  n = {len(x)}",
                 fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "syn1_residual_vs_halflife.pdf")
    plt.close(fig)
    print(f"Wrote {OUT_DIR/'syn1_halflife_from_mpn_dist.pdf'} and "
          f"{OUT_DIR/'syn1_residual_vs_halflife.pdf'}")

    # ---- Multiple regression: log10(iPM) ~ log10(TPM) + log10(t1/2) ----
    import statsmodels.api as sm

    mreg = merged.copy()
    unclipped = mreg[mreg["halflife_clipped"] != True].copy()

    def fit_mlr(df, label):
        X = sm.add_constant(df[["log10_TPM", "log10_halflife"]].values)
        y = df["log10_iPM"].values
        full = sm.OLS(y, X).fit()
        # Reduced model (TPM only) for partial-r / ΔR²
        Xr = sm.add_constant(df[["log10_TPM"]].values)
        red = sm.OLS(y, Xr).fit()
        # Partial correlation of log10(t1/2) with log10(iPM) given log10(TPM)
        rx = sm.OLS(df["log10_halflife"].values, Xr).fit().resid
        ry = red.resid
        pr, pp = pearsonr(rx, ry)
        b1, b2 = full.params[1], full.params[2]
        se1, se2 = full.bse[1], full.bse[2]
        p1, p2 = full.pvalues[1], full.pvalues[2]
        print(f"\n[{label}]  n={len(df)}")
        print(f"  log10(iPM) = {full.params[0]:.3f} + {b1:.3f}*log10(TPM) + {b2:.3f}*log10(t1/2)")
        print(f"    β(TPM)   = {b1:.3f} ± {se1:.3f}   p={p1:.2e}")
        print(f"    β(t1/2)  = {b2:.3f} ± {se2:.3f}   p={p2:.2e}   (theory: ≈ 1)")
        print(f"    R²(full) = {full.rsquared:.3f}   R²(TPM only) = {red.rsquared:.3f}   "
              f"ΔR² = {full.rsquared - red.rsquared:.3f}")
        print(f"    partial-r(log10 t1/2 | log10 TPM) = {pr:.3f}  (p={pp:.2e})")
        return full, red, pr

    fit_mlr(mreg,      "ALL half-lives (incl. clipped)")
    full, red, pr = fit_mlr(unclipped, "UNCLIPPED half-lives only")

    # Partial-residual plot: residualize both axes against log10(TPM), then scatter
    Xr = sm.add_constant(unclipped[["log10_TPM"]].values)
    rx = sm.OLS(unclipped["log10_halflife"].values, Xr).fit().resid
    ry = sm.OLS(unclipped["log10_iPM"].values, Xr).fit().resid

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    ax.scatter(rx, ry, s=14, alpha=0.6, color="#4C72B0", edgecolor="none")
    m, b = np.polyfit(rx, ry, 1)
    xs = np.linspace(rx.min(), rx.max(), 100)
    ax.plot(xs, m * xs + b, color="crimson", lw=1.5,
            label=f"slope = {m:.3f}")
    ax.axhline(0, color="gray", lw=0.7, ls=":")
    ax.axvline(0, color="gray", lw=0.7, ls=":")
    ax.set_xlabel("log₁₀(t½) | log₁₀(TPM)   (residual)")
    ax.set_ylabel("log₁₀(iPM) | log₁₀(TPM)   (residual)")
    ax.set_title(f"Partial-residual plot (unclipped, n={len(unclipped)})\n"
                 f"partial-r = {pr:.3f}")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "syn1_partial_residual_halflife.pdf")
    plt.close(fig)
    print(f"Wrote {OUT_DIR/'syn1_partial_residual_halflife.pdf'}")


if __name__ == "__main__":
    main()
