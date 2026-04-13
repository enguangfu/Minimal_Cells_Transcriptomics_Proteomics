#!/usr/bin/env python
# coding: utf-8
"""
Transcription–Translation Correlation for Syn1
===============================================
Correlates Syn1 proteomics (iPM, intensity-based protein molecules) with
transcriptomics (Illumina / PacBio sense TPM) at the gene level.

Key findings:
  - Pearson r ~ 0.6 between Illumina TPM and iPM.
  - Pearson r ~ 0.5 between PacBio TPM and iPM.
  - Cytosolic-only subset: r ~ 0.67 (slightly higher).

iPM (intensity-based Protein Molecules) is derived from iBAQ and scales
more linearly with true protein copy number than raw iBAQ.
"""

# =============================================================================
# 1. Imports
# =============================================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

# =============================================================================
# 2. Paths
# =============================================================================

HOME_DIR = ".."
PROTEOME_FOLDER   = HOME_DIR + "/Proteomics"
TRANSCRIPTOME_CSV = "./syn1_Illumina_PacBio_TPM_profiles.csv"
PROTEOME_CSV      = PROTEOME_FOLDER + "/syn1_proteomics_localization_2026.csv"

# Which transcriptome column to correlate against iPM:
#   "avg_sense_TPM"   → merged Illumina TPM
#   "PacBio_sense_TPM" → PacBio Iso-Seq TPM
TRANSCRIPTOME_COL = "avg_sense_TPM"

# =============================================================================
# 3. Load data
# =============================================================================

# --- 3a. Transcriptomics (Illumina + PacBio TPM profiles) ---
genes = pd.read_csv(TRANSCRIPTOME_CSV)
print(f"Transcriptomics loaded: {len(genes)} genes")

# --- 3b. Proteomics (iPM + localization) ---
prot = pd.read_csv(PROTEOME_CSV)
print(f"Proteomics loaded: {len(prot)} entries")
print(f"  Columns: {prot.columns.tolist()}")
print(f"  Localization categories: {sorted(prot['ptn_localization'].dropna().unique())}")

# =============================================================================
# 4. Merge proteomics onto genes dataframe
# =============================================================================

# Map iPM replicates, mean, CV, and localization onto genes by locus_tag
iPM_cols = ['iPM_rep1', 'iPM_rep2', 'iPM_rep3', 'iPM_mean', 'iPM_CV', 'ptn_localization']
iPM_map  = prot.set_index('locus_tag')[iPM_cols]

for col in iPM_cols:
    genes[col] = genes['locus_tag'].map(iPM_map[col])

# Coverage summary
n_mrna    = (genes['rna_type'] == 'mRNA').sum()
n_detected = genes.loc[genes['rna_type'] == 'mRNA', 'iPM_mean'].notna().sum()
print(f"\nProtein detection coverage (mRNA genes only):")
print(f"  Total mRNA genes : {n_mrna}")
print(f"  With iPM data    : {n_detected} ({n_detected / n_mrna * 100:.1f}%)")
print(f"  Without iPM data : {n_mrna - n_detected}")

# Save merged table
genes.to_csv("syn1_genes_transcriptomics_proteomics.csv", index=False)
print("\nSaved: syn1_genes_transcriptomics_proteomics.csv")

# =============================================================================
# 5. Correlation: iPM vs TPM
# =============================================================================

# Filter to genes with both iPM and TPM > 0 (all localization categories)
loc_categories = set(genes['ptn_localization'].dropna().unique())
mask = (
    (genes['iPM_mean'] > 0) &
    (genes[TRANSCRIPTOME_COL] > 0) &
    genes['ptn_localization'].isin(loc_categories)
)

corr_df = genes.loc[mask].copy()
corr_df['log10_iPM'] = np.log10(corr_df['iPM_mean'])
corr_df['log10_TPM'] = np.log10(corr_df[TRANSCRIPTOME_COL])

rho, p_spearman = spearmanr(corr_df['log10_iPM'], corr_df['log10_TPM'])
r,   p_pearson  = pearsonr(corr_df['log10_iPM'],  corr_df['log10_TPM'])

print(f"\nCorrelation: iPM vs {TRANSCRIPTOME_COL}")
print(f"  Localization filter : {', '.join(sorted(loc_categories))}")
print(f"  N genes             : {len(corr_df)}")
print(f"  Spearman ρ          : {rho:.4f}  (p = {p_spearman:.2e})")
print(f"  Pearson  r          : {r:.4f}  (p = {p_pearson:.2e})")

# =============================================================================
# 6. Scatter plot — log10(iPM) vs log10(TPM)
# =============================================================================
# Two panels side by side (same 6×6 per panel):
#   Left  — all detected proteins
#   Right — cytosolic proteins only (ptn_localization == 'cytoplasmic')

def _scatter_panel(ax, df, rho, r, p_rho, p_r):
    """Draw one correlation scatter panel with OLS fit and stats box."""
    z      = np.polyfit(df['log10_TPM'], df['log10_iPM'], 1)
    x_line = np.linspace(df['log10_TPM'].min(), df['log10_TPM'].max(), 100)

    ax.scatter(df['log10_TPM'], df['log10_iPM'],
               alpha=0.55, s=25, color='#4C8BB5',
               edgecolors='white', linewidths=0.3, zorder=2)
    ax.plot(x_line, np.polyval(z, x_line),
            color='#E05A5A', linewidth=1.8, alpha=0.9, zorder=3)

    ax.set_xlabel('log$_{10}$(TPM)  —  Transcriptomics', fontsize=12)
    ax.set_ylabel('log$_{10}$(iPM)  —  Proteomics', fontsize=12)

    # Stats annotation box
    ax.text(0.04, 0.97,
            f'Spearman ρ = {rho:.3f}  (p = {p_rho:.1e})\n'
            f'Pearson  r = {r:.3f}  (p = {p_r:.1e})\n'
            f'n = {len(df)} genes',
            transform=ax.transAxes, fontsize=9.5, verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.45', facecolor='#FFF8E7',
                      edgecolor='#CCBBAA', alpha=0.85))

    # Annotate top 5 positive and negative residual outliers
    # Label with the numeric suffix of locus_tag (e.g. MMSYN1_0042 → 0042)
    df = df.copy()
    df['residual'] = df['log10_iPM'] - np.polyval(z, df['log10_TPM'])
    outliers = pd.concat([df.nlargest(5, 'residual'), df.nsmallest(5, 'residual')])
    for _, row in outliers.iterrows():
        locus_num = row['locus_tag'].split('_')[-1]
        ax.annotate(locus_num, (row['log10_TPM'], row['log10_iPM']),
                    fontsize=7, color='#333333', alpha=0.85,
                    xytext=(5, 4), textcoords='offset points')

    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(True, alpha=0.25, linestyle='--')
    return df  # return with residual column attached


# --- All proteins ---
rho_all, p_rho_all = spearmanr(corr_df['log10_iPM'], corr_df['log10_TPM'])
r_all,   p_r_all   = pearsonr( corr_df['log10_iPM'], corr_df['log10_TPM'])

# --- Cytosolic proteins only ---
cyto_df = corr_df[corr_df['ptn_localization'] == 'cytoplasmic'].copy()
rho_cyto, p_rho_cyto = spearmanr(cyto_df['log10_iPM'], cyto_df['log10_TPM'])
r_cyto,   p_r_cyto   = pearsonr( cyto_df['log10_iPM'], cyto_df['log10_TPM'])

print(f"\nCorrelation — all proteins    : Spearman ρ={rho_all:.4f}, Pearson r={r_all:.4f}, n={len(corr_df)}")
print(f"Correlation — cytosolic only  : Spearman ρ={rho_cyto:.4f}, Pearson r={r_cyto:.4f}, n={len(cyto_df)}")

# All proteins
fig, ax = plt.subplots(figsize=(6, 6))
corr_df = _scatter_panel(ax, corr_df, rho_all, r_all, p_rho_all, p_r_all)
plt.tight_layout()
plt.savefig('./protein_vs_mRNA_correlation_all.pdf', dpi=300, bbox_inches='tight')
# plt.show()

# Cytosolic proteins only
fig, ax = plt.subplots(figsize=(6, 6))
cyto_df = _scatter_panel(ax, cyto_df, rho_cyto, r_cyto, p_rho_cyto, p_r_cyto)
plt.tight_layout()
plt.savefig('./protein_vs_mRNA_correlation_cytosolic.pdf', dpi=300, bbox_inches='tight')
# plt.show()

# =============================================================================
# 7. Residual analysis — genes deviating most from the mRNA–protein fit
# =============================================================================
# Positive residual: more protein than expected from mRNA
#   → high translational efficiency, or stable/long-lived protein
# Negative residual: less protein than expected from mRNA
#   → poor translation, rapid protein degradation, or low proteomics coverage

cols_show = ['locus_tag', 'gene_name', 'gene_product', 'log10_TPM', 'log10_iPM', 'residual']

print("\n=== Top 15: MORE protein than expected from mRNA ===")
print("(high translational efficiency or stable protein)\n")
print(corr_df.nlargest(15, 'residual')[cols_show].reset_index(drop=True).to_string())

print("\n=== Top 15: LESS protein than expected from mRNA ===")
print("(poor translation, rapid degradation, or low proteomics coverage)\n")
print(corr_df.nsmallest(15, 'residual')[cols_show].reset_index(drop=True).to_string())

# =============================================================================
# 8. Protein copy number distribution — linear scale
# =============================================================================
# ptn_copy_number: absolute protein molecules per cell estimated from iPM.
# Excludes Non_proteins (tRNAs, rRNAs, etc.) and missing/zero values.

LOC_CATEGORIES = ['cytoplasmic', 'membrane', 'lipoprotein', 'extracellular']

plot_df = prot[
    prot['ptn_localization'].isin(LOC_CATEGORIES) &
    prot['ptn_copy_number'].notna() &
    (prot['ptn_copy_number'] > 0)
].copy()

median_copy = plot_df['ptn_copy_number'].median()

fig, ax = plt.subplots(figsize=(6, 3))
ax.hist(plot_df['ptn_copy_number'], bins=50, color='steelblue',
        edgecolor='white', linewidth=0.3)
ax.axvline(median_copy, color='red', linewidth=1.5, linestyle='--',
           label=f'Median = {median_copy:.1f}')
ax.set_xlabel('Protein copy number per cell', fontsize=12)
ax.set_ylabel('Number of unique proteins', fontsize=12)
# ax.set_title('Syn1 Protein Copy Number Distribution (linear)', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('./protein_copy_number_linear.pdf', dpi=300, bbox_inches='tight')
# plt.show()

# =============================================================================
# 9. Protein copy number distribution — log10 scale
# =============================================================================

median_log10 = np.log10(median_copy)

fig, ax = plt.subplots(figsize=(6, 3))
ax.hist(np.log10(plot_df['ptn_copy_number']), bins=50, color='steelblue',
        edgecolor='white', linewidth=0.3)
ax.axvline(median_log10, color='red', linewidth=1.5, linestyle='--',
           label=f'Median = {median_copy:.1f} ({median_log10:.2f} in log$_{{10}}$)')
ax.set_xlabel('log$_{10}$(protein copy number per cell)', fontsize=12)
ax.set_ylabel('Number of uniqueproteins', fontsize=12)
# ax.set_title('Syn1 Protein Copy Number Distribution (log$_{10}$)', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('./protein_copy_number_log10.pdf', dpi=300, bbox_inches='tight')
# plt.show()

# =============================================================================
# 10. Copy number statistics by subcellular localization
# =============================================================================

stats_rows = []
for loc in LOC_CATEGORIES:
    subset = plot_df.loc[plot_df['ptn_localization'] == loc, 'ptn_copy_number']
    if subset.empty:
        continue
    stats_rows.append({
        'localization' : loc,
        'n'            : len(subset),
        'min'          : subset.min(),
        'Q1'           : subset.quantile(0.25),
        'median'       : subset.median(),
        'mean'         : subset.mean(),
        'Q3'           : subset.quantile(0.75),
        'max'          : subset.max(),
        'CV'           : subset.std() / subset.mean(),
    })

stats_df = pd.DataFrame(stats_rows)
print("\nProtein copy number statistics by localization (molecules per cell):")
print(stats_df.to_string(index=False, float_format=lambda x: f'{x:.2f}'))

stats_df.to_csv('./protein_copy_number_stats_by_localization.csv', index=False, float_format='%.4f')
print("\nSaved: protein_copy_number_stats_by_localization.csv")
