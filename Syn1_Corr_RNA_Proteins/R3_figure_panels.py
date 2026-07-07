#!/usr/bin/env python
# coding: utf-8
"""
R3 figure panels (transcriptome-proteome correlation in Syn1).

Regenerates the Syn1 R3 panels at final print sizes, reading the already-computed
analysis tables rather than recomputing. Born-at-size, OUTPUT.md default fonts
(7 pt labels / 6 pt ticks). Every panel carries a bold blue "Syn1.0" tag (the
Fig 5/6 organism convention); the syn3A twin (Corr_RNA_Protein_Syn3A.py) is red.

MAIN figure (correlation.pdf), syn1 row:
  b  7/3 x 7/3  Illumina TPM vs iPM (log10), FILLED circles by localization
                <- residual_analysis/gene_CAI_omics_merged.csv
  a  7/3 x 7/6  copy-number distribution by localization (4 colors)
                <- Syn1_Syn3A_Proteomics/syn1_proteomics_localization_2026.csv
  g  7/3 x 7/6  intrinsic half-life distribution (Mpn-transferred), blue bars
                <- residual_analysis/syn1_ptn_degradation_from_mpn.csv

SI figure (si-correlation.pdf), syn1 row:
  d  7/4 x 7/4  proteome residual vs log10(TIR)
  e  7/4 x 7/4  proteome residual vs CAI
  h  7/4 x 7/4  proteome residual vs log10(half-life)
  f  7/4 x 7/4  model R (baseline vs +CAI) for all and cytosolic proteins

The PacBio-vs-Illumina correlation and the length/abundance-bias panels moved to
the R0 platform-comparison figure (RNAseq_Comparison/, si-rnaseq.pdf).

Output: R3_panels/panel_{a,b,d,e,f,g,h}.pdf + syn1_omics.xlsx + R3_figure.txt
Run from Syn1_Corr_RNA_Proteins/.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, gaussian_kde
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    'font.size': 7,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Nimbus Sans', 'Helvetica', 'DejaVu Sans'],
    'axes.titlesize': 7,
    'axes.labelsize': 7,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

HALF, QUART = 7 / 2, 7 / 4          # 3.5 in, 1.75 in (SI residual panels)
W    = 7 / 3                        # 2.333 in -- main-panel width
CORR = (W, W)                       # 7/3 x 7/3  (TPM vs iPM)
DIST = (W, 7 / 6)                   # 7/3 x 7/6  (copy-number + half-life)
OUT = 'R3_panels'
GENE_TPM = '../Syn1_Transcriptomics/Gene_TPM/syn1_Illumina_PacBio_TPM_profiles.csv'
os.makedirs(OUT, exist_ok=True)

# Four distinct, consistent localization colors (panels a, b)
LOC_ORDER  = ['cytoplasmic', 'lipoprotein', 'membrane', 'extracellular']
LOC_LABEL  = {'cytoplasmic': 'Cytoplasmic', 'lipoprotein': 'Lipoprotein',
              'membrane': 'Membrane', 'extracellular': 'Extracellular'}
# Okabe-Ito colorblind-safe palette, four maximally distinct hues (panels a, b)
LOC_COLORS = {'cytoplasmic': '#0072B2',    # blue
              'lipoprotein': '#009E73',    # green
              'membrane': '#D55E00',       # vermillion
              'extracellular': '#CC79A7'}  # reddish purple

# Organism identity colors (match Fig 5/6): Syn1 blue, Syn3A red. Used only as a
# bold left-title tag on every panel and as the half-life bar fill -- never for
# the localization-coloured points.
SYN1_COL, SYN3A_COL = '#3182bd', '#c0392b'
def org_tag(ax, which='syn1'):
    name, col = ('Syn1.0', SYN1_COL) if which == 'syn1' else ('Syn3A', SYN3A_COL)
    ax.set_title(name, loc='left', color=col, fontweight='bold', fontsize=7)

log = []
def say(s):
    print(s); log.append(s)

def fit_line(ax, x, y, color='black'):
    """OLS line over the x-range."""
    b, a = np.polyfit(x, y, 1)
    xs = np.array([x.min(), x.max()])
    ax.plot(xs, b * xs + a, color=color, lw=1.0, zorder=5)

def r2_lstsq(X, y):
    """R^2 of OLS y ~ [1, X]."""
    A = np.column_stack([np.ones(len(y))] + [X[:, i] for i in range(X.shape[1])])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot

# ====================================================================== a
prot = pd.read_csv('../Syn1_Syn3A_Proteomics/syn1_proteomics_localization_2026.csv')
prot = prot[(prot['ptn_copy_number'] > 0) & prot['ptn_localization'].isin(LOC_ORDER)]
fig, ax = plt.subplots(figsize=DIST, constrained_layout=True)
allv = np.log10(prot['ptn_copy_number'].values)
lo, hi = allv.min(), allv.max()
xgrid = np.linspace(lo, hi, 300)
bw = (hi - lo) / 25.0                       # reference bin width: KDE -> protein counts
for loc in LOC_ORDER:
    v = np.log10(prot.loc[prot['ptn_localization'] == loc, 'ptn_copy_number'].values)
    if len(v) < 5:
        continue
    n = len(v); med = float(np.median(v)); med_cn = 10 ** med
    kde = gaussian_kde(v)
    ycurve = kde(xgrid) * n * bw            # area under curve = protein count for this loc
    ax.plot(xgrid, ycurve, color=LOC_COLORS[loc], lw=1.2,
            label=f"{LOC_LABEL[loc]} ({med_cn:.0f}, n={n})")
    ax.fill_between(xgrid, ycurve, color=LOC_COLORS[loc], alpha=0.13)
    ax.vlines(med, 0, float(kde(med)[0]) * n * bw, color=LOC_COLORS[loc], lw=0.9, ls='--')  # median
    say(f"a) {loc}: n={n} median={med_cn:.1f} copies")
ax.set_xlabel('Copies per cell ($\\log_{10}$)')
ax.set_ylabel('Proteins')
ax.set_xlim(lo, hi)
# compact 5 pt legend in the empty top-left corner; title flags the parenthetical
# leading number as the median copies/cell (also drawn as the dashed vertical lines)
leg = ax.legend(frameon=False, handlelength=0.8, labelspacing=0.2, borderpad=0.2,
                fontsize=5, loc='upper left',
                title='(median copies, n=unique proteins)', title_fontsize=5)
leg._legend_box.align = 'left'
ax.spines[['top', 'right']].set_visible(False)
org_tag(ax, 'syn1')
fig.savefig(f'{OUT}/panel_a_copynumber_by_localization.pdf', dpi=300); plt.close(fig)

# ====================================================================== b
cai = pd.read_csv('residual_analysis/gene_CAI_omics_merged.csv')
fig, ax = plt.subplots(figsize=CORR, constrained_layout=True)
for loc in LOC_ORDER:
    s = cai[cai['ptn_localization'] == loc]
    ax.scatter(s['log10_TPM'], s['log10_iPM'], s=9, alpha=0.7,          # filled = syn1
               c=LOC_COLORS[loc], edgecolors='none', label=f"{LOC_LABEL[loc]} (n={len(s)})")
r_all = pearsonr(cai['log10_TPM'], cai['log10_iPM'])[0]
cyto = cai[cai['ptn_localization'] == 'cytoplasmic']
r_cyto = pearsonr(cyto['log10_TPM'], cyto['log10_iPM'])[0]
fit_line(ax, cai['log10_TPM'].values, cai['log10_iPM'].values)
ax.text(0.04, 0.96, f"all $r$ = {r_all:.2f} ($n$ = {len(cai)})\n"
                    f"cytoplasmic $r$ = {r_cyto:.2f} ($n$ = {len(cyto)})",
        transform=ax.transAxes, va='top', ha='left', fontsize=6)
ax.set_xlabel('mRNA Illumina TPM ($\\log_{10}$)')
ax.set_ylabel('Protein iPM ($\\log_{10}$)')
ax.legend(frameon=False, handlelength=1.0, labelspacing=0.25, loc='lower right')
ax.spines[['top', 'right']].set_visible(False)
org_tag(ax, 'syn1')
fig.savefig(f'{OUT}/panel_b_TPM_vs_iPM.pdf', dpi=300); plt.close(fig)
say(f"b) all r={r_all:.3f} (n={len(cai)}); cytoplasmic r={r_cyto:.3f} (n={len(cyto)}); R^2 NOT shown")

# TPM profiles retained only for the syn1_omics.xlsx workbook below. The
# PacBio-vs-Illumina correlation and the length/abundance-bias panels moved to
# the R0 platform-comparison figure (RNAseq_Comparison/, si-rnaseq.pdf).
tpm = pd.read_csv(GENE_TPM)

# ====================================================================== d
tir = pd.read_csv('residual_analysis/gene_TIR_omics_merged.csv')
r_d = pearsonr(tir['log10_TIR'], tir['residual'])[0]
fig, ax = plt.subplots(figsize=(QUART, QUART), constrained_layout=True)
ax.scatter(tir['log10_TIR'], tir['residual'], s=5, alpha=0.5, c=SYN1_COL, edgecolors='none')
fit_line(ax, tir['log10_TIR'].values, tir['residual'].values, color='black')
ax.axhline(0, color='black', lw=0.5, ls=':')
ax.text(0.04, 0.96, f"$r$ = {r_d:.2f}\n$n$ = {len(tir)}", transform=ax.transAxes, va='top', fontsize=6)
ax.set_xlabel('TIR ($\\log_{10}$)')
ax.set_ylabel('Proteome residual')
ax.spines[['top', 'right']].set_visible(False)
org_tag(ax, 'syn1')
fig.savefig(f'{OUT}/panel_d_TIR_vs_residual.pdf', dpi=300); plt.close(fig)
say(f"d) residual vs log10(TIR) r={r_d:.3f} (n={len(tir)})")

# ====================================================================== e
r_e = pearsonr(cai['CAI'], cai['proteome_residual'])[0]
fig, ax = plt.subplots(figsize=(QUART, QUART), constrained_layout=True)
ax.scatter(cai['CAI'], cai['proteome_residual'], s=5, alpha=0.5, c=SYN1_COL, edgecolors='none')
fit_line(ax, cai['CAI'].values, cai['proteome_residual'].values, color='black')
ax.axhline(0, color='black', lw=0.5, ls=':')
ax.text(0.04, 0.96, f"$r$ = {r_e:.2f}\n$n$ = {len(cai)}", transform=ax.transAxes, va='top', fontsize=6)
ax.set_xlabel('CAI')
ax.set_ylabel('Proteome residual')
ax.spines[['top', 'right']].set_visible(False)
org_tag(ax, 'syn1')
fig.savefig(f'{OUT}/panel_e_CAI_vs_residual.pdf', dpi=300); plt.close(fig)
say(f"e) residual vs CAI r={r_e:.3f} (n={len(cai)})")

# ====================================================================== f
def r2_pair(df):
    y = df['log10_iPM'].values
    base = r2_lstsq(df[['log10_TPM']].values, y)
    full = r2_lstsq(df[['log10_TPM', 'CAI']].values, y)
    return base, full
b_all, f_all = r2_pair(cai)
b_cyt, f_cyt = r2_pair(cyto)
R = np.sqrt   # multiple correlation coefficient R = sqrt(R^2)
base_R = [R(b_all), R(b_cyt)]
cai_R  = [R(f_all), R(f_cyt)]
fig, ax = plt.subplots(figsize=(QUART, QUART), constrained_layout=True)
xpos = np.array([0, 1])
w = 0.38
bars1 = ax.bar(xpos - w/2, base_R, w, color='#bbbbbb', label='TPM only')
bars2 = ax.bar(xpos + w/2, cai_R,  w, color=SYN1_COL, label='+ CAI')
for bars in (bars1, bars2):
    for bar in bars:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.012,
                f"{bar.get_height():.2f}", ha='center', va='bottom', fontsize=5)
ax.set_xticks(xpos); ax.set_xticklabels([f'All\n(n={len(cai)})', f'Cytosolic\n(n={len(cyto)})'])
ax.set_ylabel('$R$')
ax.set_ylim(0, 0.95)
ax.legend(frameon=False, handlelength=1.0, labelspacing=0.25, loc='upper left')
ax.spines[['top', 'right']].set_visible(False)
org_tag(ax, 'syn1')
fig.savefig(f'{OUT}/panel_f_R_improvement.pdf', dpi=300); plt.close(fig)
say(f"f) R all {R(b_all):.3f}->{R(f_all):.3f}; cytosolic {R(b_cyt):.3f}->{R(f_cyt):.3f} "
    f"(R^2 all {b_all:.3f}->{f_all:.3f}, cyto {b_cyt:.3f}->{f_cyt:.3f})")

# ====================================================================== g
deg = pd.read_csv('residual_analysis/syn1_ptn_degradation_from_mpn.csv')
hl = deg['halflife_h_syn1'].dropna()
g_med, g_min, g_n = hl.median(), hl.min(), len(hl)
fig, ax = plt.subplots(figsize=DIST, constrained_layout=True)
ax.hist(np.log10(hl), bins=20, color=SYN1_COL, edgecolor='white', linewidth=0.3)
ax.axvline(np.log10(g_med), color='black', lw=1.0, ls='--', label=f'median {g_med:.0f} h')
ax.axvline(np.log10(g_min), color='#888888', lw=1.0, ls=':', label=f'shortest {g_min:.1f} h')
ax.set_xlabel('Half-life, h ($\\log_{10}$)')
ax.set_ylabel('Proteins')
ax.legend(frameon=False, handlelength=1.2, labelspacing=0.25, loc='upper right',
          title=f'n = {g_n}', title_fontsize=6)
ax.spines[['top', 'right']].set_visible(False)
org_tag(ax, 'syn1')
fig.savefig(f'{OUT}/panel_g_halflife_distribution.pdf', dpi=300); plt.close(fig)
say(f"g) half-life median={g_med:.1f} h, shortest={g_min:.1f} h, n={g_n}")

# ====================================================================== h
hjoin = deg.merge(cai[['locus_tag', 'proteome_residual']],
                  left_on='syn1_locus_tag', right_on='locus_tag', how='inner').dropna(
                  subset=['halflife_h_syn1', 'proteome_residual'])
xh = np.log10(hjoin['halflife_h_syn1']); yh = hjoin['proteome_residual']
r_h = pearsonr(xh, yh)[0]
fig, ax = plt.subplots(figsize=(QUART, QUART), constrained_layout=True)
ax.scatter(xh, yh, s=5, alpha=0.5, c=SYN1_COL, edgecolors='none')
fit_line(ax, xh.values, yh.values, color='black')
ax.axhline(0, color='black', lw=0.5, ls=':')
ax.text(0.04, 0.96, f"$r$ = {r_h:.2f}\n$n$ = {len(hjoin)}", transform=ax.transAxes, va='top', fontsize=6)
ax.set_xlabel('Half-life, h ($\\log_{10}$)')
ax.set_ylabel('Proteome residual')
ax.spines[['top', 'right']].set_visible(False)
org_tag(ax, 'syn1')
fig.savefig(f'{OUT}/panel_h_halflife_vs_residual.pdf', dpi=300); plt.close(fig)
say(f"h) residual vs log10(half-life) r={r_h:.3f} (n={len(hjoin)})")

# ====================================================================== workbook
# syn1_omics.xlsx -- one row per gene (all 911), gene metadata + every R3 measurement.
# Base = proteomics table (all genes + localization/iPM/copy number); left-join the
# Illumina/PacBio TPM, OSTIR TIR, CAI, and Mpn-transferred half-life by locus tag.
base   = pd.read_csv('../Syn1_Syn3A_Proteomics/syn1_proteomics_localization_2026.csv')
tir_t  = pd.read_csv('residual_analysis/gene_TIR.csv')[['locus_tag', 'TIR_weighted']]
cai_t  = pd.read_csv('residual_analysis/gene_CAI.csv')[['locus_tag', 'CAI']]
deg_t  = deg[['syn1_locus_tag', 'halflife_h_syn1']].rename(columns={'syn1_locus_tag': 'locus_tag'})

wb = (base[['locus_tag', 'gene_name', 'rna_type', 'gene_product',
            'ptn_localization', 'iPM_mean', 'ptn_copy_number']]
      .merge(tpm[['locus_tag', 'avg_sense_TPM', 'PacBio_sense_TPM']], on='locus_tag', how='left')
      .merge(tir_t, on='locus_tag', how='left')
      .merge(cai_t, on='locus_tag', how='left')
      .merge(deg_t, on='locus_tag', how='left'))

# Prefer the curated syn3A annotation for gene name + product where a syn3A
# ortholog exists (numeric locus suffix preserved: MMSYN1_NNNN <-> JCVISYN3A_NNNN);
# fall back to the syn1 annotation for genes deleted in syn3A.
syn3a = pd.read_excel('../Syn1_Syn3A_Proteomics/syn3A_proteome_annotated.xlsx',
                      sheet_name='Proteome')[['Locus Tag', 'Gene Name', 'Gene Product']]
syn3a['suffix'] = syn3a['Locus Tag'].str.extract(r'_(\d+)$')
syn3a = syn3a.dropna(subset=['suffix']).drop_duplicates('suffix').set_index('suffix')
suf = wb['locus_tag'].str.extract(r'_(\d+)$')[0]
gn3, gp3 = suf.map(syn3a['Gene Name']), suf.map(syn3a['Gene Product'])
n_upd = int(gn3.notna().sum())
wb['gene_name']    = gn3.fillna(wb['gene_name'])
wb['gene_product'] = gp3.fillna(wb['gene_product'])

copy_up = wb['ptn_copy_number'].apply(lambda x: pd.NA if pd.isna(x) else int(np.ceil(x))).astype('Int64')
wb_out = pd.DataFrame({
    'locusTag':              wb['locus_tag'],
    'gene_name':             wb['gene_name'],
    'rna_type':              wb['rna_type'],
    'gene_product':          wb['gene_product'],
    'protein_localization':  wb['ptn_localization'],
    'TPM_illumina':          wb['avg_sense_TPM'].round(2),
    'TPM_PacBio':            wb['PacBio_sense_TPM'].round(2),
    'iPM_mean':              wb['iPM_mean'].round(2),
    'protein_copy_number':   copy_up,
    'TIR':                   wb['TIR_weighted'].round(1),
    'CAI':                   wb['CAI'].round(3),
    'protein_halflife_h':    wb['halflife_h_syn1'].round(1),
})
wb_out.to_excel('syn1_omics.xlsx', index=False, sheet_name='syn1_genes')
say(f"\nsyn1_omics.xlsx: {len(wb_out)} genes "
    f"(TPM_illumina {wb_out['TPM_illumina'].notna().sum()}, TPM_PacBio {wb_out['TPM_PacBio'].notna().sum()}, "
    f"iPM {wb_out['iPM_mean'].notna().sum()}, copy {wb_out['protein_copy_number'].notna().sum()}, "
    f"TIR {wb_out['TIR'].notna().sum()}, CAI {wb_out['CAI'].notna().sum()}, "
    f"half-life {wb_out['protein_halflife_h'].notna().sum()})")

with open(f'{OUT}/R3_figure.txt', 'w') as fh:
    fh.write("R3 FIGURE PANELS (transcriptome-proteome correlation)\n")
    fh.write("=" * 60 + "\n")
    fh.write("MAIN: b 7/3 x 7/3; a,g 7/3 x 7/6.  SI: d,e,f,h 7/4 x 7/4. Default fonts.\n\n")
    fh.write("\n".join(log) + "\n")
print(f"\nSaved 7 panels (3 main + 4 SI) + syn1_omics.xlsx + R3_figure.txt to {OUT}/")
