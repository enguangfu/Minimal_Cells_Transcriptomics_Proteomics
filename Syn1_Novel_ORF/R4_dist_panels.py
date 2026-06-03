#!/usr/bin/env python
# coding: utf-8
"""
R4 distribution panels (b, c, f, g) for the novel-transcription figure.

Reads already-computed tables (no recompute of the heavy pipelines) and emits
each panel born-at-size per OUTPUT.md (Arial, 7/6/5 pt, pdf.fonttype 42).

Panel -> size (in) -> content -> source   (panel letters per reordered MANUSCRIPT.md)
  b  7/4 x 7/4  isoform LENGTH ridgeline by antisense case  <- isoform_antisense_categories.xlsx
  c  7/4 x 7/4  isoform CLUSTER-READ ridgeline by case      <- isoform_antisense_categories.xlsx
  e  7/4 x 7/4  intergenic-coverage distribution, all isoforms
                <- ../Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv
  f  7/4 x 7/4  5'/3' UTR length distribution (canonical operons)
                <- ../Syn1_Operon/operons.candidate_blocks.tsv + GFF (UTR computed
                   inline, replicating Operon_Annotation.py's canonical-operon formula)

Output: R4_panels/panel_{b,c,e,f}.pdf  + R4_dist_panels.txt
Run from Syn1_Novel_ORF/.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    'font.size': 7, 'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Nimbus Sans', 'Helvetica', 'DejaVu Sans'],
    'axes.titlesize': 7, 'axes.labelsize': 7,
    'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 6,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

QUART = 7 / 4
OUT = 'R4_panels'
os.makedirs(OUT, exist_ok=True)
CAT_XLSX  = 'isoform_antisense_categories.xlsx'
ISO_TSV   = '../Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv'
OPERON_TSV = '../Syn1_Operon/operons.candidate_blocks.tsv'
GFF       = '../Genomes_Input/syn1.genes.gff3'
MIN_READS = 10

# Three antisense cases (Okabe-Ito), ordered by abundance
CASE_ORDER  = ['spurious_prom', 'read_through', 'embedded']
CASE_LABEL  = {'spurious_prom': 'Spurious promoter', 'read_through': 'Read-through', 'embedded': 'Embedded'}
CASE_COLOR  = {'spurious_prom': '#0072B2', 'read_through': '#D55E00', 'embedded': '#009E73'}
UTR_COLOR   = {"5'": '#E69F00', "3'": '#56B4E9'}

log = []
def say(s):
    print(s); log.append(s)


def ridgeline(ax, values_by_case, xgrid, x_is_log10=False, scale=0.7,
              med_fmt=lambda m: f"{m:.0f}"):
    """Stacked KDE ridges, one per case (top->bottom = CASE_ORDER), shared x.
    Each ridge carries a rug of the actual points + a dashed median line, and the
    median VALUE is printed at the right (no case names -- panel a carries those --
    and no x ticks). Small-n cases (embedded, n=4) keep their rug so they are not
    over-read. xgrid is in plotting units (log10 of the value when x_is_log10)."""
    n = len(CASE_ORDER)
    for i, case in enumerate(CASE_ORDER):
        v = np.asarray(values_by_case[case], float)
        v = v[v > 0]
        base = (n - 1 - i)                       # first case = top baseline
        col = CASE_COLOR[case]
        vt = np.log10(v) if x_is_log10 else v
        if len(v) >= 3 and np.ptp(vt) > 0:
            dens = gaussian_kde(vt)(xgrid)
            dens = dens / dens.max() * scale
        else:
            dens = np.zeros_like(xgrid)
        ax.hlines(base, xgrid[0], xgrid[-1], color='0.75', lw=0.4, zorder=0)
        ax.fill_between(xgrid, base, base + dens, color=col, alpha=0.50, lw=0)
        ax.plot(xgrid, base + dens, color=col, lw=1.0)
        ax.plot(vt, np.full(len(vt), base), '|', color=col, ms=4, mew=0.6)   # rug
        med = float(np.median(v))
        med_plot = np.log10(med) if x_is_log10 else med
        ax.plot([med_plot, med_plot], [base, base + scale], color=col, lw=0.8, ls='--', alpha=0.85)
        ax.text(xgrid[-1], base + 0.30, med_fmt(med), ha='right', va='center',   # median on the right
                fontsize=5, fontweight='bold', color=col)
    ax.text(0.995, 0.995, 'Median', transform=ax.transAxes, ha='right', va='top',
            fontsize=5, color='0.3')
    ax.set_yticks([])                                # keep x ticks (the value scale)
    ax.set_ylim(-0.1, n - 1 + scale + 0.20)
    ax.spines[['top', 'right', 'left']].set_visible(False)


# ====================================================================== b, c (ridgelines)
cat = pd.read_excel(CAT_XLSX)
vals_len     = {c: cat.loc[cat['antisense_category'] == c, 'isoform_len_bp'].dropna().values / 1000.0
                for c in CASE_ORDER}
vals_cluster = {c: cat.loc[cat['antisense_category'] == c, 'cluster_reads'].dropna().values
                for c in CASE_ORDER}

# --- b: isoform length ridgeline (kb, shared linear x; no ticks, median on right)
allv = np.concatenate([vals_len[c] for c in CASE_ORDER])
xg = np.linspace(allv.min() * 0.85, allv.max() * 1.03, 200)
fig, ax = plt.subplots(figsize=(QUART, QUART), constrained_layout=True)
ridgeline(ax, vals_len, xg, x_is_log10=False, med_fmt=lambda m: f"{m:.1f} kb")
ax.set_xlabel('Isoform length (kb)')
fig.savefig(f'{OUT}/panel_b_length_ridge.pdf', dpi=300); plt.close(fig)
for c in CASE_ORDER:
    say(f"b) {c}: n={len(vals_len[c])} median_len={np.median(vals_len[c])*1000:.0f} bp")

# --- c: cluster-read ridgeline (log10 reads, shared x; no ticks, median on right)
allr = np.concatenate([vals_cluster[c] for c in CASE_ORDER])
xg = np.linspace(np.log10(allr.min()) - 0.05, np.log10(allr.max()) + 0.05, 200)
fig, ax = plt.subplots(figsize=(QUART, QUART), constrained_layout=True)
ridgeline(ax, vals_cluster, xg, x_is_log10=True, med_fmt=lambda m: f"{m:.0f}")
ax.set_xticks([1, 2, 3, 4]); ax.set_xticklabels(['10', '100', '1k', '10k'])
ax.set_xlabel('Isoform cluster reads')
fig.savefig(f'{OUT}/panel_c_reads_ridge.pdf', dpi=300); plt.close(fig)
for c in CASE_ORDER:
    say(f"c) {c}: n={len(vals_cluster[c])} median_cluster_reads={np.median(vals_cluster[c]):.0f} max={vals_cluster[c].max():.0f}")


# ====================================================================== f
# Canonical-operon 5'/3' UTR lengths. Replicates Operon_Annotation.py:
#   canonical = TSS and TTS both intergenic (not inside any same-strand gene body)
#   + strand: 5'UTR = first_gene.start0 - op.start0 ; 3'UTR = op.end0 - last_gene.end0
#   - strand: 5'UTR = op.end0 - last_gene.end0      ; 3'UTR = first_gene.start0 - op.start0
def load_gff(path):
    rows = []
    for line in open(path):
        if not line.strip() or line.startswith('#'):
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) != 9 or p[2] != 'gene':
            continue
        attrs = dict(kv.split('=', 1) for kv in p[8].split(';') if '=' in kv)
        rows.append({'locus_tag': attrs.get('locus_tag', ''), 'chrom': p[0],
                     'start0': int(p[3]) - 1, 'end0': int(p[4]), 'strand': p[6]})
    return pd.DataFrame(rows)

genes = load_gff(GFF).drop_duplicates('locus_tag')
gene_by_locus = genes.set_index('locus_tag')
# same-strand interval lists for intergenic test
strand_intervals = {s: g[['start0', 'end0']].values for s, g in genes.groupby('strand')}

def is_intergenic(pos, strand):
    iv = strand_intervals.get(strand)
    if iv is None:
        return True
    return not bool(((iv[:, 0] <= pos) & (pos < iv[:, 1])).any())

op = pd.read_csv(OPERON_TSV, sep='\t')
utr5, utr3 = [], []
n_canonical = 0
for _, o in op.iterrows():
    loci = [x for x in str(o['sense_gene_loci']).split(',') if x and x in gene_by_locus.index]
    if not loci:
        continue
    if not (is_intergenic(int(o['tss']), o['strand']) and is_intergenic(int(o['tts']), o['strand'])):
        continue
    n_canonical += 1
    sub = gene_by_locus.loc[loci]
    first = sub.iloc[int(np.argmin(sub['start0'].values))]   # lowest-coord gene
    last  = sub.iloc[int(np.argmax(sub['end0'].values))]     # highest-coord gene
    s0, e0 = int(o['start0']), int(o['end0'])
    if o['strand'] == '+':
        u5, u3 = int(first['start0']) - s0, e0 - int(last['end0'])
    else:
        u5, u3 = e0 - int(last['end0']), int(first['start0']) - s0
    utr5.append(max(0, u5)); utr3.append(max(0, u3))

utr5, utr3 = np.array(utr5), np.array(utr3)
fig, ax = plt.subplots(figsize=(QUART, QUART), constrained_layout=True)
bins = np.logspace(0, np.log10(max(utr5.max(), utr3.max()) + 1), 26)
for lab, arr in [("5'", utr5), ("3'", utr3)]:
    a = arr[arr > 0]
    ax.hist(a, bins=bins, histtype='stepfilled', color=UTR_COLOR[lab], alpha=0.30,
            edgecolor=UTR_COLOR[lab], linewidth=1.0)
    med = float(np.median(arr))
    ax.axvline(med, color=UTR_COLOR[lab], lw=0.9, ls='--')
    ax.plot([], [], color=UTR_COLOR[lab], lw=1.2, label=f"{lab} UTR\n(median {med:.0f} nt)")
    say(f"f) {lab}UTR: n={len(arr)} median={med:.0f} nt max={arr.max():.0f}")
say(f"f) canonical operons used: {n_canonical}")
ax.set_xscale('log')
ax.set_xlabel('UTR length (nt)')
ax.set_ylabel('Operons')
ax.legend(frameon=False, handlelength=1.1, labelspacing=0.25, loc='upper right', fontsize=5)
ax.spines[['top', 'right']].set_visible(False)
fig.savefig(f'{OUT}/panel_f_utr_distribution.pdf', dpi=300); plt.close(fig)


# ====================================================================== g
iso = pd.read_csv(ISO_TSV, sep='\t')
iso = iso[iso['n_reads'] > MIN_READS]
fi = iso['frac_intergenic'].dropna().values * 100.0          # percent
med_fi = float(np.median(fi))
fig, ax = plt.subplots(figsize=(QUART, QUART), constrained_layout=True)
ax.hist(fi, bins=np.linspace(0, 100, 51), color='#1f77b4', alpha=0.85,
        edgecolor='white', linewidth=0.3)
ax.set_yscale('log')
ax.axvline(med_fi, color='black', lw=1.0, ls='--', label=f'median {med_fi:.1f}%')
ax.set_xlim(-1, 101)
ax.set_xlabel('Intergenic coverage (%)')
ax.set_ylabel('Isoforms')
ax.legend(frameon=False, handlelength=1.1, loc='upper right', fontsize=5)
ax.spines[['top', 'right']].set_visible(False)
fig.savefig(f'{OUT}/panel_e_intergenic_coverage.pdf', dpi=300); plt.close(fig)
say(f"e) all isoforms (reads>{MIN_READS}) n={len(fi)} median_intergenic={med_fi:.1f}%")

with open(f'{OUT}/R4_dist_panels.txt', 'w') as fh:
    fh.write("R4 DISTRIBUTION PANELS (b, c, e, f)\n" + "=" * 50 + "\n")
    fh.write("Sizes 7/4 x 7/4 in. Default OUTPUT.md fonts.\n\n")
    fh.write("\n".join(log) + "\n")
print(f"\nSaved panels b/c/f/g + R4_dist_panels.txt to {OUT}/")
