#!/usr/bin/env python
# coding: utf-8
"""
R3 syn3A panels — transcriptome-proteome correlation repeated on JCVI-syn3A.

Companion to Syn1_Corr_RNA_Proteins/R3_figure_panels.py. Builds the syn3A side of
the SI figure (si-correlation.pdf, panels f-j), born-at-size, OUTPUT.md default
fonts (7 pt labels / 6 pt ticks). Same analysis as syn1, with three organism-
specific differences:

  * TIR  — regenerated FROM SCRATCH at the gene level. syn3A has no PacBio isoform
           set, so there is no read-weighting over covering isoforms; OSTIR runs
           on a fixed 30-nt window around each annotated start codon.
  * CAI  — recomputed on syn3A (genetic code 4, TGA=Trp) with syn3A's OWN top-20%-
           by-iPM reference set (the reference gene set differs from syn1).
  * t1/2 — the intrinsic Mpn half-lives are REUSED by locus-tag suffix (syn3A
           proteins are a subset of syn1's; MMSYN1_NNNN <-> JCVISYN3A_NNNN), then
           re-scaled by syn3A Lon/FtsH abundance and cell volume (no new homology
           search).

Panel -> size (in) -> content:
  f  7/2 x 7/4  copy-number distribution by localization (4 merged classes)
  g  7/2 x 7/2  Illumina sense TPM vs iPM (log10), colored by localization
  h  7/4 x 7/4  proteome residual vs log10(TIR)   [gene-level OSTIR]
  i  7/4 x 7/4  proteome residual vs CAI           [syn3A reference set]
  j  7/4 x 7/4  intrinsic half-life distribution   [Mpn t1/2 re-scaled to syn3A]
  k  7/4 x 7/4  proteome residual vs log10(half-life)
  l  7/4 x 7/4  model Pearson R (baseline vs +CAI) for all + cytoplasmic

Run from Syn3A_Corr_RNA_Proteins/ in the RNAseq conda env (needs ostir):
  conda run -n RNAseq python Corr_RNA_Protein_Syn3A.py
Output: R3_panels_syn3A/panel_{f..l}.pdf + syn3A_genes_transcriptomics_proteomics.csv
        + R3_panels_syn3A/R3_syn3A.txt (+ gene_TIR_syn3A.csv cache)
"""

import os
import warnings
from collections import defaultdict
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

HALF, QUART = 7 / 2, 7 / 4          # 3.5 in, 1.75 in
OUT = 'R3_panels_syn3A'
os.makedirs(OUT, exist_ok=True)

PROT_CSV = '../Syn1_Syn3A_Proteomics/syn3a_proteomics_summary_2026.csv'
TPM_TSV  = '../Syn3A_Transcriptomics/Gene_TPM/syn3A_TPM_Illumina_ONT.tsv'
GENOME_FA = '../Genomes_Input/syn3A_genome.fasta'
DEG_CSV  = '../Syn1_Corr_RNA_Proteins/residual_analysis/syn1_ptn_degradation_from_mpn.csv'
TIR_CACHE = f'{OUT}/gene_TIR_syn3A.csv'

GENOME_LEN = 543_379   # JCVI-syn3A (CP016816.2), circular

# --- localization: collapse syn3A's finer scheme to syn1's four classes --------
LOC_MAP = {
    'cytoplasm': 'cytoplasmic',
    'peripheral membrane': 'cytoplasmic',
    'unidentified': 'cytoplasmic',
    'trans-membrane': 'membrane',
    'lipoprotein': 'lipoprotein',
    'extracellular': 'extracellular',
}
LOC_ORDER = ['cytoplasmic', 'lipoprotein', 'membrane', 'extracellular']
LOC_LABEL = {'cytoplasmic': 'Cytoplasmic', 'lipoprotein': 'Lipoprotein',
             'membrane': 'Membrane', 'extracellular': 'Extracellular'}
LOC_COLORS = {'cytoplasmic': '#0072B2', 'lipoprotein': '#009E73',
              'membrane': '#D55E00', 'extracellular': '#CC79A7'}

# --- half-life re-scaling constants (mirror Translation_Residual_L3_degradation.py) ---
Lon_Mpn_Num, FtsH_Mpn_Num, Mpn_Volume = 122, 689, 0.05      # Maier 2011 (Mpn)
Lon_Syn3A_Num  = 518.44   # JCVISYN3A_0394, copy_number_2026
FtsH_Syn3A_Num = 260.20   # JCVISYN3A_0039, copy_number_2026
Syn3A_Volume   = (4.0 / 3.0) * np.pi * (0.20 ** 3)          # fL, r = 200 nm (Breuer 2019)
Cyto_HL_factor = (Lon_Mpn_Num / Mpn_Volume) / (Lon_Syn3A_Num / Syn3A_Volume)
Mem_HL_factor  = (FtsH_Mpn_Num / Mpn_Volume) / (FtsH_Syn3A_Num / Syn3A_Volume)

log = []
def say(s):
    print(s); log.append(s)

def fit_line(ax, x, y, color='crimson'):
    b, a = np.polyfit(x, y, 1)
    xs = np.array([x.min(), x.max()])
    ax.plot(xs, b * xs + a, color=color, lw=1.0, zorder=5)

def r2_lstsq(X, y):
    A = np.column_stack([np.ones(len(y))] + [X[:, i] for i in range(X.shape[1])])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    return 1 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)

# =============================================================== load & merge
prot = pd.read_csv(PROT_CSV)
prot['loc4'] = prot['localization'].map(LOC_MAP)
tpm = pd.read_csv(TPM_TSV, sep='\t')

df = prot.merge(
    tpm[['locus_tag', 'start0', 'end0', 'strand', 'gene_len', 'Illumina_sense_TPM']],
    on='locus_tag', how='left')
say(f"syn3A proteins: {len(prot)}; merged with TPM coords: {df['start0'].notna().sum()}")
say(f"localization merge -> " + ", ".join(
    f"{k}:{int((df['loc4'] == k).sum())}" for k in LOC_ORDER))

# =============================================================== f  copy-number
cn = df[(df['copy_number_2026'] > 0) & df['loc4'].isin(LOC_ORDER)]
fig, ax = plt.subplots(figsize=(HALF, QUART), constrained_layout=True)
allv = np.log10(cn['copy_number_2026'].values)
lo, hi = allv.min(), allv.max()
xgrid = np.linspace(lo, hi, 300)
bw = (hi - lo) / 25.0
for loc in LOC_ORDER:
    v = np.log10(cn.loc[cn['loc4'] == loc, 'copy_number_2026'].values)
    if len(v) < 5:
        say(f"f) {loc}: n={len(v)} (<5, curve skipped)")
        continue
    n = len(v); med = float(np.median(v)); med_cn = 10 ** med
    kde = gaussian_kde(v)
    ycurve = kde(xgrid) * n * bw
    ax.plot(xgrid, ycurve, color=LOC_COLORS[loc], lw=1.2,
            label=f"{LOC_LABEL[loc]} ({med_cn:.0f}, n={n})")
    ax.fill_between(xgrid, ycurve, color=LOC_COLORS[loc], alpha=0.13)
    ax.vlines(med, 0, float(kde(med)[0]) * n * bw, color=LOC_COLORS[loc], lw=0.9, ls='--')
    say(f"f) {loc}: n={n} median={med_cn:.1f} copies")
ax.set_xlabel('Copies per cell ($\\log_{10}$)')
ax.set_ylabel('Proteins')
ax.set_xlim(lo, hi)
leg = ax.legend(frameon=False, handlelength=0.8, labelspacing=0.2, borderpad=0.2,
                fontsize=5, loc='upper left',
                title='(median copies, n=unique proteins)', title_fontsize=5)
leg._legend_box.align = 'left'
ax.spines[['top', 'right']].set_visible(False)
fig.savefig(f'{OUT}/panel_f_copynumber_by_localization.pdf', dpi=300); plt.close(fig)

# =============================================================== g  TPM vs iPM
gd = df[(df['iPM_mean'] > 0) & (df['Illumina_sense_TPM'] > 0) & df['loc4'].isin(LOC_ORDER)].copy()
gd['log10_TPM'] = np.log10(gd['Illumina_sense_TPM'])
gd['log10_iPM'] = np.log10(gd['iPM_mean'])
fig, ax = plt.subplots(figsize=(HALF, HALF), constrained_layout=True)
for loc in LOC_ORDER:
    s = gd[gd['loc4'] == loc]
    ax.scatter(s['log10_TPM'], s['log10_iPM'], s=8, alpha=0.7,
               c=LOC_COLORS[loc], edgecolors='none', label=f"{LOC_LABEL[loc]} (n={len(s)})")
r_all = pearsonr(gd['log10_TPM'], gd['log10_iPM'])[0]
cyto = gd[gd['loc4'] == 'cytoplasmic']
r_cyto = pearsonr(cyto['log10_TPM'], cyto['log10_iPM'])[0]
fit_line(ax, gd['log10_TPM'].values, gd['log10_iPM'].values, color='black')
ax.text(0.04, 0.96, f"all $r$ = {r_all:.2f} ($n$ = {len(gd)})\n"
                    f"cytoplasmic $r$ = {r_cyto:.2f} ($n$ = {len(cyto)})",
        transform=ax.transAxes, va='top', ha='left', fontsize=6)
ax.set_xlabel('Illumina sense TPM ($\\log_{10}$)')
ax.set_ylabel('Protein iPM ($\\log_{10}$)')
ax.legend(frameon=False, handlelength=1.0, labelspacing=0.25, loc='lower right')
ax.spines[['top', 'right']].set_visible(False)
fig.savefig(f'{OUT}/panel_g_TPM_vs_iPM.pdf', dpi=300); plt.close(fig)
say(f"g) all r={r_all:.3f} (n={len(gd)}); cytoplasmic r={r_cyto:.3f} (n={len(cyto)})")

# =============================================================== genome + CDS
STOP_CODONS = {"TAA", "TAG"}
GENETIC_CODE_4 = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"W","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def load_single_fasta(path):
    chunks = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                chunks.append(line.strip().upper())
    return "".join(chunks)

genome = load_single_fasta(GENOME_FA)
assert len(genome) == GENOME_LEN, f"genome length {len(genome)} != {GENOME_LEN}"
_COMP = str.maketrans("ACGTN", "TGCAN")
def revcomp(s):
    return s.translate(_COMP)[::-1]

def _fetch_circular(g_start, g_end):
    g_start %= GENOME_LEN
    g_end = g_end % GENOME_LEN if g_end % GENOME_LEN != 0 else GENOME_LEN
    return genome[g_start:g_end] if g_start < g_end else genome[g_start:] + genome[:g_end]

def extract_cds(start0, end0, strand):
    seq = _fetch_circular(int(start0), int(end0))
    return revcomp(seq) if strand == "-" else seq

# protein-coding genes with valid in-frame CDS (coords from the TPM table)
cds_seqs = {}
n_bad = 0
for _, r in df.iterrows():
    if pd.isna(r['start0']):
        continue
    seq = extract_cds(r['start0'], r['end0'], r['strand'])
    if len(seq) < 6 or len(seq) % 3 != 0:
        n_bad += 1
        continue
    cds_seqs[r['locus_tag']] = seq
say(f"\nCDSs extracted: {len(cds_seqs)} (skipped non-multiple-of-3: {n_bad})")

# =============================================================== i  CAI (syn3A)
ipm_lookup = df.set_index('locus_tag')['iPM_mean'].to_dict()
prot_ranked = (df[df['locus_tag'].isin(cds_seqs) & (df['iPM_mean'] > 0)]
               .sort_values('iPM_mean', ascending=False))
n_ref = max(1, int(round(0.20 * len(prot_ranked))))
ref_loci = set(prot_ranked['locus_tag'].iloc[:n_ref])
say(f"CAI reference set (top 20% by iPM): {n_ref}/{len(prot_ranked)} genes")

SYNONYMOUS = defaultdict(list)
for codon, aa in GENETIC_CODE_4.items():
    if aa != "*":
        SYNONYMOUS[aa].append(codon)

ref_codon_counts = defaultdict(int)
for lt in ref_loci:
    codons = [cds_seqs[lt][i:i+3] for i in range(0, len(cds_seqs[lt]), 3)]
    if codons and codons[-1] in STOP_CODONS:
        codons = codons[:-1]
    for c in codons[1:]:
        if "N" in c or c in STOP_CODONS or GENETIC_CODE_4.get(c) is None:
            continue
        ref_codon_counts[c] += 1

w = {}
for aa, codons in SYNONYMOUS.items():
    counts = np.array([ref_codon_counts[c] for c in codons], dtype=float)
    if len(codons) == 1 or aa in ("M", "W") or counts.max() == 0:
        for c in codons:
            w[c] = 1.0
        continue
    rel = counts / counts.max()
    rel = np.where(rel == 0, 0.01, rel)
    for c, rv in zip(codons, rel):
        w[c] = float(rv)

def compute_cai(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    if codons and codons[-1] in STOP_CODONS:
        codons = codons[:-1]
    logs = []
    for c in codons[1:]:
        aa = GENETIC_CODE_4.get(c)
        if "N" in c or c in STOP_CODONS or aa is None or aa in ("M", "W"):
            continue
        logs.append(np.log(w[c]))
    return float(np.exp(np.mean(logs))) if logs else np.nan

cai_map = {lt: compute_cai(seq) for lt, seq in cds_seqs.items()}
df['CAI'] = df['locus_tag'].map(cai_map)

# =============================================================== h  TIR (OSTIR, gene-level, from scratch)
ANTI_SD, UTR_WINDOW, CDS_WINDOW = "ACCUCCUUU", 30, 30

def extract_initiation_seq(start0, end0, strand):
    """Fixed 30-nt window around the annotated start codon (circular)."""
    if strand == "+":
        seq = _fetch_circular(start0 - UTR_WINDOW, start0 + CDS_WINDOW)
    else:
        seq = revcomp(_fetch_circular(end0 - CDS_WINDOW, end0 + UTR_WINDOW))
    return seq, UTR_WINDOW   # start codon at 0-based index UTR_WINDOW

if os.path.exists(TIR_CACHE):
    tir_df = pd.read_csv(TIR_CACHE)
    say(f"\nTIR: loaded cache {TIR_CACHE} ({len(tir_df)} genes)")
else:
    from ostir import run_ostir
    say("\nTIR: running gene-level OSTIR on syn3A start windows ...")
    rows = []
    coord = df.dropna(subset=['start0']).set_index('locus_tag')
    for lt in cds_seqs:
        r = coord.loc[lt]
        seq, utr_len = extract_initiation_seq(int(r['start0']), int(r['end0']), r['strand'])
        start_1b = utr_len + 1
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                out = run_ostir(seq, aSD=ANTI_SD, threads=8)
        except Exception as e:
            say(f"  OSTIR error {lt}: {e}")
            continue
        if not out:
            continue
        best = next((h for h in out if h.get("start_position") == start_1b), None)
        if best is None:
            best = min(out, key=lambda x: abs(x.get("start_position", 0) - start_1b))
        rows.append({"locus_tag": lt, "TIR": float(best.get("expression", np.nan)),
                     "start_position": best.get("start_position", np.nan)})
    tir_df = pd.DataFrame(rows)
    tir_df.to_csv(TIR_CACHE, index=False)
    say(f"  TIR computed for {tir_df['TIR'].notna().sum()}/{len(cds_seqs)} genes -> {TIR_CACHE}")

df = df.merge(tir_df[['locus_tag', 'TIR']], on='locus_tag', how='left')

# =============================================================== residual helper
def residual_set(col, positive=True):
    """genes with iPM/TPM>0 and a valid `col`; return copy with residual of log10(iPM)~log10(TPM)."""
    s = df[(df['iPM_mean'] > 0) & (df['Illumina_sense_TPM'] > 0) & df[col].notna()].copy()
    if positive:
        s = s[s[col] > 0]
    s['log10_TPM'] = np.log10(s['Illumina_sense_TPM'])
    s['log10_iPM'] = np.log10(s['iPM_mean'])
    sl, ic = np.polyfit(s['log10_TPM'], s['log10_iPM'], 1)
    s['residual'] = s['log10_iPM'] - (sl * s['log10_TPM'] + ic)
    return s

# --- i  CAI vs residual
ci = residual_set('CAI', positive=False)
r_cai = pearsonr(ci['CAI'], ci['residual'])[0]
base_cai = r2_lstsq(ci[['log10_TPM']].values, ci['log10_iPM'].values)
full_cai = r2_lstsq(ci[['log10_TPM', 'CAI']].values, ci['log10_iPM'].values)
fig, ax = plt.subplots(figsize=(QUART, QUART), constrained_layout=True)
ax.scatter(ci['CAI'], ci['residual'], s=5, alpha=0.5, c='#555555', edgecolors='none')
fit_line(ax, ci['CAI'].values, ci['residual'].values)
ax.axhline(0, color='black', lw=0.5, ls=':')
ax.text(0.04, 0.96, f"$r$ = {r_cai:.2f}\n$n$ = {len(ci)}", transform=ax.transAxes, va='top', fontsize=6)
ax.set_xlabel('CAI')
ax.set_ylabel('Proteome residual')
ax.spines[['top', 'right']].set_visible(False)
fig.savefig(f'{OUT}/panel_i_CAI_vs_residual.pdf', dpi=300); plt.close(fig)
say(f"i) CAI vs residual r={r_cai:.3f} (n={len(ci)}); R2 {base_cai:.3f}->{full_cai:.3f} "
    f"(dR2={full_cai-base_cai:+.3f})")

# --- h  TIR vs residual
ti = residual_set('TIR', positive=True)
ti['log10_TIR'] = np.log10(ti['TIR'])
r_tir = pearsonr(ti['log10_TIR'], ti['residual'])[0]
base_tir = r2_lstsq(ti[['log10_TPM']].values, ti['log10_iPM'].values)
full_tir = r2_lstsq(ti[['log10_TPM', 'log10_TIR']].values, ti['log10_iPM'].values)
fig, ax = plt.subplots(figsize=(QUART, QUART), constrained_layout=True)
ax.scatter(ti['log10_TIR'], ti['residual'], s=5, alpha=0.5, c='#555555', edgecolors='none')
fit_line(ax, ti['log10_TIR'].values, ti['residual'].values)
ax.axhline(0, color='black', lw=0.5, ls=':')
ax.text(0.04, 0.96, f"$r$ = {r_tir:.2f}\n$n$ = {len(ti)}", transform=ax.transAxes, va='top', fontsize=6)
ax.set_xlabel('TIR ($\\log_{10}$)')
ax.set_ylabel('Proteome residual')
ax.spines[['top', 'right']].set_visible(False)
fig.savefig(f'{OUT}/panel_h_TIR_vs_residual.pdf', dpi=300); plt.close(fig)
say(f"h) TIR vs residual r={r_tir:.3f} (n={len(ti)}); R2 {base_tir:.3f}->{full_tir:.3f} "
    f"(dR2={full_tir-base_tir:+.3f})")

# --- R with/without CAI, all + cytoplasmic (parallel to syn1 panel d; reported only)
def R_pair(d):
    y = d['log10_iPM'].values
    return np.sqrt(r2_lstsq(d[['log10_TPM']].values, y)), \
           np.sqrt(r2_lstsq(d[['log10_TPM', 'CAI']].values, y))
cyt = ci[ci['loc4'] == 'cytoplasmic']
Rb_all, Rc_all = R_pair(ci)
Rb_cyt, Rc_cyt = R_pair(cyt)
say(f"   R all {Rb_all:.3f}->{Rc_all:.3f} (n={len(ci)}); cytoplasmic {Rb_cyt:.3f}->{Rc_cyt:.3f} (n={len(cyt)})")

# --- l  Pearson R with/without CAI, all + cytoplasmic (mirrors syn1 panel d/f)
fig, ax = plt.subplots(figsize=(QUART, QUART), constrained_layout=True)
xpos = np.array([0, 1]); wbar = 0.38
b1 = ax.bar(xpos - wbar/2, [Rb_all, Rb_cyt], wbar, color='#bbbbbb', label='TPM only')
b2 = ax.bar(xpos + wbar/2, [Rc_all, Rc_cyt], wbar, color='#0072B2', label='+ CAI')
for bars in (b1, b2):
    for bar in bars:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.012,
                f"{bar.get_height():.2f}", ha='center', va='bottom', fontsize=5)
ax.set_xticks(xpos); ax.set_xticklabels([f'All\n(n={len(ci)})', f'Cytosolic\n(n={len(cyt)})'])
ax.set_ylabel('$R$'); ax.set_ylim(0, 0.95)
ax.legend(frameon=False, handlelength=1.0, labelspacing=0.25, loc='upper left')
ax.spines[['top', 'right']].set_visible(False)
fig.savefig(f'{OUT}/panel_l_R_improvement.pdf', dpi=300); plt.close(fig)
say(f"l) R barplot: all {Rb_all:.3f}->{Rc_all:.3f}, cytosolic {Rb_cyt:.3f}->{Rc_cyt:.3f}")

# =============================================================== j  half-life (reuse + re-scale)
deg = pd.read_csv(DEG_CSV)
deg['suffix'] = deg['syn1_locus_tag'].str.extract(r'_(\d+)$')
deg = deg.dropna(subset=['suffix']).drop_duplicates('suffix')
hl = df[['locus_tag', 'loc4']].copy()
hl['suffix'] = hl['locus_tag'].str.extract(r'_(\d+)$')
hl = hl.merge(deg[['suffix', 'halflife_h']], on='suffix', how='inner').dropna(subset=['halflife_h'])

def _factor(loc):
    if loc == 'cytoplasmic':
        return Cyto_HL_factor
    if loc in ('membrane', 'lipoprotein'):
        return Mem_HL_factor
    return np.nan
hl['factor'] = hl['loc4'].map(_factor)
hl['halflife_h_syn3A'] = hl['halflife_h'] * hl['factor']
hl = hl.dropna(subset=['halflife_h_syn3A'])
df = df.merge(hl[['locus_tag', 'halflife_h_syn3A']], on='locus_tag', how='left')
say(f"\nHalf-life: {len(hl)} genes mapped by suffix; Cyto_factor={Cyto_HL_factor:.3f} "
    f"Mem_factor={Mem_HL_factor:.3f} (Syn3A_V={Syn3A_Volume:.4f} fL)")

hv = hl['halflife_h_syn3A']
DOUBLING_H = 105.0 / 60.0            # syn3A cell cycle, 105 min = 1.75 h
j_med, j_min = float(hv.median()), float(hv.min())
pct_below = float((hv < DOUBLING_H).mean() * 100)   # proteins turned over faster than dilution
fig, ax = plt.subplots(figsize=(QUART, QUART), constrained_layout=True)
ax.hist(np.log10(hv), bins=20, color='#9467bd', edgecolor='white', linewidth=0.3)
ax.axvline(np.log10(j_med), color='crimson', lw=1.0, ls='--', label=f'median {j_med:.1f} h')
ax.axvline(np.log10(DOUBLING_H), color='black', lw=1.0, ls='-',
           label=f'doubling 105 min\n(shortest {j_min:.1f} h,\n{pct_below:.0f}% below)')
ax.set_xlabel('Half-life, h ($\\log_{10}$)')
ax.set_ylabel('Proteins')
ax.legend(frameon=False, handlelength=1.2, labelspacing=0.25, loc='upper right',
          title=f'n = {len(hv)}', title_fontsize=6)
ax.spines[['top', 'right']].set_visible(False)
fig.savefig(f'{OUT}/panel_j_halflife_distribution.pdf', dpi=300); plt.close(fig)
say(f"j) half-life median={j_med:.1f} h, shortest={j_min:.1f} h, n={len(hv)}; "
    f"{pct_below:.1f}% below the 105-min doubling time")

# --- k  half-life vs residual (mirrors syn1 panel h)
hk = residual_set('halflife_h_syn3A', positive=True)
hk['log10_hl'] = np.log10(hk['halflife_h_syn3A'])
r_k = pearsonr(hk['log10_hl'], hk['residual'])[0]
fig, ax = plt.subplots(figsize=(QUART, QUART), constrained_layout=True)
ax.scatter(hk['log10_hl'], hk['residual'], s=5, alpha=0.5, c='#555555', edgecolors='none')
fit_line(ax, hk['log10_hl'].values, hk['residual'].values)
ax.axhline(0, color='black', lw=0.5, ls=':')
ax.text(0.04, 0.96, f"$r$ = {r_k:.2f}\n$n$ = {len(hk)}", transform=ax.transAxes, va='top', fontsize=6)
ax.set_xlabel('Half-life, h ($\\log_{10}$)')
ax.set_ylabel('Proteome residual')
ax.spines[['top', 'right']].set_visible(False)
fig.savefig(f'{OUT}/panel_k_halflife_vs_residual.pdf', dpi=300); plt.close(fig)
say(f"k) half-life vs residual r={r_k:.3f} (n={len(hk)})")

# =============================================================== combined table
out = df.copy()
cols = ['locus_tag', 'gene_name', 'gene_product', 'loc4', 'Illumina_sense_TPM',
        'iPM_mean', 'copy_number_2026', 'CAI', 'TIR', 'halflife_h_syn3A']
out = out[cols].rename(columns={'loc4': 'protein_localization',
                                'Illumina_sense_TPM': 'TPM_illumina',
                                'copy_number_2026': 'protein_copy_number',
                                'halflife_h_syn3A': 'protein_halflife_h'})
for c in ['TPM_illumina', 'iPM_mean', 'CAI', 'TIR', 'protein_halflife_h']:
    out[c] = out[c].round(3)
out.to_csv('syn3A_genes_transcriptomics_proteomics.csv', index=False)
say(f"\nsyn3A_genes_transcriptomics_proteomics.csv: {len(out)} genes "
    f"(TPM {out['TPM_illumina'].notna().sum()}, iPM {out['iPM_mean'].notna().sum()}, "
    f"CAI {out['CAI'].notna().sum()}, TIR {out['TIR'].notna().sum()}, "
    f"half-life {out['protein_halflife_h'].notna().sum()})")

with open(f'{OUT}/R3_syn3A.txt', 'w') as fh:
    fh.write("R3 syn3A PANELS (transcriptome-proteome correlation repeated on JCVI-syn3A)\n")
    fh.write("=" * 68 + "\n")
    fh.write("Sizes (in): f 7/2 x 7/4; g 7/2 x 7/2; h,i,j,k,l 7/4 x 7/4. Default fonts.\n\n")
    fh.write("\n".join(log) + "\n")
print(f"\nSaved 7 panels + R3_syn3A.txt to {OUT}/")
