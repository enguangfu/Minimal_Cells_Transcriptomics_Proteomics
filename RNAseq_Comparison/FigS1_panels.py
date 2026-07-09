#!/usr/bin/env python3
"""Fig S1 panels: cross-platform / cross-run RNA-seq quantification comparison.

Reads the pre-computed tables from compute_platform_TPM.py (never recomputes).
Organism colours (match Fig 5/6 + the R3 panels): Syn1 = blue (#3182bd),
Syn3A = red (#c0392b). The syn1 panels (a-c) are blue, the syn3A panel (d) red;
b/c/d also carry a bold organism tag. Identity/reference lines are neutral grey.

Panels (born at final print size, OUTPUT.md fonts, Arial, pdf.fonttype 42):
  a (7/2 x 7/2) : Syn1 4x4 scatter-plot matrix of mRNA TPM (order along both
                  axes: Illumina, PacBio, ONT run 1, ONT run 2); lower triangle =
                  log10-log10 blue scatter with identity line, diagonal = log10 TPM
                  histogram (unlabelled), upper triangle = blank. Every subpanel
                  carries [0,2,4] ticks on both axes.
  b (7/6 x 7/2) : Syn1 mapped-read length distribution, 1x4 stacked (Illumina,
                  PacBio, ONT run 1, ONT run 2) on a shared x-axis; blue histograms,
                  median marked. From readlen_syn1.tsv (compute_readlen_syn1.py).
  c (7/4 x 7/4) : Syn1 length bias, log2(PacBio / Illumina) vs gene length,
                  blue scatter + linear fit (reproduces the old Fig S1 panel).
  d (7/4 x 7/4) : Syn1 abundance bias, log2(PacBio / Illumina) vs MA-plot mean
                  abundance, blue scatter + linear fit.
  e (7/4 x 7/4) : Syn3A ONT vs Illumina TPM (log10-log10), red scatter, Pearson r.

Outputs the four individual panels to FigS1_panels/ (born-at-size). The combined
Fig S1 is assembled manually in Illustrator from these panels.
Run in the RNAseq conda env.
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, linregress

mpl.rcParams.update({
    'font.size': 7, 'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
    'axes.titlesize': 7, 'axes.labelsize': 7,
    'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 6,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
})

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "FigS1_panels")
FIGDIR = os.path.abspath(os.path.join(HERE, "..", "Manuscript", "figures"))
os.makedirs(OUT, exist_ok=True)
# Arial (this build) lacks Unicode subscript glyphs and mathtext embeds Computer
# Modern (breaks the pure-Arial rule), so use plain "log10" / "log2".
SUB10, SUB2 = "10", "2"
C_PB, C_ONT, C_GREY = "#3b6fb0", "#e08214", "#7a7a7a"
# Organism identity colors (match Fig 5/6 and the R3 panels): Syn1 blue, Syn3A red.
SYN1_COL, SYN3A_COL = "#3182bd", "#c0392b"
THR = 0.5
LIM = (-1.2, 5.2)


def org_tag(ax, which="syn1"):
    name, col = ("Syn1.0", SYN1_COL) if which == "syn1" else ("Syn3A", SYN3A_COL)
    ax.set_title(name, loc="left", color=col, fontweight="bold", fontsize=6.5)

s1 = pd.read_csv(os.path.join(HERE, "platform_TPM_syn1.tsv"), sep="\t")
s3 = pd.read_csv(os.path.join(HERE, "platform_TPM_syn3A.tsv"), sep="\t")
rl = pd.read_csv(os.path.join(HERE, "readlen_syn1.tsv"), sep="\t")  # compute_readlen_syn1.py


def l10(v):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log10(np.asarray(v, float))


def binned_median(x, y, nbins=8):
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    idx = np.array_split(np.arange(len(xs)), nbins)
    cx = [np.median(xs[k]) for k in idx if len(k)]
    cy = [np.median(ys[k]) for k in idx if len(k)]
    return np.array(cx), np.array(cy)


PLAT = [("Illumina", "Illumina_TPM"), ("PacBio", "PacBio_TPM"),
        ("ONT1", "ONT1_TPM"), ("ONT2", "ONT2_TPM")]


def draw_matrix(axes):
    """Lower triangle = log10-log10 scatter with identity line; diagonal = log10
    TPM histogram (no platform label); upper triangle = blank. Every subpanel
    carries [0,2,4] ticks on both axes (diagonal y is a count, so y-ticks off)."""
    n = len(PLAT)
    for i in range(n):
        for j in range(n):
            ax = axes[i][j]
            ax.tick_params(length=1.5, pad=1)
            if i == j:                                   # diagonal: histogram, no label
                v = l10(s1[PLAT[i][1]].values)
                v = v[np.isfinite(v) & (v > LIM[0])]
                ax.hist(v, bins=22, range=LIM, color=SYN1_COL, edgecolor="none")
                ax.set_xlim(LIM)
                ax.set_xticks([0, 2, 4]); ax.set_xticklabels(["0", "2", "4"], fontsize=4.5)
                ax.set_yticks([])                        # y is a count, not TPM
            elif i > j:                                  # lower: scatter, ticks on both axes
                a, b = s1[PLAT[j][1]].values, s1[PLAT[i][1]].values
                m = (a > THR) & (b > THR)
                ax.scatter(l10(a[m]), l10(b[m]), s=1.2, alpha=0.35, c=SYN1_COL, edgecolors="none")
                ax.plot(LIM, LIM, ls=":", lw=0.5, c="#888888")
                ax.set_xlim(LIM); ax.set_ylim(LIM)
                ax.set_xticks([0, 2, 4]); ax.set_xticklabels(["0", "2", "4"], fontsize=4.5)
                ax.set_yticks([0, 2, 4]); ax.set_yticklabels(["0", "2", "4"], fontsize=4.5)
            else:                                        # upper triangle: blank, no ticks
                ax.set_xlim(LIM); ax.set_ylim(LIM)
                ax.set_xticks([]); ax.set_yticks([])


def draw_bias(ax, mode):
    """PacBio-vs-Illumina bias only (repeats the old Fig S1 style): gray scatter,
    crimson linear fit, r / P annotation. mode='length' -> x = gene length;
    mode='abundance' -> x = MA-plot mean 0.5*(log2 PacBio + log2 Illumina)."""
    m = ((s1.Illumina_TPM > THR) & (s1.PacBio_TPM > THR)).values
    I, P = s1.Illumina_TPM.values[m], s1.PacBio_TPM.values[m]
    y = np.log2(P / I)
    if mode == "length":
        x = s1.gene_len.values[m].astype(float)
        xlabel = "Gene length (bp)"
    else:
        x = 0.5 * (np.log2(P) + np.log2(I))
        xlabel = f"Mean abundance (log{SUB2})"
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    sl, ic, r, p, _ = linregress(x, y)
    ax.scatter(x, y, s=5, alpha=0.5, c=SYN1_COL, edgecolors="none")
    ax.axhline(0, color="black", lw=0.5, ls=":")
    xs = np.array([x.min(), x.max()])
    ax.plot(xs, sl * xs + ic, color="black", lw=1.0)
    ax.text(0.04, 0.96, f"r = {r:.2f}\nP = {p:.1g}", transform=ax.transAxes, va="top", fontsize=5.5)
    ax.set_xlabel(xlabel, fontsize=6.5)
    ax.set_ylabel(f"PacBio/Illumina TPM (log{SUB2})", fontsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2)


# unit per library: cDNA (Illumina, PacBio) measured in bp; direct-RNA ONT in nt
RL_ORDER = [("Illumina", "Illumina", "bp"), ("PacBio", "PacBio", "bp"),
            ("ONT1", "ONT run 1", "nt"), ("ONT2", "ONT run 2", "nt")]
RL_XMAX, RL_BIN = 4000, 50   # shared x-axis; PacBio's ~1% tail beyond 4 kb is off-screen


def draw_readlen(axes):
    """Mapped-read length distribution of the four Syn1 platforms (top to bottom:
    Illumina, PacBio, ONT run 1, ONT run 2), one histogram per row on a SHARED
    x-axis so the read-length regimes line up. y = fraction of that platform's reads
    per 50 bp bin (own scale); median marked with a dashed line."""
    edges = np.arange(0, RL_XMAX + RL_BIN, RL_BIN)
    ctr = 0.5 * (edges[:-1] + edges[1:])
    for ax, (key, label, unit) in zip(axes, RL_ORDER):
        sub = rl[rl.platform == key].sort_values("length")
        L, C = sub.length.values.astype(float), sub["count"].values.astype(float)
        tot = C.sum()
        h, _ = np.histogram(L, bins=edges, weights=C)
        ax.bar(ctr, h / tot, width=RL_BIN, color=SYN1_COL, edgecolor="none")
        med = L[np.searchsorted(np.cumsum(C) / tot, 0.5)]
        ax.axvline(med, color="black", lw=0.6, ls="--")
        ax.text(0.97, 0.90, f"{label}\nmedian {int(med)} {unit}", transform=ax.transAxes,
                ha="right", va="top", fontsize=5, color=SYN1_COL)
        ax.set_xlim(0, RL_XMAX)
        ax.set_ylim(0, (h / tot).max() * 1.30)
        ax.set_yticks([0, round((h / tot).max(), 2)])
        ax.tick_params(length=1.5, pad=1)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)
    axes[-1].set_xticks([0, 1000, 2000, 3000, 4000])
    axes[-1].set_xticklabels(["0", "1k", "2k", "3k", "4k"])
    axes[-1].set_xlabel("Read length (bp/nt)", fontsize=6.5)


def draw_syn3a(ax):
    a, b = s3.Illumina_TPM.values, s3.ONT_TPM.values
    m = (a > THR) & (b > THR)
    x, y = l10(a[m]), l10(b[m])
    r = pearsonr(x, y)[0]
    ax.scatter(x, y, s=5, alpha=0.5, c=SYN3A_COL, edgecolors="none")  # red = syn3A
    ax.plot(LIM, LIM, ls=":", lw=0.5, c="#888888")
    ax.set_xlim(LIM); ax.set_ylim(LIM)
    ax.set_xlabel(f"Illumina TPM (log{SUB10})", fontsize=6.5)
    ax.set_ylabel(f"ONT TPM (log{SUB10})", fontsize=6.5)
    ax.text(0.05, 0.93, f"r = {r:.2f}\nn = {int(m.sum())}", transform=ax.transAxes,
            ha="left", va="top", fontsize=5.5)
    ax.tick_params(length=2)


# ------------------------- individual panels -------------------------
fig, axes = plt.subplots(4, 4, figsize=(7 / 2, 7 / 2), constrained_layout=True)
draw_matrix(axes)
fig.supxlabel(f"mRNA TPM (log{SUB10})", fontsize=6)
fig.supylabel(f"mRNA TPM (log{SUB10})", fontsize=6)
fig.savefig(os.path.join(OUT, "panel_a_syn1_TPM_matrix.pdf"), dpi=300); plt.close(fig)

fig, axes = plt.subplots(4, 1, figsize=(7 / 6, 7 / 2), constrained_layout=True, sharex=True)
draw_readlen(axes)
fig.supylabel("Fraction of reads", fontsize=6)
fig.savefig(os.path.join(OUT, "panel_b_syn1_readlen.pdf"), dpi=300); plt.close(fig)

for mode, fn in [("length", "panel_c_syn1_length_bias.pdf"),
                 ("abundance", "panel_d_syn1_abundance_bias.pdf")]:
    fig, ax = plt.subplots(figsize=(7 / 4, 7 / 4), constrained_layout=True)
    draw_bias(ax, mode)
    org_tag(ax, "syn1")
    fig.savefig(os.path.join(OUT, fn), dpi=300); plt.close(fig)

fig, ax = plt.subplots(figsize=(7 / 4, 7 / 4), constrained_layout=True)
draw_syn3a(ax)
org_tag(ax, "syn3a")
fig.savefig(os.path.join(OUT, "panel_e_syn3A_ONT_vs_Illumina.pdf"), dpi=300); plt.close(fig)

print("wrote 4 individual panels to", OUT, "(combined Fig S1 assembled manually in Illustrator)")
