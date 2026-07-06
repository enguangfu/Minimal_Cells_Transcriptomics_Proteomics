#!/usr/bin/env python3
"""Fig S1 panels: cross-platform / cross-run RNA-seq quantification comparison.

Reads the pre-computed tables from compute_platform_TPM.py (never recomputes).
Panels (born at final print size, OUTPUT.md fonts, Arial, pdf.fonttype 42):
  a (7/2 x 7/2) : Syn1 4x4 scatter-plot matrix of mRNA TPM (order along both
                  axes: Illumina, PacBio, ONT run 1, ONT run 2); lower triangle =
                  log10-log10 scatter with identity line, diagonal = log10 TPM
                  histogram (unlabelled), upper triangle = blank. Every subpanel
                  carries [0,2,4] ticks on both axes.
  b (7/4 x 7/4) : Syn1 length bias, log2(PacBio / Illumina) vs gene length,
                  gray scatter + linear fit (reproduces the old Fig S1 panel).
  c (7/4 x 7/4) : Syn1 abundance bias, log2(PacBio / Illumina) vs MA-plot mean
                  abundance, gray scatter + linear fit.
  d (7/4 x 7/4) : Syn3A ONT vs Illumina sense TPM (log10-log10), Pearson r.

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
THR = 0.5
LIM = (-1.2, 5.2)

s1 = pd.read_csv(os.path.join(HERE, "platform_TPM_syn1.tsv"), sep="\t")
s3 = pd.read_csv(os.path.join(HERE, "platform_TPM_syn3A.tsv"), sep="\t")


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
                ax.hist(v, bins=22, range=LIM, color=C_GREY, edgecolor="none")
                ax.set_xlim(LIM)
                ax.set_xticks([0, 2, 4]); ax.set_xticklabels(["0", "2", "4"], fontsize=4.5)
                ax.set_yticks([])                        # y is a count, not TPM
            elif i > j:                                  # lower: scatter, ticks on both axes
                a, b = s1[PLAT[j][1]].values, s1[PLAT[i][1]].values
                m = (a > THR) & (b > THR)
                ax.scatter(l10(a[m]), l10(b[m]), s=1.2, alpha=0.35, c=C_GREY, edgecolors="none")
                ax.plot(LIM, LIM, ls=":", lw=0.5, c="#c0392b")
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
    ax.scatter(x, y, s=5, alpha=0.5, c="#555555", edgecolors="none")
    ax.axhline(0, color="black", lw=0.5, ls=":")
    xs = np.array([x.min(), x.max()])
    ax.plot(xs, sl * xs + ic, color="crimson", lw=1.0)
    ax.text(0.04, 0.96, f"r = {r:.2f}\nP = {p:.1g}", transform=ax.transAxes, va="top", fontsize=5.5)
    ax.set_xlabel(xlabel, fontsize=6.5)
    ax.set_ylabel(f"PacBio/Illumina TPM (log{SUB2})", fontsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2)


def draw_syn3a(ax):
    a, b = s3.Illumina_TPM.values, s3.ONT_TPM.values
    m = (a > THR) & (b > THR)
    x, y = l10(a[m]), l10(b[m])
    r = pearsonr(x, y)[0]
    ax.scatter(x, y, s=2, alpha=0.5, c=C_GREY, edgecolors="none")
    ax.plot(LIM, LIM, ls=":", lw=0.5, c="#c0392b")
    ax.set_xlim(LIM); ax.set_ylim(LIM)
    ax.set_xlabel(f"Illumina sense TPM (log{SUB10})", fontsize=6.5)
    ax.set_ylabel(f"ONT sense TPM (log{SUB10})", fontsize=6.5)
    ax.text(0.05, 0.93, f"r = {r:.2f}\nn = {int(m.sum())}", transform=ax.transAxes,
            ha="left", va="top", fontsize=5.5)
    ax.tick_params(length=2)


# ------------------------- individual panels -------------------------
fig, axes = plt.subplots(4, 4, figsize=(7 / 2, 7 / 2), constrained_layout=True)
draw_matrix(axes)
fig.supxlabel(f"mRNA TPM (log{SUB10})", fontsize=6)
fig.supylabel(f"mRNA TPM (log{SUB10})", fontsize=6)
fig.savefig(os.path.join(OUT, "panel_a_syn1_TPM_matrix.pdf"), dpi=300); plt.close(fig)

for mode, fn in [("length", "panel_b_syn1_length_bias.pdf"),
                 ("abundance", "panel_c_syn1_abundance_bias.pdf")]:
    fig, ax = plt.subplots(figsize=(7 / 4, 7 / 4), constrained_layout=True)
    draw_bias(ax, mode)
    fig.savefig(os.path.join(OUT, fn), dpi=300); plt.close(fig)

fig, ax = plt.subplots(figsize=(7 / 4, 7 / 4), constrained_layout=True)
draw_syn3a(ax)
fig.savefig(os.path.join(OUT, "panel_d_syn3A_ONT_vs_Illumina.pdf"), dpi=300); plt.close(fig)

print("wrote 4 individual panels to", OUT, "(combined Fig S1 assembled manually in Illustrator)")
