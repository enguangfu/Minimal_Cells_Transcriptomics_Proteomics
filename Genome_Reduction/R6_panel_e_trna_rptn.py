#!/usr/bin/env python
"""
R6 panel e -- Syn3A tRNA-operon -> rPtn-operon junction (new neighbour from the deletion).
Co-expression test (ONT spanning/bridging + Illumina gap depth, the 06/07 method): the two
operons are SPLIT -- 0/3084 ONT reads read through and the true inter-operon middle is silent
in BOTH platforms (Illumina mean 27 = 1.2% of flanking; ONT ~0). The panel shows BOTH depth
tracks -- Illumina (uniform, the continuity readout) and ONT direct-RNA (3'-biased) -- each on
its own scale; the silent inter-operon gap is shaded in both. The ONT read-through count is in
R6_stats.txt.

Broken x-axis (minus strand, drawn 5'->3' = high genome coord on the LEFT); the visual break
between the two segments is left empty (add the // in Illustrator):
  left  : the four tRNA ORFs (Thr/Val/Glu/Asn) -> silent gap -> rPtn 5' (rpsJ ...)
  right : rPtn 3' (secY/0652 ...)

Run in the RNAseq conda env (kept consistent with 06/07):
  /home/enguang/anaconda3/envs/RNAseq/bin/python R6_panel_e_trna_rptn.py
"""
import os
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

mpl.rcParams.update({"font.size": 7, "font.family": "sans-serif", "font.sans-serif": ["Arial"],
                     "axes.labelsize": 6, "xtick.labelsize": 5, "ytick.labelsize": 5,
                     "pdf.fonttype": 42, "ps.fonttype": 42})

GR    = os.path.dirname(os.path.abspath(__file__))
ILLBG = os.path.join(GR, "../Syn3A_Transcriptomics/Illumina/Illumina_Processing/depth_bedgraph/syn3A_rep1.minus.bedGraph")
ONTBG = os.path.join(GR, "../Syn3A_Transcriptomics/ONT/ONT_Processing/depth_bedgraph/syn3A.ONT.rep1.minus.bedGraph")
OUT   = os.path.join(GR, "R6_panels/R6e_trna_rptn_syn3A.pdf")
CHROM = "CP016816.2"
GRAY, TEAL, RED = "#7a7a7a", "#2c9e8f", "#c0392b"
ILLC, ILLL = "#9ecae1", "#3182bd"     # Illumina depth fill / line (blue)
ONTC, ONTL = "#dadaeb", "#756bb1"     # ONT depth fill / line (purple)

def depth_minus(path, lo, hi):
    d = np.zeros(hi - lo)
    for ln in open(path):
        f = ln.split()
        if f[0] != CHROM:
            continue
        s, e = int(f[1]), int(f[2])
        if e < lo or s > hi:
            continue
        d[max(s, lo) - lo:min(e, hi) - lo] = float(f[3])
    return np.arange(lo, hi), d

# Syn3A coords (minus strand): rPtn operon 408864-419717; tRNAs 420560-420889
L = (417850, 421000)   # left  : 4 tRNA ORFs, silent gap, rPtn 5' (rpsJ/0672, rplC/0671, rplD/0670)
R = (408700, 411650)   # right : rPtn 3' (secY/0652, rplO/0653, rpsE/0654)
TRNAS = [(420560, 420635), (420648, 420723), (420731, 420806), (420814, 420889)]  # Thr Val Glu Asn
genesL = [(419409, 419717, "rpsJ", GRAY), (418662, 419333, "", GRAY), (418023, 418649, "", GRAY)]
genesR = [(408864, 410312, "secY", GRAY), (410312, 410749, "", GRAY), (410768, 411532, "", GRAY)]
GAP = (419784, 420350)   # true silent inter-operon middle (depth ~0 in both platforms)

fig = plt.figure(figsize=(7 / 3, 7 / 3), constrained_layout=True)
gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.5, 1.5], width_ratios=[1.7, 1.0],
                      hspace=0.10, wspace=0.06)
axg = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
axi = [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]    # Illumina depth
axo = [fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1])]    # ONT depth

def garrow(ax, s, e, fc, y=0, h=1):
    hl = min((e - s) * 0.5, 40)
    ax.add_patch(plt.Polygon([(e, y), (s + hl, y), (s, y + h / 2), (s + hl, y + h), (e, y + h)],
                             closed=True, facecolor=fc, edgecolor="none"))

imax = omax = 0
for ci, (lo, hi) in enumerate((L, R)):
    g = axg[ci]
    if ci == 0:                                   # the four tRNA ORFs (left segment only)
        for s, e in TRNAS:
            garrow(g, s, e, TEAL)
        g.text((TRNAS[0][0] + TRNAS[-1][1]) / 2, 1.45, "tRNAs\n(Asn Glu Val Thr)",
               ha="center", va="bottom", fontsize=5, color=TEAL, linespacing=0.9)
    for s, e, lab, c in (genesL if ci == 0 else genesR):
        garrow(g, s, e, c)
        if lab:
            g.text((s + e) / 2, 1.4, lab, ha="center", va="bottom", fontsize=5, color=c)
    g.set_xlim(hi, lo); g.set_ylim(-0.3, 2.4); g.axis("off")        # inverted x (high coord left)
    xi, di = depth_minus(ILLBG, lo, hi); imax = max(imax, di.max())
    axi[ci].fill_between(xi, di, color=ILLC, lw=0); axi[ci].plot(xi, di, color=ILLL, lw=0.4)
    xo, do = depth_minus(ONTBG, lo, hi); omax = max(omax, do.max())
    axo[ci].fill_between(xo, do, color=ONTC, lw=0); axo[ci].plot(xo, do, color=ONTL, lw=0.4)
    axi[ci].set_xlim(hi, lo); axo[ci].set_xlim(hi, lo)

# silent inter-operon gap: shaded in both depth tracks, labelled once (on the Illumina row)
for a in (axi[0], axo[0]):
    a.axvspan(*GAP, color=RED, alpha=0.08)
axi[0].text(sum(GAP) / 2, imax * 0.55, "silent\ngap", ha="center", va="center",
            fontsize=4.5, color=RED, linespacing=0.9)

for row, mx in ((axi, imax), (axo, omax)):
    for a in row:
        a.set_ylim(0, mx * 1.05)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
        a.tick_params(length=2, pad=1)
    row[1].spines["left"].set_visible(False); row[1].set_yticks([])

axi[0].set_yticks([0, 20000]); axi[0].set_yticklabels(["0", "20k"])
axi[0].set_ylabel("Illumina (−)", fontsize=6)
axo[0].set_yticks([0, 1000]); axo[0].set_yticklabels(["0", "1k"])
axo[0].set_ylabel("ONT (−)", fontsize=6)

# x ticks/labels only on the bottom (ONT) row; the Illumina row hides its x axis
for a in axi:
    a.tick_params(axis="x", labelbottom=False, bottom=False)
axo[0].set_xticks([420500, 419500, 418500]); axo[0].set_xticklabels(["420.5", "419.5", "418.5"])
axo[1].set_xticks([411000, 409500]); axo[1].set_xticklabels(["411", "409.5"])
fig.text(0.5, -0.02, "Syn3A genome position (kb)", fontsize=6, ha="center")

fig.savefig(OUT, dpi=300, bbox_inches="tight")
print(f"wrote {OUT}  (Illumina max {int(imax)}, ONT max {int(omax)})")
