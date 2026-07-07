#!/usr/bin/env python
"""FigS3 -- the intergenic transcription unit in the syn1 0884/0885 region, Syn1 vs Syn3A.

Born-at-size (OUTPUT.md: Arial 5-7 pt, pdf.fonttype 42, constrained_layout). Absolute Syn1 genome
coordinates on x. One gene track + six depth tracks: Syn1 PacBio, the two Syn1 ONT runs, Syn1
Illumina (replicate-weighted average), then Syn3A ONT and Syn3A Illumina mapped through the single
retained block onto the same Syn1 axis (the deleted block has no Syn3A sequence and is left blank).
Each track is normalised to that library's own BOTH-STRAND genome-mean coverage (x mean; = average
sense-gene depth, matching Figs 2/4/5/6 and Table S1). + strand up (blue), - strand down (orange).

Story: a discrete ~180 bp + strand transcript occupies the \locus{0884}(+)/\locus{0885}(-) gap
without overlapping either gene -- a candidate unannotated transcription unit, seen strongly by
short-read Illumina and one ONT run but under-sampled by the size-selected PacBio library. The
reduction deleted the flanking conserved-hypothetical gene 0884 (+ upstream), while the intergenic
transcript itself was retained but fell ~8-fold in Syn3A (Illumina 21.8x -> 2.8x).

Output: R4_panels/FigS3_intergenic_0884_0885.pdf
"""
import os
import subprocess
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Patch

mpl.rcParams.update({
    "font.size": 7, "font.family": "sans-serif", "font.sans-serif": ["Arial"],
    "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 5,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
OUT  = os.path.join(HERE, "R4_panels")
os.makedirs(OUT, exist_ok=True)
D = lambda p: os.path.join(ROOT, p)
S1_CH, S1L = "CP002027.1", 1_078_809
S3_CH, S3L = "CP016816.2", 543_379

PB  = "Syn1_Transcriptomics/PacBio/PacBio_Processing/depth_bedgraph/syn1.PacBio.FLNC.HQ"
ON1 = "Syn1_Transcriptomics/ONT/ONT_Processing/depth_bedgraph/syn1.ONT.rep1"
ON2 = "Syn1_Transcriptomics/ONT/ONT_Processing/depth_bedgraph/syn1.ONT.rep2"
ILD = "Syn1_Transcriptomics/Illumina/Illumina_Processing/depth_bedgraph"
S3ONT = "Syn3A_Transcriptomics/ONT/ONT_Processing/depth_bedgraph/syn3A.ONT.rep1"
S3ILL = "Syn3A_Transcriptomics/Illumina/Illumina_Processing/depth_bedgraph/syn3A_rep1"
ILLU = [(f"{ILD}/SRR35996296", 0.25), (f"{ILD}/SRR35996297", 0.25), (f"{ILD}/SRR35996298", 0.5)]

WIN = (1042500, 1046000)                 # absolute Syn1 coords
IG  = (1043413, 1043673)                 # 0884(+) end .. 0885(-) start = intergenic gap
DEL = (1040770, 1043410)                 # Syn1 block deleted in Syn3A (aln/raw/syn1_deleted_regions.bed)
S1B, S2B = 1043410, 533601               # retained block start: Syn1 -> Syn3A (no inversion)
GENES = [(1042654, 1043413, "+", "0884"),
         (1043672, 1045562, "-", "0885"),
         (1045709, 1047287, "-", "0886")]
C_PLUS, C_MINUS, C_GENE = "#3182bd", "#e08214", "#cfcfcf"

# track: (label, kind, payload)  kind "s1" -> list[(prefix,weight)] on Syn1 axis; "s3" -> prefix mapped
TRACKS = [
    ("Syn1 PacBio",    "s1", [(PB, 1.0)]),
    ("Syn1 ONT 1",     "s1", [(ON1, 1.0)]),
    ("Syn1 ONT 2",     "s1", [(ON2, 1.0)]),
    ("Syn1 Illumina",  "s1", ILLU),
    ("Syn3A ONT",      "s3", S3ONT),
    ("Syn3A Illumina", "s3", S3ILL),
]


def _read(path, ch, a, b):
    arr = np.zeros(b - a)
    out = subprocess.run(["awk", "-F", "\t", f'$1=="{ch}" && $3>{a} && $2<{b}', D(path)],
                         capture_output=True, text=True).stdout
    for ln in out.splitlines():
        _, s, e, v = ln.split("\t")
        s, e = int(s), int(e)
        arr[max(s, a) - a:min(e, b) - a] = float(v)
    return arr


def _gsum(path, ch):
    o = subprocess.run(["awk", "-F", "\t", '{s+=($3-$2)*$4} END{print s+0}', D(path)],
                       capture_output=True, text=True).stdout.strip()
    return float(o) if o else 0.0


def slice_track(kind, payload, strand, a, b):
    """Per-base depth on the Syn1 x-axis (Syn3A mapped through the retained block; NaN where deleted)."""
    if kind == "s1":
        cov = np.zeros(b - a)
        for pfx, wt in payload:
            cov += wt * _read(f"{pfx}.{strand}.bedGraph", S1_CH, a, b)
        return cov
    cov = np.full(b - a, np.nan)                     # Syn3A: blank in the deleted part
    lo = max(S1B, a)
    if lo < b:
        q0 = S2B + (lo - S1B)
        cov[lo - a:b - a] = _read(f"{payload}.{strand}.bedGraph", S3_CH, q0, q0 + (b - lo))
    return cov


def track_mean(kind, payload):
    if kind == "s1":
        return sum(wt * (_gsum(f"{pfx}.plus.bedGraph", S1_CH) + _gsum(f"{pfx}.minus.bedGraph", S1_CH))
                   for pfx, wt in payload) / S1L
    return (_gsum(f"{payload}.plus.bedGraph", S3_CH) + _gsum(f"{payload}.minus.bedGraph", S3_CH)) / S3L


w0, w1 = WIN
xg = np.arange(w0, w1)
_dl, _dr = max(DEL[0], w0), min(DEL[1], w1)
fig, axes = plt.subplots(len(TRACKS) + 1, 1, figsize=(7, 4.0),
                         height_ratios=[0.5] + [1.0] * len(TRACKS),
                         constrained_layout=True, sharex=True)

# ---- gene track ----
axG = axes[0]
axG.axvspan(_dl, _dr, color="#e8736a", alpha=0.22, lw=0, zorder=0)
axG.text((_dl + _dr) / 2, 0.45, "deleted in Syn3A", ha="center", va="bottom",
         fontsize=5, color="#c0392b", clip_on=False)
for s0, e0, st, lab in GENES:
    a, b = max(s0, w0), min(e0, w1)
    if a >= b:
        continue
    if st == "+":
        axG.add_patch(FancyArrow(a, 0, b - a, 0, width=0.5, head_width=0.5, head_length=min(120, b - a),
                                 length_includes_head=True, color=C_GENE, ec="#8a8a8a", lw=0.4))
    else:
        axG.add_patch(FancyArrow(b, 0, a - b, 0, width=0.5, head_width=0.5, head_length=min(120, b - a),
                                 length_includes_head=True, color=C_GENE, ec="#8a8a8a", lw=0.4))
    axG.text((a + b) / 2, 0, f"{lab}{'+' if st=='+' else '−'}", ha="center", va="center", fontsize=6)
axG.set_ylim(-0.5, 1.0); axG.set_yticks([]); axG.axis("off")

# ---- depth tracks (mirror: + up, - down) ----
report = []
for ax, (name, kind, payload) in zip(axes[1:], TRACKS):
    m = track_mean(kind, payload)
    up = slice_track(kind, payload, "plus", w0, w1) / m
    dn = slice_track(kind, payload, "minus", w0, w1) / m
    ax.axvspan(_dl, _dr, color="#e8736a", alpha=0.22, lw=0, zorder=0)      # deletion block
    ax.axvspan(IG[0], IG[1], color="#fff3d6", zorder=0)                    # intergenic unit
    ax.fill_between(xg, 0, up, color=C_PLUS, lw=0, zorder=2)
    ax.fill_between(xg, 0, -dn, color=C_MINUS, lw=0, zorder=2)
    ax.axhline(0, color="#888888", lw=0.4, zorder=3)
    fu, fd = up[np.isfinite(up)], dn[np.isfinite(dn)]
    top = max(1.0, fu.max() if fu.size else 1.0, fd.max() if fd.size else 1.0)
    ax.set_ylim(-top * 1.08, top * 1.08)
    ax.set_yticks([-round(top), 0, round(top)] if top >= 1.5 else [-round(top, 1), 0, round(top, 1)])
    ax.set_yticklabels([f"{abs(t):g}×" for t in ax.get_yticks()])
    ax.set_ylabel(name, fontsize=6, rotation=0, ha="right", va="center")
    ax.set_xlim(w0, w1)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ig_plus = np.nanmean(slice_track(kind, payload, "plus", IG[0], IG[1]) / m)
    report.append(f"  {name:15s} mean={m:7.0f}  intergenic + = {ig_plus:.2f}x")

axes[-1].set_xlabel("Syn1 genome position (bp)", fontsize=7)
axes[-1].ticklabel_format(axis="x", style="plain")
axes[1].legend(handles=[Patch(color=C_PLUS, label="+ strand"), Patch(color=C_MINUS, label="− strand")],
               frameon=False, fontsize=5, loc="upper right", handlelength=1.0, labelspacing=0.2,
               ncol=2, columnspacing=0.8)

out = os.path.join(OUT, "FigS3_intergenic_0884_0885.pdf")
fig.savefig(out, dpi=300)
plt.close(fig)
print("intergenic (+ strand) mean depth over 1043413-1043673, x both-strand genome-mean:")
print("\n".join(report))
print("wrote", out)
