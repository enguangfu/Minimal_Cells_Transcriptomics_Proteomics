#!/usr/bin/env python
"""Investigate the DCW / mra (division and cell wall) operon across all RNA-seq data.

Operon MMSYN1_0527..0520 (== JCVISYN3A_0527..0520), all on the MINUS strand, retained
and syntenic in Syn3A (same 6,026 bp span). Genes 5'->3':
  0527 (hyp) - rpmF/0526 - mraZ/0525 - mraW/0524 - ftsA/0523 - ftsZ/0522 - sepF/0521 - 0520 (hyp)
(ftsA and sepF are annotated as conserved-hypothetical in Syn1 but named in Syn3A.)

x-axis = relative genome position from the 5' end of 0527 (nt), 5'->3' left to right
(rel = anchor - genomic_pos, anchor = 0527 end0; each organism anchored at its OWN 0527).
For every track the SENSE (minus-strand) depth is drawn upward and the ANTISENSE (plus)
depth downward, each normalised to that library's own BOTH-STRAND genome-mean coverage
(x mean; = average sense-gene depth, matching Figs 2/4/5/6). Syn1 = PacBio, the two ONT
runs, and replicate-weighted Illumina; Syn3A = ONT and Illumina.

Output: dcw_operon/dcw_operon_rnaseq.pdf + dcw_operon_stats.txt
Run from Genome_Reduction/ in the RNAseq env.
"""
import os
import sys
import subprocess
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, FancyArrowPatch, Patch

mpl.rcParams.update({
    "font.size": 7, "font.family": "sans-serif", "font.sans-serif": ["Arial"],
    "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 5,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(HERE, "dcw_operon")
os.makedirs(OUT, exist_ok=True)
D = lambda p: os.path.join(ROOT, p)
S1_CH, S1L = "CP002027.1", 1_078_809
S3_CH, S3L = "CP016816.2", 543_379

PB  = "Syn1_Transcriptomics/PacBio/PacBio_Processing/depth_bedgraph/syn1.PacBio.FLNC.HQ"
ON1 = "Syn1_Transcriptomics/ONT/ONT_Processing/depth_bedgraph/syn1.ONT.rep1"
ON2 = "Syn1_Transcriptomics/ONT/ONT_Processing/depth_bedgraph/syn1.ONT.rep2"
ILD = "Syn1_Transcriptomics/Illumina/Illumina_Processing/depth_bedgraph"
ILLU = [(f"{ILD}/SRR35996296", 0.25), (f"{ILD}/SRR35996297", 0.25), (f"{ILD}/SRR35996298", 0.5)]
S3ONT = "Syn3A_Transcriptomics/ONT/ONT_Processing/depth_bedgraph/syn3A.ONT.rep1"
S3ILL = "Syn3A_Transcriptomics/Illumina/Illumina_Processing/depth_bedgraph/syn3A_rep1"
ISO_TSV = "Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv"

ANC1, ANC3 = 628640, 336007            # 5' end (end0, minus strand) of 0527 in Syn1 / Syn3A
REL_LO, REL_HI = -300, 6350            # rel window: a little 5' of 0527 through past 0520
SYN1_COL, SYN3A_COL = "#3182bd", "#c0392b"
C_ANTI, C_GENE = "#e08214", "#cfcfcf"

# (label, kind, prefix-weight list, anchor, chrom, genome_len). All six depth tracks share one
# log scale (x mean); PacBio is size-selected, so its absolute level is low but its shape is shown.
TRACKS = [
    ("Syn1 PacBio",    "s1", [(PB, 1.0)], ANC1, S1_CH, S1L),
    ("Syn1 ONT 1",     "s1", [(ON1, 1.0)], ANC1, S1_CH, S1L),
    ("Syn1 ONT 2",     "s1", [(ON2, 1.0)], ANC1, S1_CH, S1L),
    ("Syn1 Illumina",  "s1", ILLU,         ANC1, S1_CH, S1L),
    ("Syn3A ONT",      "s3", [(S3ONT, 1.0)], ANC3, S3_CH, S3L),
    ("Syn3A Illumina", "s3", [(S3ILL, 1.0)], ANC3, S3_CH, S3L),
]

# Syn1 genes over the window, with the Syn3A functional name for the two Syn1-hypotheticals.
SYN3A_NAME = {"0523": "ftsA", "0521": "sepF"}


def load_genes():
    rows = []
    for line in open(D("Genomes_Input/syn1.genes.gff3")):
        if line.startswith("#") or "\tgene\t" not in line:
            continue
        p = line.rstrip("\n").split("\t")
        a = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
        rows.append((int(p[3]) - 1, int(p[4]), p[6], a.get("locus_tag", "")[-4:],
                     a.get("Name") or a.get("gene") or ""))
    return rows


def _read(pwl, strand, ch, g_lo, g_hi):
    arr = np.zeros(g_hi - g_lo)
    for pfx, wt in pwl:
        out = subprocess.run(["awk", "-F", "\t", f'$1=="{ch}" && $3>{g_lo} && $2<{g_hi}',
                              D(f"{pfx}.{strand}.bedGraph")], capture_output=True, text=True).stdout
        for ln in out.splitlines():
            _, s, e, v = ln.split("\t")
            s, e = int(s), int(e)
            arr[max(s, g_lo) - g_lo:min(e, g_hi) - g_lo] += wt * float(v)
    return arr


def _gsum(pwl, strand, ch):
    tot = 0.0
    for pfx, wt in pwl:
        o = subprocess.run(["awk", "-F", "\t", '{s+=($3-$2)*$4} END{print s+0}',
                            D(f"{pfx}.{strand}.bedGraph")], capture_output=True, text=True).stdout.strip()
        tot += wt * (float(o) if o else 0.0)
    return tot


def track_mean(pwl, ch, glen):
    return (_gsum(pwl, "plus", ch) + _gsum(pwl, "minus", ch)) / glen


def rel_arrays(pwl, anchor, ch, glen):
    """Sense (minus) and antisense (plus) depth on the rel axis (ascending), x mean."""
    g_lo, g_hi = anchor - REL_HI, anchor - REL_LO
    m = track_mean(pwl, ch, glen)
    sense = _read(pwl, "minus", ch, g_lo, g_hi)
    anti = _read(pwl, "plus", ch, g_lo, g_hi)
    rel = anchor - np.arange(g_lo, g_hi)
    order = np.argsort(rel)
    return rel[order], sense[order] / m, anti[order] / m, m


# ----------------------------------------------------------------- assemble
genes = [g for g in load_genes() if g[0] < ANC1 - REL_LO and g[1] > ANC1 - REL_HI]
xref = None
depth = {}
for name, kind, pwl, anchor, ch, glen in TRACKS:
    xg, sense, anti, m = rel_arrays(pwl, anchor, ch, glen)
    depth[name] = (xg, sense, anti, m)
    xref = xg

# PacBio isoforms (minus strand) spanning any part of the operon window (Syn1 rel axis)
iso = pd.read_csv(D(ISO_TSV), sep="\t")
iso = iso[(iso.chrom == S1_CH) & (iso.strand == "-") &
          (iso.start0 < ANC1 - REL_LO) & (iso.end0 > ANC1 - REL_HI) & (iso.n_reads >= 10)].copy()
iso["r5"] = ANC1 - iso.end0      # 5' end (minus = end0) -> small rel
iso["r3"] = ANC1 - iso.start0    # 3' end -> large rel
iso = iso.sort_values("n_reads", ascending=False).head(40)


def pack(intervals, gap=60):
    order = sorted(range(len(intervals)), key=lambda i: intervals[i][0])
    ends, rows = [], [0] * len(intervals)
    for i in order:
        a, b = intervals[i]
        for r, e in enumerate(ends):
            if a >= e + gap:
                ends[r] = b; rows[i] = r; break
        else:
            ends.append(b); rows[i] = len(ends) - 1
    return rows


nrow = 1 + 1 + len(TRACKS)                      # genes + isoforms + 5 depth
fig, axes = plt.subplots(nrow, 1, figsize=(7, 5.2), sharex=True,
                         height_ratios=[0.7, 1.3] + [1.0] * len(TRACKS), constrained_layout=True)

# ---- gene track (pheT/0528 is a flanking gene, not drawn; labels staggered) ----
axG = axes[0]
op_genes = [g for g in sorted(genes, key=lambda g: ANC1 - g[1]) if g[3] != "0528"]
for i, (s0, e0, st, num, nm) in enumerate(op_genes):
    r_a, r_b = ANC1 - e0, ANC1 - s0          # minus strand: 5'=e0 (small rel), 3'=s0 (large rel)
    lab = f"{SYN3A_NAME[num]}/{num}" if num in SYN3A_NAME else (f"{nm}/{num}" if nm else num)
    axG.add_patch(FancyArrow(r_a, 0, r_b - r_a, 0, width=0.34, head_width=0.34,
                             head_length=min(90, r_b - r_a), length_includes_head=True,
                             color=C_GENE, ec="#8a8a8a", lw=0.4))
    axG.text((r_a + r_b) / 2, 0.46 if i % 2 == 0 else -0.46, lab,
             ha="center", va="center", fontsize=5.0)
axG.set_ylim(-0.8, 0.8); axG.axis("off")

# ---- PacBio isoform track ----
axI = axes[1]
ints = list(zip(iso.r5.astype(int), iso.r3.astype(int)))
rows = pack(ints)
nmax = max(1, iso.n_reads.max())
for (xa, xb), ri, (_, r) in zip(ints, rows, iso.iterrows()):
    lw = float(np.clip(0.3 + 0.7 * np.log10(max(1, r.n_reads)), 0.5, 2.8))
    axI.add_patch(FancyArrowPatch((xa, ri), (xb, ri), arrowstyle="-|>", mutation_scale=5,
                                  lw=lw, color=SYN1_COL, alpha=0.8, shrinkA=0, shrinkB=0))
axI.set_ylim(-1, (max(rows) if rows else 0) + 1.5)
axI.set_ylabel("Syn1 PacBio\nisoforms", fontsize=5.5, rotation=0, ha="right", va="center", color=SYN1_COL)
axI.set_yticks([]); axI.spines[["top", "right", "left"]].set_visible(False)

# rel marks: promoters (green dotted) and TransTermHP terminators (red dotted). Text labels are
# left off for Illustrator. TERM 225 after rpmF; internal mraZ TSS (ONT-defined, rel ~800);
# TERM 223/224 after 0520. (label kept in the tuple only as a code-side reference)
GREEN = "#2ca02c"
TSS1, TSS2, TSS_MRAZ = ANC1 - 628239, ANC1 - 625215, ANC1 - 627840
MARKS = [(TSS1, GREEN, "P (0527)"),
         (TSS_MRAZ, GREEN, "P (mraZ)"),
         (TSS2, GREEN, "P (ftsZ)"),
         (ANC1 - 627898, "#c0392b", "term (TERM 225)"),
         (ANC1 - 622520, "#c0392b", "term (TERM 223/224)")]

# ---- depth tracks (sense only; shared LOG scale, x mean, across all five) ----
FLOOR, TOPY = 0.05, 50.0
report = []
for ax, (name, kind, pwl, anchor, ch, glen) in zip(axes[2:], TRACKS):
    xg, sense, anti, m = depth[name]
    col = SYN1_COL if kind == "s1" else SYN3A_COL
    ax.set_yscale("log")
    ax.fill_between(xg, FLOOR, np.clip(sense, FLOOR, None), color=col, lw=0, zorder=2)
    ax.set_ylim(FLOOR, TOPY)
    ax.set_yticks([0.1, 1, 10]); ax.set_yticklabels(["0.1×", "1×", "10×"])
    ax.set_ylabel(name, fontsize=6, rotation=0, ha="right", va="center", color=col)
    ax.set_xlim(REL_LO, REL_HI)
    ax.spines[["top", "right"]].set_visible(False)
    report.append((name, m))

# genes deleted in JCVI-syn3.0 but restored in Syn3A: 0527 and the 0522-0520 block.
# Red band drawn ONLY on the gene-arrow track.
DEL30 = [(ANC1 - 628640, ANC1 - 628121),          # 0527
         (ANC1 - 625201, ANC1 - 622613)]          # 0522 (ftsZ) .. 0520
for a, b in DEL30:
    axG.axvspan(a, b, facecolor="#e8736a", alpha=0.20, lw=0, zorder=0)
    axG.text((a + b) / 2, 1.05, "deleted in Syn3.0 but added back to Syn3A",
             ha="center", va="bottom", fontsize=4.0, color="#c0392b", clip_on=False)
axG.set_ylim(-0.8, 1.5)

# TSS (green) / terminator (red) guide lines across every track; text added later in Illustrator
for ax in axes:
    for xr, cc, _ in MARKS:
        ax.axvline(xr, color=cc, lw=0.6, ls=(0, (2, 2)), alpha=0.85, zorder=3)

axes[-1].set_xlabel("Relative genome position from the 5′ end of 0527 (nt)", fontsize=7)
axes[-1].ticklabel_format(axis="x", style="plain")

out_pdf = os.path.join(OUT, "dcw_operon_rnaseq.pdf")
fig.savefig(out_pdf, dpi=300)
plt.close(fig)

# ----------------------------------------------------------------- per-gene stats
op = [g for g in genes if 520 <= int(g[3]) <= 527]
op = sorted(op, key=lambda g: ANC1 - g[1])       # 5'->3' order (0527 first)
lines = ["DCW / mra operon (MMSYN1_0527..0520), minus strand", "=" * 60,
         f"anchor (5' of 0527): Syn1 {ANC1}, Syn3A {ANC3};  rel = anchor - genomic_pos", ""]
lines.append("Both-strand genome-mean depth per library (normaliser):")
for name, m in report:
    lines.append(f"  {name:15s} mean = {m:8.1f}")
hdr = f"\n{'gene':16s}" + "".join(f"{n.split()[0][:4]+n.split()[1][:3]:>10s}" for n, _ in report)
lines.append("Per-gene mean SENSE depth (x both-strand genome mean):")
lines.append(hdr)
for s0, e0, st, num, nm in op:
    r_a, r_b = ANC1 - e0, ANC1 - s0
    lab = (SYN3A_NAME.get(num) or nm or "hyp") + "/" + num
    row = f"{lab:16s}"
    for name, m in report:
        xg, sense, anti, _ = depth[name]
        mask = (xg >= r_a) & (xg < r_b)
        row += f"{np.nanmean(sense[mask]):10.2f}"
    lines.append(row)

# paired reduction table context (relTPM + fold change) if available
paired = D("Genome_Reduction/Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv")
if os.path.isfile(paired):
    dfp = pd.read_csv(paired, sep="\t")
    dfp["num"] = dfp.locus_syn1.str[-4:]
    sub = dfp[dfp.num.isin([g[3] for g in op])]
    lines += ["", "Reduction paired table (mean-normalised rel units; TPM=Illumina):",
              f"{'gene':10s}{'relTPM_s1':>11s}{'relTPM_s3':>11s}{'TPM_FC':>9s}{'relIPM_s1':>11s}{'relIPM_s3':>11s}{'iPM_FC':>9s}"]
    for _, r in sub.sort_values("num", ascending=False).iterrows():
        lines.append(f"{r.num:10s}{r.get('relTPM_syn1', np.nan):11.2f}{r.get('relTPM_syn3a', np.nan):11.2f}"
                     f"{r.get('TPM_fold_change', np.nan):9.2f}{r.get('relIPM_syn1', np.nan):11.2f}"
                     f"{r.get('relIPM_syn3a', np.nan):11.2f}{r.get('iPM_fold_change', np.nan):9.2f}")

# ---- promoter (-10) scan at the two block TSSs + TransTermHP terminators in the window ----
sys.path.insert(0, D("Syn1_Operon"))
try:
    import promoter_motif as pm
    lines += ["", "Promoter (-10) scan at the block 5' ends (minus strand):"]
    for tag, tss in [("OP_00267  0527/rpmF block         TSS 628239 (inside 0527, rel 401)", 628239),
                     ("internal mraZ-mraW-ftsA promoter  TSS ~627840 (rel ~800, ONT-defined)", 627840),
                     ("OP_00266  ftsZ block              TSS 625215 (rel 3425)", 625215)]:
        r = pm.scan_minus10(tss, S1_CH, "-")
        lines.append(f"  {tag}")
        lines.append(f"      -10 6-mer {r['minus10_6mer']} (mm{r['mm6']}); 9-mer {r['minus10_9mer']} "
                     f"(mm{r['mm9']}); tier {r['motif_tier']}")
    lines.append("  note: mraZ TSS taken from the least-truncated ONT 5' ends (rel ~800, both organisms);")
    lines.append("        the dominant ONT 5' pile-up at rel ~823 is a downstream 5' truncation hotspot.")
except Exception as e:                                    # noqa
    lines.append(f"  [promoter scan skipped: {e}]")

tt = []
for ln in open(D("Syn1_Operon/syn1_TransTermHP.txt")):
    if ln.strip().startswith("TERM"):
        f = ln.split()
        a, b = int(f[2]), int(f[4])
        if min(a, b) < ANC1 - REL_LO and max(a, b) > ANC1 - REL_HI:
            tt.append(ln.strip())
lines += ["", "TransTermHP intrinsic terminators in the window (rel = anchor - coord):"]
lines += (["  " + t for t in tt] or ["  (none predicted)"])

txt = "\n".join(lines)
with open(os.path.join(OUT, "dcw_operon_stats.txt"), "w") as fh:
    fh.write(txt + "\n")
print(txt)
print("\nwrote", out_pdf)
