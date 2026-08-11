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

The three Syn1 Illumina replicates are combined exactly as Gene_TPM/Gene_Transcriptomics.py
combines them (two-step: technical reps SRR35996296+97 -> sample_95, then sample_95 and
sample_enr/SRR35996298 equally, i.e. weights 0.25/0.25/0.5), and - as there - each replicate
is normalised to its own library size BEFORE averaging. See rel_arrays().

Two figure versions, selected by the single CLI argument (default `pacbio`). They differ
ONLY in how the transcript-level (isoform) evidence is displayed - the six depth tracks are
identical in both:

  pacbio   ONE isoform panel, CLUSTER display, Syn1 PacBio only.
           Rows are isoform CLUSTERS from the canonical pipeline table (ISO_TSV: raw
           (5',3') tuples merged by complete-linkage at eps=10 bp), pre-filtered to
           n_reads >= 10, then the top 40 by read count, packed onto as few rows as fit.
           Line width encodes log10(n_reads). So one row = one recurrent transcript form,
           and rare forms are absent by construction.
           -> dcw_operon/dcw_operon_rnaseq.pdf + dcw_operon_stats.txt

  stack    FOUR isoform panels, READ-STACK display, one per long-read run (Syn1 PacBio,
           Syn1 ONT 1, Syn1 ONT 2, Syn3A ONT), each placed directly ABOVE its own depth
           track. Rows are individual primary minus-strand alignments - no clustering, no
           abundance filter - sorted by 5' end and evenly subsampled to N_STACK rows, so
           every row is one molecule and rare long transcripts survive. All four runs get
           the identical treatment (same filters, same sort, same row budget), which is
           what makes PacBio vs ONT comparable at ftsA/0523 and 0520.
           (`ont` is accepted as a legacy alias for `stack`.)
           -> dcw_operon/dcw_operon_rnaseq_readstack.pdf + dcw_operon_stats_readstack.txt

Cluster display answers "which transcript forms recur?"; read-stack answers "what did each
library actually sample?". Illumina has no isoform panel in either: 51-nt paired-end
fragments are not transcript-length observations.

Run from Genome_Reduction/ in the RNAseq env (needs samtools on PATH for `stack`).
"""
import os
import re
import sys
import subprocess
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import FancyArrow, FancyArrowPatch, Patch

mpl.rcParams.update({
    "font.size": 7, "font.family": "sans-serif", "font.sans-serif": ["Arial"],
    "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 5,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})

MODE = (sys.argv[1] if len(sys.argv) > 1 else "pacbio").lower()
MODE = "stack" if MODE == "ont" else MODE            # legacy alias
if MODE not in ("pacbio", "stack"):
    sys.exit(f"usage: {os.path.basename(__file__)} [pacbio|stack] [abs|linear|log]")
# Depth y-scale. `pacbio` is frozen on log (it is the published panel). `stack` defaults to
# ABS - raw per-base depth, each run on its own 0..in-window-max axis, with that library's
# genome-average depth drawn as a reference line. No normaliser, nothing to misread: a panel
# says how deep the library actually is here, and the reference line says whether that is
# high or low FOR THAT LIBRARY. The two normalised variants remain available:
#   linear  x genome mean, shared 0..4x axis   (fair across libraries, but PacBio goes flat)
#   log     x genome mean, shared log axis     (shows shape, but visually rescues low coverage:
#           ONT's 0.25x at ftsA fills half a panel while PacBio's 0.065x hugs the floor)
SCALE = (sys.argv[2] if len(sys.argv) > 2 else ("abs" if MODE == "stack" else "log")).lower()
if MODE == "pacbio":
    SCALE = "log"
if SCALE not in ("log", "linear", "abs"):
    sys.exit(f"usage: {os.path.basename(__file__)} [pacbio|stack] [abs|linear|log]")

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
# Syn1 Illumina replicate combination. For the NORMALISED displays we reproduce
# Gene_TPM/Gene_Transcriptomics.py exactly: two-step averaging (tech reps 296+97 -> sample_95,
# then sample_95 vs sample_enr/298 equally = weights .25/.25/.5) applied AFTER each replicate is
# put on its own library-size scale. For the ABSOLUTE display there is nothing to normalise, so
# the three replicates are pooled (summed) - the depth you would get by merging the three BAMs.
ILLU_W = [(f"{ILD}/SRR35996296", 0.25), (f"{ILD}/SRR35996297", 0.25), (f"{ILD}/SRR35996298", 0.5)]
ILLU_POOL = [(pfx, 1.0) for pfx, _ in ILLU_W]
ILLU = ILLU_POOL if SCALE == "abs" else ILLU_W
S3ONT = "Syn3A_Transcriptomics/ONT/ONT_Processing/depth_bedgraph/syn3A.ONT.rep1"
S3ILL = "Syn3A_Transcriptomics/Illumina/Illumina_Processing/depth_bedgraph/syn3A_rep1"
ISO_TSV = "Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv"
# Alignment file behind each read-stack panel, keyed by the depth-track label it pairs with.
# Illumina is deliberately absent: 51-nt paired-end fragments are not transcript observations.
READ_BAM = {
    "Syn1 PacBio": "Syn1_Transcriptomics/PacBio/PacBio_Processing/syn1.PacBio.FLNC.sorted.HQ.bam",
    "Syn1 ONT 1":  "Syn1_Transcriptomics/ONT/ONT_Processing/syn1.ONT.rep1.sorted.bam",
    "Syn1 ONT 2":  "Syn1_Transcriptomics/ONT/ONT_Processing/syn1.ONT.rep2.sorted.bam",
    "Syn3A ONT":   "Syn3A_Transcriptomics/ONT/ONT_Processing/syn3A.ONT.rep1.sorted.bam",
}

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

# Long-read runs that get a read-stack panel, in TRACKS order; each sits above its own depth track.
STACK_TRACKS = [t for t in TRACKS if t[0] in READ_BAM]

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
    """Sense (minus) and antisense (plus) depth on the rel axis (ascending), x mean.

    Replicates are combined the way Gene_TPM/Gene_Transcriptomics.py combines them: each
    replicate is normalised to ITS OWN library size first (per-sample TPM there, depth /
    own both-strand genome mean here), and only then weight-averaged. Doing it the other
    way round - averaging raw depth and dividing by one pooled mean - silently re-weights
    by library size, and SRR35996298 is 2.2x deeper than SRR35996296/97, so it would pull
    the merged Syn1 Illumina track well beyond its intended 0.5 share.
    The returned `m` is the weight-averaged per-replicate genome mean (what "1x" refers to);
    for the five single-library tracks this is identical to the old pooled behaviour. Under
    SCALE == "abs" nothing is divided: arrays stay in raw per-base depth and `m` is the pooled
    genome-average depth of the same weighted set, drawn as the reference line.
    """
    g_lo, g_hi = anchor - REL_HI, anchor - REL_LO
    sense = np.zeros(g_hi - g_lo)
    anti = np.zeros(g_hi - g_lo)
    m = 0.0
    for pfx, wt in pwl:
        one = [(pfx, 1.0)]
        m_i = track_mean(one, ch, glen)
        s_i = _read(one, "minus", ch, g_lo, g_hi)
        a_i = _read(one, "plus", ch, g_lo, g_hi)
        if SCALE != "abs":
            s_i, a_i = s_i / m_i, a_i / m_i
        sense += wt * s_i
        anti += wt * a_i
        m += wt * m_i
    rel = anchor - np.arange(g_lo, g_hi)
    order = np.argsort(rel)
    return rel[order], sense[order], anti[order], m


N_STACK = 190   # rows drawn per read-stack panel (even subsample of the 5'-sorted read list)


def read_spans(bam, ch, anchor):
    """(rel5', rel3') per primary MINUS-strand alignment overlapping the window, sorted by 5'.

    Identical treatment for all four long-read runs, which is what makes them comparable.
    PacBio FLNC and ONT direct-RNA reads are both 5'->3' sense, so BAM strand == transcript
    strand (same premise as the two 0*_sequencing_depth.sh scripts): -f 0x10 keeps the
    minus-strand = sense reads of this operon, -F 0x904 drops unmapped/secondary/supplementary
    (a no-op on the PacBio FLNC BAM, which is primary-only). Reference span comes from the
    CIGAR with soft clips excluded, so one row is exactly one of the alignments the depth
    track below it is built from.
    """
    g_lo, g_hi = anchor - REL_HI, anchor - REL_LO
    out = subprocess.run(["samtools", "view", "-F", "0x904", "-f", "0x10", D(bam),
                          f"{ch}:{g_lo + 1}-{g_hi}"], capture_output=True, text=True).stdout
    spans = []
    for ln in out.splitlines():
        f = ln.split("\t")
        s0 = int(f[3]) - 1
        e0 = s0 + sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]) if op in "MDN=X")
        spans.append((anchor - e0, anchor - s0))       # minus strand: 5' = e0 -> small rel
    return np.array(sorted(spans), dtype=float).reshape(-1, 2)


def stack_rows(sp, n=N_STACK):
    """Evenly subsample the 5'-sorted read list down to <=n rows, preserving the landscape."""
    if len(sp) <= n:
        return sp
    return sp[np.unique(np.linspace(0, len(sp) - 1, n).round().astype(int))]


# ----------------------------------------------------------------- assemble
genes = [g for g in load_genes() if g[0] < ANC1 - REL_LO and g[1] > ANC1 - REL_HI]
xref = None
depth = {}
for name, kind, pwl, anchor, ch, glen in TRACKS:
    xg, sense, anti, m = rel_arrays(pwl, anchor, ch, glen)
    depth[name] = (xg, sense, anti, m)
    xref = xg

iso, spans = None, {}
if MODE == "pacbio":
    # PacBio isoform CLUSTERS (minus strand) spanning any part of the window (Syn1 rel axis)
    iso = pd.read_csv(D(ISO_TSV), sep="\t")
    iso = iso[(iso.chrom == S1_CH) & (iso.strand == "-") &
              (iso.start0 < ANC1 - REL_LO) & (iso.end0 > ANC1 - REL_HI) & (iso.n_reads >= 10)].copy()
    iso["r5"] = ANC1 - iso.end0      # 5' end (minus = end0) -> small rel
    iso["r3"] = ANC1 - iso.start0    # 3' end -> large rel
    iso = iso.sort_values("n_reads", ascending=False).head(40)
else:
    for nm, _kind, _pwl, anchor, ch, _glen in STACK_TRACKS:
        spans[nm] = read_spans(READ_BAM[nm], ch, anchor)


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


# Row plan. pacbio: gene, one cluster panel, then the six depth tracks.
# stack:  gene, then for each depth track its read stack (if it has one) immediately above it.
if MODE == "pacbio":
    ROWS = [("gene", None), ("iso", None)] + [("depth", t) for t in TRACKS]
    HGT = [1.6, 1.3] + [1.0] * len(TRACKS)
    fig_h = 5.9
else:
    ROWS, HGT = [("gene", None)], [1.6]
    for i, t in enumerate(TRACKS):
        if i:
            ROWS.append(("gap", None)); HGT.append(0.45)   # blank row between runs
        if t[0] in READ_BAM:
            ROWS.append(("stack", t)); HGT.append(1.25)
        ROWS.append(("depth", t)); HGT.append(1.0)
    fig_h = 9.0
fig, axes = plt.subplots(len(ROWS), 1, figsize=(7, fig_h), sharex=True,
                         height_ratios=HGT, constrained_layout=True)
if MODE == "stack":
    # near-zero row gap so each read stack sits flush on its own depth track; the spacer
    # rows above then carry all of the visual separation BETWEEN runs.
    fig.get_layout_engine().set(h_pad=0.02, hspace=0.005)
    for ax, (rk, _t) in zip(axes, ROWS):
        if rk == "gap":
            ax.axis("off")
        if rk in ("gene", "stack", "gap"):
            ax.tick_params(axis="x", length=0)   # ticks belong to the bottom axis only

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

# ---- isoform panel(s): PacBio clusters (packed arrows) or ONT read stacks ----
if MODE == "pacbio":
    axI = axes[1]
    ints = list(zip(iso.r5.astype(int), iso.r3.astype(int)))
    rows = pack(ints)
    for (xa, xb), ri, (_, r) in zip(ints, rows, iso.iterrows()):
        lw = float(np.clip(0.3 + 0.7 * np.log10(max(1, r.n_reads)), 0.5, 2.8))
        axI.add_patch(FancyArrowPatch((xa, ri), (xb, ri), arrowstyle="-|>", mutation_scale=5,
                                      lw=lw, color=SYN1_COL, alpha=0.8, shrinkA=0, shrinkB=0))
    axI.set_ylim(-1, (max(rows) if rows else 0) + 1.5)
    axI.set_ylabel("Syn1 PacBio\nisoforms", fontsize=5.5, rotation=0, ha="right", va="center",
                   color=SYN1_COL)
    axI.set_yticks([]); axI.spines[["top", "right", "left"]].set_visible(False)
else:
    # one thin line per read, rows sorted by 5' end and evenly subsampled to N_STACK.
    # y is inverted so the stack reads top->bottom in 5'->3' order (IGV-like). Every panel
    # gets the same row budget so the four runs are directly comparable by eye.
    for ax, (rk, t) in zip(axes, ROWS):
        if rk != "stack":
            continue
        name, kind = t[0], t[1]
        sp = spans[name]
        shown = stack_rows(sp)
        col = SYN1_COL if kind == "s1" else SYN3A_COL
        ax.add_collection(LineCollection([[(a, i), (b, i)] for i, (a, b) in enumerate(shown)],
                                         colors=col, linewidths=0.35, alpha=0.75, zorder=2))
        ax.set_ylim(len(shown) + 1, -2)
        lab = f"{name} reads\nn={len(sp):,}"
        if len(shown) < len(sp):
            lab += f" ({len(shown)} drawn)"
        ax.set_ylabel(lab, fontsize=5.5, rotation=0, ha="right", va="center", color=col)
        ax.set_yticks([]); ax.spines[["top", "right", "bottom", "left"]].set_visible(False)

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

# ---- depth tracks (sense only; one shared scale, x mean, across all six) ----
FLOOR, TOPY = 0.05, 50.0     # log-scale floor / ceiling
# Shared LINEAR ceiling. 4x keeps the region of interest (ftsA/0520 span 0.01-2.4x) legible
# instead of squashing it into the bottom third; the few higher peaks clip and are labelled.
YMAX_LIN = 4.0
report = []
for ax, (rk, t) in zip(axes, ROWS):
    if rk != "depth":
        continue
    name, kind = t[0], t[1]
    xg, sense, anti, m = depth[name]
    col = SYN1_COL if kind == "s1" else SYN3A_COL
    if SCALE == "abs":
        # raw per-base depth, own axis per run: 0 .. in-window max. Nothing is normalised, so
        # the panel cannot understate a library the way a shared normalised axis can.
        top = float(max(sense.max(), 1.0)) * 1.08
        ax.fill_between(xg, 0, sense, color=col, lw=0, zorder=2)
        ax.set_ylim(0, top)
        ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=3, integer=True))
        # that library's genome-average depth: the "is this high or low for me" reference.
        if m <= top:
            ax.axhline(m, color="#555555", lw=0.6, ls=(0, (3, 2)), zorder=4)
            ax.annotate(f"genome avg {m:,.0f}", xy=(REL_HI, m), xytext=(-2, 1.5),
                        textcoords="offset points", ha="right", va="bottom", fontsize=4.5,
                        color="#555555", zorder=6,
                        bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.8))
        else:   # PacBio: the whole operon sits below its own library average
            ax.annotate(f"genome avg {m:,.0f} (above panel, max here {sense.max():,.0f})",
                        xy=(REL_HI, top), xytext=(-2, -1.5), textcoords="offset points",
                        ha="right", va="top", fontsize=4.5, color="#555555", zorder=6,
                        bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.8))
    elif SCALE == "log":
        ax.set_yscale("log")
        ax.fill_between(xg, FLOOR, np.clip(sense, FLOOR, None), color=col, lw=0, zorder=2)
        ax.set_ylim(FLOOR, TOPY)
        ax.set_yticks([0.1, 1, 10]); ax.set_yticklabels(["0.1×", "1×", "10×"])
    else:
        ax.fill_between(xg, 0, np.clip(sense, 0, YMAX_LIN), color=col, lw=0, zorder=2)
        ax.set_ylim(0, YMAX_LIN * 1.02)
        ax.set_yticks([0, 1, 2, 4]); ax.set_yticklabels(["0", "1×", "2×", "4×"])
        ax.axhline(1.0, color="#555555", lw=0.5, ls=(0, (3, 3)), zorder=4)   # library average
        if (sense > YMAX_LIN).any():      # keep the clipped peak honest by naming its height
            i = int(np.argmax(sense))
            ax.annotate(f"peak {sense[i]:.0f}×", xy=(xg[i], YMAX_LIN), xytext=(3, -0.5),
                        textcoords="offset points", ha="left", va="top", fontsize=5, color=col,
                        zorder=6, bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.85))
    if MODE == "stack":
        # Mark bases with genuinely ZERO coverage, so "low" can never be misread as "absent".
        # Only Syn1 Illumina has any here (23.4% of ftsA); PacBio's minimum over the whole
        # window is 169x raw, i.e. never zero despite sitting near the log floor.
        lo_, hi_ = ax.get_ylim()
        top = lo_ * 1.7 if SCALE == "log" else lo_ + 0.055 * (hi_ - lo_)
        ax.fill_between(xg, lo_, top, where=(sense <= 0), color="#111111", lw=0, zorder=5)
        # 1x is that library's own both-strand genome-mean raw depth -> absolute depth stays
        # recoverable from the panel, which the normalised axis alone hides.
        ylab = name if SCALE == "abs" else f"{name}\n1× = {m:,.0f}"
        ax.set_ylabel(ylab, fontsize=6, rotation=0, ha="right", va="center", color=col)
    else:
        ax.set_ylabel(name, fontsize=6, rotation=0, ha="right", va="center", color=col)
    ax.set_xlim(REL_LO, REL_HI)
    ax.spines[["top", "right"]].set_visible(False)
    report.append((name, m))

# genes deleted in JCVI-syn3.0 but restored in Syn3A: 0527 and the 0522-0520 block.
# Red band drawn ONLY on the gene-arrow track.
DEL30 = [(ANC1 - 628640, ANC1 - 628121),          # 0527
         (ANC1 - 625201, ANC1 - 622613)]          # 0522 (ftsZ) .. 0520
for a, b in DEL30:                                 # red band on the gene-arrow track only; label added in Illustrator
    axG.axvspan(a, b, facecolor="#e8736a", alpha=0.20, lw=0, zorder=0)
axG.set_ylim(-1.2, 1.8)                             # extra head/foot room for hand-added markers

# TSS (green) / terminator (red) guide lines across every track; text added later in Illustrator
for ax, (rk, _t) in zip(axes, ROWS):          # ROWS has no "gap" rows in pacbio mode
    if rk == "gap":
        continue
    for xr, cc, _ in MARKS:
        ax.axvline(xr, color=cc, lw=0.6, ls=(0, (2, 2)), alpha=0.85, zorder=3)

xlab = "Relative genome position from the 5′ end of 0527 (nt)"
if SCALE == "abs":
    xlab += "        y = sense-strand read depth (raw ×)"
axes[-1].set_xlabel(xlab, fontsize=7)
axes[-1].ticklabel_format(axis="x", style="plain")

sfx = "" if MODE == "pacbio" else "_readstack_" + {"abs": "abs", "linear": "rel", "log": "log"}[SCALE]
out_pdf = os.path.join(OUT, f"dcw_operon_rnaseq{sfx}.pdf")
fig.savefig(out_pdf, dpi=300)
plt.close(fig)

# ----------------------------------------------------------------- per-gene stats
op = [g for g in genes if 520 <= int(g[3]) <= 527]
op = sorted(op, key=lambda g: ANC1 - g[1])       # 5'->3' order (0527 first)
lines = ["DCW / mra operon (MMSYN1_0527..0520), minus strand", "=" * 60,
         f"anchor (5' of 0527): Syn1 {ANC1}, Syn3A {ANC3};  rel = anchor - genomic_pos", ""]
lines.append("Both-strand genome-mean depth per library "
             + ("(reference line in the figure):" if SCALE == "abs" else "(normaliser):"))
for name, m in report:
    lines.append(f"  {name:15s} mean = {m:8.1f}")
hdr = f"\n{'gene':16s}" + "".join(f"{n.split()[0][:4]+n.split()[1][:3]:>10s}" for n, _ in report)
lines.append("Per-gene mean SENSE depth "
             + ("(RAW per-base depth):" if SCALE == "abs" else "(x both-strand genome mean):"))
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

if MODE == "stack":
    lines += ["", "Read stacks - primary MINUS-strand (= sense) alignments in the window.",
              "Syn3A is syntenic over this operon (same 6,026 bp span), so the same rel window applies.",
              f"{'run':12s}{'reads':>8s}{'med_len':>9s}{'p90_len':>9s}{'rows_drawn':>12s}"]
    for rn, *_ in STACK_TRACKS:
        sp = spans[rn]
        L = sp[:, 1] - sp[:, 0]
        lines.append(f"{rn:12s}{len(sp):8d}{np.median(L):9.0f}{np.percentile(L, 90):9.0f}"
                     f"{len(stack_rows(sp)):12d}")
    lines += ["", "Reads on each gene, per run  (overlap / covering >=50% of the gene / full span):",
              f"{'gene':12s}" + "".join(f"{rn:>20s}" for rn, *_ in STACK_TRACKS)]
    for s0, e0, st, num, gnm in op:
        r_a, r_b = ANC1 - e0, ANC1 - s0
        row = f"{(SYN3A_NAME.get(num) or gnm or 'hyp') + '/' + num:12s}"
        for rn, *_ in STACK_TRACKS:
            sp = spans[rn]
            ov = int(((sp[:, 0] < r_b) & (sp[:, 1] > r_a)).sum())
            half = int((np.minimum(sp[:, 1], r_b) - np.maximum(sp[:, 0], r_a) >= 0.5 * (r_b - r_a)).sum())
            full = int(((sp[:, 0] <= r_a) & (sp[:, 1] >= r_b)).sum())
            row += f"{ov:>10d} /{half:4d} /{full:<4d}"
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
with open(os.path.join(OUT, f"dcw_operon_stats{sfx}.txt"), "w") as fh:
    fh.write(txt + "\n")
print(txt)
print("\nwrote", out_pdf)
