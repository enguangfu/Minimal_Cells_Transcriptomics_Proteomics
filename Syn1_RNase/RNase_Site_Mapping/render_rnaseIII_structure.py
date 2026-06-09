#!/usr/bin/env python
"""
render_rnaseIII_structure.py
============================

Draw the local RNA secondary structure at each homology-anchored RNase III cut.
RNase III is a dsRNA endonuclease: when Taggart detected TWO 5' ends for one gene
they are usually the two staggered cuts on the two strands of ONE stem (the 2-nt
3' overhang), so for a 2-cut gene we fold ONE window spanning both cuts and -- if
they base-pair across that fold -- join them with a dashed connector (the dsRNA
double cut).  For a 1-cut gene we fold a symmetric window around the single cut.

Reads the anchored-cleavage TSV written by map_bsub_rnase_to_syn1.py
(syn1_projected_cut_1b per site) and re-folds locally; nucleotides are drawn as
A/G/C/U letters with the cut base(s) in bold red.

Usage:
    render_rnaseIII_structure.py                 -> atpA only (default)
    render_rnaseIII_structure.py all [flank] [fs]-> all 18 homology-hit genes + overview grid
    render_rnaseIII_structure.py MMSYN1_0157     -> one gene by Syn1 locus tag

Individual PDFs: 7/4 x 7/4 in, transparent (drop onto panel f / SI).  Overview
grid: line-style contact sheet of all genes.  Pure Arial; primes are Unicode.
Run in the RNAseq env (needs ViennaRNA `RNA`).
"""
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO
import RNA

HERE = Path(__file__).resolve().parent
PROJECT = Path("/data/enguang/Transcriptomics/Minimal_Cells_Transcriptomics_Proteomics")
CUTS_TSV = HERE / "output" / "rnaseIII" / "rnaseIII_syn1_predicted_cleavage_pairs.tsv"   # structure-derived cut (rank-1 per gene)
SYN1_FASTA = PROJECT / "Genomes_Input" / "syn1_genome.fasta"
OUTDIR = HERE / "output" / "rnaseIII"
STEMDIR = OUTDIR / "stems"

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})

CUT_COL = "#c0392b"      # cut nucleotides (red, matches panel-f cut lines)
BB_COL = "#888888"       # backbone
BP_COL = "#cccccc"       # base-pair rungs
BASE_TXT = "#333333"     # nucleotide letters

LET_FS = 3.2             # nucleotide-letter fontsize for individual PDFs
FLANK_EACH_SIDE = 50     # nt on EACH side of the cleavage-site pair, RNAfold default parameters --
                         #   the exact Taggart/Li recipe (intervening RNA + 50 nt either side; ViennaRNA).
                         #   At this width the two cuts fold into separate local stems within one connected
                         #   structure (not a single duplex), drawn without a connector -- the honest atpA view.
DEFAULT_FLANK = 10       # (legacy) tight flank that folds the single conserved atpA stem
SINGLE_HALF = 45         # (legacy) half-window for a single-cut gene
PAIR_TOL = 4             # nt: two cuts are one duplex if a cut's partner lands within this of the other cut

_GENOME = None


def genome():
    global _GENOME
    if _GENOME is None:
        _GENOME = str(next(SeqIO.parse(str(SYN1_FASTA), "fasta")).seq).upper()
    return _GENOME


def revcomp(s):
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def fold_gene_window(cuts, strand, flank=None):
    """Fold the intervening RNA between the cut(s) plus FLANK_EACH_SIDE nt on either
    side, with RNAfold default parameters -- the exact Taggart/Li recipe.  Returns
    seq, struct, mfe, {genomic_cut: fold_idx}, (window_start, window_end)."""
    g = genome()
    flank = FLANK_EACH_SIDE if flank is None else int(flank)
    gs, ge = max(1, min(cuts) - flank), min(len(g), max(cuts) + flank)
    frag = g[gs - 1:ge]
    rna = (frag if strand == "+" else revcomp(frag)).replace("T", "U")
    struct, mfe = RNA.fold(rna)   # default parameters (37 C, dangles=2, lonely pairs allowed)
    idx = {c: (c - gs if strand == "+" else ge - c) for c in cuts}
    return rna, struct, mfe, idx, (gs, ge)


def _layout(struct):
    co = RNA.naview_xy_coordinates(struct)
    n = len(struct)
    return np.array([co[i].X for i in range(n)]), np.array([co[i].Y for i in range(n)])


def _pairs0(struct):
    pt = RNA.ptable(struct)
    return [(i - 1, pt[i] - 1) for i in range(1, pt[0] + 1) if pt[i] > i], pt


def _align_horizontal(xs, ys):
    """Rotate coordinates so the structure's principal (long-stem) axis is horizontal,
    i.e. parallel to the figure width."""
    pts = np.column_stack([xs, ys]).astype(float)
    pts -= pts.mean(0)
    _, _, vt = np.linalg.svd(pts, full_matrices=False)      # vt[0] = principal axis
    theta = -np.arctan2(vt[0, 1], vt[0, 0])                 # rotate that axis onto x
    rot = pts @ np.array([[np.cos(theta), -np.sin(theta)],
                          [np.sin(theta),  np.cos(theta)]]).T
    return rot[:, 0], rot[:, 1]


def draw_stem(ax, rna, struct, cut_idxs, letters=True, let_fs=LET_FS,
              backbone=True, rotate=False, letter_rot=0):
    """Render one folded window onto ax; mark the cut base(s) red, and connect a
    base-pairing cut pair with a dashed line (the dsRNA double cut).
    backbone=False drops the chain trace (it blurs the letters); rotate=True lays
    the long stem along the width; letter_rot rotates the A/G/C/U glyphs."""
    xs, ys = _layout(struct)
    if rotate:
        xs, ys = _align_horizontal(xs, ys)
    pairs, pt = _pairs0(struct)
    ci = list(cut_idxs)
    ci_set = set(ci)

    for i, j in pairs:
        ax.plot([xs[i], xs[j]], [ys[i], ys[j]], color=BP_COL, lw=0.4, zorder=1)
    if backbone:
        ax.plot(xs, ys, color=BB_COL, lw=0.5, zorder=2)

    if letters:
        for k, ch in enumerate(rna):
            if k in ci_set:
                continue
            ax.text(xs[k], ys[k], ch, fontsize=let_fs, color=BASE_TXT,
                    ha="center", va="center", rotation=letter_rot, zorder=3)
    else:
        ax.scatter(xs, ys, s=2.5, color="#e6e6e6", edgecolors=BB_COL, linewidths=0.1, zorder=3)

    ax.annotate("5′", (xs[0], ys[0]), fontsize=5, color="#666", ha="center", va="center",
                xytext=(-4, 0), textcoords="offset points")
    ax.annotate("3′", (xs[-1], ys[-1]), fontsize=5, color="#666", ha="center", va="center",
                xytext=(4, 0), textcoords="offset points")

    # dashed connector only when the two cuts truly base-pair across one stem
    paired = False
    if len(ci) == 2:
        a, b = ci
        pa = pt[a + 1] - 1 if pt[a + 1] else None
        pb = pt[b + 1] - 1 if pt[b + 1] else None
        paired = (pa is not None and abs(pa - b) <= PAIR_TOL) or (pb is not None and abs(pb - a) <= PAIR_TOL)
        if paired:
            ax.plot([xs[a], xs[b]], [ys[a], ys[b]], color=CUT_COL, lw=0.9, ls=(0, (2, 1.5)), zorder=4)

    for k in ci:
        if letters:
            ax.text(xs[k], ys[k], rna[k], fontsize=let_fs + 1.2, color=CUT_COL,
                    ha="center", va="center", fontweight="bold", rotation=letter_rot, zorder=5)
        else:
            ax.scatter([xs[k]], [ys[k]], s=12, color=CUT_COL, edgecolors="white", linewidths=0.25, zorder=5)

    ax.set_aspect("equal")
    ax.axis("off")
    ax.margins(0.14)
    return paired


def get_cuts(df, locus):
    """The structure-derived rank-1 RNase III cut pair for a gene (genomic cuts)."""
    sub = df[(df["syn1_locus_tag"] == locus) & (df["rank_within_gene"] == 1)]
    if sub.empty:
        sub = df[df["syn1_locus_tag"] == locus]
    r = sub.iloc[0]
    c2 = r.get("genomic_cut2", np.nan)
    cuts = sorted([int(r["genomic_cut1"]), int(c2)]) if pd.notna(c2) else [int(r["genomic_cut1"])]
    strand = str(r["syn1_strand"])
    gene = r["syn1_gene"]
    if not (isinstance(gene, str) and gene.strip()):
        gene = str(r["bsub_gene"])
    return cuts, strand, (gene or locus)


def render_one(df, locus, flank, let_fs):
    STEMDIR.mkdir(parents=True, exist_ok=True)
    cuts, strand, gene = get_cuts(df, locus)
    rna, struct, mfe, idx, _ = fold_gene_window(cuts, strand, flank)
    fig = plt.figure(figsize=(7 / 4, 7 / 4), constrained_layout=True)
    ax = fig.add_subplot(111)
    # backbone off (it blurs the letters); long stem laid along the width; letters kept horizontal
    paired = draw_stem(ax, rna, struct, [idx[c] for c in cuts], letters=True, let_fs=let_fs,
                       backbone=False, rotate=True, letter_rot=0)
    out = STEMDIR / f"R2_{locus}_{gene}_rnaseIII_stem.pdf"
    # crop the PDF box tight to the drawn structure (no blank margins) so it can be
    # placed and scaled freely in Illustrator -- this is a free-floating graphic, not an
    # assembled born-at-size panel, so the tight bbox is intentional here.
    fig.savefig(out, dpi=300, transparent=True, bbox_inches="tight", pad_inches=0.01)
    plt.close(fig)
    tag = "duplex pair" if (len(cuts) == 2 and paired) else f"{len(cuts)} cut(s)"
    print(f"  {locus:14s} {gene:8s} {tag:12s} MFE {mfe:6.1f}  cuts={cuts}  -> {out.name}")
    return cuts, mfe, paired


def render_grid(df, loci, flank):
    cols = 4
    rows = math.ceil(len(loci) / cols)
    fig = plt.figure(figsize=(7, 7 / cols * rows), constrained_layout=True)
    for n, locus in enumerate(loci):
        ax = fig.add_subplot(rows, cols, n + 1)
        cuts, strand, gene = get_cuts(df, locus)
        rna, struct, mfe, idx, _ = fold_gene_window(cuts, strand, flank)
        paired = draw_stem(ax, rna, struct, [idx[c] for c in cuts], letters=False)
        note = "duplex" if (len(cuts) == 2 and paired) else f"{len(cuts)} cut"
        ax.set_title(f"{gene} {locus.split('_')[-1]}  ·  {note}  ·  MFE {mfe:.0f}",
                     fontsize=4.5, color="#333", pad=1)
    out = OUTDIR / "rnaseIII_all_stems_overview.pdf"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"\nOverview grid: {out}")


def main():
    args = sys.argv[1:]
    flank = int(args[1]) if len(args) > 1 and args[1].isdigit() else FLANK_EACH_SIDE
    let_fs = float(args[2]) if len(args) > 2 else LET_FS

    df = pd.read_csv(CUTS_TSV, sep="\t")
    if args and args[0] == "all":
        loci = list(dict.fromkeys(df["syn1_locus_tag"]))   # preserve TSV order
        print(f"Rendering {len(loci)} homology-hit genes -> {STEMDIR}/")
        for locus in loci:
            render_one(df, locus, flank, let_fs)
        render_grid(df, loci, flank)
    else:
        locus = args[0] if args else "MMSYN1_0792"
        render_one(df, locus, flank, let_fs)


if __name__ == "__main__":
    main()
