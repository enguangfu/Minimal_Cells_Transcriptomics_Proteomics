#!/usr/bin/env python3
"""
05_deletion_operon.py

Single-bp overlap of Syn1 -> Syn3A deletion intervals against Syn1 operons.

Per-hit classification — every overlapping deletion is classified on its own,
then aggregated for the operon. Per-hit cases:

    fully_deleted         - this deletion covers both operon boundaries
    5'_truncation_gene    - cuts the 5' boundary AND reaches into a sense-gene CDS
    5'_truncation_UTR     - cuts the 5' boundary but stops in the 5' UTR
    3'_truncation_gene    - cuts the 3' boundary AND reaches into a sense-gene CDS
    3'_truncation_UTR     - cuts the 3' boundary but stops in the 3' UTR
    intra_truncated       - lies strictly inside the operon, neither boundary touched

5'/3' is strand-aware: + strand 5' = start0, - strand 5' = end0.

Operon-level aggregate:
    intact   - no overlapping deletion
    <case>   - single hit
    multi:<case_a>+<case_b>+...  - >=2 hits, unique cases sorted

Coordinates are 0-based half-open [start0, end0).
"""

from __future__ import annotations

import io
import re
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from PyPDF2 import PdfReader, PdfWriter

ROOT     = Path(__file__).resolve().parent.parent
OPERONS  = ROOT / "Operon_Annotation_Visualization" / "operons.candidate_blocks.tsv"
GFF      = ROOT / "Genomes_Input" / "syn1.genes.gff3"
DEL_BED  = Path(__file__).resolve().parent / "aln" / "raw" / "syn1_deleted_regions.bed"
OUT_DIR  = Path(__file__).resolve().parent / "aln" / "analysis"
OUT_TSV  = OUT_DIR / "operon_deletion_classification.tsv"
OPERON_PLOTS_DIR = ROOT / "Operon_Annotation_Visualization" / "operon_plots"


# ---------------------------------------------------------------- loaders

def load_deletions(bed: Path) -> pd.DataFrame:
    df = pd.read_csv(bed, sep="\t")
    df = df.rename(columns={c: c.lstrip("#") for c in df.columns})
    return df[["chrom", "start0", "end"]].rename(columns={"end": "end0"})


def load_genes(gff: Path) -> pd.DataFrame:
    """Return one row per locus_tag with 0-based half-open coordinates."""
    rows = []
    pat_locus = re.compile(r"locus_tag=([^;]+)")
    pat_name  = re.compile(r"(?:Name|gene)=([^;]+)")
    with gff.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom, start1, end1, strand, attr = f[0], int(f[3]), int(f[4]), f[6], f[8]
            m = pat_locus.search(attr)
            if not m:
                continue
            n = pat_name.search(attr)
            rows.append({
                "locus_tag": m.group(1),
                "name": n.group(1) if n else "",
                "chrom": chrom,
                "strand": strand,
                "start0": start1 - 1,
                "end0":   end1,
            })
    df = pd.DataFrame(rows).drop_duplicates(subset="locus_tag").set_index("locus_tag")
    return df


# ---------------------------------------------------------------- per-hit classify

def classify_hit(del_s: int, del_e: int,
                 op_s: int, op_e: int, strand: str,
                 first_gene_s: int | None, first_gene_e: int | None,
                 last_gene_s:  int | None, last_gene_e:  int | None) -> str:
    """Classify a single deletion vs one operon. Coords half-open."""
    left  = del_s <= op_s and del_e > op_s     # covers operon's left boundary
    right = del_s <  op_e and del_e >= op_e    # covers operon's right boundary

    if left and right:
        return "fully_deleted"

    if not left and not right:
        return "intra_truncated"

    # Decide which side (5' vs 3') the boundary belongs to:
    if strand == "+":
        side = "5" if left else "3"
    else:
        side = "5" if right else "3"

    # Determine if the deletion reaches a sense-gene CDS or stays in UTR.
    # The "first gene" is the 5'-most sense gene; "last gene" is the 3'-most.
    # On + strand: first gene has smallest start0; last gene has largest end0.
    # On - strand: first gene has largest end0;   last gene has smallest start0.
    in_gene = False
    if side == "5":
        # 5' UTR boundary side. Does the deletion reach the first gene's CDS?
        if first_gene_s is not None:
            if strand == "+":
                # 5' UTR = [op_s, first_gene_s); CDS reached if del_e > first_gene_s
                in_gene = del_e > first_gene_s
            else:
                # 5' UTR = [first_gene_e, op_e); CDS reached if del_s < first_gene_e
                in_gene = del_s < first_gene_e
    else:  # side == "3"
        if last_gene_s is not None:
            if strand == "+":
                # 3' UTR = [last_gene_e, op_e); CDS reached if del_s < last_gene_e
                in_gene = del_s < last_gene_e
            else:
                # 3' UTR = [op_s, last_gene_s); CDS reached if del_e > last_gene_s
                in_gene = del_e > last_gene_s

    suffix = "gene" if in_gene else "UTR"
    return f"{side}'_truncation_{suffix}"


def operon_first_last(op, genes: pd.DataFrame) -> tuple[tuple|None, tuple|None]:
    """Return ((first_s,first_e), (last_s,last_e)) sense genes, strand-aware.
    Returns (None, None) for either side when no sense gene is annotated."""
    raw = op.sense_gene_loci if isinstance(op.sense_gene_loci, str) else ""
    seen: set = set()
    loci = []
    for t in raw.split(","):
        t = t.strip()
        if t and t not in seen and t in genes.index:
            seen.add(t)
            loci.append(t)
    sense = genes.loc[loci] if loci else pd.DataFrame()
    if not sense.empty:
        # Operon membership in candidate_blocks.tsv is loose — keep only sense
        # genes whose CDS span is fully contained within [op.start0, op.end0).
        sense = sense[(sense.start0 >= op.start0) & (sense.end0 <= op.end0)]
    if sense.empty:
        return None, None
    if op.strand == "+":
        first = sense.loc[sense.start0.idxmin()]
        last  = sense.loc[sense.end0.idxmax()]
    else:
        first = sense.loc[sense.end0.idxmax()]
        last  = sense.loc[sense.start0.idxmin()]
    return (int(first.start0), int(first.end0)), (int(last.start0), int(last.end0))


# ---------------------------------------------------------------- main analysis

def analyse(operons: pd.DataFrame, dels: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for op in operons.itertuples():
        d_chr = dels[dels.chrom == op.chrom]
        hit = d_chr[(d_chr.start0 < op.end0) & (d_chr.end0 > op.start0)].copy()

        first_g, last_g = operon_first_last(op, genes)
        fgs, fge = first_g if first_g else (None, None)
        lgs, lge = last_g  if last_g  else (None, None)

        per_hit_cases: list[str] = []
        overlap_strs:  list[str] = []
        for h in hit.itertuples():
            c = classify_hit(int(h.start0), int(h.end0),
                             int(op.start0), int(op.end0), op.strand,
                             fgs, fge, lgs, lge)
            per_hit_cases.append(c)
            overlap_strs.append(f"{h.chrom}:{h.start0}-{h.end0}")

        if not per_hit_cases:
            agg = "intact"
        elif len(per_hit_cases) == 1:
            agg = per_hit_cases[0]
        else:
            agg = "multi:" + "+".join(sorted(set(per_hit_cases)))

        rows.append({
            "operon_id": op.operon_id,
            "chrom": op.chrom,
            "strand": op.strand,
            "start0": op.start0,
            "end0": op.end0,
            "length": op.end0 - op.start0,
            "sense_gene_names": op.sense_gene_names,
            "deletion_class": agg,
            "n_overlapping_deletions": len(per_hit_cases),
            "per_hit_classes": ";".join(per_hit_cases),
            "overlapping_deletions": ";".join(overlap_strs),
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------- visualization

CASE_COLOR = {
    "fully_deleted":         "#444444",
    "5'_truncation_gene":    "#d62728",
    "5'_truncation_UTR":     "#ff9896",
    "3'_truncation_gene":    "#1f77b4",
    "3'_truncation_UTR":     "#aec7e8",
    "intra_truncated":       "#9467bd",
}


CATEGORY_FILES = {
    "5'_truncation_gene":  "operon_deletion__5p_trunc_gene.pdf",
    "5'_truncation_UTR":   "operon_deletion__5p_trunc_UTR.pdf",
    "3'_truncation_gene":  "operon_deletion__3p_trunc_gene.pdf",
    "3'_truncation_UTR":   "operon_deletion__3p_trunc_UTR.pdf",
    "intra_truncated":     "operon_deletion__intra_truncated.pdf",
    "multi":               "operon_deletion__multi_hit.pdf",
}


def _bucket(deletion_class: str) -> str | None:
    if deletion_class.startswith("multi:"):
        return "multi"
    if deletion_class in CATEGORY_FILES:
        return deletion_class
    return None  # intact / fully_deleted -> no panel


PANEL_FIGSIZE = (19.9, 4.0)  # inches; width matches OP_*_wdepth.pdf for side-by-side stacking


def _draw_panel(writer: PdfWriter, r, op, dels: pd.DataFrame,
                genes: pd.DataFrame, append_wdepth: bool = True) -> None:
    op_s, op_e = int(op.start0), int(op.end0)
    pad = max(200, int(0.05 * (op_e - op_s)))
    xlim = (op_s - pad, op_e + pad)

    fig, ax = plt.subplots(figsize=PANEL_FIGSIZE)

    # y bands
    Y_OPERON   = (0.78, 0.15)   # (bottom, height)
    Y_SENSE    = 0.62           # arrow center
    Y_SENSE_LB = 0.53
    Y_DEL      = (0.32, 0.15)
    Y_ANTI     = 0.15           # arrow center for antisense
    Y_ANTI_LB  = 0.06

    # operon span
    ax.add_patch(mpatches.Rectangle(
        (op_s, Y_OPERON[0]), op_e - op_s, Y_OPERON[1],
        facecolor="#cccccc", edgecolor="black", lw=0.6, zorder=1))
    ax.text((op_s + op_e) / 2, Y_OPERON[0] + Y_OPERON[1] + 0.05,
            f"{r.operon_id}  ({op.strand})  {r.deletion_class}",
            ha="center", va="bottom", fontsize=18)

    def draw_gene(g, y, color):
        # Draw arrow in the gene's own transcription direction (data coords).
        if g.strand == "+":
            x0, dx = int(g.start0), int(g.end0) - int(g.start0)
        else:
            x0, dx = int(g.end0), int(g.start0) - int(g.end0)
        ax.add_patch(mpatches.FancyArrow(
            x0, y, dx, 0,
            width=0.06, head_width=0.10,
            head_length=min(120, abs(dx) * 0.25),
            length_includes_head=True,
            facecolor=color, edgecolor="black", lw=0.4, zorder=2))

    def gene_subset(loci_field: str | float):
        raw = loci_field if isinstance(loci_field, str) else ""
        out = []
        for lt in [t.strip() for t in raw.split(",") if t.strip()]:
            if lt not in genes.index:
                continue
            g = genes.loc[lt]
            if int(g.start0) < op_s or int(g.end0) > op_e:
                continue
            out.append((lt, g))
        return out

    # sense genes (own strand == operon strand)
    for lt, g in gene_subset(op.sense_gene_loci):
        draw_gene(g, Y_SENSE, "#2ca02c")
        lbl = g["name"] or lt
        ax.text((int(g.start0) + int(g.end0)) / 2, Y_SENSE_LB,
                f"{lbl}\n{int(g.start0)}-{int(g.end0)}",
                ha="center", va="top", fontsize=12)

    # antisense genes (opposite strand, fully contained)
    for lt, g in gene_subset(op.antisense_gene_loci):
        draw_gene(g, Y_ANTI, "#ff7f0e")
        lbl = g["name"] or lt
        ax.text((int(g.start0) + int(g.end0)) / 2, Y_ANTI_LB,
                f"{lbl}\n{int(g.start0)}-{int(g.end0)}",
                ha="center", va="top", fontsize=12)

    # deletions
    d_chr = dels[(dels.chrom == op.chrom)
                 & (dels.start0 < op_e) & (dels.end0 > op_s)]
    per_hit = (r.per_hit_classes or "").split(";")
    for (h, case) in zip(d_chr.itertuples(), per_hit):
        color = CASE_COLOR.get(case, "#888888")
        ax.add_patch(mpatches.Rectangle(
            (h.start0, Y_DEL[0]), h.end0 - h.start0, Y_DEL[1],
            facecolor=color, edgecolor="black", lw=0.5, alpha=0.85, zorder=3))
        ax.text((h.start0 + h.end0) / 2, Y_DEL[0] + Y_DEL[1] / 2,
                f"{case}\n{h.start0}-{h.end0}",
                ha="center", va="center", fontsize=12,
                color="white" if case == "fully_deleted" else "black")

    # operon boundary ticks
    for x in (op_s, op_e):
        ax.axvline(x, color="black", lw=0.4, ls=":")
        ax.text(x, 1.02, f"{x}", ha="center", va="bottom", fontsize=14)

    ax.set_xlim(*xlim)
    if op.strand == "-":
        # transcription left-to-right: high coord on the left, low on the right
        ax.invert_xaxis()
    ax.set_ylim(0, 1.10)
    ax.set_yticks([])
    ax.set_xlabel("Syn1 genome coordinate (bp)" + ("  [- strand: inverted]" if op.strand == "-" else ""))
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)

    handles = [mpatches.Patch(color=c, label=k) for k, c in CASE_COLOR.items()]
    handles.append(mpatches.Patch(color="#2ca02c", label="sense gene"))
    handles.append(mpatches.Patch(color="#ff7f0e", label="antisense gene"))
    handles.append(mpatches.Patch(color="#cccccc", label="operon span"))
    ax.legend(handles=handles, loc="upper right", fontsize=12, ncol=3,
              frameon=False, bbox_to_anchor=(1.0, -0.13))

    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf")
    plt.close(fig)
    buf.seek(0)
    for page in PdfReader(buf).pages:
        writer.add_page(page)

    # Append the matching operon depth plot, if present.
    if append_wdepth:
        op_plot = OPERON_PLOTS_DIR / f"{r.operon_id}_wdepth.pdf"
        if op_plot.exists():
            for page in PdfReader(str(op_plot)).pages:
                writer.add_page(page)


def visualize(report: pd.DataFrame, operons: pd.DataFrame,
              dels: pd.DataFrame, genes: pd.DataFrame, out_dir: Path) -> dict:
    """Write two PDFs per category: panel-only and panel+wdepth.
    Returns {category: n_panels}."""
    op_idx = operons.set_index("operon_id")
    out_dir.mkdir(parents=True, exist_ok=True)

    buckets: dict[str, list] = {k: [] for k in CATEGORY_FILES}
    for r in report.itertuples():
        b = _bucket(r.deletion_class)
        if b is not None:
            buckets[b].append(r)

    counts: dict[str, int] = {}
    for cat, rows in buckets.items():
        if not rows:
            counts[cat] = 0
            continue
        rows = sorted(rows, key=lambda r: (r.deletion_class, r.operon_id))

        base = CATEGORY_FILES[cat]
        out_plain  = out_dir / base
        out_wdepth = out_dir / base.replace(".pdf", "_wdepth.pdf")

        for out_pdf, with_depth in [(out_plain, False), (out_wdepth, True)]:
            writer = PdfWriter()
            for r in rows:
                _draw_panel(writer, r, op_idx.loc[r.operon_id], dels, genes,
                            append_wdepth=with_depth)
            with open(out_pdf, "wb") as fh:
                writer.write(fh)

        counts[cat] = len(rows)
    return counts


# ---------------------------------------------------------------- entry

def main() -> None:
    operons = pd.read_csv(OPERONS, sep="\t")
    dels    = load_deletions(DEL_BED)
    genes   = load_genes(GFF)

    report = analyse(operons, dels, genes)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(OUT_TSV, sep="\t", index=False)

    # console summary
    print(f"Operons classified: {len(report)}")
    print("\nAggregate deletion_class:")
    for k, v in report.deletion_class.value_counts().items():
        print(f"  {k:<48s} {v:5d}")

    print("\nPer-hit case totals (across all hits):")
    all_cases = report.per_hit_classes.fillna("").str.split(";").explode()
    all_cases = all_cases[all_cases != ""]
    for k, v in all_cases.value_counts().items():
        print(f"  {k:<32s} {v:5d}")

    counts = visualize(report, operons, dels, genes, OUT_DIR)
    print(f"\nWrote: {OUT_TSV}")
    for cat, n in counts.items():
        base = CATEGORY_FILES[cat]
        print(f"  {OUT_DIR / base}  ({n} panels)")
        print(f"  {OUT_DIR / base.replace('.pdf', '_wdepth.pdf')}  ({n} panels + wdepth)")


if __name__ == "__main__":
    main()
