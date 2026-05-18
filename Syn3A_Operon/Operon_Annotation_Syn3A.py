#!/usr/bin/env python3
"""
Operon_Annotation_Syn3A.py

For Syn3A operons derived from ONT isoform clusters in
`operons.candidate_blocks.tsv`, quantify how well isoform-defined operons
contain the full ORFs of their sense genes.

ONT cDNA reads are heavily 5'/3'-degraded compared with PacBio FLNC, so even
when an operon is correctly identified at the gene level, its longest isoform
may not reach both ORF ends. This script restricts attention to operons whose
`segmentation_type` is `isoform_operon` or `isoform_operon_merged` (i.e. those
defined by isoform evidence rather than by gene-rescue fallbacks) and asks:

    * What fraction of those operons cover every sense gene's ORF in full?
    * Among the rest, how much sequence is missing at the 5' and 3' boundaries?

Outputs (in this folder):
    operon_orf_coverage.tsv      per-operon classification + per-gene shortfalls
    operon_orf_coverage.txt      short numeric summary
"""

from __future__ import annotations

import re
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

OPERONS_TSV = HERE / "operons.candidate_blocks.tsv"
GFF_PATH    = ROOT / "Genomes_Input" / "syn3a_genome.gff3"
OUT_TSV     = HERE / "operon_orf_coverage.tsv"
OUT_TXT     = HERE / "operon_orf_coverage.txt"
OUT_ATTR_TSV = HERE / "gene_operon_attribution.tsv"
OUT_PDF_TMPL = HERE / "operon_multiplicity_{n}.pdf"

TARGET_TYPES = {"isoform_operon", "isoform_operon_merged"}


# ----------------------------------------------------------- GFF gene table

_PAT_LOCUS = re.compile(r"locus_tag=([^;]+)")
_PAT_NAME  = re.compile(r"(?:Name|gene)=([^;]+)")


def load_genes(gff: Path) -> pd.DataFrame:
    rows = []
    with gff.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
                continue
            chrom, start1, end1, strand, attr = (
                f[0], int(f[3]), int(f[4]), f[6], f[8])
            m = _PAT_LOCUS.search(attr)
            if not m:
                continue
            n = _PAT_NAME.search(attr)
            rows.append({
                "locus_tag": m.group(1),
                "name":      n.group(1) if n else "",
                "chrom":     chrom,
                "strand":    strand,
                "start0":    start1 - 1,
                "end0":      end1,
            })
    return (pd.DataFrame(rows)
            .drop_duplicates(subset="locus_tag")
            .set_index("locus_tag"))


# ----------------------------------------------------------- coverage logic

def analyse_operon(op, genes: pd.DataFrame) -> dict:
    """For one operon row, return per-operon coverage stats and a per-gene
    shortfall string."""
    op_s, op_e = int(op.start0), int(op.end0)
    raw = op.sense_gene_loci if isinstance(op.sense_gene_loci, str) else ""
    loci = [t.strip() for t in raw.split(",") if t.strip()]

    n_sense = 0
    n_fully_covered = 0
    shortfalls: list[str] = []
    worst_5p = 0
    worst_3p = 0
    total_5p = 0
    total_3p = 0
    max_missing_frac = 0.0

    for lt in loci:
        if lt not in genes.index:
            continue
        g = genes.loc[lt]
        if str(g.strand) != str(op.strand):
            continue  # only sense genes contribute
        n_sense += 1
        g_s, g_e = int(g.start0), int(g.end0)
        miss_left  = max(0, op_s - g_s)        # ORF extends before operon
        miss_right = max(0, g_e - op_e)        # ORF extends past operon
        if op.strand == "+":
            miss_5p, miss_3p = miss_left, miss_right
        else:
            miss_5p, miss_3p = miss_right, miss_left
        missing = miss_left + miss_right
        if missing == 0:
            n_fully_covered += 1
        else:
            shortfalls.append(f"{lt}(5p:{miss_5p},3p:{miss_3p})")
            total_5p += miss_5p
            total_3p += miss_3p
            worst_5p = max(worst_5p, miss_5p)
            worst_3p = max(worst_3p, miss_3p)
            orf_len = g_e - g_s
            if orf_len > 0:
                max_missing_frac = max(max_missing_frac, missing / orf_len)

    if n_sense == 0:
        verdict = "no_sense_gene"
    elif n_fully_covered == n_sense:
        verdict = "fully_covers_ORFs"
    else:
        verdict = "truncates_ORFs"

    return {
        "n_sense_genes":         n_sense,
        "n_genes_fully_covered": n_fully_covered,
        "orf_coverage_verdict":  verdict,
        "n_genes_truncated":     len(shortfalls),
        "worst_5p_missing_bps":  worst_5p,
        "worst_3p_missing_bps":  worst_3p,
        "total_5p_missing_bps":  total_5p,
        "total_3p_missing_bps":  total_3p,
        "max_missing_orf_frac":  round(max_missing_frac, 4),
        "gene_shortfalls":       ";".join(shortfalls),
    }


# ----------------------------------------------------------- gene attribution

SEG_ORDER = [
    "isoform_operon",
    "isoform_operon_merged",
    "isoform_gene_combined",
    "rescue_single",
    "rescue_multiple",
]


def _split(field) -> list[str]:
    if not isinstance(field, str) or not field:
        return []
    return [t.strip() for t in field.split(",") if t.strip()]


def gene_attribution(operons: pd.DataFrame, genes: pd.DataFrame,
                     out_tsv: Path) -> list[str]:
    """Per gene, find which operon(s) attribute it (sense or antisense),
    grouped by segmentation_type. Write a per-gene table and return a
    list of summary lines for the .txt report."""
    n_loci = len(genes)
    # Per-gene records: which operons (sense / antisense) include this locus
    sense_op: dict[str, list[str]]      = {lt: [] for lt in genes.index}
    antisense_op: dict[str, list[str]]  = {lt: [] for lt in genes.index}
    sense_seg: dict[str, list[str]]     = {lt: [] for lt in genes.index}
    antisense_seg: dict[str, list[str]] = {lt: [] for lt in genes.index}

    for op in operons.itertuples():
        for lt in _split(op.sense_gene_loci):
            if lt in sense_op:
                sense_op[lt].append(op.operon_id)
                sense_seg[lt].append(op.segmentation_type)
        for lt in _split(op.antisense_gene_loci):
            if lt in antisense_op:
                antisense_op[lt].append(op.operon_id)
                antisense_seg[lt].append(op.segmentation_type)

    rows = []
    for lt in genes.index:
        g = genes.loc[lt]
        sense_types = sense_seg[lt]
        anti_types  = antisense_seg[lt]
        if sense_types:
            role = "sense"
        elif anti_types:
            role = "antisense_only"
        else:
            role = "uncovered"
        # multiset of segmentation types across the gene's sense operons
        # (NOT deduplicated — so e.g. a gene in two isoform_operon operons is
        # 'isoform_operon|isoform_operon').
        sense_combo = "|".join(sorted(sense_types)) if sense_types else ""
        rows.append({
            "locus_tag":      lt,
            "name":           g["name"],
            "chrom":          g.chrom,
            "strand":         g.strand,
            "start0":         int(g.start0),
            "end0":           int(g.end0),
            "role":           role,
            "n_sense_operons":     len(sense_op[lt]),
            "n_antisense_operons": len(antisense_op[lt]),
            "sense_segmentation_combo":      sense_combo,
            "sense_operon_ids":      ",".join(sense_op[lt]),
            "antisense_operon_ids":  ",".join(antisense_op[lt]),
            "sense_segmentation_types":      ",".join(sense_types),
            "antisense_segmentation_types":  ",".join(anti_types),
        })
    attr_df = pd.DataFrame(rows)
    attr_df.to_csv(out_tsv, sep="\t", index=False)

    # ---- assemble narrative summary ----
    lines: list[str] = []
    lines.append("=" * 64)
    lines.append("GENE ATTRIBUTION BY OPERON SEGMENTATION TYPE")
    lines.append("=" * 64)
    lines.append(f"Total syn3A loci (gene + pseudogene) : {n_loci}")
    lines.append("")

    # 1. Operon-level: counts + sense / antisense gene volume per segmentation type
    lines.append("Per-segmentation-type, summed across operons:")
    header = (f"  {'segmentation_type':<24s} {'n_op':>5s} "
              f"{'sense_genes':>11s} {'antisense_genes':>16s}")
    lines.append(header)
    lines.append("  " + "-" * (len(header) - 2))
    for st in SEG_ORDER + sorted(set(operons.segmentation_type) - set(SEG_ORDER)):
        sub = operons[operons.segmentation_type == st]
        if sub.empty:
            continue
        n_op = len(sub)
        n_s  = sum(len(_split(v)) for v in sub.sense_gene_loci)
        n_a  = sum(len(_split(v)) for v in sub.antisense_gene_loci) \
               if "antisense_gene_loci" in sub.columns else 0
        lines.append(f"  {st:<24s} {n_op:>5d} {n_s:>11d} {n_a:>16d}")

    # 2. Gene-level role tally (each gene counted once).
    lines.append("")
    lines.append("Gene role tally (each locus counted once):")
    role_tab = attr_df.role.value_counts().reindex(
        ["sense", "antisense_only", "uncovered"], fill_value=0)
    for k, v in role_tab.items():
        pct = 100.0 * v / n_loci
        lines.append(f"  {k:<16s} : {int(v):5d}  ({pct:5.1f}%)")

    # 3. Multiplicity: how many sense operons each gene sits in.
    sense_only = attr_df[attr_df.role == "sense"]
    n_sense = len(sense_only)
    lines.append("")
    lines.append("Sense-operon multiplicity (#sense operons attributing the same gene):")
    mult_tab = (sense_only.n_sense_operons
                .value_counts().sort_index())
    for k, v in mult_tab.items():
        pct = 100.0 * int(v) / max(1, n_sense)
        lines.append(f"  {int(k):>2d} operon(s) : {int(v):5d}  ({pct:5.1f}%)")

    # 4. Per-multiplicity breakdown: for each #sense-operons bucket, list every
    # combination of segmentation types (as an ordered multiset) and the
    # number of genes with that exact combination.
    lines.append("")
    lines.append("Per-multiplicity sense_segmentation combos:")
    for n_ops in sorted(sense_only.n_sense_operons.unique()):
        sub = sense_only[sense_only.n_sense_operons == n_ops]
        lines.append(f"  {n_ops} sense operon(s) — {len(sub)} genes:")
        combo_tab = sub.sense_segmentation_combo.value_counts()
        for combo, n in combo_tab.items():
            pct_of_group = 100.0 * int(n) / len(sub)
            lines.append(f"    {combo:<60s} : {int(n):5d}  ({pct_of_group:5.1f}%)")

    # 5. Antisense-operon multiplicity sanity check
    n_multi_anti = int((attr_df.n_antisense_operons >  1).sum())
    lines.append("")
    lines.append(f"Genes in >1 antisense operon : {n_multi_anti}")

    lines.append("")
    lines.append(f"Wrote: {out_tsv}")
    return lines


# ----------------------------------------------------------- visualization

PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
           "#8c564b", "#e377c2"]


def _plot_one_focal(ax, focal_locus: str, focal_gene, focal_ops: pd.DataFrame,
                    genes: pd.DataFrame) -> None:
    """Render a single panel: gene track at y=0, operons stacked above."""
    n = len(focal_ops)
    op_y = list(range(1, n + 1))

    view_s = int(focal_ops.start0.min()) - 200
    view_e = int(focal_ops.end0.max())   + 200

    # ---- gene track at y=0 ----
    region = genes[(genes.chrom == focal_gene.chrom) &
                   (genes.end0   > view_s) &
                   (genes.start0 < view_e)]
    GENE_H = 0.30
    for lt, g in region.iterrows():
        is_focal = (lt == focal_locus)
        color = "#C44E52" if is_focal else "#888888"
        ax.broken_barh([(int(g.start0), int(g.end0 - g.start0))],
                       (-GENE_H / 2, GENE_H),
                       facecolors=color, alpha=0.85, zorder=2)
        label = g["name"] or lt
        if is_focal:
            label = f"{label}  ({lt})"
        ax.text((int(g.start0) + int(g.end0)) / 2, -GENE_H / 2 - 0.05,
                label, fontsize=6, ha="center", va="top",
                color="#C44E52" if is_focal else "#444444",
                fontweight="bold" if is_focal else "normal")
        # tiny strand caret on each gene
        marker = "▶" if g.strand == "+" else "◀"
        ax.text((int(g.start0) + int(g.end0)) / 2, 0,
                marker, fontsize=7, ha="center", va="center", color="white")

    # ---- operon tracks ----
    OP_H = 0.40
    for i, (_, op) in enumerate(focal_ops.iterrows()):
        y = op_y[i]
        color = PALETTE[i % len(PALETTE)]
        op_s, op_e = int(op.start0), int(op.end0)
        ax.broken_barh([(op_s, op_e - op_s)], (y - OP_H / 2, OP_H),
                       facecolors=color, alpha=0.70,
                       edgecolor="black", linewidth=0.5, zorder=2)
        ax.text((op_s + op_e) / 2, y,
                f"{op.operon_id}  {op.segmentation_type}  ({op_e - op_s} bp)",
                fontsize=7, ha="center", va="center",
                color="white", fontweight="bold")

    ax.set_xlim(view_s, view_e)
    ax.set_ylim(-0.9, n + 0.8)
    ax.set_yticks([0] + op_y)
    ax.set_yticklabels(["Genes"] + [f"Op {i+1}" for i in range(n)], fontsize=7)
    ax.set_xlabel("Genomic coordinate (bp)", fontsize=8)
    ax.tick_params(axis="x", labelsize=7)

    focal_name = focal_gene["name"] or focal_locus
    ax.set_title(
        f"{focal_locus} ({focal_name})  strand {focal_gene.strand}  —  "
        f"{n} sense operons; ORF {int(focal_gene.start0)}-{int(focal_gene.end0)}",
        fontsize=8, pad=3,
    )


def plot_multiplicity_panels(operons: pd.DataFrame, genes: pd.DataFrame,
                             attr_df: pd.DataFrame) -> dict:
    """For each multiplicity n in 2..5, write one PDF with one panel per gene
    showing the operons that attribute it."""
    op_idx = operons.set_index("operon_id")
    written: dict = {}
    for n in (2, 3, 4, 5):
        sub = attr_df[(attr_df.role == "sense") & (attr_df.n_sense_operons == n)]
        if sub.empty:
            written[n] = 0
            continue
        # sort by chrom / start for stable order
        sub = sub.sort_values(["chrom", "start0"]).reset_index(drop=True)
        out_pdf = Path(str(OUT_PDF_TMPL).format(n=n))
        with PdfPages(out_pdf) as pdf:
            for r in sub.itertuples():
                op_ids = [x for x in r.sense_operon_ids.split(",") if x]
                focal_ops = op_idx.loc[op_ids].sort_values("start0").reset_index()
                focal_gene = genes.loc[r.locus_tag]
                fig, ax = plt.subplots(figsize=(11, 1.6 + 0.55 * n))
                _plot_one_focal(ax, r.locus_tag, focal_gene, focal_ops, genes)
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
        written[n] = len(sub)
        print(f"  wrote {out_pdf}  ({len(sub)} panels)")
    return written


# ----------------------------------------------------------- main

def main() -> None:
    operons = pd.read_csv(OPERONS_TSV, sep="\t")
    genes   = load_genes(GFF_PATH)

    sub = operons[operons.segmentation_type.isin(TARGET_TYPES)] \
                 .copy().reset_index(drop=True)

    rows = []
    for op in sub.itertuples():
        stats = analyse_operon(op, genes)
        rows.append({
            "operon_id":         op.operon_id,
            "chrom":             op.chrom,
            "strand":            op.strand,
            "start0":            int(op.start0),
            "end0":              int(op.end0),
            "length":            int(op.end0) - int(op.start0),
            "segmentation_type": op.segmentation_type,
            "sense_gene_loci":   op.sense_gene_loci,
            **stats,
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT_TSV, sep="\t", index=False)

    # --------- summary ----------
    lines = []
    lines.append("=" * 64)
    lines.append("ORF-COVERAGE SUMMARY  (isoform_operon + isoform_operon_merged)")
    lines.append("=" * 64)
    lines.append(f"isoform-defined operons examined : {len(out)}")
    if len(out) == 0:
        OUT_TXT.write_text("\n".join(lines) + "\n")
        print("\n".join(lines))
        return

    by_type = out.segmentation_type.value_counts().to_dict()
    for t in ("isoform_operon", "isoform_operon_merged"):
        lines.append(f"  {t:<24s} : {by_type.get(t, 0)}")

    counts = out.orf_coverage_verdict.value_counts()
    n_full  = int(counts.get("fully_covers_ORFs", 0))
    n_trunc = int(counts.get("truncates_ORFs", 0))
    n_none  = int(counts.get("no_sense_gene", 0))
    total = len(out)
    lines.append("")
    lines.append("Per-operon verdict:")
    lines.append(f"  fully_covers_ORFs : {n_full:5d}  ({n_full/total:.1%})")
    lines.append(f"  truncates_ORFs    : {n_trunc:5d}  ({n_trunc/total:.1%})")
    lines.append(f"  no_sense_gene     : {n_none:5d}")

    total_genes = int(out.n_sense_genes.sum())
    full_genes  = int(out.n_genes_fully_covered.sum())
    trunc_genes = total_genes - full_genes
    if total_genes > 0:
        lines.append("")
        lines.append("Per-sense-gene verdict (across these operons):")
        lines.append(f"  total sense genes        : {total_genes}")
        lines.append(f"  ORF fully inside operon  : {full_genes}  ({full_genes/total_genes:.1%})")
        lines.append(f"  ORF partially truncated  : {trunc_genes}  ({trunc_genes/total_genes:.1%})")

    if n_trunc > 0:
        trunc_rows = out[out.orf_coverage_verdict == "truncates_ORFs"]
        m5 = trunc_rows.worst_5p_missing_bps
        m3 = trunc_rows.worst_3p_missing_bps
        lines.append("")
        lines.append("Among truncating operons, worst-gene bp missing:")
        lines.append("  5' side  median={:>5}  mean={:>6.1f}  p90={:>5}  max={:>5}".format(
            int(m5.median()), float(m5.mean()), int(m5.quantile(0.9)), int(m5.max())))
        lines.append("  3' side  median={:>5}  mean={:>6.1f}  p90={:>5}  max={:>5}".format(
            int(m3.median()), float(m3.mean()), int(m3.quantile(0.9)), int(m3.max())))

        thresholds = [10, 30, 100, 300]
        lines.append("")
        lines.append("Operons with worst-gene shortfall exceeding threshold:")
        lines.append(f"  {'bp':>5}  {'5p-side':>9}  {'3p-side':>9}")
        for t in thresholds:
            n5 = int((m5 >= t).sum())
            n3 = int((m3 >= t).sum())
            lines.append(f"  {t:>5}  {n5:>9}  {n3:>9}")

    lines.append("")
    lines.append(f"Wrote: {OUT_TSV}")

    # Gene attribution across all operon types (uses the full operons table).
    lines.append("")
    lines.extend(gene_attribution(operons, genes, OUT_ATTR_TSV))

    # Multiplicity diagrams (one PDF per multiplicity level).
    attr_df = pd.read_csv(OUT_ATTR_TSV, sep="\t")
    lines.append("")
    lines.append("Multi-operon diagrams (one panel per focal gene):")
    written = plot_multiplicity_panels(operons, genes, attr_df)
    for n, k in written.items():
        path = str(OUT_PDF_TMPL).format(n=n)
        lines.append(f"  multiplicity {n}: {k} panels -> {path}")

    lines.append("")
    lines.append(f"Wrote: {OUT_TXT}")
    text = "\n".join(lines) + "\n"
    OUT_TXT.write_text(text)
    print(text)


if __name__ == "__main__":
    main()
