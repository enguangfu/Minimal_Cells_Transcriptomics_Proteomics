#!/usr/bin/env python3
"""
07_operon_change.py

Operon-context reshaping by the JCVI-syn1 → syn3A genome reduction, tested
against ONT direct-RNA reads on syn3A using a single primitive: count of
ONT reads that fully span a given (gene_a, gene_b) pair on the operon strand.

Because ONT direct-RNA reads on syn3A are heavily fragmented (read length
distribution: 49-2858 nt, mean 383 nt) the de-novo operon segmentation in
`Syn3A_Operon/` over-splits real polycistrons into RNase-cleavage products
(525 syn3A "operons" covering 496 genes at MIN_READS=5). Rather than
trying to recover operon boundaries, this script uses the syn1 operon set
as the hypothesis library and tests each hypothesis directly against the
syn3A reads — a single ONT read only needs to span ONE adjacent gene pair
of the operon to provide evidence for that pair, so we sidestep the
fragmentation problem.

Three biological questions, one shared primitive:

  Q1 + Q2  Are syn1 operons preserved as transcription units in syn3A?
           For every syn1 operon, for every consecutive pair of RETAINED
           sense genes in transcription order, count ONT reads in syn3A
           that span the pair. Pairs are flagged as `deletion_collapsed`
           when the corresponding syn1 distance was much larger than the
           syn3A distance (i.e. a deletion sat between them in syn1).
           Per-pair and per-operon verdicts ("preserved" / "split").

  Q3       Have new gene proximities created new operons?
           For every consecutive same-strand gene pair in syn3A whose
           syn1 ancestors were NOT consecutive in syn1 (i.e. one or more
           syn1 genes between them were deleted), count ONT reads in
           syn3A that span the pair. Flag candidates as new transcription
           units.

Outputs (under Genome_Reduction/operon_change/):
    Q1Q2_pair_preservation.tsv    one row per consecutive retained pair
    Q1Q2_operon_preservation.tsv  one row per syn1 operon (aggregated)
    Q3_new_pair_candidates.tsv    one row per newly-adjacent syn3A pair
    operon_change_summary.txt     narrative numeric summary
"""

from __future__ import annotations

import re
from pathlib import Path
import pandas as pd
import pysam


# ---------- paths ----------
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

SYN1_GFF       = ROOT / "Genomes_Input" / "syn1.genes.gff3"
SYN3A_GFF      = ROOT / "Genomes_Input" / "syn3a_genome.gff3"
ONT_BAM        = ROOT / "ONT_Processing" / "syn3A.ONT.rep1.sorted.bam"
OPERON_TSV     = HERE / "aln" / "analysis" / "operon_deletion_classification.tsv"
SYN1_OPERONS_TSV = ROOT / "Syn1_Operon" / "operons.candidate_blocks.tsv"

OUT_DIR        = HERE / "operon_change"
OUT_Q1Q2_PAIRS = OUT_DIR / "Q1Q2_pair_preservation.tsv"
OUT_Q1Q2_OPS   = OUT_DIR / "Q1Q2_operon_preservation.tsv"
OUT_Q3         = OUT_DIR / "Q3_new_pair_candidates.tsv"
OUT_SUMMARY    = OUT_DIR / "operon_change_summary.txt"
OUT_COMPARE_DIR = OUT_DIR / "comparison_plots"

# Thresholds
MIN_SPAN_FOR_PRESERVED  = 1     # ≥ this many strict spanning reads = pair preserved (strict)
MIN_BRIDGE_FOR_PRESERVED = 1    # ≥ this many bridging reads = pair preserved (loose)
BRIDGE_MIN_OVERLAP_BP    = 50   # min bp the read must cover inside each CDS to count as bridging
# Operon-level verdict: an operon is called preserved only when EVERY consecutive
# retained-gene pair passes the metric. strict > loose > split.


# ---------- GFF loaders ----------

_PAT_LOCUS = re.compile(r"locus_tag=([^;]+)")
_PAT_NAME  = re.compile(r"(?:Name|gene)=([^;]+)")


def load_gff_genes(gff: Path) -> pd.DataFrame:
    rows = []
    with gff.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
                continue
            chrom, start1, end1, strand, attr = f[0], int(f[3]), int(f[4]), f[6], f[8]
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
    df = pd.DataFrame(rows).drop_duplicates(subset="locus_tag")
    return df


def locus_suffix(lt: str) -> str:
    return lt.rsplit("_", 1)[-1]


# ---------- the one shared primitive ----------

def count_spanning_reads(bam: pysam.AlignmentFile,
                         chrom: str, span_s: int, span_e: int,
                         strand: str) -> int:
    """STRICT: count primary ONT alignments on `strand` whose reference span
    fully encloses [span_s, span_e). Secondary / supplementary / unmapped
    excluded."""
    n = 0
    for read in bam.fetch(chrom, span_s, span_e):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        rs = "-" if read.is_reverse else "+"
        if rs != strand:
            continue
        if read.reference_start <= span_s and read.reference_end >= span_e:
            n += 1
    return n


def count_bridging_reads(bam: pysam.AlignmentFile,
                         chrom: str,
                         a_s: int, a_e: int, b_s: int, b_e: int,
                         strand: str,
                         min_overlap: int = BRIDGE_MIN_OVERLAP_BP) -> int:
    """LOOSE: count primary ONT alignments on `strand` that:
      * overlap CDS_A by >= min_overlap bp, AND
      * overlap CDS_B by >= min_overlap bp, AND
      * fully cover the intergenic region between A and B.

    Same biology as `count_spanning_reads` (one transcript across the
    A-intergenic-B boundary) but doesn't require the read to reach the
    far ends of either CDS — appropriate for short ONT direct-RNA reads.

    Inputs may have A, B in either order — function normalises by genomic
    coordinate (left gene, right gene)."""
    # Normalise so L = left gene by coordinate, R = right gene
    if a_s <= b_s:
        L_s, L_e, R_s, R_e = a_s, a_e, b_s, b_e
    else:
        L_s, L_e, R_s, R_e = b_s, b_e, a_s, a_e
    fetch_s = min(L_s, R_s)
    fetch_e = max(L_e, R_e)
    n = 0
    for read in bam.fetch(chrom, fetch_s, fetch_e):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        rs = "-" if read.is_reverse else "+"
        if rs != strand:
            continue
        r_s = read.reference_start
        r_e = read.reference_end
        # overlap with each CDS
        ov_L = max(0, min(r_e, L_e) - max(r_s, L_s))
        ov_R = max(0, min(r_e, R_e) - max(r_s, R_s))
        if ov_L < min_overlap or ov_R < min_overlap:
            continue
        # full coverage of the intergenic region between L and R
        # intergenic = [L_e, R_s) when L_e <= R_s; if they overlap, no intergenic check needed
        if L_e <= R_s:
            if r_s > L_e or r_e < R_s:
                continue
        n += 1
    return n


# ---------- Q1 + Q2  ----------

def analyse_syn1_operons(operons_tsv: Path,
                         syn1_genes: pd.DataFrame,
                         syn3a_genes: pd.DataFrame,
                         bam: pysam.AlignmentFile) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per syn1 operon, walk consecutive RETAINED genes in transcription order
    and test each pair for spanning reads in syn3A.

    Returns (pair_df, operon_df).
    """
    df = pd.read_csv(operons_tsv, sep="\t")
    # We need: operon_id, strand, sense_gene_locusNums (transcription order),
    # retained_genes (transcription order), gene_deletion_pattern.

    syn1_by_lt  = {r.locus_tag: r for r in syn1_genes.itertuples()}
    syn3a_by_lt = {r.locus_tag: r for r in syn3a_genes.itertuples()}

    pair_rows  = []
    op_rows    = []
    for op in df.itertuples():
        gdp = str(op.gene_deletion_pattern) if isinstance(op.gene_deletion_pattern, str) else ""
        # Skip operons with no retained sense gene to test (all_deleted, fully_deleted, no_sense_gene)
        retained_str = "" if not isinstance(op.retained_genes, str) else op.retained_genes
        sense_order  = "" if not isinstance(op.sense_gene_locusNums, str) else op.sense_gene_locusNums
        retained_set = {x.strip() for x in retained_str.split(",") if x.strip()}
        if len(retained_set) < 2:
            # need at least 2 retained genes to form a pair
            op_rows.append({
                "operon_id": op.operon_id,
                "strand":    op.strand,
                "gene_deletion_pattern": gdp,
                "n_retained_sense_genes": len(retained_set),
                "n_pairs":               0,
                "n_pairs_with_spanning": 0,
                "frac_pairs_spanned":    float("nan"),
                "operon_verdict":        "untestable_<2_retained",
            })
            continue

        # Retain syn1 transcription order, keeping only retained loci.
        ordered_locusNums = [x.strip() for x in sense_order.split(",") if x.strip()]
        ordered_retained = [ln for ln in ordered_locusNums if ln in retained_set]

        n_pairs = 0
        n_spanned = 0
        n_bridged = 0
        for i in range(len(ordered_retained) - 1):
            ln_a, ln_b = ordered_retained[i], ordered_retained[i + 1]
            syn1_a_lt = f"MMSYN1_{ln_a}"
            syn1_b_lt = f"MMSYN1_{ln_b}"
            syn3a_a_lt = f"JCVISYN3A_{ln_a}"
            syn3a_b_lt = f"JCVISYN3A_{ln_b}"

            ga1 = syn1_by_lt.get(syn1_a_lt)
            gb1 = syn1_by_lt.get(syn1_b_lt)
            ga3 = syn3a_by_lt.get(syn3a_a_lt)
            gb3 = syn3a_by_lt.get(syn3a_b_lt)
            if ga3 is None or gb3 is None:
                continue
            # syn3A pair span on the operon strand
            span_s = min(int(ga3.start0), int(gb3.start0))
            span_e = max(int(ga3.end0),   int(gb3.end0))
            chrom3 = str(ga3.chrom)
            strand = str(ga3.strand)

            n_pairs += 1
            n_span = count_spanning_reads(bam, chrom3, span_s, span_e, strand)
            n_bridge = count_bridging_reads(bam, chrom3,
                                            int(ga3.start0), int(ga3.end0),
                                            int(gb3.start0), int(gb3.end0),
                                            strand)
            spanned  = n_span   >= MIN_SPAN_FOR_PRESERVED
            bridged  = n_bridge >= MIN_BRIDGE_FOR_PRESERVED
            if spanned:
                n_spanned += 1
            if bridged:
                n_bridged += 1

            # Distance comparison: syn1 vs syn3A intergenic
            if ga1 is not None and gb1 is not None:
                syn1_dist = max(0, max(int(ga1.start0), int(gb1.start0))
                                  - min(int(ga1.end0),   int(gb1.end0)))
            else:
                syn1_dist = -1
            syn3a_dist = max(0, max(int(ga3.start0), int(gb3.start0))
                              - min(int(ga3.end0),   int(gb3.end0)))
            deletion_collapsed = (syn1_dist >= 0 and syn1_dist > syn3a_dist + 50)

            pair_rows.append({
                "operon_id":         op.operon_id,
                "strand":            strand,
                "gene_a_locusNum":   ln_a,
                "gene_b_locusNum":   ln_b,
                "syn3A_chrom":       chrom3,
                "syn3A_span_s":      span_s,
                "syn3A_span_e":      span_e,
                "syn1_intergenic_bp":  syn1_dist,
                "syn3A_intergenic_bp": syn3a_dist,
                "deletion_collapsed":  deletion_collapsed,
                "n_spanning_reads":  n_span,
                "n_bridging_reads":  n_bridge,
                "pair_preserved_strict": spanned,
                "pair_preserved_loose":  bridged,
            })

        frac_strict = (n_spanned / n_pairs) if n_pairs else float("nan")
        frac_loose  = (n_bridged / n_pairs) if n_pairs else float("nan")
        if n_pairs == 0:
            verdict = "untestable_no_pairs"
        elif n_spanned == n_pairs:
            verdict = "preserved_strict"
        elif n_bridged == n_pairs:
            verdict = "preserved_loose"
        else:
            verdict = "split"

        op_rows.append({
            "operon_id":              op.operon_id,
            "strand":                 op.strand,
            "gene_deletion_pattern":  gdp,
            "n_retained_sense_genes": len(retained_set),
            "n_pairs":                n_pairs,
            "n_pairs_spanning":       n_spanned,
            "n_pairs_bridging":       n_bridged,
            "frac_pairs_spanning":    round(frac_strict, 4) if n_pairs else float("nan"),
            "frac_pairs_bridging":    round(frac_loose,  4) if n_pairs else float("nan"),
            "operon_verdict":         verdict,
        })

    return pd.DataFrame(pair_rows), pd.DataFrame(op_rows)


# ---------- Q3 ----------

def analyse_new_pairs(syn1_genes: pd.DataFrame,
                      syn3a_genes: pd.DataFrame,
                      bam: pysam.AlignmentFile,
                      max_intergenic_bp: int | None = None) -> pd.DataFrame:
    """For each consecutive same-strand gene pair in syn3A whose syn1
    ancestors are NOT consecutive (≥1 syn1 gene of any strand was deleted
    between them), test whether the new pair is co-transcribed in syn3A.

    `max_intergenic_bp` optionally filters extremely-distant syn3A pairs;
    set to None to keep all novel-proximity pairs regardless of distance.
    The `syn3A_intergenic_bp` column is reported in the output so you can
    post-filter."""
    syn1_by_suf = {locus_suffix(r.locus_tag): r for r in syn1_genes.itertuples()}
    # Sort syn1 genes once for "between" lookup.
    syn1_sorted = syn1_genes.sort_values(["chrom", "start0"]).reset_index(drop=True)

    rows = []
    for strand in ("+", "-"):
        sub = syn3a_genes[syn3a_genes.strand == strand] \
                       .sort_values("start0").reset_index(drop=True)
        # Transcription order
        if strand == "-":
            sub = sub.iloc[::-1].reset_index(drop=True)
        for i in range(len(sub) - 1):
            ga3, gb3 = sub.iloc[i], sub.iloc[i + 1]
            ln_a, ln_b = locus_suffix(ga3.locus_tag), locus_suffix(gb3.locus_tag)
            chrom3 = str(ga3.chrom)
            # syn3A intergenic in transcription direction
            if strand == "+":
                inter_s = int(ga3.end0)
                inter_e = int(gb3.start0)
            else:
                inter_s = int(gb3.end0)
                inter_e = int(ga3.start0)
            syn3a_dist = max(0, inter_e - inter_s)
            if max_intergenic_bp is not None and syn3a_dist > max_intergenic_bp:
                continue

            # Are the syn1 ancestors annotated? If either has no MMSYN1 mate,
            # the pair is by definition new (newly-introduced gene).
            ga1 = syn1_by_suf.get(ln_a)
            gb1 = syn1_by_suf.get(ln_b)
            if ga1 is None or gb1 is None:
                novel_reason = "no_syn1_ancestor"
                n_intervening = -1
            else:
                # Count syn1 genes (any strand) strictly between the two ancestors.
                lo = min(int(ga1.end0), int(gb1.end0))
                hi = max(int(ga1.start0), int(gb1.start0))
                if hi <= lo:
                    n_intervening = 0
                else:
                    inter = syn1_sorted[
                        (syn1_sorted.chrom == ga1.chrom) &
                        (syn1_sorted.start0 >= lo) &
                        (syn1_sorted.end0   <= hi) &
                        (syn1_sorted.locus_tag != ga1.locus_tag) &
                        (syn1_sorted.locus_tag != gb1.locus_tag)
                    ]
                    n_intervening = int(len(inter))
                if n_intervening == 0:
                    continue  # pair was already adjacent in syn1 — not novel
                novel_reason = "deletion_brought_together"

            span_s = min(int(ga3.start0), int(gb3.start0))
            span_e = max(int(ga3.end0),   int(gb3.end0))
            n_span = count_spanning_reads(bam, chrom3, span_s, span_e, strand)
            n_bridge = count_bridging_reads(bam, chrom3,
                                            int(ga3.start0), int(ga3.end0),
                                            int(gb3.start0), int(gb3.end0),
                                            strand)

            rows.append({
                "syn3A_gene_a":      ga3.locus_tag,
                "syn3A_gene_b":      gb3.locus_tag,
                "strand":            strand,
                "syn3A_chrom":       chrom3,
                "syn3A_intergenic_bp": syn3a_dist,
                "n_intervening_syn1_genes": n_intervening,
                "novel_reason":      novel_reason,
                "syn3A_span_s":      span_s,
                "syn3A_span_e":      span_e,
                "n_spanning_reads":  n_span,
                "n_bridging_reads":  n_bridge,
                "co_transcribed_strict": n_span   >= MIN_SPAN_FOR_PRESERVED,
                "co_transcribed_loose":  n_bridge >= MIN_BRIDGE_FOR_PRESERVED,
            })
    return pd.DataFrame(rows)


# ---------- comparison plots (calls Operon_Comparison_Syn1_Syn3A.plot_one_operon_comparison) ----------

def plot_comparison_for_flagged_categories(op_df: pd.DataFrame,
                                           out_dir: Path) -> dict:
    """For each of the four flagged categories, render side-by-side
    syn1-vs-syn3A operon comparison PDFs into a category subfolder.

    Categories
      split_intact     : gene_deletion_pattern == 'intact' AND operon_verdict == 'split'
      leading_deleted  : gene_deletion_pattern == 'leading_deleted'
      intra_deleted    : gene_deletion_pattern == 'intra_deleted'
      lagging_deleted  : gene_deletion_pattern == 'lagging_deleted'
    """
    import sys
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    import Operon_Comparison_Syn1_Syn3A as comp  # noqa: E402

    syn1_ops = pd.read_csv(SYN1_OPERONS_TSV, sep="\t")
    syn1_ops_by_id = {r.operon_id: r for r in syn1_ops.itertuples()}

    cats = {
        "split_intact":    op_df[(op_df.gene_deletion_pattern == "intact") &
                                 (op_df.operon_verdict == "split")],
        "leading_deleted": op_df[op_df.gene_deletion_pattern == "leading_deleted"],
        "intra_deleted":   op_df[op_df.gene_deletion_pattern == "intra_deleted"],
        "lagging_deleted": op_df[op_df.gene_deletion_pattern == "lagging_deleted"],
    }
    counts: dict = {}
    for cat, sub in cats.items():
        cat_dir = out_dir / cat
        cat_dir.mkdir(parents=True, exist_ok=True)
        ok = 0
        for op_id in sub.operon_id.tolist():
            op_row = syn1_ops_by_id.get(op_id)
            if op_row is None:
                continue
            # convert namedtuple to dict for the plotter (uses .get with defaults)
            op_dict = op_row._asdict() if hasattr(op_row, "_asdict") else dict(op_row.__dict__)
            save = cat_dir / f"{op_id}.pdf"
            try:
                wrote = comp.plot_one_operon_comparison(op_dict, str(save), PLOT_DEPTH=True)
                if wrote:
                    ok += 1
            except Exception as e:
                print(f"  WARN: {op_id} in {cat} failed to plot: {e}")
        counts[cat] = ok
        print(f"  {cat:<18s}: {ok}/{len(sub)} plotted -> {cat_dir}")
    return counts


# ---------- main ----------

def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading GFFs and BAM ...")
    syn1_genes  = load_gff_genes(SYN1_GFF)
    syn3a_genes = load_gff_genes(SYN3A_GFF)
    bam = pysam.AlignmentFile(str(ONT_BAM), "rb")

    print(f"  syn1 genes : {len(syn1_genes)}")
    print(f"  syn3A loci : {len(syn3a_genes)}")

    print("\nQ1+Q2: Walking syn1 operons, testing consecutive retained pairs ...")
    pair_df, op_df = analyse_syn1_operons(OPERON_TSV, syn1_genes, syn3a_genes, bam)
    pair_df.to_csv(OUT_Q1Q2_PAIRS, sep="\t", index=False)
    op_df.to_csv(OUT_Q1Q2_OPS, sep="\t", index=False)

    print(f"\nQ3: Scanning newly-adjacent same-strand syn3A pairs ...")
    q3_df = analyse_new_pairs(syn1_genes, syn3a_genes, bam)
    q3_df.to_csv(OUT_Q3, sep="\t", index=False)
    bam.close()

    # ---- narrative summary ----
    lines = []
    lines.append("=" * 64)
    lines.append("OPERON-CONTEXT RESHAPING (syn1 -> syn3A) — ONT spanning evidence")
    lines.append("=" * 64)

    # ---- Per-pattern gene-attribution counts (corrected from earlier) ----
    lines.append("")
    lines.append("Gene-attribution counts per syn1 gene_deletion_pattern")
    lines.append("(each cell is the sum of gene-attribution slots across operons of that pattern)")
    lines.append("-" * 64)
    pat_df = pd.read_csv(OPERON_TSV, sep="\t")
    def _ncsv(s: str | float) -> int:
        if not isinstance(s, str) or not s.strip():
            return 0
        return len([x for x in s.split(",") if x.strip()])
    def _npartial(s: str | float) -> int:
        # entries are 'NUM(frac, Xbps),NUM(frac, Xbps),...' — count '(' to avoid
        # splitting on the comma inside the parentheses.
        if not isinstance(s, str) or not s.strip():
            return 0
        return s.count("(")
    pat_df["n_total_sense"] = pat_df["sense_gene_locusNums"].fillna("").map(_ncsv)
    pat_df["n_retained"]    = pat_df["retained_genes"].fillna("").map(_ncsv)
    pat_df["n_partial"]     = pat_df["partially_deleted_genes"].fillna("").map(_npartial)
    pat_df["n_fully"]       = pat_df["fully_deleted_genes"].fillna("").map(_ncsv)
    pat_agg = pat_df.groupby("gene_deletion_pattern").agg(
        n_operons        = ("operon_id",       "size"),
        sense_genes_syn1 = ("n_total_sense",   "sum"),
        fully_retained   = ("n_retained",      "sum"),
        partially_deleted= ("n_partial",       "sum"),
        fully_deleted    = ("n_fully",         "sum"),
        mean_n_sense     = ("n_total_sense",   "mean"),
        mean_n_retained  = ("n_retained",      "mean"),
    ).round(2)
    pat_order = ["intact", "leading_deleted", "intra_deleted",
                 "lagging_deleted", "all_deleted", "fully_deleted",
                 "no_sense_gene"]
    pat_agg = pat_agg.reindex([k for k in pat_order if k in pat_agg.index])
    lines.append(pat_agg.to_string())
    tot_sense   = int(pat_df["n_total_sense"].sum())
    tot_retained = int(pat_df["n_retained"].sum())
    tot_partial  = int(pat_df["n_partial"].sum())
    tot_fully    = int(pat_df["n_fully"].sum())
    lines.append("")
    lines.append(f"  Totals: {tot_sense} sense-attributions; "
                 f"{tot_retained} fully retained, "
                 f"{tot_partial} partially deleted, "
                 f"{tot_fully} fully deleted.")

    # Q1+Q2 — by gene_deletion_pattern
    lines.append("")
    lines.append("Q1+Q2  Syn1 operons preserved as transcription units in syn3A")
    lines.append("-" * 64)
    lines.append(f"  strict (spanning): read covers entire CDS_A + intergenic + CDS_B")
    lines.append(f"  loose (bridging) : read enters each CDS by >={BRIDGE_MIN_OVERLAP_BP} bp and covers the intergenic region")
    lines.append(f"  operon verdict   : preserved_strict if EVERY pair has a spanning read;")
    lines.append(f"                     preserved_loose if EVERY pair has a bridging read (but not strict);")
    lines.append(f"                     split otherwise.")

    testable = op_df[op_df.n_pairs > 0]
    untestable = op_df[op_df.n_pairs == 0]
    lines.append(f"  testable operons (>=2 retained genes)        : {len(testable)}")
    lines.append(f"  untestable operons (<2 retained or no pairs) : {len(untestable)}")
    if len(testable):
        n_strict = int((testable.operon_verdict == "preserved_strict").sum())
        n_loose  = int((testable.operon_verdict == "preserved_loose").sum())
        n_split  = int((testable.operon_verdict == "split").sum())
        lines.append(f"  preserved_strict : {n_strict:5d}  ({n_strict/len(testable):.1%})  (full-span reads)")
        lines.append(f"  preserved_loose  : {n_loose:5d}  ({n_loose/len(testable):.1%})  (only bridging reads pass)")
        lines.append(f"  split            : {n_split:5d}  ({n_split/len(testable):.1%})")

    # By gene_deletion_pattern
    if len(testable):
        lines.append("")
        lines.append("  Breakdown by syn1 gene_deletion_pattern:")
        tab = (testable
               .groupby("gene_deletion_pattern")["operon_verdict"]
               .value_counts()
               .unstack(fill_value=0))
        lines.append(tab.to_string())

    # Pair-level summary, especially deletion-collapsed pairs (the Q2 core test)
    if len(pair_df):
        lines.append("")
        lines.append("  Pair-level totals:")
        lines.append(f"    total pairs tested          : {len(pair_df)}")
        n_strict = int(pair_df.pair_preserved_strict.sum())
        n_loose  = int(pair_df.pair_preserved_loose.sum())
        lines.append(f"    spanning (strict)           : {n_strict:5d}  ({n_strict/len(pair_df):.1%})")
        lines.append(f"    bridging (loose)            : {n_loose:5d}  ({n_loose/len(pair_df):.1%})")
        n_dc = int(pair_df.deletion_collapsed.sum())
        lines.append(f"    deletion_collapsed pairs    : {n_dc}  (pairs that lost intervening syn1 distance)")
        if n_dc:
            sub = pair_df[pair_df.deletion_collapsed]
            n_dc_strict = int(sub.pair_preserved_strict.sum())
            n_dc_loose  = int(sub.pair_preserved_loose.sum())
            lines.append(f"      spanning preserved        : {n_dc_strict}  ({n_dc_strict/n_dc:.1%})")
            lines.append(f"      bridging preserved        : {n_dc_loose}  ({n_dc_loose/n_dc:.1%})")

    # Q3
    lines.append("")
    lines.append("Q3  Newly-adjacent same-strand syn3A pairs (deletion-induced)")
    lines.append("-" * 64)
    if len(q3_df) == 0:
        lines.append("  No newly-adjacent same-strand pairs found.")
    else:
        lines.append(f"  pairs evaluated                  : {len(q3_df)}")
        n_co_strict = int(q3_df.co_transcribed_strict.sum())
        n_co_loose  = int(q3_df.co_transcribed_loose.sum())
        lines.append(f"  co-transcribed (spanning, strict): {n_co_strict}  ({n_co_strict/len(q3_df):.1%})")
        lines.append(f"  co-transcribed (bridging, loose) : {n_co_loose}  ({n_co_loose/len(q3_df):.1%})")
        co = q3_df[q3_df.co_transcribed_loose]
        if len(co):
            d = co.syn3A_intergenic_bp
            lines.append(f"  intergenic_bp of co-transcribed pairs (loose): "
                         f"median={int(d.median())}  mean={d.mean():.1f}  "
                         f"p90={int(d.quantile(0.9))}  max={int(d.max())}")
        lines.append("  Top 10 newly-co-transcribed pairs by bridging-read count:")
        top = q3_df.sort_values("n_bridging_reads", ascending=False).head(10)
        for r in top.itertuples():
            lines.append(f"    {r.syn3A_gene_a} -> {r.syn3A_gene_b}  "
                         f"({r.strand}, {r.syn3A_intergenic_bp} bp apart, "
                         f"{r.n_intervening_syn1_genes} syn1 genes lost; "
                         f"n_span={r.n_spanning_reads}, n_bridge={r.n_bridging_reads})")

    lines.append("")
    lines.append(f"Wrote: {OUT_Q1Q2_PAIRS}")
    lines.append(f"Wrote: {OUT_Q1Q2_OPS}")
    lines.append(f"Wrote: {OUT_Q3}")
    lines.append(f"Wrote: {OUT_SUMMARY}")

    # ---- side-by-side syn1 vs syn3A comparison plots for the 4 flagged classes ----
    print("\nRendering syn1-vs-syn3A comparison plots for flagged categories ...")
    counts = plot_comparison_for_flagged_categories(op_df, OUT_COMPARE_DIR)
    lines.append("")
    lines.append("Side-by-side comparison plots written:")
    for cat, n in counts.items():
        lines.append(f"  {cat:<18s} : {n} PDFs  ->  {OUT_COMPARE_DIR / cat}")

    text = "\n".join(lines) + "\n"
    OUT_SUMMARY.write_text(text)
    print(text)


if __name__ == "__main__":
    main()
