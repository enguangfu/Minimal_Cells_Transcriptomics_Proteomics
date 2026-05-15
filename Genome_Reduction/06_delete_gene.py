#!/usr/bin/env python3
"""
06_delete_gene.py

For every syn1 gene retained in syn3A (no bp of its ORF overlaps any deletion),
find the same-strand transcription-direction upstream and downstream neighbor
in *both* syn1 and syn3A, and report how many contiguous unaltered base pairs
of context sit upstream and downstream in the syn1 frame.

Both genomes are treated as circular: neighbor lookup and unaltered-bp scans
wrap around the chromosome end.

Outputs (in aln/analysis/):
    retained_gene_context.tsv           one row per retained syn1 gene
    retained_gene_context_summary.txt   short narrative + retained-set check
"""

from __future__ import annotations

import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT     = Path(__file__).resolve().parent.parent
SYN1_GFF = ROOT / "Genomes_Input" / "syn1.genes.gff3"
SYN3A_GFF = ROOT / "Genomes_Input" / "syn3a_genome.gff3"
DEL_BED  = Path(__file__).resolve().parent / "aln" / "raw" / "syn1_deleted_regions.bed"
OUT_DIR  = Path(__file__).resolve().parent / "aln" / "analysis"
OPERON_TSV = OUT_DIR / "operon_deletion_classification.tsv"
RNA_CSV    = Path(__file__).resolve().parent / "syn1_vs_syn3a_RNA_protein.csv"
OUT_TSV  = OUT_DIR / "retained_gene_context.tsv"
OUT_TXT  = OUT_DIR / "retained_gene_context_summary.txt"
OUT_PDF  = OUT_DIR / "operon_change_TPM_fold_boxplot.pdf"

SYN1_LEN  = 1_078_809   # CP002027.1, circular
SYN3A_LEN =   543_379   # CP016816.2, circular


# ---------------------------------------------------------------- loaders

_PAT_LOCUS   = re.compile(r"locus_tag=([^;]+)")
_PAT_NAME    = re.compile(r"(?:Name|gene)=([^;]+)")
_PAT_PRODUCT = re.compile(r"product=([^;]+)")

# Feature types we treat as a "gene record" for context analysis.
_GENE_FEATURES = {"gene", "pseudogene"}


def load_gff_genes(gff: Path) -> pd.DataFrame:
    """Return one row per gene/pseudogene record with 0-based half-open coords."""
    rows = []
    with gff.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in _GENE_FEATURES:
                continue
            chrom, ftype, start1, end1, strand, attr = (
                f[0], f[2], int(f[3]), int(f[4]), f[6], f[8])
            m = _PAT_LOCUS.search(attr)
            if not m:
                continue
            n = _PAT_NAME.search(attr)
            p = _PAT_PRODUCT.search(attr)
            rows.append({
                "locus_tag": m.group(1),
                "name":      n.group(1) if n else "",
                "chrom":     chrom,
                "strand":    strand,
                "start0":    start1 - 1,
                "end0":      end1,
                "feature":   ftype,
                "is_pseudo": ftype == "pseudogene" or "pseudo=true" in attr,
                "product":   p.group(1) if p else "",
                "raw":       line.rstrip("\n"),
            })
    return pd.DataFrame(rows).drop_duplicates(subset="locus_tag")


def load_deletions(bed: Path) -> pd.DataFrame:
    df = pd.read_csv(bed, sep="\t")
    df = df.rename(columns={c: c.lstrip("#") for c in df.columns})
    return df[["chrom", "start0", "end"]].rename(columns={"end": "end0"})


# ---------------------------------------------------------------- retained set

def compute_deleted_bp(genes: pd.DataFrame, dels: pd.DataFrame) -> pd.Series:
    """Return a Series indexed like `genes` with bp of each ORF overlapping any
    deletion."""
    out = []
    for g in genes.itertuples():
        d = dels[dels.chrom == g.chrom]
        bp = 0
        for r in d.itertuples():
            ov = min(int(r.end0), int(g.end0)) - max(int(r.start0), int(g.start0))
            if ov > 0:
                bp += ov
        out.append(bp)
    return pd.Series(out, index=genes.index)


# ---------------------------------------------------------------- lift syn1 -> syn3A

def build_shift_table(dels: pd.DataFrame) -> list[tuple[int, int]]:
    """Cumulative deletion length seen *before* a given syn1 coordinate.
    Returns a list of (start0, cum_deleted_bp_before_this_start) sorted by
    start0. To lift a syn1 coordinate p to syn3A: subtract the cum length of
    deletions that end at or before p.

    Assumes deletions are non-overlapping on a single chromosome (true for our
    BED file)."""
    d = dels.sort_values("start0").reset_index(drop=True)
    pairs = []
    cum = 0
    for r in d.itertuples():
        pairs.append((int(r.start0), cum))
        cum += int(r.end0) - int(r.start0)
    pairs.append((10**12, cum))  # sentinel past genome end
    return pairs


def lift_syn1_to_syn3a(p: int, shift_table: list[tuple[int, int]]) -> int | None:
    """Lift a syn1 coordinate `p` to the corresponding syn3A coordinate, or
    None if `p` falls inside a deletion."""
    # cum = sum of deleted bp whose start0 < p
    cum = 0
    for s, c in shift_table:
        if s <= p:
            cum = c + max(0, p - s) if False else c  # we'll recompute below
        else:
            break
    # Recompute cleanly: walk through deletions.
    return _lift(p, shift_table)


def _lift(p: int, shift_table: list[tuple[int, int]]) -> int | None:
    # Reconstruct deletions from shift_table.
    # shift_table[i] = (start_i, cum_before_i); deletion length_i = cum_{i+1} - cum_i
    cum = 0
    for i in range(len(shift_table) - 1):
        s_i, c_i = shift_table[i]
        s_next, c_next = shift_table[i + 1]
        del_len = c_next - c_i
        e_i = s_i + del_len
        if p < s_i:
            return p - c_i
        if p < e_i:
            return None  # inside the deletion
        # otherwise continue past this deletion
    # past all deletions
    return p - shift_table[-1][1]


# ---------------------------------------------------------------- neighbors (circular, same-strand)

def _neighbor_pairs(sub: pd.DataFrame, chrom_len: int) -> dict:
    """For each gene in `sub` (any strand or filtered), return its physical
    left (ccw / decreasing coord) and right (cw / increasing coord) neighbor
    with circular wrap. Distance is the intergenic bp gap between CDS spans."""
    sub = sub.sort_values("start0").reset_index(drop=True)
    n = len(sub)
    out: dict = {}
    if n == 0:
        return out
    for i, g in enumerate(sub.itertuples()):
        li, ri = (i - 1) % n, (i + 1) % n
        left  = sub.iloc[li]
        right = sub.iloc[ri]
        if n == 1:
            ld = rd = chrom_len - int(g.end0) + int(g.start0)
        else:
            ld = (chrom_len - int(left.end0) + int(g.start0)
                  if li == n - 1 else int(g.start0) - int(left.end0))
            rd = (chrom_len - int(g.end0) + int(right.start0)
                  if ri == 0 else int(right.start0) - int(g.end0))
        out[g.locus_tag] = {
            "ccw_locus":  left.locus_tag,
            "ccw_dist":   int(ld),
            "ccw_strand": left.strand,
            "cw_locus":   right.locus_tag,
            "cw_dist":    int(rd),
            "cw_strand":  right.strand,
        }
    return out


def same_strand_neighbors(genes: pd.DataFrame, chrom_len: int) -> dict:
    out: dict = {}
    for strand in ("+", "-"):
        out.update(_neighbor_pairs(genes[genes.strand == strand], chrom_len))
    return out


def any_strand_neighbors(genes: pd.DataFrame, chrom_len: int) -> dict:
    return _neighbor_pairs(genes, chrom_len)


def transcription_neighbors(g_strand: str, nb: dict) -> tuple[str, int, str, int]:
    """Map ccw / cw to transcription upstream / downstream.
    + strand: upstream = ccw, downstream = cw
    - strand: upstream = cw,  downstream = ccw"""
    if g_strand == "+":
        return nb["ccw_locus"], nb["ccw_dist"], nb["cw_locus"], nb["cw_dist"]
    return nb["cw_locus"], nb["cw_dist"], nb["ccw_locus"], nb["ccw_dist"]


# ---------------------------------------------------------------- unaltered bp in syn1 frame

def unaltered_bp(g_start: int, g_end: int,
                 del_sorted: pd.DataFrame, chrom_len: int) -> tuple[int, int]:
    """Return (unaltered_ccw_bps, unaltered_cw_bps) for a gene in the syn1
    frame, with circular wrap. Strand-agnostic: ccw = decreasing coord side,
    cw = increasing coord side."""
    starts = del_sorted.start0.to_numpy()
    ends   = del_sorted.end0.to_numpy()
    n = len(del_sorted)

    # physical-left wall: largest end0 <= g_start; if none, wrap to largest end0
    left_wall_end = None
    # iterate deletions whose end0 <= g_start
    candidates = ends[ends <= g_start]
    if candidates.size:
        left_wall_end = int(candidates.max())
        left_unalt = g_start - left_wall_end
    else:
        # wrap: take the deletion with the largest end0 in the chromosome
        if n:
            left_wall_end = int(ends.max())
            left_unalt = (chrom_len - left_wall_end) + g_start
        else:
            left_unalt = chrom_len  # no deletions at all

    # physical-right wall: smallest start0 >= g_end
    candidates = starts[starts >= g_end]
    if candidates.size:
        right_wall_start = int(candidates.min())
        right_unalt = right_wall_start - g_end
    else:
        if n:
            right_wall_start = int(starts.min())
            right_unalt = (chrom_len - g_end) + right_wall_start
        else:
            right_unalt = chrom_len

    return int(left_unalt), int(right_unalt)


# ---------------------------------------------------------------- operon change

def annotate_operon_change(out: pd.DataFrame, operon_tsv: Path) -> list[str]:
    """For each retained gene in `out`, look up the operon it belongs to
    (via sense_gene_locusNums) and assign an operon_change tag.

    Rules (only the first matching rule applies; otherwise empty):
        leading_gene_deleted : operon gene_deletion_pattern == 'leading_deleted'
        promoter_deleted     : operon per_hit_classes == "5'_truncation_UTR"
                               AND operon gene_deletion_pattern == 'intact'
    """
    op = pd.read_csv(operon_tsv, sep="\t", dtype=str)
    # Map locusNum -> (operon_id, gene_deletion_pattern, per_hit_classes)
    lookup: dict[str, tuple[str, str, str]] = {}
    for r in op.itertuples():
        raw = r.sense_gene_locusNums if isinstance(r.sense_gene_locusNums, str) else ""
        for ln in raw.split(","):
            ln = ln.strip()
            if ln:
                lookup[ln] = (r.operon_id,
                              r.gene_deletion_pattern or "",
                              r.per_hit_classes or "")

    tags: list[str] = []
    for r in out.itertuples():
        ln = r.locus_tag.rsplit("_", 1)[-1]
        info = lookup.get(ln)
        if info is None:
            tags.append("")
            continue
        _, gdp, phc = info
        if gdp == "leading_deleted":
            tags.append("leading_gene_deleted")
        elif phc == "5'_truncation_UTR" and gdp == "intact":
            tags.append("promoter_deleted")
        else:
            tags.append("")
    return tags


# ---------------------------------------------------------------- boxplot

def boxplot_TPM_fold_by_operon_change(out: pd.DataFrame, rna_csv: Path,
                                      out_pdf: Path) -> str:
    """Box-plot of TPM_fold_change_Illumina (syn3A / syn1) for retained
    mRNA-coding genes grouped by operon_change tag (leading_gene_deleted,
    promoter_deleted, others). Returns a short text summary."""
    rna = pd.read_csv(rna_csv)
    rna = rna[rna.rna_type == "mRNA"].copy()
    rna["TPM_fold_change_Illumina"] = pd.to_numeric(
        rna["TPM_fold_change_Illumina"], errors="coerce")
    rna = rna.dropna(subset=["TPM_fold_change_Illumina"])

    merged = out.merge(rna, left_on="locus_tag", right_on="locus_syn1", how="inner")
    merged["category"] = merged.operon_change.fillna("").replace("", "others")
    merged.loc[~merged.category.isin(
        ["leading_gene_deleted", "promoter_deleted"]), "category"] = "others"

    cats   = ["leading_gene_deleted", "promoter_deleted", "others"]
    colors = ["#d62728", "#ff7f0e", "#7f7f7f"]
    data   = [merged.loc[merged.category == c, "TPM_fold_change_Illumina"].to_numpy()
              for c in cats]

    fig, ax = plt.subplots(figsize=(6.0, 4.5))
    bp = ax.boxplot(data, tick_labels=cats, patch_artist=True,
                    widths=0.55, showfliers=False)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.55); patch.set_edgecolor("black")
    for med in bp["medians"]:
        med.set_color("black"); med.set_linewidth(1.6)

    rng = np.random.default_rng(7)
    for i, vals in enumerate(data, start=1):
        x = rng.normal(loc=i, scale=0.06, size=len(vals))
        ax.scatter(x, vals, s=10, color=colors[i - 1],
                   edgecolor="black", linewidth=0.3, alpha=0.7, zorder=3)

    ax.set_yscale("log")
    ax.axhline(1.0, color="black", lw=0.5, ls=":")
    ax.set_ylabel("TPM fold change (Illumina, syn3A / syn1)")
    ax.set_xlabel("operon_change (retained mRNA genes)")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ymin = ax.get_ylim()[0]
    for i, vals in enumerate(data, start=1):
        if len(vals):
            ax.text(i, ymin, f"n={len(vals)}\nmed={np.median(vals):.2f}",
                    ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf)
    plt.close(fig)

    lines = ["", "TPM_fold_change_Illumina (Illumina) by operon_change:",
             f"  {'group':<24s} {'n':>4s} {'median':>8s} {'mean':>7s} {'min':>7s} {'max':>7s}"]
    for c, vals in zip(cats, data):
        if len(vals) == 0:
            lines.append(f"  {c:<24s} {0:>4d} {'NA':>8s} {'NA':>7s} {'NA':>7s} {'NA':>7s}")
        else:
            lines.append(f"  {c:<24s} {len(vals):>4d} {np.median(vals):>8.3f} "
                         f"{np.mean(vals):>7.3f} {np.min(vals):>7.3f} {np.max(vals):>7.3f}")
    return "\n".join(lines)


# ---------------------------------------------------------------- main

def main() -> None:
    syn1  = load_gff_genes(SYN1_GFF)
    syn3a = load_gff_genes(SYN3A_GFF)
    dels  = load_deletions(DEL_BED)

    # ---- retained set
    syn1 = syn1.reset_index(drop=True)
    syn1["deleted_bp"] = compute_deleted_bp(syn1, dels)
    syn1["orf_len"]    = syn1.end0 - syn1.start0
    retained = syn1[syn1.deleted_bp == 0].copy()

    # ---- cross-check against syn3A annotation by locus-tag numeric suffix.
    # syn3A preserves syn1 numbering: MMSYN1_0025 <-> JCVISYN3A_0025.
    def _suffix(lt: str) -> str:
        return lt.rsplit("_", 1)[-1]

    syn3a_by_suffix = {_suffix(r.locus_tag): r for r in syn3a.itertuples()}

    syn3a_locus_col, syn3a_name_col = [], []
    syn3a_s_col, syn3a_e_col = [], []
    matched, unmatched = 0, 0
    unmatched_loci = []
    for r in retained.itertuples():
        m = syn3a_by_suffix.get(_suffix(r.locus_tag))
        if m is not None:
            syn3a_locus_col.append(m.locus_tag)
            syn3a_name_col.append(m.name)
            syn3a_s_col.append(int(m.start0))
            syn3a_e_col.append(int(m.end0))
            matched += 1
        else:
            syn3a_locus_col.append("")
            syn3a_name_col.append("")
            syn3a_s_col.append(None)
            syn3a_e_col.append(None)
            unmatched += 1
            unmatched_loci.append(r.locus_tag)
    retained["syn3A_locus_tag"]  = syn3a_locus_col
    retained["syn3A_name"]       = syn3a_name_col
    retained["syn3A_start0"]     = syn3a_s_col
    retained["syn3A_end0"]       = syn3a_e_col

    # ---- neighbors: same-strand (transcription context) + any-strand (genomic).
    syn1_tx   = same_strand_neighbors(syn1,  SYN1_LEN)
    syn3a_tx  = same_strand_neighbors(syn3a, SYN3A_LEN)
    syn1_geo  = any_strand_neighbors(syn1,  SYN1_LEN)
    syn3a_geo = any_strand_neighbors(syn3a, SYN3A_LEN)

    # ---- unaltered bp (in syn1 frame, strand-agnostic)
    del_sorted = dels[dels.chrom == "CP002027.1"].sort_values("start0").reset_index(drop=True)

    rows = []
    for r in retained.itertuples():
        # transcription-direction (same-strand) neighbors
        nb = syn1_tx.get(r.locus_tag)
        if nb:
            s1_up_loc, s1_up_dist, s1_dn_loc, s1_dn_dist = transcription_neighbors(r.strand, nb)
        else:
            s1_up_loc = s1_dn_loc = ""
            s1_up_dist = s1_dn_dist = 0

        s3_up_loc = s3_dn_loc = ""
        s3_up_dist = s3_dn_dist = 0
        if r.syn3A_locus_tag and r.syn3A_locus_tag in syn3a_tx:
            s3_up_loc, s3_up_dist, s3_dn_loc, s3_dn_dist = \
                transcription_neighbors(r.strand, syn3a_tx[r.syn3A_locus_tag])

        # genomic (any-strand, ccw/cw) neighbors
        g1 = syn1_geo.get(r.locus_tag, {})
        g3 = syn3a_geo.get(r.syn3A_locus_tag, {}) if r.syn3A_locus_tag else {}

        unalt_ccw, unalt_cw = unaltered_bp(int(r.start0), int(r.end0),
                                           del_sorted, SYN1_LEN)

        rows.append({
            "locus_tag":  r.locus_tag,
            "name":       r.name,
            "chrom":      r.chrom,
            "start0":     int(r.start0),
            "end0":       int(r.end0),
            "strand":     r.strand,
            "syn3A_locus_tag":  r.syn3A_locus_tag,
            "syn3A_name":       r.syn3A_name,
            "syn3A_start0":     r.syn3A_start0,
            "syn3A_end0":       r.syn3A_end0,

            "syn1_upstream_locus":   s1_up_loc,
            "syn1_upstream_dist":    s1_up_dist,
            "syn1_downstream_locus": s1_dn_loc,
            "syn1_downstream_dist":  s1_dn_dist,
            "syn3A_upstream_locus":   s3_up_loc,
            "syn3A_upstream_dist":    s3_up_dist,
            "syn3A_downstream_locus": s3_dn_loc,
            "syn3A_downstream_dist":  s3_dn_dist,

            "syn1_cw_neighbor":   g1.get("cw_locus", ""),
            "syn1_cw_dist":       g1.get("cw_dist", 0),
            "syn1_cw_strand":     g1.get("cw_strand", ""),
            "syn1_ccw_neighbor":  g1.get("ccw_locus", ""),
            "syn1_ccw_dist":      g1.get("ccw_dist", 0),
            "syn1_ccw_strand":    g1.get("ccw_strand", ""),
            "syn3A_cw_neighbor":  g3.get("cw_locus", ""),
            "syn3A_cw_dist":      g3.get("cw_dist", 0),
            "syn3A_cw_strand":    g3.get("cw_strand", ""),
            "syn3A_ccw_neighbor": g3.get("ccw_locus", ""),
            "syn3A_ccw_dist":     g3.get("ccw_dist", 0),
            "syn3A_ccw_strand":   g3.get("ccw_strand", ""),

            "unaltered_ccw_bps": unalt_ccw,
            "unaltered_cw_bps":  unalt_cw,
        })

    out = pd.DataFrame(rows)

    # cw / ccw context_changed: True when the syn3A neighbor (mapped back to
    # its MMSYN1 ancestor by locus-tag suffix) differs from the syn1 neighbor.
    def _changed(syn3a_nb_loc: str, syn1_nb_loc: str) -> object:
        if not syn3a_nb_loc:
            return ""
        suf = syn3a_nb_loc.rsplit("_", 1)[-1]
        ancestor = f"MMSYN1_{suf}"
        return ancestor != syn1_nb_loc
    out["cw_context_changed"]  = [_changed(a, b) for a, b in
                                  zip(out.syn3A_cw_neighbor, out.syn1_cw_neighbor)]
    out["ccw_context_changed"] = [_changed(a, b) for a, b in
                                  zip(out.syn3A_ccw_neighbor, out.syn1_ccw_neighbor)]

    # ---- operon_change: classify by operon-level deletion patterns.
    out["operon_change"] = annotate_operon_change(out, OPERON_TSV)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)

    boxplot_summary = boxplot_TPM_fold_by_operon_change(out, RNA_CSV, OUT_PDF)

    # ---- summary
    lines = []
    lines.append("=" * 64)
    lines.append("RETAINED-GENE CONTEXT SUMMARY")
    lines.append("=" * 64)
    lines.append(f"syn1 genes total              : {len(syn1):5d}")
    lines.append(f"syn1 genes retained (no dele) : {len(retained):5d}")
    lines.append(f"syn3A annotated genes         : {len(syn3a):5d}")
    lines.append("")
    lines.append("Cross-check: syn1 retained <-> syn3A annotation (match by locus_tag suffix)")
    lines.append(f"  matched   : {matched}")
    lines.append(f"  unmatched : {unmatched}")
    # Reverse: syn3A genes not corresponding to any retained syn1 gene
    matched_syn3a = set(out[out.syn3A_locus_tag != ""].syn3A_locus_tag)
    syn3a_extra = [r.locus_tag for r in syn3a.itertuples() if r.locus_tag not in matched_syn3a]
    lines.append(f"  syn3A loci with no syn1 retained mate : {len(syn3a_extra)}")

    # Abnormal cases: full annotation lines from both GFFs side by side.
    syn1_by_lt  = {r.locus_tag: r for r in syn1.itertuples()}
    syn3a_by_lt = {r.locus_tag: r for r in syn3a.itertuples()}
    suf2syn3a = {_suffix(lt): lt for lt in syn3a_by_lt}
    suf2syn1  = {_suffix(lt): lt for lt in syn1_by_lt}

    def _annot_lines(header: str, loci_syn1: list[str], loci_syn3a: list[str]):
        lines.append("")
        lines.append("-" * 64)
        lines.append(header)
        lines.append("-" * 64)
        seen = set()
        for lt in loci_syn1:
            if lt in seen:
                continue
            seen.add(lt)
            lines.append(f"\n[{lt}]")
            r1 = syn1_by_lt.get(lt)
            lines.append("  syn1 : " + (r1.raw if r1 else "<not in syn1 gff>"))
            mate = syn3a_by_lt.get(suf2syn3a.get(_suffix(lt), ""))
            lines.append("  syn3A: " + (mate.raw if mate else "<no syn3A record>"))
        for lt in loci_syn3a:
            mate_lt = suf2syn1.get(_suffix(lt), "")
            if mate_lt in seen:
                continue
            seen.add(lt)
            lines.append(f"\n[{lt}]")
            mate1 = syn1_by_lt.get(mate_lt)
            lines.append("  syn1 : " + (mate1.raw if mate1 else "<no syn1 record>"))
            r3 = syn3a_by_lt.get(lt)
            lines.append("  syn3A: " + (r3.raw if r3 else "<not in syn3a gff>"))

    _annot_lines(
        "ABNORMAL CASES (retained syn1 not matched, or syn3A without syn1 mate)",
        unmatched_loci, syn3a_extra,
    )
    lines.append("")
    n_cw  = int(out.cw_context_changed.apply(lambda v: v is True).sum())
    n_ccw = int(out.ccw_context_changed.apply(lambda v: v is True).sum())
    lines.append(f"cw  neighbor changed in syn3A : {n_cw} retained genes")
    lines.append(f"ccw neighbor changed in syn3A : {n_ccw} retained genes")
    lines.append("")
    lines.append("operon_change tags:")
    oc = out.operon_change.value_counts(dropna=False)
    for k, v in oc.items():
        label = k if k else "(none)"
        lines.append(f"  {label:<24s} {v:5d}")
    lines.append(boxplot_summary)
    txt = "\n".join(lines) + "\n"

    OUT_TXT.write_text(txt)
    print(txt)
    print(f"Wrote: {OUT_TSV}")
    print(f"Wrote: {OUT_TXT}")


if __name__ == "__main__":
    main()
