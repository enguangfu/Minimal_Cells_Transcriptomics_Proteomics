"""
Degradation Asymmetry Analysis (per-isoform, genomic-context based)
====================================================================
Quantify the directionality of mRNA degradation in Syn1 by labelling
each PacBio FLNC isoform endpoint as intragenic or intergenic.

Logic:
  A canonical, intact transcript starts at a promoter (5' end in an
  intergenic region) and ends at a terminator (3' end in an intergenic
  region). Endpoints that fall *inside* a gene body cannot be explained
  by promoter/terminator usage and must come from RNA processing or
  exonucleolytic erosion:

    - 5' end intragenic  →  evidence of 5'→3' exoribonuclease activity
                            (RNase J1/J2)
    - 3' end intragenic  →  evidence of 3'→5' exoribonuclease activity
                            (RNase R / YhaM)

  Each isoform is then placed into one of four categories:

    ┌──────────────────────────┬─────────────┬─────────────┐
    │  Category                │ 5' end      │ 3' end      │
    │──────────────────────────┼─────────────┼─────────────│
    │  canonical               │ intergenic  │ intergenic  │
    │  5p_intragenic_only      │ intragenic  │ intergenic  │
    │  3p_intragenic_only      │ intergenic  │ intragenic  │
    │  both_intragenic         │ intragenic  │ intragenic  │
    └──────────────────────────┴─────────────┴─────────────┘

  The ratio
        N(5'-intragenic isoforms) / N(3'-intragenic isoforms)
  (with both single- and read-weighted forms) directly reports the
  asymmetry between 5'→3' and 3'→5' exoribonuclease activities.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import os

os.makedirs("RNase", exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════
MOTHER_FOLDER = ".."
ISOFORMS_TSV  = MOTHER_FOLDER + "/isoform_annotation/isoform_clusters_annotated.tsv"
GFF3_FILE     = MOTHER_FOLDER + "/Genomes_Input/syn1.genes.gff3"
OUT_FOLDER    = "./RNase"

MIN_READS = 10    # only keep well-supported isoforms


# ═══════════════════════════════════════════════════════════════════════════════
# Step 1 — Load gene intervals → per-strand intragenic mask
# ═══════════════════════════════════════════════════════════════════════════════
def load_intragenic_mask(gff3_path: str) -> dict:
    """Return {'+': bool array, '-': bool array} marking intragenic
    positions on each strand, indexed in 0-based genomic coordinates.
    Only protein-coding genes (rna_type=mRNA) are included."""
    gene_rows = []
    chrom_len = 0
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[2] != "gene":
                continue
            attrs = f[8] if len(f) >= 9 else ""
            if "rna_type=mRNA" not in attrs:
                continue
            s1, e1, strand = int(f[3]), int(f[4]), f[6]
            gene_rows.append((s1 - 1, e1, strand))   # → 0-based half-open
            if e1 > chrom_len:
                chrom_len = e1
    mask = {"+": np.zeros(chrom_len, dtype=bool),
            "-": np.zeros(chrom_len, dtype=bool)}
    for s0, e0, strand in gene_rows:
        if strand in mask:
            mask[strand][s0:e0] = True
    print(f"Loaded {len(gene_rows)} protein-coding gene intervals from GFF3 "
          f"(chrom length {chrom_len:,})")
    return mask

INTRAGENIC = load_intragenic_mask(GFF3_FILE)


# ═══════════════════════════════════════════════════════════════════════════════
# Step 2 — Load isoforms and label endpoint context
# ═══════════════════════════════════════════════════════════════════════════════
df = pd.read_csv(ISOFORMS_TSV, sep="\t")
df_iso = df[df["n_reads"] >= MIN_READS].copy().reset_index(drop=True)
print(f"Isoforms loaded (n_reads >= {MIN_READS}): {len(df_iso)}")


def endpoint_positions(row):
    """Return (5'-end pos0, 3'-end pos0) for an isoform."""
    if row["strand"] == "+":
        return int(row["start0"]), int(row["end0"]) - 1
    else:
        return int(row["end0"]) - 1, int(row["start0"])


def is_intragenic(pos0: int, strand: str) -> bool:
    mask = INTRAGENIC[strand]
    return bool(0 <= pos0 < len(mask) and mask[pos0])


pos5_list, pos3_list, in5_list, in3_list = [], [], [], []
for _, row in df_iso.iterrows():
    p5, p3 = endpoint_positions(row)
    pos5_list.append(p5)
    pos3_list.append(p3)
    in5_list.append(is_intragenic(p5, row["strand"]))
    in3_list.append(is_intragenic(p3, row["strand"]))

df_iso["pos5p_0"]       = pos5_list
df_iso["pos3p_0"]       = pos3_list
df_iso["intragenic_5p"] = in5_list
df_iso["intragenic_3p"] = in3_list


def categorise(row) -> str:
    if row["intragenic_5p"] and row["intragenic_3p"]:
        return "both_intragenic"
    if row["intragenic_5p"]:
        return "5p_intragenic_only"
    if row["intragenic_3p"]:
        return "3p_intragenic_only"
    return "canonical"

df_iso["category"] = df_iso.apply(categorise, axis=1)


# ═══════════════════════════════════════════════════════════════════════════════
# Step 3 — Summary statistics
# ═══════════════════════════════════════════════════════════════════════════════
CATS = ["canonical", "5p_intragenic_only", "3p_intragenic_only", "both_intragenic"]
LABELS = ["Canonical\n(5' & 3' intergenic)",
          "5' intragenic only",
          "3' intragenic only",
          "Both ends eroded)"]
COLORS = ["#969696", "#2166AC", "#B2182B", "#762A83"]

cat_counts = df_iso["category"].value_counts()
cat_reads  = df_iso.groupby("category")["n_reads"].sum()
total_iso  = len(df_iso)
total_reads = int(df_iso["n_reads"].sum())

print("\n" + "="*70)
print("ENDPOINT-CONTEXT CATEGORY SUMMARY")
print("="*70)
print(f"{'Category':<24} {'Isoforms':>10} {'%':>8} {'Reads':>14} {'%':>8}")
print("-"*70)
for cat in CATS:
    n_iso  = int(cat_counts.get(cat, 0))
    n_read = int(cat_reads.get(cat, 0))
    print(f"{cat:<24} {n_iso:>10,} {n_iso/total_iso*100:>7.1f}% "
          f"{n_read:>14,} {n_read/total_reads*100:>7.1f}%")
print("-"*70)
print(f"{'TOTAL':<24} {total_iso:>10,} {'100.0%':>8} "
      f"{total_reads:>14,} {'100.0%':>8}")

# Overall (any) intragenic-end counts: union over the two single-end groups
n_5p_intra = int((df_iso["intragenic_5p"]).sum())
n_3p_intra = int((df_iso["intragenic_3p"]).sum())
r_5p_intra = int(df_iso.loc[df_iso["intragenic_5p"], "n_reads"].sum())
r_3p_intra = int(df_iso.loc[df_iso["intragenic_3p"], "n_reads"].sum())

print("\nIntragenic-endpoint totals (counting any isoform with that end intragenic):")
print(f"  5' intragenic: {n_5p_intra:,} isoforms / {r_5p_intra:,} reads")
print(f"  3' intragenic: {n_3p_intra:,} isoforms / {r_3p_intra:,} reads")

ratio_iso  = n_5p_intra / n_3p_intra if n_3p_intra > 0 else float("inf")
ratio_read = r_5p_intra / r_3p_intra if r_3p_intra > 0 else float("inf")
print(f"\n5'→3' / 3'→5' ratio (by unique isoform kinds):  {ratio_iso:.2f}")
print(f"5'→3' / 3'→5' ratio (by unique isoform counts): {ratio_read:.2f}")

# Per-strand breakdown
print("\nPer-strand breakdown:")
for strand in ["+", "-"]:
    sub = df_iso[df_iso["strand"] == strand]
    n5 = int(sub["intragenic_5p"].sum())
    n3 = int(sub["intragenic_3p"].sum())
    r5 = int(sub.loc[sub["intragenic_5p"], "n_reads"].sum())
    r3 = int(sub.loc[sub["intragenic_3p"], "n_reads"].sum())
    rstr_i = f"{n5/n3:.2f}" if n3 > 0 else "inf"
    rstr_r = f"{r5/r3:.2f}" if r3 > 0 else "inf"
    print(f"  {strand} strand: 5'-intra={n5} ({r5:,} reads), "
          f"3'-intra={n3} ({r3:,} reads), "
          f"ratio_iso={rstr_i}, ratio_read={rstr_r}")


# ═══════════════════════════════════════════════════════════════════════════════
# Step 4 — Figures (one PDF per panel)
# ═══════════════════════════════════════════════════════════════════════════════

# --- Panel A: Category bar chart (unique-isoform vs read-weighted) ---
figA, ax = plt.subplots(figsize=(8, 5))
iso_vals  = [int(cat_counts.get(c, 0)) for c in CATS]
read_vals = [int(cat_reads.get(c, 0))  for c in CATS]
iso_pcts  = [v / total_iso  * 100 for v in iso_vals]
read_pcts = [v / total_reads * 100 for v in read_vals]

x = np.arange(len(CATS))
w = 0.35
bars1 = ax.bar(x - w/2, iso_pcts,  w, color=COLORS, alpha=0.85,
               edgecolor="white", label="By unique isoform kinds")
bars2 = ax.bar(x + w/2, read_pcts, w, color=COLORS, alpha=0.50,
               edgecolor="white", hatch="//", label="By unique isoform counts")
ax.set_xticks(x)
ax.set_xticklabels(LABELS, fontsize=12, ha="center")
ax.set_ylabel("Percentage (%)", fontsize=12)
# ax.set_title("A. Isoform endpoint-context categories",
            #  fontsize=11, fontweight="bold")
ax.legend(fontsize=12, loc="upper left")
ax.yaxis.set_major_formatter(mticker.PercentFormatter(decimals=0))
for bar, val in zip(bars1, iso_vals):
    if val > 0:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f"n={val:,}", ha="center", va="bottom", fontsize=7.5)
for bar, val in zip(bars2, read_vals):
    if val > 0:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f"{round(val/1000)}k", ha="center", va="bottom", fontsize=7.5)
figA.tight_layout()
figA.savefig(f"{OUT_FOLDER}/RNase_asymmetry_isoform_categories.pdf",
             dpi=300, bbox_inches="tight")
plt.close(figA)


# --- Panel B: 5'→3' vs 3'→5' exo activity bar chart + ratio ---
figB, ax = plt.subplots(figsize=(7, 5))
groups = ["5'→3' exo\n(5' intragenic)", "3'→5' exo\n(3' intragenic)"]
iso_counts_grp  = [n_5p_intra, n_3p_intra]
read_counts_grp = [r_5p_intra, r_3p_intra]
group_colors    = ["#2166AC", "#B2182B"]

x = np.arange(len(groups))
w = 0.35
b1 = ax.bar(x - w/2, iso_counts_grp,  w, color=group_colors, alpha=0.85,
            edgecolor="white", label="By unique isoform kinds")
ax2 = ax.twinx()
b2 = ax2.bar(x + w/2, read_counts_grp, w, color=group_colors, alpha=0.45,
             edgecolor="white", hatch="//", label="By unique isoform counts")

ax.set_xticks(x)
ax.set_xticklabels(groups, fontsize=10)
ax.set_ylabel("Unique isoform kinds", fontsize=10, color="#222222")
ax2.set_ylabel("Unique isoform counts (reads)", fontsize=10, color="#222222")
# ax.set_title("B. Exoribonuclease activity asymmetry\n"
#              f"ratio (kinds) = {ratio_iso:.2f},  ratio (counts) = {ratio_read:.2f}",
#              fontsize=11, fontweight="bold")

for bar, val in zip(b1, iso_counts_grp):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
            f"{val:,}", ha="center", va="bottom", fontsize=8.5)
for bar, val in zip(b2, read_counts_grp):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
             f"{val:,}", ha="center", va="bottom", fontsize=8.5)

# combined legend
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=8.5, loc="upper right")

figB.tight_layout()
figB.savefig(f"{OUT_FOLDER}/RNase_asymmetry_exoribo_categories.pdf",
             dpi=300, bbox_inches="tight")
plt.close(figB)

print("\nFigures saved:")
print(f"  {OUT_FOLDER}/RNase_asymmetry_isoform_categories.pdf")
print(f"  {OUT_FOLDER}/RNase_asymmetry_exoribo_categories.pdf")


# ═══════════════════════════════════════════════════════════════════════════════
# Step 5 — Save annotated isoform table
# ═══════════════════════════════════════════════════════════════════════════════
out_tsv = f"{OUT_FOLDER}/isoform_endpoint_context.tsv"
keep_cols = ["isoform_id", "chrom", "strand", "start0", "end0", "n_reads",
             "pos5p_0", "pos3p_0", "intragenic_5p", "intragenic_3p", "category"]
df_iso[[c for c in keep_cols if c in df_iso.columns]].to_csv(
    out_tsv, sep="\t", index=False)
print(f"Saved: {out_tsv} ({len(df_iso)} rows)")


# ═══════════════════════════════════════════════════════════════════════════════
# Step 6 — Estimate ORFs with start codon but no stop codon
# ═══════════════════════════════════════════════════════════════════════════════
# For every filtered isoform we walk over the canonical (GFF3) genes on the
# same strand. Whenever the gene's *start codon* lies fully inside the isoform
# we count this as one observed ORF (one start/stop pair). We then check
# whether the matching stop codon is also inside the isoform — if not, the
# 3' end of the read is truncated relative to the canonical CDS.
#
# Convention: GFF3 gene coordinates already include the stop codon.
#   + strand:  start codon = [start0,   start0+3)
#              stop  codon = [end0-3,   end0)
#   − strand:  start codon = [end0-3,   end0)
#              stop  codon = [start0,   start0+3)

def load_gene_intervals(gff3_path: str) -> pd.DataFrame:
    """Load gene intervals, keeping only protein-coding genes (rna_type=mRNA)."""
    rows = []
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[2] != "gene":
                continue
            attrs = f[8] if len(f) >= 9 else ""
            if "rna_type=mRNA" not in attrs:
                continue
            s1, e1, strand = int(f[3]), int(f[4]), f[6]
            locus = ""
            for kv in attrs.split(";"):
                if kv.startswith("locus_tag="):
                    locus = kv.split("=", 1)[1]
                    break
            rows.append({
                "locus_tag": locus,
                "strand":    strand,
                "start0":    s1 - 1,
                "end0":      e1,
            })
    return pd.DataFrame(rows)


genes_df = load_gene_intervals(GFF3_FILE)
print(f"\nLoaded {len(genes_df)} gene intervals for ORF analysis")

# Split by strand for fast per-isoform lookup
genes_by_strand = {
    s: g.sort_values("start0").reset_index(drop=True)
    for s, g in genes_df.groupby("strand")
}

iso_orf_rows = []
for _, iso in df_iso.iterrows():
    strand = iso["strand"]
    iso_s, iso_e = int(iso["start0"]), int(iso["end0"])
    g = genes_by_strand.get(strand)
    if g is None:
        continue
    overlap = g[(g["start0"] < iso_e) & (g["end0"] > iso_s)]

    # Walk overlapping genes in transcription order along the isoform
    # (5'→3'). For + strand that's ascending start0, for − strand it's
    # descending end0. We then collect every gene whose start codon is
    # fully contained inside the isoform — each such (gene, isoform) pair
    # is one ORF observation.
    if strand == "+":
        ordered = overlap.sort_values("start0")
    else:
        ordered = overlap.sort_values("end0", ascending=False)

    orf_loci, has_stop_flags = [], []
    for _, gene in ordered.iterrows():
        gs, ge = int(gene["start0"]), int(gene["end0"])
        if strand == "+":
            start_codon = (gs, gs + 3)
            stop_codon  = (ge - 3, ge)
        else:
            start_codon = (ge - 3, ge)
            stop_codon  = (gs, gs + 3)

        start_in = (start_codon[0] >= iso_s) and (start_codon[1] <= iso_e)
        if not start_in:
            continue
        stop_in = (stop_codon[0] >= iso_s) and (stop_codon[1] <= iso_e)
        orf_loci.append(str(gene["locus_tag"]))
        has_stop_flags.append(bool(stop_in))

    if not orf_loci:
        continue

    # Compact text view of the ORFs along the isoform, e.g.
    #   "MMSYN1_0001✓ MMSYN1_0002✓ MMSYN1_0003✗"
    orfs_str = " ".join(
        f"{lt}{'✓' if hs else '✗'}"
        for lt, hs in zip(orf_loci, has_stop_flags)
    )
    iso_orf_rows.append({
        "orfs":       orfs_str,
        "isoform_id": iso["isoform_id"],
        "n_reads":    int(iso["n_reads"]),
        "strand":     strand,
        "n_orfs":     len(orf_loci),
        "n_with_stop":    int(sum(has_stop_flags)),
        "n_without_stop": int(len(orf_loci) - sum(has_stop_flags)),
    })

# Per-isoform table — the format the user asked for
orf_df = pd.DataFrame(iso_orf_rows,
                      columns=["orfs", "isoform_id", "n_reads", "strand",
                               "n_orfs", "n_with_stop", "n_without_stop"])

# Per-ORF totals derived directly from the isoform-level table.
# Each ORF in a polycistronic isoform is counted independently:
# n_orfs gives the multiplicity, so summing n_with_stop / n_without_stop
# yields the per-ORF totals exactly.
n_total     = int(orf_df["n_orfs"].sum())
n_no_stop   = int(orf_df["n_without_stop"].sum())
n_isoforms  = len(orf_df)
n_polycis   = int((orf_df["n_orfs"] >= 2).sum())

print(f"ORFs (start codon contained in an isoform): {n_total:,}")
print(f"  from {n_isoforms:,} isoforms ({n_polycis:,} carry ≥2 ORFs)")

if n_total > 0:
    pct_no_stop = n_no_stop / n_total * 100

    # Read-weighted: each ORF carries its isoform's read count
    r_total   = int((orf_df["n_orfs"]        * orf_df["n_reads"]).sum())
    r_no_stop = int((orf_df["n_without_stop"] * orf_df["n_reads"]).sum())
    pct_no_stop_r = r_no_stop / r_total * 100

    print("\n" + "="*70)
    print("ORFs WITH START CODON BUT NO STOP CODON")
    print("="*70)
    print("  (each ORF in a polycistronic isoform is counted independently)")
    print(f"  By unique ORF observations:")
    print(f"    no-stop / total = {n_no_stop:,} / {n_total:,} "
          f"= {pct_no_stop:.2f}%")
    print(f"  By read-weighted observations:")
    print(f"    no-stop / total = {r_no_stop:,} / {r_total:,} "
          f"= {pct_no_stop_r:.2f}%")

    print(f"\n  ORFs-per-isoform distribution: "
          f"mean={orf_df['n_orfs'].mean():.2f}, max={int(orf_df['n_orfs'].max())}")

    print("\n  Per-strand:")
    for strand in ["+", "-"]:
        sub = orf_df[orf_df["strand"] == strand]
        if len(sub) == 0:
            continue
        s_total = int(sub["n_orfs"].sum())
        s_ns    = int(sub["n_without_stop"].sum())
        print(f"    {strand} strand: no-stop {s_ns:,} / {s_total:,} "
              f"({s_ns/s_total*100:.2f}%)")

    orf_out = f"{OUT_FOLDER}/orf_start_stop_observations.tsv"
    orf_df.to_csv(orf_out, sep="\t", index=False)
    print(f"\nSaved: {orf_out} ({len(orf_df)} rows)")

print("\nDone.")
