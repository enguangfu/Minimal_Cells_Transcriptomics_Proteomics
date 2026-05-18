#!/usr/bin/env python
# coding: utf-8

# ## Operon Segmentation — Isoform-based
# 
# Segment the *Syn3A* minimal cell genome (~543 kb, 493 genes + 3 pseudogenes) into operons directly from ONT cDNA long-read RNA-seq.
#
# **Biological rationale:** In bacteria, genes are co-transcribed in polycistronic operons from a shared promoter. A single full-length RNA isoform spanning multiple genes is therefore the natural evidence for operon membership. ONT cDNA reads (after polyA trimming and length-based QC upstream) capture near-full-length transcripts from TSS to TTS, making them suitable for operon discovery without relying on depth drop-offs or intergenic gap heuristics.
# 
# **Key insight:** Under containment-only clustering, each operon is defined by its **longest isoform** — the full-length transcript spanning the true TSS to TTS. All shorter co-clustered isoforms are nested subsets arising from RNase degradation or premature termination and do not change the boundary.
# 
# **Complications addressed:**
# - *RNase activity* — bacteria actively degrade mRNA from the 3′ end and internally; this fragments a single operon transcript into multiple shorter isoforms that can appear as separate clusters.
# - *Low-abundance transcripts* — essential genes under low expression may produce too few reads to pass the isoform threshold; these are rescued using BAM-level read evidence.
# - *rRNA operons* — rRNA transcripts are too abundant and structured for isoform clustering; these two operons (each encoding 5S–23S–16S) are annotated directly from the genome.
# 
# **Single tunable parameter:** `MIN_READS` — minimum read support per isoform cluster.
# 
# **Input:** `isoform_clusters_annotated.tsv` — 80 K clustered isoforms with genomic coordinates, strand, and read support.
# 
# **Output:** `operons.candidate_blocks.tsv` — full operon table with gene annotations and segmentation type.

# ## Step 1 — Load isoforms
# 
# Filter the ~80 K clustered isoforms to those with at least `MIN_READS` supporting reads. Lower thresholds retain more isoforms but risk including sequencing noise and spurious clusters; higher thresholds miss lowly-expressed operons. At `MIN_READS = 50` the dominant transcription units are well-supported while RNase-generated fragments (which typically accumulate fewer reads than full-length transcripts) are partially suppressed.

# In[1]:


from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────────
MOTHER_FOLDER = ".."
ISOFORMS_TSV  = MOTHER_FOLDER + "/Syn3A_Transcriptomics/Isoform_Cluster/isoform_clusters_annotated.tsv"
OUT_FOLDER    = "."
Path(OUT_FOLDER).mkdir(parents=True, exist_ok=True)

# ── Single parameter ───────────────────────────────────────────────────────
MIN_READS = 5   # minimum read support per isoform — only tunable parameter

# ── Load ───────────────────────────────────────────────────────────────────
df = pd.read_csv(ISOFORMS_TSV, sep="\t")
df_iso = df[df["n_reads"] >= MIN_READS].copy()

print(f"Total isoforms loaded: {len(df)}")
print(f"Isoforms with n_reads >= {MIN_READS}: {len(df_iso)}")
print(f"\nStrand breakdown:")
print(df_iso["strand"].value_counts().to_string())


# ## Step 2 — Containment clustering → initial operons
# 
# **Biological basis:** Each operon produces one primary full-length transcript (longest isoform) and many degradation/truncation products. All products of the same operon are nested subsets of the full-length transcript. Two isoforms from *different* operons may overlap at their boundaries but one will not be fully contained within the other.
# 
# **Clustering rule:** Isoforms i and j join the same cluster iff one is fully contained within the other (within `BOUNDARY_TOL` bp tolerance to absorb soft-clipping). Partially overlapping isoforms — characteristic of adjacent operons whose 3′ and 5′ ends interleave — are **never** merged here.
# 
# **Mathematical guarantee:** The cluster boundary always equals the longest isoform. If any member extended beyond the longest isoform's boundaries, that member would itself be longer — a contradiction.
# 
# **Parameters:**
# - `MIN_READS` — isoform read-support threshold (Step 1)
# - `BOUNDARY_TOL` — boundary tolerance in bp; absorbs minor coordinate variation from soft-clipping at read ends
# 
# **Implementation:** Union-Find with path compression on containment pairs; sweep-line reduces candidate comparisons from O(n²) to O(n log n).

# In[2]:


BOUNDARY_TOL = 10   # bp — isoforms within this many bp of containment are treated as contained

def cluster_isoforms(isoforms: pd.DataFrame,
                     tol: int = BOUNDARY_TOL) -> list[dict]:
    """
    Cluster isoforms into operons by containment (with tolerance).

    Two isoforms i and j are in the same cluster iff one is within `tol` bp
    of being fully contained within the other:
        start_i >= start_j - tol  AND  end_i <= end_j + tol   (i inside j)
        OR vice-versa.

    The cluster boundary = the outermost (longest) isoform.
    Uses union-find + sweep-line for efficiency.
    """
    if isoforms.empty:
        return []

    iso    = isoforms.reset_index(drop=True)
    starts = iso["start0"].astype(int).tolist()
    ends   = iso["end0"].astype(int).tolist()
    n      = len(iso)

    # ── Union-Find with path compression ───────────────────────────────────
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        parent[find(x)] = find(y)

    # ── Sweep: for each isoform i, check candidates j that start near i ────
    order = sorted(range(n), key=lambda i: starts[i])
    for idx, i in enumerate(order):
        si, ei = starts[i], ends[i]
        for j in order[idx + 1:]:
            sj = starts[j]
            if sj >= ei + tol:   # j starts well after i ends — no containment possible
                break
            ej = ends[j]
            # near-containment: i inside j (with tol), or j inside i (with tol)
            i_in_j = (si >= sj - tol) and (ei <= ej + tol)
            j_in_i = (sj >= si - tol) and (ej <= ei + tol)
            if i_in_j or j_in_i:
                union(i, j)

    # ── Collect connected components → one block each ──────────────────────
    components = defaultdict(list)
    for i in range(n):
        components[find(i)].append(i)

    chrom  = iso.loc[0, "chrom"]
    strand = iso.loc[0, "strand"]
    blocks = []
    for indices in components.values():
        members = [iso.loc[k, "isoform_id"] for k in indices]
        reads   = [int(iso.loc[k, "n_reads"])  for k in indices]
        s0 = min(starts[k] for k in indices)
        e0 = max(ends[k]   for k in indices)
        blocks.append({
            "chrom": chrom, "strand": strand,
            "start0": s0, "end0": e0,
            "members": members, "reads": reads,
        })
    return sorted(blocks, key=lambda b: b["start0"])


def blocks_to_df(blocks: list[dict]) -> pd.DataFrame:
    rows = [{
        "chrom":         b["chrom"],
        "strand":        b["strand"],
        "start0":        b["start0"],
        "end0":          b["end0"],
        "n_isoforms":    len(b["members"]),
        "n_reads_total": int(np.sum(b["reads"])),
        "member_ids":    ",".join(b["members"]),
    } for b in blocks]
    return pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["chrom","strand","start0","end0","n_isoforms","n_reads_total","member_ids"])


# ── Cluster per strand ──────────────────────────────────────────────────────
operons_plus  = cluster_isoforms(df_iso[df_iso["strand"] == "+"])
operons_minus = cluster_isoforms(df_iso[df_iso["strand"] == "-"])

GENOME_LEN = 543_379     # CP016816.2 Syn3A genome length (bp)
def genome_coverage_pct(blocks, genome_len=GENOME_LEN):
    covered = sum(b["end0"] - b["start0"] for b in blocks)
    return covered, covered / genome_len * 100

cov_plus_bp,  cov_plus_pct  = genome_coverage_pct(operons_plus)
cov_minus_bp, cov_minus_pct = genome_coverage_pct(operons_minus)
print(f"BOUNDARY_TOL = {BOUNDARY_TOL} bp")
print(f"Operons — plus: {len(operons_plus)}, minus: {len(operons_minus)}, total: {len(operons_plus)+len(operons_minus)}")
print(f"Genome coverage — plus:  {cov_plus_bp:>8,} bp  ({cov_plus_pct:.1f}%)")
print(f"                  minus: {cov_minus_bp:>8,} bp  ({cov_minus_pct:.1f}%)")

# ── Build final DataFrame ───────────────────────────────────────────────────
all_blocks  = sorted(operons_plus + operons_minus, key=lambda b: (b["chrom"], b["start0"]))
operons_all = blocks_to_df(all_blocks)
operons_all.insert(0, "operon_id", [f"OP3A_{i+1:05d}" for i in range(len(operons_all))])
operons_all["length"] = operons_all["end0"] - operons_all["start0"]

# TSS/TTS: directional 5\'→3\' coordinates
# + strand: TSS = start0 (low), TTS = end0 (high)
# − strand: TSS = end0  (high), TTS = start0 (low)
operons_all["tss"] = np.where(operons_all["strand"] == "+",
                               operons_all["start0"], operons_all["end0"])
operons_all["tts"] = np.where(operons_all["strand"] == "+",
                               operons_all["end0"],  operons_all["start0"])

print(f"\nLength distribution (bp):")
print(operons_all["length"].describe().dropna().astype(int).to_string())

print(operons_all.drop(columns="member_ids").head(20))


# ## Step 3 — Gene annotation
# 
# Map the 496 annotated Syn3A loci (493 genes + 3 pseudogenes, from GFF3) onto the isoform-derived operons. Each operon is annotated with:
# - **Sense genes** — annotated genes on the *same strand* that overlap the operon coordinates; these are the genes the operon is transcribing.
# - **Antisense genes** — annotated genes on the *opposite strand* that overlap the same genomic window; these reflect the dense gene packing of the minimal genome, where antisense transcription or convergent/divergent gene pairs are common.
# 
# Running gene annotation here — before the overlap and coverage analyses — makes gene context immediately available for diagnosing same-strand operon boundary conflicts in Step 4.

# In[3]:


import re

GFF3_FILE = MOTHER_FOLDER + "/Genomes_Input/syn3a_genome.gff3"

# ── Parse GFF3 (1-based closed → 0-based half-open) ────────────────────────
gene_records = []
with open(GFF3_FILE) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        # Syn3A annotates a handful of pseudogenes (0051, 0546, 0602) only
        # under the `pseudogene` feature type — include them alongside `gene`.
        if len(parts) < 9 or parts[2] not in ("gene", "pseudogene"):
            continue
        chrom  = parts[0]
        start1 = int(parts[3])
        end1   = int(parts[4])
        strand = parts[6]
        attrs  = parts[8]
        locus  = re.search(r"locus_tag=([^;]+)", attrs)
        name   = re.search(r"Name=([^;]+)",      attrs)
        prod   = re.search(r"product=([^;]+)",   attrs)
        gene_records.append({
            "chrom":     chrom,
            "strand":    strand,
            "start0":    start1 - 1,
            "end0":      end1,
            "locus_tag": locus.group(1) if locus else "",
            "gene_name": name.group(1)  if name  else "",
            "product":   prod.group(1)  if prod  else "",
        })

genes_df = pd.DataFrame(gene_records)
print(f"Genes loaded: {len(genes_df)}")
print(genes_df["strand"].value_counts().to_string())


# ── Annotate each operon with sense AND antisense overlapping genes ─────────
def annotate_operons_with_genes(operons: pd.DataFrame,
                                genes: pd.DataFrame) -> pd.DataFrame:
    """
    For each operon add:
      sense_gene_count / sense_gene_loci / sense_gene_names
          — genes on the SAME strand that overlap the operon
      antisense_gene_count / antisense_gene_loci / antisense_gene_names
          — genes on the OPPOSITE strand that overlap the operon
    """
    ops = operons.copy()
    opp = {"+": "-", "-": "+"}
    s_counts, s_loci, s_names = [], [], []
    a_counts, a_loci, a_names = [], [], []

    for _, op in ops.iterrows():
        overlaps  = (genes["start0"] < op["end0"]) & (genes["end0"] > op["start0"])
        sense     = genes[overlaps & (genes["strand"] == op["strand"])]
        antisense = genes[overlaps & (genes["strand"] == opp[op["strand"]])]
        s_counts.append(len(sense));  s_loci.append(",".join(sense["locus_tag"]));  s_names.append(",".join(sense["gene_name"]))
        a_counts.append(len(antisense)); a_loci.append(",".join(antisense["locus_tag"])); a_names.append(",".join(antisense["gene_name"]))

    ops["sense_gene_count"]     = s_counts
    ops["sense_gene_loci"]      = s_loci
    ops["sense_gene_names"]     = s_names
    ops["antisense_gene_count"] = a_counts
    ops["antisense_gene_loci"]  = a_loci
    ops["antisense_gene_names"] = a_names
    ops["gene_count"] = ops["sense_gene_count"] + ops["antisense_gene_count"]
    ops["gene_loci"]  = ops["sense_gene_loci"]
    ops["gene_names"] = ops["sense_gene_names"]
    return ops

operons_annotated = annotate_operons_with_genes(operons_all, genes_df)

print(f"\nOperons with sense genes:")
print(f"  0 sense genes:   {(operons_annotated['sense_gene_count'] == 0).sum()}")
print(f"  1 sense gene:    {(operons_annotated['sense_gene_count'] == 1).sum()}")
print(f"  2 sense genes:   {(operons_annotated['sense_gene_count'] == 2).sum()}")
print(f"  3+ sense genes:  {(operons_annotated['sense_gene_count'] >= 3).sum()}")
print(f"\nOperons with >=1 antisense gene: {(operons_annotated['antisense_gene_count'] > 0).sum()}")


# In[4]:


operons_annotated.head(20)


# In[5]:


operons_annotated.sort_values("gene_count", ascending=False).head(20)


# ## Step 4 — Operon overlap analysis
# 
# **Same-strand overlaps** should be zero by construction from Step 2 (containment clustering guarantees non-overlapping clusters on each strand). Any same-strand overlap therefore indicates a pair of isoforms that partially overlap without containment — typically caused by RNase cleavage producing a 5′ fragment whose 3′ end falls inside the next operon's 5′ region. These pairs are resolved by the merge rule in Step 5a:
# 
# - *Merge* if the downstream operon's TSS lies **inside a gene body** → the apparent TSS is an RNase cleavage site, not a genuine promoter.
# - *Separate* if the TSS is **intergenic or within 10 nt upstream of a gene start** → a genuine promoter with a short bacterial 5′ UTR.
# 
# **Cross-strand (antisense) overlaps** are biologically expected in the compact Syn3A genome, where convergent and divergent gene pairs and naturally antisense transcription are common. These are reported for information but not resolved here.

# In[26]:


def find_overlapping_pairs(df_a: pd.DataFrame, df_b: pd.DataFrame,
                           same_strand: bool) -> pd.DataFrame:
    """
    Sweep-line overlap detection between two operon DataFrames.
    Returns one row per overlapping pair with overlap_bp and overlap fractions.
    If same_strand=True, df_a and df_b should be the same strand subset.
    """
    a = df_a.sort_values("start0").reset_index(drop=True)
    b = df_b.sort_values("start0").reset_index(drop=True)

    rows = []
    bi = 0
    for _, ra in a.iterrows():
        # advance b pointer past intervals that end before ra starts
        while bi < len(b) and b.loc[bi, "end0"] <= ra["start0"]:
            bi += 1
        j = bi
        while j < len(b) and b.loc[j, "start0"] < ra["end0"]:
            rb = b.loc[j]
            # skip self-comparison
            if same_strand and ra["operon_id"] == rb["operon_id"]:
                j += 1
                continue
            ovl = min(ra["end0"], rb["end0"]) - max(ra["start0"], rb["start0"])
            if ovl > 0:
                len_a = ra["end0"] - ra["start0"]
                len_b = rb["end0"] - rb["start0"]
                rows.append({
                    "operon_a":    ra["operon_id"],
                    "strand_a":    ra["strand"],
                    "start_a":     ra["start0"],
                    "end_a":       ra["end0"],
                    "operon_b":    rb["operon_id"],
                    "strand_b":    rb["strand"],
                    "start_b":     rb["start0"],
                    "end_b":       rb["end0"],
                    "overlap_bp":  ovl,
                    "frac_of_a":   round(ovl / len_a, 3),
                    "frac_of_b":   round(ovl / len_b, 3),
                    "relationship": (
                        "a_contains_b" if ra["start0"] <= rb["start0"] and ra["end0"] >= rb["end0"]
                        else "b_contains_a" if rb["start0"] <= ra["start0"] and rb["end0"] >= ra["end0"]
                        else "partial"
                    ),
                })
            j += 1

    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=[
        "operon_a","strand_a","start_a","end_a",
        "operon_b","strand_b","start_b","end_b",
        "overlap_bp","frac_of_a","frac_of_b","relationship"])


# In[27]:


# ── 1. Same-strand overlap (should be zero by construction) ─────────────────
# operons_annotated carries sense_gene_loci so we can show gene context for
# each conflicting pair — useful for deciding which to merge in the RNase step.
plus_ops  = operons_annotated[operons_annotated["strand"] == "+"].copy()
minus_ops = operons_annotated[operons_annotated["strand"] == "-"].copy()

ss_plus  = find_overlapping_pairs(plus_ops,  plus_ops,  same_strand=True)
ss_minus = find_overlapping_pairs(minus_ops, minus_ops, same_strand=True)
ss_plus  = ss_plus[ss_plus["operon_a"]  < ss_plus["operon_b"]]
ss_minus = ss_minus[ss_minus["operon_a"] < ss_minus["operon_b"]]
same_strand_overlaps = pd.concat([ss_plus, ss_minus], ignore_index=True)

# Annotate each pair with sense gene loci from operons_annotated
op_gene_map = operons_annotated.set_index("operon_id")["sense_gene_loci"].to_dict()
if len(same_strand_overlaps):
    same_strand_overlaps["genes_a"] = same_strand_overlaps["operon_a"].map(op_gene_map)
    same_strand_overlaps["genes_b"] = same_strand_overlaps["operon_b"].map(op_gene_map)

print(f"Same-strand overlapping operon pairs: {len(same_strand_overlaps)}")
if len(same_strand_overlaps) > 0:
    print("  -> Caused by isoforms that partially overlap without containment.")
    print("  -> These pairs are candidates for merging in the RNase-rescue step.")
    # print(same_strand_overlaps[["operon_a","operon_b","strand_a",
    #                                "start_a","end_a","start_b","end_b",
    #                                "overlap_bp","frac_of_a","frac_of_b",
    #                                "relationship","genes_a","genes_b"]])


# In[28]:


# ── Pairs where genes_a and genes_b share at least one gene locus ────────────
# A shared gene means both operons claim the same gene on the same strand —
# a direct conflict that should be resolved by merging the two operons.

if len(same_strand_overlaps):
    def shared_loci(row):
        a = set(str(row["genes_a"]).split(",")) - {"", "nan"}
        b = set(str(row["genes_b"]).split(",")) - {"", "nan"}
        shared = a & b
        return ",".join(sorted(shared)) if shared else ""

    same_strand_overlaps["shared_genes"] = same_strand_overlaps.apply(shared_loci, axis=1)
    conflicts = same_strand_overlaps[same_strand_overlaps["shared_genes"] != ""].copy()

    print(f"Pairs with shared (conflicting) genes: {len(conflicts)} / {len(same_strand_overlaps)} total same-strand overlaps")

    if len(conflicts):
        print(conflicts[["operon_a","operon_b","strand_a",
                            "start_a","end_a","start_b","end_b",
                            "overlap_bp","frac_of_a","frac_of_b",
                            "relationship","genes_a","genes_b","shared_genes"]]
                .reset_index(drop=True))

        # Save to TSV for inspection
        OUT_SS = OUT_FOLDER + "/same_strand_gene_conflicts.tsv"
        conflicts.to_csv(OUT_SS, sep="\t", index=False)
        print(f"\nSaved: {OUT_SS}")
    else:
        print("No shared genes across conflicting pairs — overlaps are between adjacent distinct operons.")
        print(same_strand_overlaps[["operon_a","operon_b","strand_a",
                                       "start_a","end_a","start_b","end_b",
                                       "overlap_bp","relationship",
                                       "genes_a","genes_b"]].reset_index(drop=True))
else:
    print("No same-strand overlaps to analyse.")


# ## Step 5a — Merge possible RNase-Processed Operons

# In[29]:


# ══════════════════════════════════════════════════════════════════════════════
# Merge same-strand overlapping operons using transcription-direction logic
#
UTR_MAX_NT = 10   # nt upstream of a gene start counted as a valid promoter position

# Rule:
#   SEPARATE (genuine promoter) if EITHER:
#     a) TSS is intergenic — not inside any same-strand gene body, OR
#     b) TSS is within UTR_MAX_NT bp upstream of the nearest downstream
#        same-strand gene start (i.e. inside the 5' UTR window).
#   MERGE (RNase artifact) only when BOTH conditions fail:
#        TSS is inside a gene body AND further than UTR_MAX_NT from any gene start.
# ══════════════════════════════════════════════════════════════════════════════

def get_downstream_tss(row):
    """
    Return (downstream_operon_id, downstream_tss, upstream_operon_id).
    Downstream = further along the direction of transcription:
      + strand → larger start0 coordinate
      − strand → smaller end0 coordinate
    """
    strand = row["strand_a"]
    if strand == "+":
        if row["start_a"] >= row["start_b"]:
            return row["operon_a"], int(row["start_a"]), row["operon_b"]
        else:
            return row["operon_b"], int(row["start_b"]), row["operon_a"]
    else:
        if row["end_a"] <= row["end_b"]:
            return row["operon_a"], int(row["end_a"]), row["operon_b"]
        else:
            return row["operon_b"], int(row["end_b"]), row["operon_a"]


def tss_is_separate(tss: int, strand: str,
                    genes_on_strand: pd.DataFrame,
                    utr_max: int = UTR_MAX_NT):
    """
    Return (separate: bool, reason_str).
    Separate (genuine promoter) if EITHER:
      a) TSS is intergenic — not inside any same-strand gene body, OR
      b) TSS is within utr_max bp upstream of the nearest downstream gene start.
    Merge (RNase artifact) only when BOTH fail.
    """
    # Condition a: intergenic
    inside = genes_on_strand[
        (genes_on_strand["start0"] <= tss) &
        (genes_on_strand["end0"]   >  tss)
    ]
    if inside.empty:
        return True, f"intergenic (TSS {tss} not inside any gene)"

    # Condition b: within utr_max nt upstream of a downstream gene start
    # For + strand: downstream gene start > tss, distance = start - tss
    # For − strand: downstream gene end  < tss, distance = tss - end
    if strand == "+":
        nearby = genes_on_strand[
            (genes_on_strand["start0"] > tss) &
            (genes_on_strand["start0"] - tss <= utr_max)
        ]
    else:
        nearby = genes_on_strand[
            (genes_on_strand["end0"] < tss) &
            (tss - genes_on_strand["end0"] <= utr_max)
        ]
    if not nearby.empty:
        locus = nearby.iloc[0]["locus_tag"]
        dist  = (int(nearby.iloc[0]["start0"]) - tss
                 if strand == "+" else
                 tss - int(nearby.iloc[0]["end0"]))
        return True, (f"UTR window: TSS {tss} is {dist} nt upstream of "
                      f"{locus} (within {utr_max} nt)")

    loci = ",".join(inside["locus_tag"].tolist())
    return False, f"TSS {tss} inside gene(s) {loci} and >{utr_max} nt from any gene start"


# ── Evaluate each same-strand overlap pair ───────────────────────────────────
genes_plus  = genes_df[genes_df["strand"] == "+"].copy()
genes_minus = genes_df[genes_df["strand"] == "-"].copy()

merge_decisions = []
for _, row in same_strand_overlaps.iterrows():
    strand = row["strand_a"]
    down_id, tss, up_id = get_downstream_tss(row)
    g_strand = genes_plus if strand == "+" else genes_minus
    is_intergenic, reason = tss_is_separate(tss, strand, g_strand)
    merge_decisions.append({
        "operon_a":       row["operon_a"],
        "operon_b":       row["operon_b"],
        "strand":         strand,
        "upstream_op":    up_id,
        "downstream_op":  down_id,
        "downstream_tss": tss,
        "overlap_bp":     row["overlap_bp"],
        "genes_a":        row.get("genes_a", ""),
        "genes_b":        row.get("genes_b", ""),
        "shared_genes":   row.get("shared_genes", ""),
        "merge":          not is_intergenic,
        "reason":         reason,
    })

decisions_df = pd.DataFrame(merge_decisions)
n_merge    = decisions_df["merge"].sum()
n_separate = len(decisions_df) - n_merge

print(f"Same-strand overlap pairs: {len(decisions_df)}")
print(f"  Separate (intergenic OR within {UTR_MAX_NT} nt of gene start): {n_separate}")
print(f"  Merge    (inside gene body AND >{UTR_MAX_NT} nt from gene start): {n_merge}")
print(decisions_df[["operon_a","operon_b","strand","upstream_op","downstream_op",
                       "downstream_tss","overlap_bp","genes_a","genes_b",
                       "shared_genes","merge","reason"]])

# ── Apply merges to operons_annotated ────────────────────────────────────────
to_merge = decisions_df[decisions_df["merge"]]

if len(to_merge) == 0:
    print("\nNo merges to apply.")
    operons_merged = operons_annotated.copy()
else:
    from collections import defaultdict as _dd

    op_ids = list(operons_annotated["operon_id"])
    idx_of = {oid: i for i, oid in enumerate(op_ids)}
    parent = list(range(len(op_ids)))

    def _find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def _union(x, y):
        parent[_find(x)] = _find(y)

    skipped = []
    for _, row in to_merge.iterrows():
        a_ok = row["operon_a"] in idx_of
        b_ok = row["operon_b"] in idx_of
        if a_ok and b_ok:
            _union(idx_of[row["operon_a"]], idx_of[row["operon_b"]])
        else:
            skipped.append((row["operon_a"], row["operon_b"],
                            "operon_a missing" if not a_ok else "operon_b missing"))
    if skipped:
        print(f"WARNING: {len(skipped)} merge pair(s) skipped — operon ID not found in operons_annotated:")
        for a, b, reason in skipped:
            print(f"  {a} + {b}: {reason}")
        print("  These operons likely have IDs from a previous run. Re-run from Step 2.")

    # Build clusters and report sizes for full transparency
    clusters = _dd(list)
    for i, oid in enumerate(op_ids):
        clusters[_find(i)].append(oid)
    cluster_sizes = [len(v) for v in clusters.values()]
    n_merged_clusters  = sum(1 for s in cluster_sizes if s >= 2)
    n_singleton_ops    = sum(1 for s in cluster_sizes if s == 1)
    pairs_accounted    = sum(s - 1 for s in cluster_sizes if s >= 2)
    print(f"Merge cluster breakdown:")
    from collections import Counter as _Ctr
    for size, count in sorted(_Ctr(cluster_sizes).items()):
        if size == 1:
            print(f"  size 1 (unchanged):  {count} operons")
        else:
            print(f"  size {size} cluster:       {count} × ({size} operons → 1 merged) "
                  f"= {count*(size-1)} pairs consumed")
    print(f"  Total pairs consumed by Union-Find: {pairs_accounted}")
    print(f"  Merge pairs requested: {len(to_merge) - len(skipped)}")
    if pairs_accounted != len(to_merge) - len(skipped):
        print(f"  WARNING: mismatch — {len(to_merge)-len(skipped)} pairs requested "
              f"but {pairs_accounted} accounted for in clusters")

    ops = operons_annotated.set_index("operon_id")
    merged_rows = []
    for root, members in clusters.items():
        if len(members) == 1:
            merged_rows.append(ops.loc[members[0]].to_dict())
        else:
            group  = ops.loc[members]
            s0     = int(group["start0"].min())
            e0     = int(group["end0"].max())
            strand = ops.loc[members[0]]["strand"]
            # Dedupe gene lists across the merged group, preserving order —
            # multiple pre-merge operons typically overlap the same gene, so a
            # plain join would repeat that locus once per pre-merge row.
            def _dedup_join(series):
                seen = set(); out = []
                for s in series.fillna(""):
                    for t in str(s).split(","):
                        t = t.strip()
                        if t and t not in seen:
                            seen.add(t); out.append(t)
                return ",".join(out)
            s_loci   = _dedup_join(group["sense_gene_loci"])
            s_names  = _dedup_join(group["sense_gene_names"])
            a_loci   = _dedup_join(group["antisense_gene_loci"])
            a_names  = _dedup_join(group["antisense_gene_names"])
            merged_rows.append({
                **ops.loc[members[0]].to_dict(),
                "start0":               s0,
                "end0":                 e0,
                "length":               e0 - s0,
                "tss":                  s0 if strand == "+" else e0,
                "tts":                  e0 if strand == "+" else s0,
                "n_isoforms":           int(group["n_isoforms"].sum()),
                "n_reads_total":        int(group["n_reads_total"].sum()),
                "sense_gene_count":     0 if not s_loci else s_loci.count(",") + 1,
                "sense_gene_loci":      s_loci,
                "sense_gene_names":     s_names,
                "antisense_gene_count": 0 if not a_loci else a_loci.count(",") + 1,
                "antisense_gene_loci":  a_loci,
                "antisense_gene_names": a_names,
                "member_ids":           ",".join(filter(None, group["member_ids"])),
                "segmentation_type":          "isoform_operon_merged",
            })

    operons_merged = (pd.DataFrame(merged_rows)
                      .sort_values(["chrom", "start0"])
                      .reset_index(drop=True))

    if "operon_id" in operons_merged.columns:
        operons_merged = operons_merged.drop(columns=["operon_id"])
    operons_merged.insert(0, "operon_id",
        [f"OP3A_{i+1:05d}" for i in range(len(operons_merged))])

    n_before = len(operons_annotated)
    n_after  = len(operons_merged)
    print(f"\nOperons before merge: {n_before}")
    print(f"Operons after merge:  {n_after}  (reduced by {n_before - n_after})")
    print(f"  = {n_singleton_ops} unchanged + {n_merged_clusters} merged clusters")
    print(f"  (isoform_operon_merged = number of merged OUTPUT operons, "
          f"not number of pairs consumed)")

    pm = operons_merged[operons_merged["strand"] == "+"]
    mm = operons_merged[operons_merged["strand"] == "-"]
    ss_check = pd.concat([
        find_overlapping_pairs(pm, pm, same_strand=True),
        find_overlapping_pairs(mm, mm, same_strand=True),
    ])
    ss_check = ss_check[ss_check["operon_a"] < ss_check["operon_b"]]
    print(f"Remaining same-strand overlaps after merge: {len(ss_check)}")


# In[30]:


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def plot_operon_conflict(row, genes_df, merged: bool = False, ax=None):
    """
    Draw a single same-strand operon conflict pair.

    Layout:
      Top track    — Operon A bar
      Middle track — gene bars
      Bottom track — Operon B bar
      Gold shading — overlap region

    merged=True  → frame and title in green  ("MERGED")
    merged=False → frame and title in orange ("KEPT SEPARATE")
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 2.2))

    strand   = row["strand_a"]
    s_a, e_a = int(row["start_a"]), int(row["end_a"])
    s_b, e_b = int(row["start_b"]), int(row["end_b"])
    ovl_s    = max(s_a, s_b)
    ovl_e    = min(e_a, e_b)

    shared   = set(str(row.get("shared_genes","")).split(",")) - {"","nan"}
    loci_a   = set(str(row.get("genes_a","")).split(","))      - {"","nan"}
    loci_b   = set(str(row.get("genes_b","")).split(","))      - {"","nan"}
    all_loci = loci_a | loci_b

    region_genes = genes_df[genes_df["locus_tag"].isin(all_loci)].copy()

    view_s = min(s_a, s_b)
    view_e = max(e_a, e_b)
    pad    = max((view_e - view_s) * 0.05, 50)
    view_s -= pad;  view_e += pad

    # Overlap shading
    ax.axvspan(ovl_s, ovl_e, color="gold", alpha=0.25, zorder=0)

    OPERON_Y = {"A": 1.0, "B": 0.0}
    OPERON_H = 0.25
    GENE_H   = 0.14
    COLORS   = {"A": "#4C72B0", "B": "#DD8452"}

    for label, (s, e) in [("A", (s_a, e_a)), ("B", (s_b, e_b))]:
        y = OPERON_Y[label]
        ax.broken_barh([(s, e - s)], (y - OPERON_H/2, OPERON_H),
                       facecolors=COLORS[label], alpha=0.85, zorder=2)
        mid = (s + e) / 2
        dx  = (e - s) * 0.08 * (1 if strand == "+" else -1)
        ax.annotate("", xy=(mid + dx, y), xytext=(mid, y),
                    arrowprops=dict(arrowstyle="->", color="white", lw=1.2), zorder=3)
        ax.text(s + (e - s) * 0.02, y + OPERON_H/2 + 0.04,
                f"Operon {label}: {row[f'operon_{label.lower()}']}  ({e-s:,} bp)",
                fontsize=7, va="bottom", color=COLORS[label], fontweight="bold")

    for _, gene in region_genes.iterrows():
        lt        = gene["locus_tag"]
        gs, ge    = int(gene["start0"]), int(gene["end0"])
        is_shared = lt in shared
        in_a, in_b = lt in loci_a, lt in loci_b
        color = "#C44E52" if is_shared else "#888888"
        ax.broken_barh([(gs, ge - gs)], (0.5 - GENE_H/2, GENE_H),
                       facecolors=color, alpha=0.75, zorder=2)
        gname = gene["gene_name"]
        claim = "A+B" if (in_a and in_b) else ("A" if in_a else "B")
        parts = [gname if gname else lt]
        if is_shared: parts.append("*")
        parts.append(f"[{claim}]")
        ax.text((gs + ge) / 2, 0.5 - GENE_H/2 - 0.06, " ".join(parts),
                fontsize=6, ha="center", va="top",
                color="#C44E52" if is_shared else "#444444")

    ax.set_xlim(view_s, view_e)
    ax.set_ylim(-0.45, 1.55)
    ax.set_xlabel("Genomic coordinate (bp)", fontsize=8)
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_yticklabels(["Operon B", "Genes", "Operon A"], fontsize=7)
    ax.tick_params(axis="x", labelsize=7)

    # ── Merge / separate annotation ──────────────────────────────────────────
    decision_color = "#2ca02c" if merged else "#ff7f0e"   # green / orange
    decision_label = "MERGED" if merged else "KEPT SEPARATE"
    reason         = row.get("reason", "")

    shared_str = f"  shared: {', '.join(sorted(shared))}" if shared else ""
    ax.set_title(
        f"Strand {strand}  |  overlap {int(row['overlap_bp'])} bp  "
        f"({row['relationship']}){shared_str}\n"
        f"{reason}",
        fontsize=7, pad=3
    )

    # Coloured border to instantly signal decision
    for spine in ax.spines.values():
        spine.set_edgecolor(decision_color)
        spine.set_linewidth(2.0)

    # Decision badge (top-left corner)
    ax.text(0.01, 0.97, decision_label,
            transform=ax.transAxes, fontsize=8, fontweight="bold",
            color="white", va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.2", facecolor=decision_color,
                      edgecolor="none", alpha=0.85))

    patches = [
        mpatches.Patch(color=COLORS["A"], label=f"Operon A ({row['operon_a']})"),
        mpatches.Patch(color=COLORS["B"], label=f"Operon B ({row['operon_b']})"),
        mpatches.Patch(color="#C44E52",   label="Shared gene"),
        mpatches.Patch(color="#888888",   label="Unique gene"),
        mpatches.Patch(color="gold",      alpha=0.5, label="Overlap region"),
    ]
    ax.legend(handles=patches, fontsize=6, loc="upper right",
              framealpha=0.7, ncol=2)
    return ax


# ── Plot all same-strand overlap pairs, annotated with merge decision ────────
if len(same_strand_overlaps) == 0:
    print("No same-strand overlaps to plot.")
else:
    # Join merge decision onto same_strand_overlaps via operon_a + operon_b key
    plot_df = same_strand_overlaps.merge(
        decisions_df[["operon_a","operon_b","merge","reason","downstream_tss",
                      "downstream_op","upstream_op"]],
        on=["operon_a","operon_b"], how="left"
    )
    plot_df["merge"] = plot_df["merge"].fillna(False)

    n   = len(plot_df)
    fig, axes = plt.subplots(n, 1, figsize=(12, 2.8 * n), constrained_layout=True)
    if n == 1:
        axes = [axes]

    for ax, (_, row) in zip(axes, plot_df.iterrows()):
        plot_operon_conflict(row, genes_df, merged=bool(row["merge"]), ax=ax)

    n_merged   = int(plot_df["merge"].sum())
    n_separate = n - n_merged
    fig.suptitle(
        f"Same-strand operon overlaps  "
        f"(MIN_READS={MIN_READS}, BOUNDARY_TOL={BOUNDARY_TOL})  |  "
        f"merged: {n_merged}  kept separate: {n_separate}",
        fontsize=10, y=1.01
    )

    OUT_FIG = OUT_FOLDER + "/same_strand_conflicts.pdf"
    fig.savefig(OUT_FIG, dpi=300, bbox_inches="tight")
    print(f"Saved: {OUT_FIG}")
    # plt.show()


# ## Step 5b — Merge operons sharing the same sense gene
# 
# A separate failure mode occurs when RNase activity is so extensive that no isoform spans the full gene — the surviving isoform clusters land on either side of the cleavage site, each claiming the same gene. These appear as two same-strand operons with identical `sense_gene_loci`. Their coordinate union recovers the true operon span.

# In[31]:


# ── Find same-strand operon pairs with identical non-empty sense_gene_loci ──
from collections import defaultdict as _dd

def merge_duplicate_sense_loci(operons: pd.DataFrame) -> pd.DataFrame:
    """
    Merge same-strand operons that share identical sense_gene_loci into one
    operon labelled segmentation_type="isoform_gene_combined".
    Uses Union-Find so transitive chains (A=B, B=C → merge all three) are handled.
    """
    ops = operons.copy().reset_index(drop=True)
    n   = len(ops)

    parent = list(range(n))
    def _find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def _union(x, y):
        parent[_find(x)] = _find(y)

    # Index by (strand, sense_gene_loci) — only non-empty loci
    loci_index = _dd(list)   # (strand, loci_str) → [row indices]
    for i, row in ops.iterrows():
        loci = str(row.get("sense_gene_loci", "") or "").strip()
        if loci:
            loci_index[(row["strand"], loci)].append(i)

    n_groups_merged = 0
    for (strand, loci), indices in loci_index.items():
        if len(indices) > 1:
            for idx in indices[1:]:
                _union(indices[0], idx)
            n_groups_merged += 1

    if n_groups_merged == 0:
        print("No duplicate sense_gene_loci found — nothing to merge.")
        return ops

    clusters = _dd(list)
    for i in range(n):
        clusters[_find(i)].append(i)

    merged_rows = []
    for root, members in clusters.items():
        if len(members) == 1:
            merged_rows.append(ops.iloc[members[0]].to_dict())
        else:
            group  = ops.iloc[members]
            s0     = int(group["start0"].min())
            e0     = int(group["end0"].max())
            strand = group.iloc[0]["strand"]
            base   = group.iloc[0].to_dict()
            base.update({
                "start0":               s0,
                "end0":                 e0,
                "length":               e0 - s0,
                "tss":                  s0 if strand == "+" else e0,
                "tts":                  e0 if strand == "+" else s0,
                "n_isoforms":           int(group["n_isoforms"].sum()),
                "n_reads_total":        int(group["n_reads_total"].sum()),
                "sense_gene_count":     int(group["sense_gene_count"].max()),
                "sense_gene_loci":      group.iloc[0]["sense_gene_loci"],
                "sense_gene_names":     group.iloc[0]["sense_gene_names"],
                "antisense_gene_count": int(group["antisense_gene_count"].max()),
                "antisense_gene_loci":  ",".join(filter(None,
                    group["antisense_gene_loci"].fillna("").tolist())),
                "antisense_gene_names": ",".join(filter(None,
                    group["antisense_gene_names"].fillna("").tolist())),
                "member_ids":           ",".join(filter(None,
                    group["member_ids"].fillna("").tolist())),
                "segmentation_type":    "isoform_gene_combined",
            })
            print(f"  Merged {len(members)} operons for gene(s) "
                  f"'{group.iloc[0]['sense_gene_loci']}' "
                  f"(strand {strand}, {s0}–{e0})")
            merged_rows.append(base)

    result = (pd.DataFrame(merged_rows)
              .sort_values(["chrom", "start0"])
              .reset_index(drop=True))
    if "operon_id" in result.columns:
        result = result.drop(columns=["operon_id"])
    result.insert(0, "operon_id",
        [f"OP3A_{i+1:05d}" for i in range(len(result))])
    return result


print("Merging operons with duplicate sense_gene_loci...")
operons_merged = merge_duplicate_sense_loci(operons_merged)

for stype in ["isoform_operon", "isoform_operon_merged", "isoform_gene_combined"]:
    n = (operons_merged["segmentation_type"] == stype).sum()
    print(f"  {stype:<30} {n}")


# In[32]:


# ── 2. Cross-strand overlap (antisense pairs) ───────────────────────────────
# Use operons_merged so cross-strand analysis reflects the post-merge state
plus_ops  = operons_merged[operons_merged["strand"] == "+"].copy()
minus_ops = operons_merged[operons_merged["strand"] == "-"].copy()
cross_strand = find_overlapping_pairs(plus_ops, minus_ops, same_strand=False)

print(f"\nCross-strand (antisense) overlapping operon pairs: {len(cross_strand)}")
if len(cross_strand) > 0:
    print(f"\nOverlap length distribution (bp):")
    print(cross_strand["overlap_bp"].describe().dropna().astype(int).to_string())
    print(f"\nRelationship breakdown:")
    print(cross_strand["relationship"].value_counts().to_string())
    print()
    print(cross_strand.sort_values("overlap_bp", ascending=False).head(20))


# ## Step 6 — Gene coverage
# 
# For each of the 496 annotated Syn3A loci (genes + pseudogenes), determine whether it is transcriptionally covered by the isoform-derived operons. Four coverage categories:
# 
# | Category | Meaning |
# |---|---|
# | `both` | Covered by operons on both strands (sense + antisense) |
# | `sense_only` | Covered only by a same-strand operon — normal transcription |
# | `antisense_only` | Only an antisense operon overlaps — the gene itself may be silent or very lowly expressed |
# | `uncovered` | No operon on either strand — candidate for Step 7 rescue |
# 
# Genes in the `uncovered` category are most likely lowly expressed (below `MIN_READS = 50`) rather than truly silent, since Syn3A is a minimal cell containing only essential genes.

# In[33]:


# ── Gene coverage: derived by inverting operons_annotated ──────────────────
def find_gene_coverage(genes: pd.DataFrame,
                       operons: pd.DataFrame) -> pd.DataFrame:
    """
    Invert sense/antisense loci columns in `operons` to map each gene back to
    its covering operons.  O(n_operons) — no per-gene operon scan needed.
    Returns genes_df augmented with:
      sense_operon_count / sense_covering_ops
      antisense_operon_count / antisense_covering_ops
      coverage_type: "both" | "sense_only" | "antisense_only" | "uncovered"
    """
    sense_map     = defaultdict(list)
    antisense_map = defaultdict(list)
    for _, op in operons.iterrows():
        oid = op["operon_id"]
        for locus in str(op.get("sense_gene_loci", "") or "").split(","):
            if locus:
                sense_map[locus].append(oid)
        for locus in str(op.get("antisense_gene_loci", "") or "").split(","):
            if locus:
                antisense_map[locus].append(oid)

    gdf = genes.copy()
    s_counts, s_ops, a_counts, a_ops, cov_types = [], [], [], [], []
    for _, gene in gdf.iterrows():
        lt = gene["locus_tag"]
        s  = sense_map.get(lt, [])
        a  = antisense_map.get(lt, [])
        s_counts.append(len(s));  s_ops.append(",".join(s))
        a_counts.append(len(a));  a_ops.append(",".join(a))
        if s and a:    cov_types.append("both")
        elif s:        cov_types.append("sense_only")
        elif a:        cov_types.append("antisense_only")
        else:          cov_types.append("uncovered")

    gdf["sense_operon_count"]     = s_counts
    gdf["sense_covering_ops"]     = s_ops
    gdf["antisense_operon_count"] = a_counts
    gdf["antisense_covering_ops"] = a_ops
    gdf["coverage_type"]          = cov_types
    return gdf

genes_covered = find_gene_coverage(genes_df, operons_merged)

n_both      = (genes_covered["coverage_type"] == "both").sum()
n_sense     = (genes_covered["coverage_type"] == "sense_only").sum()
n_anti_only = (genes_covered["coverage_type"] == "antisense_only").sum()
n_uncov     = (genes_covered["coverage_type"] == "uncovered").sum()
total       = len(genes_df)

print(f"Gene coverage summary (out of {total} genes):")
print(f"  Both sense + antisense:              {n_both:>4}  ({n_both/total*100:.1f}%)")
print(f"  Sense only:                          {n_sense:>4}  ({n_sense/total*100:.1f}%)")
print(f"  Antisense only:                      {n_anti_only:>4}  ({n_anti_only/total*100:.1f}%)")
print(f"  Uncovered (no operon either strand): {n_uncov:>4}  ({n_uncov/total*100:.1f}%)")

if n_anti_only > 0:
    print(f"\nGenes covered by antisense operon only:")
    print(genes_covered[genes_covered["coverage_type"] == "antisense_only"]
            [["locus_tag","gene_name","product","strand","start0","end0",
              "antisense_covering_ops"]].reset_index(drop=True))

uncovered_genes = genes_covered[genes_covered["coverage_type"] == "uncovered"].copy()
uncovered_genes["gene_length"] = uncovered_genes["end0"] - uncovered_genes["start0"]
print(f"\nUncovered genes by strand:")
print(uncovered_genes["strand"].value_counts().to_string())
print(f"\nUncovered gene length distribution (bp):")
print(uncovered_genes["gene_length"].describe().dropna().astype(int).to_string())
print(f"\nUncovered genes (first 30):")
print(uncovered_genes[["locus_tag","gene_name","product","chrom","strand",
                          "start0","end0","gene_length"]]
        .head(30).reset_index(drop=True))


# ## Step 7 — Rescue uncovered genes
# 
# Genes not covered by any isoform operon (Step 6 `uncovered`) are recovered through three complementary strategies:
# 
# **7a — BAM spanning read count:** Count ONT cDNA reads that fully span each uncovered gene's coordinates in the raw BAM file. A read must map from ≤ gene start to ≥ gene end — this is stricter than counting any overlapping read and avoids contamination from adjacent high-expression operons. Even 1–2 spanning reads confirms basal transcription activity.
# 
# **7b — Consecutive gene rescue:** When multiple consecutive uncovered genes on the same strand are separated by gaps ≤ 500 bp, they likely form a single lowly-expressed operon. Any isoform (no read-depth filter) that spans the entire gene group is used as evidence; the best-supported spanning isoform defines the operon boundaries.
# 
# **7c — rRNA operons (hard-coded):** The two rRNA operons (5S–23S–16S on the minus strand, at ~52–57 kb and ~340–345 kb in syn3A) are not captured by isoform clustering because rRNA is enormously abundant and heavily processed into mature subunits rather than read as full-length polycistronic transcripts. These operons are annotated directly from GFF3 coordinates.
# 
# All rescued operons are passed through `annotate_operons_with_genes` to populate `sense_gene_loci` / `antisense_gene_loci`, ensuring they are correctly counted in the final gene coverage summary.

# In[34]:


import pysam

BAM_FILE = MOTHER_FOLDER + "/Syn3A_Transcriptomics/ONT/ONT_Processing/syn3A.ONT.rep1.sorted.bam"

# ══════════════════════════════════════════════════════════════════════════════
# 4a — Count reads that span the ENTIRE gene region (reference_start ≤ gene_start
#      AND reference_end ≥ gene_end).  A simple bam.count() would include any
#      read that merely overlaps the region by even 1 bp — reads from adjacent
#      high-expression operons would inflate the count for low-expression neighbors.
#      Full-span counting is stricter and biologically more meaningful: it asks
#      "is there a transcript that covers this whole gene?"
# ══════════════════════════════════════════════════════════════════════════════
def count_spanning_reads(bam_path: str, chrom: str,
                         start0: int, end0: int) -> int:
    """Count reads whose mapped span fully covers [start0, end0)."""
    n = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom, start0, end0):
            if (not read.is_unmapped and
                    read.reference_start <= start0 and
                    read.reference_end   >= end0):
                n += 1
    return n

spanning_counts = []
for _, gene in uncovered_genes.iterrows():
    spanning_counts.append(
        count_spanning_reads(BAM_FILE, gene["chrom"],
                             gene["start0"], gene["end0"])
    )

uncovered_genes = uncovered_genes.copy()
uncovered_genes["spanning_read_count"] = spanning_counts

print("Spanning-read distribution for uncovered genes:")
bins   = [0, 1, 5, 10, 20, 50, 100, 99999]
labels = ["0", "1-4", "5-9", "10-19", "20-49", "50-99", "≥100"]
for lo, hi, lbl in zip(bins, bins[1:], labels):
    n = ((uncovered_genes["spanning_read_count"] >= lo) &
         (uncovered_genes["spanning_read_count"] < hi)).sum()
    print(f"  {lbl:>6} reads: {n} genes")


# In[ ]:


uncovered_genes.sort_values("spanning_read_count", ascending=False).head(20)


# In[ ]:


# ══════════════════════════════════════════════════════════════════════════════
# 4b — Consecutive-gene rescue: check for spanning isoforms (any depth)
# ══════════════════════════════════════════════════════════════════════════════
# Sort uncovered genes by chrom, strand, start
ug = uncovered_genes.sort_values(["chrom", "strand", "start0"]).reset_index(drop=True)

# Group consecutive genes on same strand (gap ≤ MAX_GENE_GAP bp)
MAX_GENE_GAP = 500   # bp between adjacent gene ends — tune if needed

def group_consecutive(df):
    """Return list of gene-index groups that are consecutive on same strand."""
    groups = []
    current = [0]
    for i in range(1, len(df)):
        prev, cur = df.iloc[i - 1], df.iloc[i]
        if (cur["strand"] == prev["strand"] and
                cur["chrom"] == prev["chrom"] and
                cur["start0"] - prev["end0"] <= MAX_GENE_GAP):
            current.append(i)
        else:
            groups.append(current)
            current = [i]
    groups.append(current)
    return groups

cons_groups = group_consecutive(ug)
multi_groups = [g for g in cons_groups if len(g) > 1]
print(f"\nConsecutive uncovered gene groups (gap ≤ {MAX_GENE_GAP} bp): {len(multi_groups)}")

# For each multi-gene group, check if any isoform spans the full region
# (using ALL isoforms from df, regardless of read depth)
df_all_strand = {s: df[df["strand"] == s].reset_index(drop=True)
                 for s in ["+", "-"]}

spanning_rescues = []
for grp in multi_groups:
    genes_in_grp = ug.iloc[grp]
    strand  = genes_in_grp.iloc[0]["strand"]
    chrom   = genes_in_grp.iloc[0]["chrom"]
    reg_s   = genes_in_grp["start0"].min()
    reg_e   = genes_in_grp["end0"].max()
    loci    = ",".join(genes_in_grp["locus_tag"].tolist())

    iso_s = df_all_strand[strand]
    # isoform must span: start <= reg_s+BOUNDARY_TOL AND end >= reg_e-BOUNDARY_TOL
    spanning = iso_s[
        (iso_s["start0"] <= reg_s + BOUNDARY_TOL) &
        (iso_s["end0"]   >= reg_e - BOUNDARY_TOL)
    ]
    if len(spanning) > 0:
        best = spanning.loc[spanning["n_reads"].idxmax()]
        spanning_rescues.append({
            "chrom":      chrom,
            "strand":     strand,
            "start0":     int(reg_s),
            "end0":       int(reg_e),
            "n_genes":    len(grp),
            "gene_loci":  loci,
            "best_isoform": best["isoform_id"],
            "n_reads":    int(best["n_reads"]),
            "segmentation_type": "rescue_multiple",
        })

span_df = pd.DataFrame(spanning_rescues) if spanning_rescues else pd.DataFrame()
print(f"Consecutive groups rescued by spanning isoform: {len(span_df)}")
if len(span_df):
    print(span_df)


# In[ ]:


# ══════════════════════════════════════════════════════════════════════════════
# 4c — Build rescue operons
# ══════════════════════════════════════════════════════════════════════════════
rescued_loci = set()
if len(span_df):
    for loci_str in span_df["gene_loci"]:
        rescued_loci.update(loci_str.split(","))

# rRNA genes are handled separately in 4d — exclude them here
rrna_loci = {"JCVISYN3A_0067","JCVISYN3A_0068","JCVISYN3A_0069",
             "JCVISYN3A_0532","JCVISYN3A_0533","JCVISYN3A_0534"}
rescued_loci.update(rrna_loci)

remaining_ug = ug[~ug["locus_tag"].isin(rescued_loci)].copy()

single_gene_ops = []
for _, gene in remaining_ug.iterrows():
    single_gene_ops.append({
        "chrom":         gene["chrom"],
        "strand":        gene["strand"],
        "start0":        int(gene["start0"]),
        "end0":          int(gene["end0"]),
        "n_isoforms":    0,
        "n_reads_total": int(gene["spanning_read_count"]),
        "member_ids":    "",
        "gene_count":    1,
        "gene_loci":     gene["locus_tag"],
        "gene_names":    gene["gene_name"],
        "length":        int(gene["end0"] - gene["start0"]),
        "segmentation_type":   "rescue_single",
    })

print(f"Single-gene BAM-rescue operons: {len(single_gene_ops)}")
print(f"  (of which 0 spanning reads:  "
      f"{sum(1 for g in single_gene_ops if g['n_reads_total'] == 0)})")


# In[ ]:


# ══════════════════════════════════════════════════════════════════════════════
# 4d — rRNA operons (hard-coded)
# Both rRNA clusters are on the minus strand.
# Transcription goes 16S → 23S → 5S (high→low coordinate).
# So TSS = end0 (high coord), TTS = start0 (low coord).
# ══════════════════════════════════════════════════════════════════════════════
RRNA_OPERONS = [
    {   # cluster 1: JCVISYN3A_0067(rrfA, 5S)–0068(rrlA, 23S)–0069(rrsA, 16S), minus strand
        # GFF coords (1-based): rrfA 52410-52518, rrlA 52580-55493, rrsA 55713-57247
        # 0-based: genomic min = 52409, genomic max = 57247
        "chrom": "CP016816.2", "strand": "-",
        "start0": 52409,   # genomic min (0-based) = TTS for minus strand
        "end0":   57247,   # genomic max           = TSS for minus strand
        "n_isoforms": 0, "n_reads_total": 0, "member_ids": "",
        "gene_count": 3,
        "gene_loci":  "JCVISYN3A_0067,JCVISYN3A_0068,JCVISYN3A_0069",
        "gene_names": "rrfA,rrlA,rrsA",
        "segmentation_type": "rescue_multiple",
    },
    {   # cluster 2: JCVISYN3A_0532(rrfB, 5S)–0533(rrlB, 23S)–0534(rrsB, 16S), minus strand
        # GFF coords (1-based): rrfB 340214-340322, rrlB 340384-343297, rrsB 343517-345051
        # 0-based: genomic min = 340213, genomic max = 345051
        "chrom": "CP016816.2", "strand": "-",
        "start0": 340213,
        "end0":   345051,
        "n_isoforms": 0, "n_reads_total": 0, "member_ids": "",
        "gene_count": 3,
        "gene_loci":  "JCVISYN3A_0532,JCVISYN3A_0533,JCVISYN3A_0534",
        "gene_names": "rrfB,rrlB,rrsB",
        "segmentation_type": "rescue_multiple",
    },
]

# Count spanning reads across each rRNA operon
for rop in RRNA_OPERONS:
    rop["n_reads_total"] = count_spanning_reads(
        BAM_FILE, rop["chrom"], rop["start0"], rop["end0"])

print("rRNA operons:")
for rop in RRNA_OPERONS:
    print(f"  {rop['gene_loci']}  {rop['strand']}  "
          f"{rop['start0']}–{rop['end0']}  {rop['n_reads_total']} spanning reads")


# In[ ]:


# ══════════════════════════════════════════════════════════════════════════════
# 4e — Merge rescued operons + add TSS/TTS directional coordinates
# ══════════════════════════════════════════════════════════════════════════════
rescue_rows = []

if len(span_df):
    for _, r in span_df.iterrows():
        rescue_rows.append({
            "chrom":         r["chrom"],
            "strand":        r["strand"],
            "start0":        int(r["start0"]),
            "end0":          int(r["end0"]),
            "n_isoforms":    1,
            "n_reads_total": int(r["n_reads"]),
            "member_ids":    r["best_isoform"],
            "gene_count":    int(r["n_genes"]),
            "gene_loci":     r["gene_loci"],
            "gene_names":    "",
            "length":        int(r["end0"] - r["start0"]),
            "segmentation_type":   "rescue_multiple",
        })

rescue_rows.extend(single_gene_ops)
rescue_rows.extend(RRNA_OPERONS)

rescue_df = pd.DataFrame(rescue_rows)
rescue_df["length"] = rescue_df["end0"] - rescue_df["start0"]

# ── Populate sense/antisense gene loci for find_gene_coverage ──────────────
# rescue_df rows already have gene_loci set to exactly the intended locus/loci
# (one gene for single_gene_bam, N genes for consecutive_spanning/rRNA).
# We must NOT let annotate_operons_with_genes overwrite gene_loci/gene_count
# because it would pick up any gene that genomically overlaps the coordinates,
# potentially adding unintended neighbours to a single_gene_bam operon.
# Instead, copy gene_loci directly into sense_gene_loci (all rescue operons
# are sense by construction) and compute antisense_gene_loci independently.
opp = {"+": "-", "-": "+"}
s_loci, s_names, s_counts = [], [], []
a_loci, a_names, a_counts = [], [], []
for _, op in rescue_df.iterrows():
    overlaps  = (genes_df["start0"] < op["end0"]) & (genes_df["end0"] > op["start0"])
    antisense = genes_df[overlaps & (genes_df["strand"] == opp[op["strand"]])]
    # sense: use the already-correct gene_loci, not a fresh scan
    s_loci.append(str(op.get("gene_loci", "") or ""))
    s_names.append(str(op.get("gene_names", "") or ""))
    s_counts.append(len([l for l in str(op.get("gene_loci","") or "").split(",") if l]))
    a_loci.append(",".join(antisense["locus_tag"].tolist()))
    a_names.append(",".join(antisense["gene_name"].tolist()))
    a_counts.append(len(antisense))
rescue_df["sense_gene_loci"]      = s_loci
rescue_df["sense_gene_names"]      = s_names
rescue_df["sense_gene_count"]      = s_counts
rescue_df["antisense_gene_loci"]   = a_loci
rescue_df["antisense_gene_names"]  = a_names
rescue_df["antisense_gene_count"]  = a_counts

# Tag segmentation_type on operons_merged:
# rows that went through ss_merge already have "isoform_operon_merged";
# rows that were never in a conflict pair have NaN — fill those now.
operons_merged = operons_merged.copy()
operons_merged["segmentation_type"] = (
    operons_merged.get("segmentation_type", pd.Series(dtype=str))
    .fillna("isoform_operon")
)

# Combine all operons — sort by genomic position (start0 always = min coord)
all_operons = pd.concat(
    [operons_merged, rescue_df], ignore_index=True
).sort_values(["chrom", "start0"]).reset_index(drop=True)

# Drop any pre-existing operon_id before re-assigning
if "operon_id" in all_operons.columns:
    all_operons = all_operons.drop(columns=["operon_id"])
all_operons.insert(0, "operon_id",
    [f"OP3A_{i+1:05d}" for i in range(len(all_operons))])

# ── Add TSS / TTS directional coordinates ───────────────────────────────────
# For + strand: TSS = start0 (5′ end, low coord), TTS = end0 (3′ end, high coord)
# For − strand: TSS = end0  (5′ end, high coord), TTS = start0 (3′ end, low coord)
all_operons["tss"] = np.where(all_operons["strand"] == "+",
                               all_operons["start0"], all_operons["end0"])
all_operons["tts"] = np.where(all_operons["strand"] == "+",
                               all_operons["end0"],  all_operons["start0"])

print(f"Final operon count: {len(all_operons)}")
for stype in all_operons["segmentation_type"].unique():
    n = (all_operons["segmentation_type"] == stype).sum()
    print(f"  {stype:<30} {n}")
n_nan = all_operons["segmentation_type"].isna().sum()
if n_nan:
    print(f"  WARNING: {n_nan} operons with no segmentation_type (NaN)")

print("\nRescue operons (non-isoform):")
print(all_operons[all_operons["segmentation_type"] != "isoform_operon"]
        [["operon_id","strand","tss","tts","length","gene_loci",
          "n_reads_total","segmentation_type"]].reset_index(drop=True))


# ## Update isoform memberships

# ## Step 8 — Save outputs
# 
# Write the final operon table and gene coverage files. The operon table (`operons.candidate_blocks.tsv`) includes all segmentation types:
# 
# | `segmentation_type` | Origin |
# |---|---|
# | `isoform_operon` | Directly from containment clustering (Step 2) |
# | `isoform_operon_merged` | Two overlapping isoform operons merged by TSS rule (Step 5a) |
# | `isoform_gene_combined` | Two isoform operons sharing a gene merged by duplicate-loci rule (Step 5b) |
# | `rescue_multiple` | Consecutive low-expression genes or rRNA operon (Step 7b/7c) |
# | `rescue_single` | Single lowly-expressed gene from BAM evidence (Step 7a) |
# 
# BED6 files (plus and minus strand separately) are written for direct loading into IGV.

# In[ ]:


# ── Final gene coverage check ───────────────────────────────────────────────
genes_covered_final = find_gene_coverage(genes_df, all_operons)

n_both      = (genes_covered_final["coverage_type"] == "both").sum()
n_sense     = (genes_covered_final["coverage_type"] == "sense_only").sum()
n_anti_only = (genes_covered_final["coverage_type"] == "antisense_only").sum()
n_uncov     = (genes_covered_final["coverage_type"] == "uncovered").sum()
total       = len(genes_df)

print(f"Final gene coverage ({total} genes):")
print(f"  Both sense + antisense:              {n_both:>4}  ({n_both/total*100:.1f}%)")
print(f"  Sense only:                          {n_sense:>4}  ({n_sense/total*100:.1f}%)")
print(f"  Antisense only:                      {n_anti_only:>4}  ({n_anti_only/total*100:.1f}%)")
print(f"  Uncovered (no operon either strand): {n_uncov:>4}  ({n_uncov/total*100:.1f}%)")

if n_uncov > 0:
    print(genes_covered_final[genes_covered_final["coverage_type"] == "uncovered"]
            [["locus_tag","gene_name","product","strand","start0","end0"]]
            .reset_index(drop=True))

# ── Save all outputs ─────────────────────────────────────────────────────────
OUT_TSV = OUT_FOLDER + "/operons.candidate_blocks.tsv"
all_operons.to_csv(OUT_TSV, sep="\t", index=False)
print(f"\nSaved: {OUT_TSV}  ({len(all_operons)} operons)")
print(f"Columns: {list(all_operons.columns)}")

OUT_GCOV = OUT_FOLDER + "/gene_operon_coverage.tsv"
genes_covered_final.to_csv(OUT_GCOV, sep="\t", index=False)
print(f"Saved: {OUT_GCOV}")

OUT_UNCOV = OUT_FOLDER + "/uncovered_genes.tsv"
still_uncov = genes_covered_final[genes_covered_final["coverage_type"] == "uncovered"]
still_uncov.to_csv(OUT_UNCOV, sep="\t", index=False)
print(f"Saved: {OUT_UNCOV}  ({len(still_uncov)} genes)")

# ── BED6 for IGV ─────────────────────────────────────────────────────────────
def to_bed6(df, score_col="n_reads_total", name_col="operon_id"):
    bed = df[["chrom", "start0", "end0", name_col, score_col, "strand"]].copy()
    bed.columns = ["chrom", "chromStart", "chromEnd", "name", "score", "strand"]
    bed["score"] = bed["score"].clip(upper=1000).astype(int)
    return bed

to_bed6(all_operons[all_operons["strand"] == "+"]).to_csv(
    OUT_FOLDER + "/operons.candidate_blocks.plus.bed",  sep="\t", index=False, header=False)
to_bed6(all_operons[all_operons["strand"] == "-"]).to_csv(
    OUT_FOLDER + "/operons.candidate_blocks.minus.bed", sep="\t", index=False, header=False)
print("BED files written for IGV inspection.")

