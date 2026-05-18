"""
Residual analysis of differential translation — Level 2: Translation Elongation Efficiency

Part of a three-level effort to explain the Pearson r ≈ 0.6 between transcriptome
(avg_sense_TPM) and proteome (iPM_mean) in JCVI-Syn1.0. Level 1 (initiation rate
via OSTIR) lives in Translation_Residual_L1_initiation.py; Level 3 (degradation)
lives in Translation_Residual_L3_degradation.py. A separate script summarises all
three levels together.

Objective
---------
Estimate relative elongation speed based on codon optimality and tRNA availability,
then test how much of the log10(protein) ~ log10(mRNA) residual (after controlling
for Level 1 initiation) is explained by elongation efficiency.

Method
------
- Quantify relative tRNA abundance from Illumina + PacBio RNA-seq data.

- CAI is a common codon optimality metric but relies on a predefined set of "highly expressed" genes.
Choose 20% top expressed proteins as reference set for CAI calculation.

- tAI (dos Reis, Savva & Wernisch 2004) is a codon-optimality metric that
  replaces CAI's reference-set codon counts with the actual tRNA pool decoding
  each codon, weighted by wobble efficiency. We use *measured* tRNA abundance
  (Illumina avg_sense_TPM) rather than gene copy number, since the 30 Syn1
  tRNAs span a large dynamic range.

- Anticodons: Syn1 GenBank lacks /anticodon qualifiers, so we run tRNAscan-SE
  (bacterial mode) on the genome and match each predicted anticodon back to the
  omics tRNA rows by coordinate overlap.

- Decoding matrix (62 sense codons × N tRNAs) uses the s-values from
  Maier et al., "Comprehensive quantitative modeling of translation
  efficiency in a genome-reduced bacterium" (Mol Syst Biol 2011, Mpn):

      Watson–Crick (A:U, U:A, G:C, C:G)       s = 0.00
      G34 : U3                                s = 0.41   (dos Reis default)
      U34 : G3                                s = 0.68
      I34 : U3                                s = 0.00   (Arg-ACG only)
      I34 : C3                                s = 0.28
      I34 : A3                                s = 0.99
      L34 (lysidine-CAU) : A3                 s = 0.89   (Ile2-CAU)
      4-codon-box extension (four-way wobble):
          U34 : U3                            s = 0.70
          U34 : C3                            s = 0.95

  The four-way wobble is applied only to U34 tRNAs whose amino acid has a
  single U34 decoder covering a full 4-codon box. For Syn1 that is:

      Ala-UGC, Pro-UGG, Val-UAC           (explicit in Maier et al.)
      Gly-UCC, Leu-UAG, Ser-UGA, Thr-UGU  (same logic — single U34 decoder
                                           in a 4-codon box under the
                                           Syn1 tRNA complement)

  In all split 2-codon boxes (Phe, Tyr, His, Gln, Asn, Lys, Asp, Glu, Cys,
  TTR-Leu, AGR-Arg, AGY-Ser, AUH-Ile, GAR-Glu, and the UGA→Trp Sup-UCA),
  U34 is assumed hypermodified: U34-U and U34-C pairings are forbidden, so
  e.g. Lys-UUU never mis-reads the Asn codons AAU/AAC, and the Mycoplasma
  Sup-UCA (cmnm⁵Um-hypermodified) reads only UGA/UGG and never UGU/UGC.

  Bacterial inosine at Arg-ACG (the one bacterial I34 case) is the only
  A34 → I34 conversion applied; all other A34 anticodons are kept as plain
  A34 and therefore read only the NNU codon via Watson–Crick.

- Absolute adaptiveness W_i = Σ_j (1 − s_ij) · t_j  where t_j = Illumina
  avg_sense_TPM of tRNA j. Relative w_i = W_i / max(W). Met (ATG) and Trp
  (TGG, TGA) are forced to w = 1 (single-codon / Myco-specific families).
  Stops are excluded.

- Per-gene tAI_g = exp(mean(log w)) over informative codons — same machinery
  as CAI, deferred to a later section.
- Also consider internal-SD–like sequences that can cause ribosome pausing.

"""

# =============================================================================
# 1. Imports
# =============================================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text
from scipy.stats import pearsonr, spearmanr

# =============================================================================
# 2. Paths
# =============================================================================

HOME_DIR   = ".."
OMICS_CSV  = "./syn1_genes_transcriptomics_proteomics.csv"
OPERON_TSV = HOME_DIR + "/Syn1_Operon/operons.candidate_blocks.tsv"
OUT_DIR    = "./residual_analysis"

GENOME_LEN = 1_078_809   # JCVI-Syn1.0 (CP002027.1), circular

os.makedirs(OUT_DIR, exist_ok=True)

# =============================================================================
# 3. Load omics table and extract tRNAs
# =============================================================================

print("Loading transcriptomics + proteomics table ...")
omics = pd.read_csv(OMICS_CSV)
trna = omics[omics["rna_type"] == "tRNA"].copy()
trna["aa"] = trna["gene_product"].str.replace("tRNA-", "", regex=False)
trna["locus_num"] = trna["locus_tag"].str.extract(r"_(\d+)$")[0]
trna["mid0"] = ((trna["start0"] + trna["end0"]) / 2).astype(int)
print(f"  tRNA genes: {len(trna)}")

# =============================================================================
# 4. Annotate operon membership (via sense_gene_loci in operons table)
# =============================================================================

print("Loading operons ...")
operons = pd.read_csv(OPERON_TSV, sep="\t")

# Build locus_tag → operon_id mapping from the sense_gene_loci field
locus_to_operon = {}
for _, op in operons.iterrows():
    loci = str(op.get("sense_gene_loci", "") or "")
    if not loci or loci == "nan":
        continue
    for lt in loci.split(","):
        lt = lt.strip()
        if lt:
            locus_to_operon[lt] = op["operon_id"]

trna["operon_id"] = trna["locus_tag"].map(locus_to_operon)

# Group consecutive tRNAs that belong to the same operon
trna_sorted = trna.sort_values("start0").reset_index(drop=True)
operon_trna_counts = (
    trna_sorted.dropna(subset=["operon_id"])
               .groupby("operon_id")
               .size()
               .sort_values(ascending=False)
)
print("\ntRNA-containing operons (>= 2 tRNAs):")
for op_id, n in operon_trna_counts.items():
    if n >= 2:
        members = trna_sorted[trna_sorted["operon_id"] == op_id]
        aa_list = ",".join(members["aa"].tolist())
        strand  = members["strand"].iloc[0]
        span    = f"{members['start0'].min():>7d}-{members['end0'].max():<7d}"
        print(f"  {op_id:<14s}  {strand}  {span}  n={n:2d}  {aa_list}")

trna_sorted.to_csv(f"{OUT_DIR}/trna_genome_annotation.csv", index=False)

# =============================================================================
# 5. Visualize tRNA distribution along the circular genome
#    Strategy adapted from Statistics_Analysis.ipynb:
#    - Spoke from centre to dot at the strand-specific radius.
#    - One merged text block per operon (all members on one block), singletons
#      get an individual block. Greedy polar placement with small angular +
#      radial nudges to avoid overlaps.
#    - Strand encoded by colour: forward = black, reverse = red. Forward labels
#      placed outside the backbone, reverse labels inside.
# =============================================================================

CIRCLE_R_FWD =  4.0   # forward-strand dot radius (outer ring)
CIRCLE_R_REV =  3.5   # reverse-strand dot radius (inner ring)
R_MIN        =  1.5   # crop inner empty area via set_rmin
R_MAX        =  4.5   # outer plot limit

COL_FWD = "black"
COL_REV = "#c0392b"

def pos_to_theta(pos):
    return 2 * np.pi * (pos / GENOME_LEN)

def circ_mean(theta_arr):
    s = np.sin(theta_arr).mean()
    c = np.cos(theta_arr).mean()
    return float(np.arctan2(s, c) % (2 * np.pi))

def ang_dist(a, b):
    return abs((a - b + np.pi) % (2 * np.pi) - np.pi)

fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"polar": True}, dpi=300)
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)

trna_sorted["theta"] = trna_sorted["mid0"].map(pos_to_theta)

# Dots only — no spokes to centre
dot_theta = trna_sorted["theta"].to_numpy()
dot_r     = np.where(trna_sorted["strand"].to_numpy() == "+",
                     CIRCLE_R_FWD, CIRCLE_R_REV)
dot_col   = np.where(trna_sorted["strand"].to_numpy() == "+",
                     COL_FWD, COL_REV)
dot_scatter = ax.scatter(dot_theta, dot_r, s=70, c=dot_col,
                         edgecolor="k", linewidths=0.5, zorder=3)

# --- Greedy polar-aware placer with stagger at multiple depths --------------
# (adjust_text chokes on polar axes — KDTree receives NaN bboxes — so we
# resolve overlaps ourselves: many angular slots × many radial depths, plus
# a leader line from the tRNA dot to wherever the label finally lands.)
placed = []  # (theta, r, ang_clear, r_clear)
dthetas = np.deg2rad([0, 2, -2, 4, -4, 6, -6, 8, -8, 10, -10, 13, -13, 16, -16])
drs_outer = [0.30, 0.45, 0.60, 0.75]            # outward space is tight (≤ R_MAX − 4.0)
drs_inner = [0.30, 0.55, 0.80, 1.05, 1.30, 1.55, 1.80]  # inward space is generous

def add_label(theta0, r0, text, n_lines, direction, color, fontsize):
    ang_clear = np.deg2rad(3.0) + np.deg2rad(0.7) * max(0, n_lines - 1)
    r_clear   = 0.18 + 0.11 * max(0, n_lines - 1)
    drs       = drs_outer if direction > 0 else drs_inner

    chosen = None
    for dr in drs:
        for dt in dthetas:
            th = (theta0 + dt) % (2 * np.pi)
            rr = r0 + direction * dr
            if direction < 0 and rr < R_MIN + 0.15:
                continue
            if direction > 0 and rr > R_MAX - 0.05:
                continue
            conflict = any(
                ang_dist(th, pth) < (ang_clear + pang) and abs(rr - pr) < (r_clear + prc)
                for (pth, pr, pang, prc) in placed
            )
            if conflict:
                continue
            chosen = (th, rr)
            break
        if chosen:
            break
    if chosen is None:
        chosen = (theta0, r0 + direction * drs[-1])

    th, rr = chosen
    # Leader line from the dot to the label if we moved angularly or deeply
    if ang_dist(th, theta0) > np.deg2rad(1.5) or abs((rr - r0) * direction) > drs[0] + 0.02:
        ax.plot([theta0, th], [r0, rr], color="#888888",
                linewidth=0.4, alpha=0.8, zorder=2)

    ax.text(th, rr, text, ha="center", va="center",
            fontsize=fontsize, color=color,
            bbox=dict(facecolor="white", edgecolor="none",
                      pad=0.4, alpha=0.80),
            zorder=4)
    placed.append((th, rr, ang_clear, r_clear))

# --- 1) Operon groups: ONE block per operon with all members --------------
multi_ops = set(operon_trna_counts[operon_trna_counts >= 2].index)
for op_id in multi_ops:
    g = trna_sorted[trna_sorted["operon_id"] == op_id].sort_values("mid0")
    strand = g["strand"].iloc[0]
    theta_center = circ_mean(g["theta"].to_numpy(dtype=float))

    lines = [f"{row['locus_num']}  {row['aa']}" for _, row in g.iterrows()]
    block = f"{op_id}\n" + "\n".join(lines)

    color     = COL_FWD if strand == "+" else COL_REV
    r0        = CIRCLE_R_FWD if strand == "+" else CIRCLE_R_REV
    direction = +1 if strand == "+" else -1
    add_label(theta_center, r0, block, n_lines=1 + len(lines),
              direction=direction, color=color, fontsize=10)

# --- 2) Non-operon / singleton-operon tRNAs: individual blocks -------------
iso_mask = ~trna_sorted["operon_id"].isin(multi_ops)
for _, row in trna_sorted[iso_mask].iterrows():
    block = f"{row['locus_num']}\n{row['aa']}"
    color     = COL_FWD if row["strand"] == "+" else COL_REV
    r0        = CIRCLE_R_FWD if row["strand"] == "+" else CIRCLE_R_REV
    direction = +1 if row["strand"] == "+" else -1
    add_label(row["theta"], r0, block, n_lines=2,
              direction=direction, color=color, fontsize=10)

# --- Backbone rings ---------------------------------------------------------
theta_bg = np.linspace(0, 2 * np.pi, 500)
ax.plot(theta_bg, np.full_like(theta_bg, CIRCLE_R_FWD), color="#cccccc",
        linewidth=0.8, zorder=0)
ax.plot(theta_bg, np.full_like(theta_bg, CIRCLE_R_REV), color="#e8c7c4",
        linewidth=0.8, zorder=0)

ax.set_rmin(R_MIN)
ax.set_rmax(R_MAX)

ax.grid(False)
ax.set_rticks([])
ax.set_yticklabels([])
ax.set_xticks([])
ax.spines["polar"].set_visible(False)

fwd_patch = plt.Line2D([0], [0], marker="o", color="w",
                       markerfacecolor=COL_FWD, markeredgecolor="k",
                       markersize=8, label="(+) strand tRNA")
rev_patch = plt.Line2D([0], [0], marker="o", color="w",
                       markerfacecolor=COL_REV, markeredgecolor="k",
                       markersize=8, label="(−) strand tRNA")
# ax.legend(handles=[fwd_patch, rev_patch],
#           loc="lower right", bbox_to_anchor=(1.08, -0.02),
#           frameon=False, fontsize=9)

# ax.set_title(
#     f"tRNA distribution on JCVI-Syn1.0 circular genome\n"
#     f"{len(trna_sorted)} tRNAs; outer ring = (+) strand, inner ring = (−) strand; "
#     f"operon members merged into one label block",
#     fontsize=11, pad=18,
# )
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/trna_genome_map.pdf", dpi=300)
plt.close(fig)
print("Saved: trna_genome_map.pdf")

# =============================================================================
# 6. Illumina vs PacBio TPM correlation for tRNAs
# =============================================================================

tpm = trna_sorted[["locus_tag", "locus_num", "aa", "operon_id",
                   "avg_sense_TPM", "PacBio_sense_TPM"]].copy()
tpm = tpm[(tpm["avg_sense_TPM"] > 0) & (tpm["PacBio_sense_TPM"] > 0)].copy()
tpm["log10_illumina"] = np.log10(tpm["avg_sense_TPM"])
tpm["log10_pacbio"]   = np.log10(tpm["PacBio_sense_TPM"])

r_p, p_p = pearsonr(tpm["log10_illumina"], tpm["log10_pacbio"])
r_s, p_s = spearmanr(tpm["avg_sense_TPM"], tpm["PacBio_sense_TPM"])
print(f"\ntRNA TPM correlation (n = {len(tpm)})")
print(f"  Pearson  (log10) r = {r_p:.3f}  p = {p_p:.2e}")
print(f"  Spearman         r = {r_s:.3f}  p = {p_s:.2e}")

fig, ax = plt.subplots(figsize=(6.5, 6))
ax.scatter(tpm["log10_illumina"], tpm["log10_pacbio"],
           s=38, alpha=0.8, color="#4C8BB5", edgecolors="white", linewidths=0.4)

texts = [
    ax.text(row["log10_illumina"], row["log10_pacbio"],
            f"{row['aa']}_{row['locus_num']}",
            fontsize=7, color="#444444")
    for _, row in tpm.iterrows()
]
adjust_text(
    texts,
    ax=ax,
    arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.5),
    expand_text=(1.2, 1.4),
    expand_points=(1.3, 1.5),
    force_text=(0.6, 0.8),
)

lo = min(tpm["log10_illumina"].min(), tpm["log10_pacbio"].min())
hi = max(tpm["log10_illumina"].max(), tpm["log10_pacbio"].max())
pad = 0.1 * (hi - lo)
ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
        color="#999999", linewidth=1, linestyle="--", label="y = x")

ax.set_xlabel("log₁₀(Illumina avg_sense_TPM)", fontsize=11)
ax.set_ylabel("log₁₀(PacBio_sense_TPM)", fontsize=11)
ax.set_title(
    f"tRNA expression: Illumina vs PacBio\n"
    f"Pearson r (log₁₀) = {r_p:.3f}   Spearman ρ = {r_s:.3f}   n = {len(tpm)}",
    fontsize=11,
)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, fontsize=9, loc="upper left")
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/trna_illumina_vs_pacbio_TPM.pdf", dpi=300)
plt.close(fig)
print("Saved: trna_illumina_vs_pacbio_TPM.pdf")

tpm.to_csv(f"{OUT_DIR}/trna_TPM_illumina_vs_pacbio.csv", index=False)

print("\nLevel 2 preliminary (tRNA map + TPM correlation) complete.")
print(f"Outputs in: {OUT_DIR}/")


# Construct the decoding matrix between 30 tRNA isoforms and 61 sense codons, accounting for the Mycoplasma genetic code and wobble rules. 
# This will be used to calculate elongation efficiency scores for each gene based on its codon usage and the measured tRNA abundances.

# ----------------------------
# Genetic code 4 (Mycoplasma: TGA = W)
# ----------------------------
STOP_CODONS = {"TAA", "TAG"}

GENETIC_CODE_4 = {
	"TTT":"F","TTC":"F","TTA":"L","TTG":"L",
	"TCT":"S","TCC":"S","TCA":"S","TCG":"S",
	"TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
	"TGT":"C","TGC":"C","TGA":"W","TGG":"W",
	"CTT":"L","CTC":"L","CTA":"L","CTG":"L",
	"CCT":"P","CCC":"P","CCA":"P","CCG":"P",
	"CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
	"CGT":"R","CGC":"R","CGA":"R","CGG":"R",
	"ATT":"I","ATC":"I","ATA":"I","ATG":"M",
	"ACT":"T","ACC":"T","ACA":"T","ACG":"T",
	"AAT":"N","AAC":"N","AAA":"K","AAG":"K",
	"AGT":"S","AGC":"S","AGA":"R","AGG":"R",
	"GTT":"V","GTC":"V","GTA":"V","GTG":"V",
	"GCT":"A","GCC":"A","GCA":"A","GCG":"A",
	"GAT":"D","GAC":"D","GAA":"E","GAG":"E",
	"GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

# =============================================================================
# 7. CAI (Codon Adaptation Index)
#    Reference set = top 20% protein-coding genes by iPM_mean.
#    Codon counts over reference CDSs → RSCU → relative adaptiveness w.
#    CAI_g = geometric mean of w over informative codons (exclude start, stop,
#    single-codon families Met/Trp, and codons with w undefined).
# =============================================================================

PROTEIN_FAA = HOME_DIR + "/Genomes_Input/syn1_proteins.faa"
GENOME_FA   = HOME_DIR + "/Genomes_Input/syn1_genome.fasta"

# --- 7.1 Load genome (single circular contig) -------------------------------
def load_single_fasta(path):
    seq_chunks = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq_chunks.append(line.strip().upper())
    return "".join(seq_chunks)

genome = load_single_fasta(GENOME_FA)
assert len(genome) == GENOME_LEN, f"genome length {len(genome)} != {GENOME_LEN}"

_COMP = str.maketrans("ACGTN", "TGCAN")
def revcomp(s):
    return s.translate(_COMP)[::-1]

def extract_cds(start0, end0, strand):
    """Extract CDS nucleotides, handling circular wraparound."""
    if end0 <= GENOME_LEN:
        seq = genome[start0:end0]
    else:
        seq = genome[start0:] + genome[:end0 - GENOME_LEN]
    if strand == "-":
        seq = revcomp(seq)
    return seq

# --- 7.2 Canonical protein-coding locus list from faa -----------------------
faa_locus_tags = set()
with open(PROTEIN_FAA) as fh:
    for line in fh:
        if line.startswith(">"):
            faa_locus_tags.add(line[1:].split("|", 1)[0].strip())
print(f"\nProtein-coding loci in faa: {len(faa_locus_tags)}")

# --- 7.3 CDS coordinates come from the omics table (already has start0/end0/strand)
from collections import defaultdict
mrna = omics[(omics["rna_type"] == "mRNA") & omics["locus_tag"].isin(faa_locus_tags)].copy()

cds_seqs = {}
n_bad_frame = 0
for _, row in mrna.iterrows():
    seq = extract_cds(int(row["start0"]), int(row["end0"]), row["strand"])
    if len(seq) < 6 or len(seq) % 3 != 0:
        n_bad_frame += 1
        continue
    cds_seqs[row["locus_tag"]] = seq
print(f"  CDSs extracted: {len(cds_seqs)} (skipped non-multiple-of-3: {n_bad_frame})")

# --- 7.4 Reference set: top 20% by iPM_mean ---------------------------------
prot = omics[omics["locus_tag"].isin(cds_seqs.keys())].copy()
prot = prot[prot["iPM_mean"].notna() & (prot["iPM_mean"] > 0)].copy()
prot = prot.sort_values("iPM_mean", ascending=False).reset_index(drop=True)
n_ref = max(1, int(round(0.20 * len(prot))))
ref_loci = prot["locus_tag"].iloc[:n_ref].tolist()
print(f"  CAI reference set (top 20% by iPM_mean): {n_ref}/{len(prot)} genes")

# Localization breakdown of the top-20% iPM reference set
if "ptn_localization" in prot.columns:
    ref_loc = prot["ptn_localization"].iloc[:n_ref].fillna("unknown")
    loc_frac = ref_loc.value_counts(normalize=True).sort_values(ascending=False)
    loc_cnt  = ref_loc.value_counts().reindex(loc_frac.index)
    print("  Localization fractions in top 20% iPM_mean:")
    for loc, frac in loc_frac.items():
        print(f"    {loc:<16s}: {frac*100:5.1f}%  (n = {int(loc_cnt[loc])})")

# --- 7.5 Count reference codons → RSCU → w ----------------------------------
#  Exclude start codon and stop codon; skip any premature stop codons.
SYNONYMOUS = defaultdict(list)
for codon, aa in GENETIC_CODE_4.items():
    if aa != "*":
        SYNONYMOUS[aa].append(codon)

ref_codon_counts = defaultdict(int)
for lt in ref_loci:
    seq = cds_seqs[lt]
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    if codons and codons[-1] in STOP_CODONS:
        codons = codons[:-1]
    codons = codons[1:]  # drop start
    for c in codons:
        if "N" in c or c in STOP_CODONS:
            continue
        if GENETIC_CODE_4.get(c) is None:
            continue
        ref_codon_counts[c] += 1

w = {}
for aa, codons in SYNONYMOUS.items():
    counts = np.array([ref_codon_counts[c] for c in codons], dtype=float)
    if len(codons) == 1 or aa in ("M", "W"):
        for c in codons:
            w[c] = 1.0
        continue
    if counts.max() == 0:
        for c in codons:
            w[c] = 1.0
        continue
    rel = counts / counts.max()
    # Avoid log(0): floor at 0.01 (standard Sharp & Li practice)
    rel = np.where(rel == 0, 0.01, rel)
    for c, rv in zip(codons, rel):
        w[c] = float(rv)

# --- 7.6 CAI per gene -------------------------------------------------------
def compute_cai(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    if codons and codons[-1] in STOP_CODONS:
        codons = codons[:-1]
    codons = codons[1:]
    logs = []
    for c in codons:
        if "N" in c or c in STOP_CODONS:
            continue
        aa = GENETIC_CODE_4.get(c)
        if aa is None or aa in ("M", "W"):
            continue
        logs.append(np.log(w[c]))
    if not logs:
        return np.nan, 0
    return float(np.exp(np.mean(logs))), len(logs)

cai_rows = []
for lt, seq in cds_seqs.items():
    cai, n_inf = compute_cai(seq)
    cai_rows.append({"locus_tag": lt, "length_codons": len(seq)//3,
                     "n_informative": n_inf, "CAI": cai,
                     "in_reference": lt in set(ref_loci)})
cai_df = pd.DataFrame(cai_rows)
cai_df.to_csv(f"{OUT_DIR}/gene_CAI.csv", index=False)
print(f"  gene_CAI.csv written ({len(cai_df)} genes)")

# --- 7.7 Sanity check: CAI vs iPM_mean and avg_sense_TPM --------------------
merged = omics.merge(cai_df[["locus_tag", "CAI", "in_reference"]],
                     on="locus_tag", how="inner")
ok = merged[merged["CAI"].notna() & (merged["iPM_mean"] > 0)
            & (merged["avg_sense_TPM"] > 0)].copy()
r_pi, p_pi = pearsonr(ok["CAI"], np.log10(ok["iPM_mean"]))
r_si, p_si = spearmanr(ok["CAI"], ok["iPM_mean"])
r_pt, p_pt = pearsonr(ok["CAI"], np.log10(ok["avg_sense_TPM"]))
print(f"\nCAI vs log10(iPM_mean):   Pearson r = {r_pi:.3f} (p={p_pi:.2e})   "
      f"Spearman ρ = {r_si:.3f} (p={p_si:.2e})   n={len(ok)}")
print(f"CAI vs log10(avg_TPM):    Pearson r = {r_pt:.3f} (p={p_pt:.2e})")
print(f"  mean CAI in reference set: {ok.loc[ok['in_reference'], 'CAI'].mean():.3f}")
print(f"  mean CAI outside:          {ok.loc[~ok['in_reference'], 'CAI'].mean():.3f}")

# =============================================================================
# 8. CAI vs transcriptome/proteome — scatter plots + proteome residual
#    No TIR / Level-1 coupling here.
# =============================================================================

ok["log10_TPM"]  = np.log10(ok["avg_sense_TPM"])
ok["log10_iPM"]  = np.log10(ok["iPM_mean"])

# Residual of log10(iPM) ~ log10(TPM): the part of proteome not explained by mRNA
slope, intercept = np.polyfit(ok["log10_TPM"], ok["log10_iPM"], 1)
ok["proteome_residual"] = ok["log10_iPM"] - (slope * ok["log10_TPM"] + intercept)
r_base, _ = pearsonr(ok["log10_TPM"], ok["log10_iPM"])
print(f"\nBaseline log10(iPM) ~ log10(TPM): Pearson r = {r_base:.3f}  "
      f"(R² = {r_base**2:.3f})  n = {len(ok)}")

r_res_p, p_res_p = pearsonr(ok["CAI"], ok["proteome_residual"])
r_res_s, p_res_s = spearmanr(ok["CAI"], ok["proteome_residual"])
print(f"CAI vs proteome residual:  Pearson r = {r_res_p:.3f} (p={p_res_p:.2e})  "
      f"Spearman ρ = {r_res_s:.3f} (p={p_res_s:.2e})")

# --- ΔR² from adding CAI to the baseline log10(iPM) ~ log10(TPM) OLS --------
y  = ok["log10_iPM"].to_numpy()
X1 = np.column_stack([np.ones(len(ok)), ok["log10_TPM"].to_numpy()])
X2 = np.column_stack([X1, ok["CAI"].to_numpy()])

def _r2(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1.0 - ss_res / ss_tot, beta

r2_base, _   = _r2(X1, y)
r2_aug,  b2  = _r2(X2, y)
delta_r2     = r2_aug - r2_base
print(f"\nBaseline   R² (log10 iPM ~ log10 TPM)          = {r2_base:.4f}")
print(f"Augmented  R² (log10 iPM ~ log10 TPM + CAI)    = {r2_aug:.4f}")
print(f"ΔR² from CAI                                    = {delta_r2:+.4f}  "
      f"({100*delta_r2/r2_base:+.1f}% of baseline)")
print(f"  Augmented coefficients: intercept={b2[0]:.3f}, "
      f"β_logTPM={b2[1]:.3f}, β_CAI={b2[2]:.3f}")

# --- Cytosolic-only: baseline R² and ΔR² from adding CAI --------------------
if "ptn_localization" in ok.columns:
    cyto = ok[ok["ptn_localization"] == "cytoplasmic"].copy()
    if len(cyto) > 3:
        y_c  = cyto["log10_iPM"].to_numpy()
        X1_c = np.column_stack([np.ones(len(cyto)), cyto["log10_TPM"].to_numpy()])
        X2_c = np.column_stack([X1_c, cyto["CAI"].to_numpy()])
        r2_base_c, _  = _r2(X1_c, y_c)
        r2_aug_c,  _  = _r2(X2_c, y_c)
        delta_r2_c    = r2_aug_c - r2_base_c
        print(f"\n[Cytosolic only, n = {len(cyto)}]")
        print(f"  Baseline  R² (log10 iPM ~ log10 TPM)        = {r2_base_c:.4f}")
        print(f"  Augmented R² (log10 iPM ~ log10 TPM + CAI)  = {r2_aug_c:.4f}")
        print(f"  ΔR² from CAI                                 = {delta_r2_c:+.4f}")

ok.to_csv(f"{OUT_DIR}/gene_CAI_omics_merged.csv", index=False)

# --- Three-panel figure: CAI vs TPM, CAI vs iPM, CAI vs residual -----------
fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))

def _scatter(ax, x, y, xlab, ylab, title, r, p):
    ax.scatter(x, y, s=12, alpha=0.55, color="#4C8BB5",
               edgecolors="none")
    m, b = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, m * xs + b, color="#c0392b", linewidth=1.2)
    ax.set_xlabel(xlab, fontsize=10)
    ax.set_ylabel(ylab, fontsize=10)
    ax.set_title(f"{title}\nr = {r:.3f}  p = {p:.1e}  n = {len(x)}",
                 fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)

_scatter(axes[0], ok["CAI"], ok["log10_TPM"],
         "CAI", "log₁₀(avg_sense_TPM)",
         "CAI vs transcriptome", r_pt, p_pt)
_scatter(axes[1], ok["CAI"], ok["log10_iPM"],
         "CAI", "log₁₀(iPM_mean)",
         "CAI vs proteome", r_pi, p_pi)
_scatter(axes[2], ok["CAI"], ok["proteome_residual"],
         "CAI", "log₁₀(iPM) − fit(log₁₀TPM)",
         "CAI vs proteome residual", r_res_p, p_res_p)

fig.tight_layout()
fig.savefig(f"{OUT_DIR}/CAI_vs_omics.pdf", dpi=300)
plt.close(fig)
print(f"Saved: CAI_vs_omics.pdf")

# Standalone CAI vs proteome residual (main-text panel)
# All proteins
fig, ax = plt.subplots(figsize=(3, 3))
_scatter(ax, ok["CAI"], ok["proteome_residual"],
         "CAI", "log₁₀(iPM) − fit(log₁₀TPM)",
         "", r_res_p, p_res_p)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/CAI_vs_residual_all.pdf", dpi=300)
plt.close(fig)
print("Saved: CAI_vs_residual_all.pdf")

# Cytosolic only (refit residual within cytosolic)
if "ptn_localization" in ok.columns:
    cyto_res = ok[ok["ptn_localization"] == "cytoplasmic"].copy()
    if len(cyto_res) > 3:
        sl_c, in_c = np.polyfit(cyto_res["log10_TPM"], cyto_res["log10_iPM"], 1)
        cyto_res["proteome_residual"] = (
            cyto_res["log10_iPM"] - (sl_c * cyto_res["log10_TPM"] + in_c)
        )
        r_res_c, p_res_c = pearsonr(cyto_res["CAI"], cyto_res["proteome_residual"])
        fig, ax = plt.subplots(figsize=(3, 3))
        _scatter(ax, cyto_res["CAI"], cyto_res["proteome_residual"],
                 "CAI", "log₁₀(iPM) − fit(log₁₀TPM)",
                 "", r_res_c, p_res_c)
        fig.tight_layout()
        fig.savefig(f"{OUT_DIR}/CAI_vs_residual_cytosolic.pdf", dpi=300)
        plt.close(fig)
        print("Saved: CAI_vs_residual_cytosolic.pdf")

# Bar chart: baseline vs +CAI R² (all proteins and cytosolic-only)
cyto_bar = ok[ok["ptn_localization"] == "cytoplasmic"].copy() \
    if "ptn_localization" in ok.columns else ok.iloc[0:0]
if len(cyto_bar) > 3:
    y_c  = cyto_bar["log10_iPM"].to_numpy()
    X1_c = np.column_stack([np.ones(len(cyto_bar)), cyto_bar["log10_TPM"].to_numpy()])
    X2_c = np.column_stack([X1_c, cyto_bar["CAI"].to_numpy()])
    r2_base_cyto, _ = _r2(X1_c, y_c)
    r2_aug_cyto,  _ = _r2(X2_c, y_c)

    fig, ax = plt.subplots(figsize=(3, 3))
    groups  = ["All proteins", "Cytosolic only"]
    x_pos   = np.arange(len(groups))
    width   = 0.38
    base_vals = [r2_base,     r2_base_cyto]
    aug_vals  = [r2_aug,      r2_aug_cyto]
    b1 = ax.bar(x_pos - width/2, base_vals, width,
                color="#9CA3AF", edgecolor="white", label="Baseline\n(log₁₀TPM)")
    b2 = ax.bar(x_pos + width/2, aug_vals,  width,
                color="#4C8BB5", edgecolor="white", label="+ CAI")
    for bars, vals in ((b1, base_vals), (b2, aug_vals)):
        for rect, v in zip(bars, vals):
            ax.text(rect.get_x() + rect.get_width()/2, v + 0.01,
                    f"{v:.3f}", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(groups)
    ax.set_ylabel("R²  (log₁₀ iPM model)", fontsize=10)
    ax.set_ylim(0, max(aug_vals) * 1.18)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(f"{OUT_DIR}/CAI_deltaR2_bar.pdf", dpi=300)
    plt.close(fig)
    print("Saved: CAI_deltaR2_bar.pdf")

# =============================================================================
# 9. tAI — build the 61×N decoding matrix and compute relative adaptiveness w
#    Per-gene tAI is deferred to a later section (Step 5 per the plan).
#
#    Step 1: parse tRNAscan-SE output for anticodons (run externally:
#            `tRNAscan-SE -B -o residual_analysis/syn1_trnascan.out
#             -f residual_analysis/syn1_trnascan.ss -q
#             ../Genomes_Input/syn1_genome.fasta`).
#    Step 2: match each tRNAscan row to an omics tRNA by coordinate overlap
#            and pull avg_sense_TPM as the abundance t_j.
#    Step 3: build a 61 × N_tRNA decoding matrix with s-values.
#    Step 4: W_i = Σ_j (1 − s_ij) t_j;   w_i = W_i / max(W_i).
# =============================================================================

TRNASCAN_OUT = f"{OUT_DIR}/syn1_trnascan.out"

# --- 9.1 Parse tRNAscan output ----------------------------------------------
ts_rows = []
with open(TRNASCAN_OUT) as fh:
    for ln in fh:
        p = ln.split("\t")
        if len(p) < 9 or not p[1].strip().isdigit():
            continue
        b = int(p[2]); e = int(p[3])
        start0, end0 = (b - 1, e) if b < e else (e - 1, b)
        ts_rows.append({
            "trnascan_id": int(p[1]),
            "aa": p[4].strip(),
            "anticodon": p[5].strip().upper().replace("U", "T"),
            "start0": start0,
            "end0":   end0,
            "score":  float(p[8]),
        })
ts = pd.DataFrame(ts_rows)
print(f"\n[tAI] tRNAscan-SE rows: {len(ts)}")

# --- 9.2 Match tRNAscan rows to omics tRNAs by coordinate overlap -----------
def _overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))

matches = []
for _, r in ts.iterrows():
    best = None
    for _, t in trna_sorted.iterrows():
        ov = _overlap(int(r["start0"]), int(r["end0"]),
                      int(t["start0"]), int(t["end0"]))
        if ov > 0 and (best is None or ov > best[0]):
            best = (ov, t)
    if best is None:
        continue
    t = best[1]
    matches.append({
        "locus_tag":  t["locus_tag"],
        "locus_num":  t["locus_num"],
        "aa":         r["aa"],
        "anticodon":  r["anticodon"],
        "tpm":        float(t["avg_sense_TPM"]),
        "strand":     t["strand"],
        "operon_id":  t["operon_id"],
        "ts_score":   r["score"],
    })
trna_tab = pd.DataFrame(matches).drop_duplicates("locus_tag").reset_index(drop=True)
print(f"[tAI] tRNAs matched to omics: {len(trna_tab)} / {len(trna_sorted)}")
skipped = set(ts.index) - {i for i, r in ts.iterrows()
                            if any(_overlap(int(r["start0"]), int(r["end0"]),
                                            int(t["start0"]), int(t["end0"])) > 0
                                   for _, t in trna_sorted.iterrows())}
print(f"[tAI] tRNAscan rows with no omics match (dropped): "
      f"{ts.loc[list(skipped), ['aa','anticodon','start0','end0','score']].to_dict('records') if skipped else 'none'}")

# --- 9.3 Decoding matrix (61 sense codons × N tRNAs) with wobble s-values ---
SENSE_CODONS = [c for c, aa in GENETIC_CODE_4.items() if aa != "*"]
assert len(SENSE_CODONS) == 62   # Mycoplasma code 4: TAA, TAG stop; TGA = Trp

# Anticodon (5'→3') N34-N35-N36 base-pairs antiparallel to codon (5'→3')
# N1-N2-N3, so anticodon position 34 pairs with codon position 3.
#   codon[0] (pos 1) ↔ anticodon[2] (pos 36)
#   codon[1] (pos 2) ↔ anticodon[1] (pos 35)
#   codon[2] (pos 3) ↔ anticodon[0] (pos 34)   <-- wobble position
# tRNAscan reports anticodons 5'→3' as DNA; we keep them as DNA letters.

WC = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}

# U34 tRNAs that get four-way wobble (U34:U3 = 0.70, U34:C3 = 0.95).
# Maier et al. explicitly cover Ala/Pro/Val in Mpn; extended to the other
# single-U34-decoder 4-codon boxes in the Syn1 tRNA complement.
FOURWAY_U34 = {
    ("Ala", "TGC"),
    ("Pro", "TGG"),
    ("Val", "TAC"),
    ("Gly", "TCC"),
    ("Leu", "TAG"),
    ("Ser", "TGA"),
    ("Thr", "TGT"),
}

def _trna_mods(aa, anticodon):
    """Modification flags for a tRNA: L=lysidine, I=inosine, 4=four-way U34."""
    m = set()
    if aa == "Ile2" and anticodon == "CAT":
        m.add("L")
    if aa == "Arg" and anticodon == "ACG":
        m.add("I")
    if (aa, anticodon) in FOURWAY_U34:
        m.add("4")
    return m

def _wobble_s(cpos3, ap34, mods):
    """Return the pairing s-value or None if the pairing is disallowed."""
    if "L" in mods:                             # L34 reads only A3
        return 0.89 if cpos3 == "A" else None
    if "I" in mods:                             # I34 (Arg-ACG)
        if cpos3 == "T": return 0.00
        if cpos3 == "C": return 0.28
        if cpos3 == "A": return 0.99
        return None                             # I:G disallowed
    if (ap34, cpos3) in WC:
        return 0.00
    if ap34 == "G" and cpos3 == "T":            # G34:U3 wobble
        return 0.41
    if ap34 == "T" and cpos3 == "G":            # U34:G3 wobble
        return 0.68
    if "4" in mods and ap34 == "T":             # four-way U34 wobble
        if cpos3 == "T": return 0.70
        if cpos3 == "C": return 0.95
    return None

def _pairs_13_23(codon, ac):
    """Watson-Crick requirement on codon positions 1–2 ↔ anticodon 36,35."""
    return (codon[0], ac[2]) in WC and (codon[1], ac[1]) in WC

trna_tab["mods"] = [_trna_mods(r["aa"], r["anticodon"])
                    for _, r in trna_tab.iterrows()]

N = len(trna_tab)
N_CODON = len(SENSE_CODONS)
S = np.full((N_CODON, N), np.nan)   # s-value matrix; NaN = no edge
for i, codon in enumerate(SENSE_CODONS):
    for j, row in trna_tab.iterrows():
        ac = row["anticodon"]
        if len(ac) != 3 or set(ac) - set("ACGT"):
            continue
        if not _pairs_13_23(codon, ac):
            continue
        s = _wobble_s(codon[2], ac[0], row["mods"])
        if s is None:
            continue
        S[i, j] = s

# --- 9.4 W_i = Σ_j (1 − s_ij) · t_j  and  w_i = W_i / max(W_i) --------------
t_vec = trna_tab["tpm"].to_numpy(dtype=float)
edges = ~np.isnan(S)
contrib = np.where(edges, (1.0 - np.nan_to_num(S, nan=0.0)) * t_vec[None, :], 0.0)
W = contrib.sum(axis=1)

# Quality check: every codon should have at least one decoder
no_decoder = [SENSE_CODONS[i] for i in range(N_CODON) if W[i] == 0]
print(f"[tAI] codons with no decoding tRNA: "
      f"{no_decoder if no_decoder else 'none'}")

W_max = W[W > 0].max()
w_vec = np.where(W > 0, W / W_max, np.nan)

# Force Met and Trp (including Myco UGA→Trp) to w = 1
for special in ("ATG", "TGG", "TGA"):
    w_vec[SENSE_CODONS.index(special)] = 1.0

# --- 9.5 Codon usage across all CDSs: unweighted and protein-weighted -------
# Unweighted fraction f_raw(c)   = Σ_g n_g(c)          / Σ_g Σ_c n_g(c)
# Protein-weighted   f_iPM(c)    = Σ_g iPM_g · n_g(c)  / Σ_g iPM_g · Σ_c n_g(c)
# These reflect, respectively, codon composition of the coded proteome and
# the translational *demand* for each codon (how often ribosomes actually
# decode it, given protein abundance).
from collections import Counter
gene_codon_counts = {}   # locus_tag -> Counter
for lt, seq in cds_seqs.items():
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    if codons and codons[-1] in STOP_CODONS:
        codons = codons[:-1]
    codons = codons[1:]  # drop start
    c = Counter(x for x in codons if x in GENETIC_CODE_4 and GENETIC_CODE_4[x] != "*")
    gene_codon_counts[lt] = c

raw_counts = Counter()
for c in gene_codon_counts.values():
    raw_counts.update(c)
raw_total = sum(raw_counts.values())

ipm_lookup = omics.set_index("locus_tag")["iPM_mean"].to_dict()
wtd_counts = Counter()
for lt, c in gene_codon_counts.items():
    ipm = ipm_lookup.get(lt, np.nan)
    if not (isinstance(ipm, (int, float)) and ipm > 0) or np.isnan(ipm):
        continue
    for codon, n in c.items():
        wtd_counts[codon] += ipm * n
wtd_total = sum(wtd_counts.values())

# Per-amino-acid-normalized w: divide W within each synonymous family by
# the family max. Asks "given your amino acid, how good is this codon?"
# rather than "how abundant is this codon's tRNA in absolute terms".
w_aa_vec = np.full_like(W, np.nan, dtype=float)
aa_of = [GENETIC_CODE_4[c] for c in SENSE_CODONS]
for aa_fam in set(aa_of):
    idx = [i for i, a in enumerate(aa_of) if a == aa_fam]
    Wf = W[idx]
    if Wf.max() > 0:
        w_aa_vec[idx] = Wf / Wf.max()
    else:
        w_aa_vec[idx] = np.nan
# Force Met/Trp = 1 (single-codon families under Myco code)
for special in ("ATG", "TGG", "TGA"):
    w_aa_vec[SENSE_CODONS.index(special)] = 1.0

tai_w = pd.DataFrame({
    "codon":      SENSE_CODONS,
    "aa":         aa_of,
    "W_abs":      W,
    "w_rel":      w_vec,
    "w_rel_aa":   w_aa_vec,
    "n_tRNAs":    edges.sum(axis=1),
    "count_raw":  [raw_counts[c]        for c in SENSE_CODONS],
    "frac_raw":   [raw_counts[c]/raw_total   for c in SENSE_CODONS],
    "frac_iPM":   [wtd_counts[c]/wtd_total   for c in SENSE_CODONS],
})
tai_w.to_csv(f"{OUT_DIR}/tAI_codon_weights.csv", index=False)
trna_tab.to_csv(f"{OUT_DIR}/tAI_tRNA_table.csv", index=False)

print("\n[tAI] codon weights (w_rel, rounded):")
for aa_order in "ACDEFGHIKLMNPQRSTVWY":
    row = tai_w[tai_w["aa"] == aa_order]
    parts = [f"{r['codon']}:{r['w_rel']:.2f}(n={r['n_tRNAs']})"
             for _, r in row.iterrows()]
    print(f"  {aa_order}  " + "  ".join(parts))
print(f"[tAI] Saved: tAI_codon_weights.csv, tAI_tRNA_table.csv")

# =============================================================================
# 10. Per-gene tAI, correlations with omics, and ΔR² over the baseline model
#     tAI_g = exp( mean_i log(w_i) )  over informative codons in gene g
#     (skip start, stops, Met/Trp/single-codon families — they are w=1 by
#      construction and carry no discriminating information; skip codons
#      whose w is undefined; for the orphan codon CGG use the geometric mean
#      of the other defined w's as per dos Reis 2004).
# =============================================================================

w_series = pd.Series(tai_w["w_rel"].values, index=tai_w["codon"].values)
defined = w_series[w_series.notna() & (w_series > 0)]
# Skip Met/Trp (w = 1 by construction, non-informative for ranking genes)
informative_w = defined.drop(labels=[c for c in ("ATG", "TGG", "TGA")
                                     if c in defined.index])
# Fallback for orphan CGG: geometric mean of informative w's
geo_fallback = float(np.exp(np.log(informative_w.values).mean()))

w_used = w_series.copy()
w_used["CGG"] = geo_fallback   # fill orphan

def compute_tai(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    if codons and codons[-1] in STOP_CODONS:
        codons = codons[:-1]
    codons = codons[1:]  # drop start
    logs = []
    for c in codons:
        if c in STOP_CODONS or "N" in c:
            continue
        aa = GENETIC_CODE_4.get(c)
        if aa is None or aa in ("M", "W"):
            continue
        wv = w_used.get(c)
        if wv is None or not np.isfinite(wv) or wv <= 0:
            continue
        logs.append(np.log(wv))
    if not logs:
        return np.nan, 0
    return float(np.exp(np.mean(logs))), len(logs)

tai_rows = []
for lt, seq in cds_seqs.items():
    tai, n_inf = compute_tai(seq)
    tai_rows.append({"locus_tag": lt, "length_codons": len(seq)//3,
                     "n_informative": n_inf, "tAI": tai})
tai_df = pd.DataFrame(tai_rows)
tai_df.to_csv(f"{OUT_DIR}/gene_tAI.csv", index=False)
print(f"\n[tAI] gene_tAI.csv written ({len(tai_df)} genes, "
      f"CGG fallback w = {geo_fallback:.3f})")

# --- 10.1 Merge with omics and run the same correlation suite as CAI -------
merged2 = omics.merge(tai_df[["locus_tag", "tAI"]], on="locus_tag", how="inner")
# Also pull CAI so we can do the joint model
merged2 = merged2.merge(cai_df[["locus_tag", "CAI"]], on="locus_tag", how="inner")
ok2 = merged2[merged2["tAI"].notna() & merged2["CAI"].notna()
              & (merged2["iPM_mean"] > 0) & (merged2["avg_sense_TPM"] > 0)].copy()
ok2["log10_TPM"] = np.log10(ok2["avg_sense_TPM"])
ok2["log10_iPM"] = np.log10(ok2["iPM_mean"])

slope2, int2 = np.polyfit(ok2["log10_TPM"], ok2["log10_iPM"], 1)
ok2["proteome_residual"] = ok2["log10_iPM"] - (slope2 * ok2["log10_TPM"] + int2)

rt_p, pt_p = pearsonr(ok2["tAI"], ok2["log10_TPM"])
ri_p, pi_p = pearsonr(ok2["tAI"], ok2["log10_iPM"])
ri_s, pi_s = spearmanr(ok2["tAI"], ok2["iPM_mean"])
rr_p, pr_p = pearsonr(ok2["tAI"], ok2["proteome_residual"])
rr_s, pr_s = spearmanr(ok2["tAI"], ok2["proteome_residual"])

rct, pct = pearsonr(ok2["tAI"], ok2["CAI"])

print(f"\ntAI vs log10(avg_TPM):    Pearson r = {rt_p:.3f} (p={pt_p:.2e})")
print(f"tAI vs log10(iPM_mean):   Pearson r = {ri_p:.3f} (p={pi_p:.2e})  "
      f"Spearman ρ = {ri_s:.3f}")
print(f"tAI vs proteome residual: Pearson r = {rr_p:.3f} (p={pr_p:.2e})  "
      f"Spearman ρ = {rr_s:.3f}")
print(f"tAI vs CAI (codon-metric agreement): r = {rct:.3f} (p={pct:.2e})")

# --- 10.2 ΔR² on baseline log10(iPM) ~ log10(TPM) ---------------------------
y  = ok2["log10_iPM"].to_numpy()
X1 = np.column_stack([np.ones(len(ok2)), ok2["log10_TPM"].to_numpy()])
X2_tai = np.column_stack([X1, ok2["tAI"].to_numpy()])
X2_cai = np.column_stack([X1, ok2["CAI"].to_numpy()])
X3     = np.column_stack([X1, ok2["CAI"].to_numpy(), ok2["tAI"].to_numpy()])

def _r2(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1.0 - ss_res / ss_tot, beta

r2_b,  _   = _r2(X1,      y)
r2_t,  bt  = _r2(X2_tai,  y)
r2_c,  bc  = _r2(X2_cai,  y)
r2_ct, bct = _r2(X3,      y)

print(f"\nn = {len(ok2)}")
print(f"Baseline   R² (log10 iPM ~ log10 TPM)                = {r2_b:.4f}")
print(f"+ tAI       R² (log10 iPM ~ log10 TPM + tAI)         = {r2_t:.4f}  "
      f"ΔR² = {r2_t - r2_b:+.4f}")
print(f"+ CAI       R² (log10 iPM ~ log10 TPM + CAI)         = {r2_c:.4f}  "
      f"ΔR² = {r2_c - r2_b:+.4f}")
print(f"+ CAI + tAI R² (log10 iPM ~ log10 TPM + CAI + tAI)   = {r2_ct:.4f}  "
      f"ΔR² = {r2_ct - r2_b:+.4f}")
print(f"    coeffs: int={bct[0]:.3f}  β_logTPM={bct[1]:.3f}  "
      f"β_CAI={bct[2]:.3f}  β_tAI={bct[3]:.3f}")
print(f"    ΔR² attributable to tAI *beyond CAI* = {r2_ct - r2_c:+.4f}")

# --- Cytosolic-only: repeat CAI/tAI co-correlation & joint ΔR² --------------
if "ptn_localization" in ok2.columns:
    cyto2 = ok2[ok2["ptn_localization"] == "cytoplasmic"].copy()
    if len(cyto2) > 4:
        y_c  = cyto2["log10_iPM"].to_numpy()
        X1_c = np.column_stack([np.ones(len(cyto2)), cyto2["log10_TPM"].to_numpy()])
        X2_t = np.column_stack([X1_c, cyto2["tAI"].to_numpy()])
        X2_c = np.column_stack([X1_c, cyto2["CAI"].to_numpy()])
        X3_c = np.column_stack([X1_c, cyto2["CAI"].to_numpy(), cyto2["tAI"].to_numpy()])
        r2_b_c, _   = _r2(X1_c, y_c)
        r2_t_c, _   = _r2(X2_t, y_c)
        r2_c_c, _   = _r2(X2_c, y_c)
        r2_ct_c, _  = _r2(X3_c, y_c)
        rct_c, pct_c = pearsonr(cyto2["tAI"], cyto2["CAI"])
        print(f"\n[Cytosolic only, n = {len(cyto2)}]  CAI/tAI co-correlation")
        print(f"  tAI vs CAI:                                  r = {rct_c:.3f} (p={pct_c:.2e})")
        print(f"  Baseline   R² (log10 iPM ~ log10 TPM)        = {r2_b_c:.4f}")
        print(f"  + tAI       R² = {r2_t_c:.4f}  ΔR² = {r2_t_c - r2_b_c:+.4f}")
        print(f"  + CAI       R² = {r2_c_c:.4f}  ΔR² = {r2_c_c - r2_b_c:+.4f}")
        print(f"  + CAI + tAI R² = {r2_ct_c:.4f}  ΔR² = {r2_ct_c - r2_b_c:+.4f}")
        print(f"    ΔR² of tAI beyond CAI = {r2_ct_c - r2_c_c:+.4f}")

ok2.to_csv(f"{OUT_DIR}/gene_tAI_omics_merged.csv", index=False)

# --- 10.3 Three-panel scatter: tAI vs TPM, iPM, residual --------------------
fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))
_scatter(axes[0], ok2["tAI"], ok2["log10_TPM"],
         "tAI", "log₁₀(avg_sense_TPM)",
         "tAI vs transcriptome", rt_p, pt_p)
_scatter(axes[1], ok2["tAI"], ok2["log10_iPM"],
         "tAI", "log₁₀(iPM_mean)",
         "tAI vs proteome", ri_p, pi_p)
_scatter(axes[2], ok2["tAI"], ok2["proteome_residual"],
         "tAI", "log₁₀(iPM) − fit(log₁₀TPM)",
         "tAI vs proteome residual", rr_p, pr_p)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/tAI_vs_omics.pdf", dpi=300)
plt.close(fig)
print(f"Saved: tAI_vs_omics.pdf")

# =============================================================================
# 11. Per-amino-acid-normalized tAI (sensitivity check)
#     Uses w_rel_aa (max within synonymous family) instead of the global
#     w = W/max(W). Asks "given your aa, is this the best codon?" rather
#     than "how abundant is your cognate tRNA in absolute terms?".
# =============================================================================

w_aa_series = pd.Series(tai_w["w_rel_aa"].values, index=tai_w["codon"].values)
defined_aa = w_aa_series[w_aa_series.notna() & (w_aa_series > 0)]
informative_aa = defined_aa.drop(labels=[c for c in ("ATG", "TGG", "TGA")
                                         if c in defined_aa.index])
geo_fallback_aa = float(np.exp(np.log(informative_aa.values).mean()))

w_used_aa = w_aa_series.copy()
w_used_aa["CGG"] = geo_fallback_aa

def compute_tai_aa(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    if codons and codons[-1] in STOP_CODONS:
        codons = codons[:-1]
    codons = codons[1:]
    logs = []
    for c in codons:
        if c in STOP_CODONS or "N" in c:
            continue
        aa = GENETIC_CODE_4.get(c)
        if aa is None or aa in ("M", "W"):
            continue
        wv = w_used_aa.get(c)
        if wv is None or not np.isfinite(wv) or wv <= 0:
            continue
        logs.append(np.log(wv))
    if not logs:
        return np.nan, 0
    return float(np.exp(np.mean(logs))), len(logs)

tai_aa_rows = []
for lt, seq in cds_seqs.items():
    tai, n_inf = compute_tai_aa(seq)
    tai_aa_rows.append({"locus_tag": lt, "tAI_aa": tai,
                        "n_informative_aa": n_inf})
tai_aa_df = pd.DataFrame(tai_aa_rows)

# Merge with existing per-gene tAI table and re-save
tai_df = tai_df.merge(tai_aa_df, on="locus_tag", how="left")
tai_df.to_csv(f"{OUT_DIR}/gene_tAI.csv", index=False)

# --- 11.1 Correlations and ΔR² for the per-aa variant ----------------------
ok3 = ok2.merge(tai_aa_df[["locus_tag", "tAI_aa"]], on="locus_tag", how="left")
ok3 = ok3[ok3["tAI_aa"].notna()].copy()

r_taa_tpm, p_taa_tpm = pearsonr(ok3["tAI_aa"], ok3["log10_TPM"])
r_taa_ipm, p_taa_ipm = pearsonr(ok3["tAI_aa"], ok3["log10_iPM"])
r_taa_res, p_taa_res = pearsonr(ok3["tAI_aa"], ok3["proteome_residual"])
r_taa_cai, _         = pearsonr(ok3["tAI_aa"], ok3["CAI"])
r_taa_tai, _         = pearsonr(ok3["tAI_aa"], ok3["tAI"])
print("\n[tAI_aa] per-amino-acid-normalized tAI  (n = {})".format(len(ok3)))
print(f"  vs log10(TPM):           r = {r_taa_tpm:.3f} (p={p_taa_tpm:.2e})")
print(f"  vs log10(iPM):           r = {r_taa_ipm:.3f} (p={p_taa_ipm:.2e})")
print(f"  vs proteome residual:    r = {r_taa_res:.3f} (p={p_taa_res:.2e})")
print(f"  vs CAI:                  r = {r_taa_cai:.3f}")
print(f"  vs global-normalized tAI: r = {r_taa_tai:.3f}")

y  = ok3["log10_iPM"].to_numpy()
X1 = np.column_stack([np.ones(len(ok3)), ok3["log10_TPM"].to_numpy()])
X2_taa = np.column_stack([X1, ok3["tAI_aa"].to_numpy()])
X3_cta = np.column_stack([X1, ok3["CAI"].to_numpy(), ok3["tAI_aa"].to_numpy()])
r2_b3,  _   = _r2(X1,      y)
r2_taa, btt = _r2(X2_taa,  y)
r2_cta, bct2 = _r2(X3_cta, y)
print(f"\n  Baseline R²                                  = {r2_b3:.4f}")
print(f"  + tAI_aa            R² = {r2_taa:.4f}  ΔR² = {r2_taa - r2_b3:+.4f}")
print(f"  + CAI + tAI_aa      R² = {r2_cta:.4f}  ΔR² = {r2_cta - r2_b3:+.4f}")
print(f"    ΔR² of tAI_aa beyond CAI = {r2_cta - r2_c:+.4f}")
print(f"    coeffs: int={bct2[0]:.3f}  β_logTPM={bct2[1]:.3f}  "
      f"β_CAI={bct2[2]:.3f}  β_tAI_aa={bct2[3]:.3f}")

# --- Cytosolic-only: repeat CAI/tAI_aa co-correlation & joint ΔR² ----------
if "ptn_localization" in ok3.columns:
    cyto3 = ok3[ok3["ptn_localization"] == "cytoplasmic"].copy()
    if len(cyto3) > 4:
        y_c  = cyto3["log10_iPM"].to_numpy()
        X1_c = np.column_stack([np.ones(len(cyto3)), cyto3["log10_TPM"].to_numpy()])
        X2_taa_c = np.column_stack([X1_c, cyto3["tAI_aa"].to_numpy()])
        X2_cai_c = np.column_stack([X1_c, cyto3["CAI"].to_numpy()])
        X3_cta_c = np.column_stack([X1_c, cyto3["CAI"].to_numpy(), cyto3["tAI_aa"].to_numpy()])
        r2_b_c,   _ = _r2(X1_c,     y_c)
        r2_taa_c, _ = _r2(X2_taa_c, y_c)
        r2_cai_c, _ = _r2(X2_cai_c, y_c)
        r2_cta_c, _ = _r2(X3_cta_c, y_c)
        r_ctaa_c, p_ctaa_c = pearsonr(cyto3["tAI_aa"], cyto3["CAI"])
        print(f"\n[Cytosolic only, n = {len(cyto3)}]  CAI/tAI_aa co-correlation")
        print(f"  tAI_aa vs CAI:                               r = {r_ctaa_c:.3f} (p={p_ctaa_c:.2e})")
        print(f"  Baseline      R²                            = {r2_b_c:.4f}")
        print(f"  + tAI_aa       R² = {r2_taa_c:.4f}  ΔR² = {r2_taa_c - r2_b_c:+.4f}")
        print(f"  + CAI          R² = {r2_cai_c:.4f}  ΔR² = {r2_cai_c - r2_b_c:+.4f}")
        print(f"  + CAI + tAI_aa R² = {r2_cta_c:.4f}  ΔR² = {r2_cta_c - r2_b_c:+.4f}")
        print(f"    ΔR² of tAI_aa beyond CAI = {r2_cta_c - r2_cai_c:+.4f}")

# Three-panel scatter for tAI_aa
fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))
_scatter(axes[0], ok3["tAI_aa"], ok3["log10_TPM"],
         "tAI (per-aa)", "log₁₀(avg_sense_TPM)",
         "tAI_aa vs transcriptome", r_taa_tpm, p_taa_tpm)
_scatter(axes[1], ok3["tAI_aa"], ok3["log10_iPM"],
         "tAI (per-aa)", "log₁₀(iPM_mean)",
         "tAI_aa vs proteome", r_taa_ipm, p_taa_ipm)
_scatter(axes[2], ok3["tAI_aa"], ok3["proteome_residual"],
         "tAI (per-aa)", "log₁₀(iPM) − fit(log₁₀TPM)",
         "tAI_aa vs proteome residual", r_taa_res, p_taa_res)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/tAI_aa_vs_omics.pdf", dpi=300)
plt.close(fig)
print("Saved: tAI_aa_vs_omics.pdf")

# Standalone tAI_aa vs proteome residual (SI panel)
fig, ax = plt.subplots(figsize=(3, 3))
_scatter(ax, ok3["tAI_aa"], ok3["proteome_residual"],
         "tAI (per-aa)", "log₁₀(iPM) − fit(log₁₀TPM)",
         "tAI_aa vs proteome residual", r_taa_res, p_taa_res)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/tAI_aa_vs_residual.pdf", dpi=300)
plt.close(fig)
print("Saved: tAI_aa_vs_residual.pdf")

# =============================================================================
# 11b. Codon-level: tRNA-derived w_rel vs proteome-weighted codon usage
#      One point per codon — does tRNA abundance track how heavily each
#      codon is actually used in the translated proteome?
# =============================================================================

codon_df = tai_w[(tai_w["w_rel"] > 0) & (tai_w["frac_iPM"] > 0)].copy()
codon_df["log10_w"]    = np.log10(codon_df["w_rel"])
codon_df["log10_fraw"] = np.log10(codon_df["frac_raw"])
codon_df["log10_fipm"] = np.log10(codon_df["frac_iPM"])

r_ipm, p_ipm = pearsonr(codon_df["log10_w"], codon_df["log10_fipm"])
s_ipm, q_ipm = spearmanr(codon_df["w_rel"],  codon_df["frac_iPM"])
r_raw, p_raw = pearsonr(codon_df["log10_w"], codon_df["log10_fraw"])
s_raw, q_raw = spearmanr(codon_df["w_rel"],  codon_df["frac_raw"])
print(f"\n[codon-level] w_rel vs frac_iPM   Pearson (log10) r = {r_ipm:+.3f} "
      f"(p={p_ipm:.2e})   Spearman ρ = {s_ipm:+.3f}")
print(f"[codon-level] w_rel vs frac_raw   Pearson (log10) r = {r_raw:+.3f} "
      f"(p={p_raw:.2e})   Spearman ρ = {s_raw:+.3f}")

fig, axes = plt.subplots(1, 2, figsize=(13, 6))
for ax, ycol, ylab, rr, pp, title in [
    (axes[0], "log10_fipm", "log₁₀(frac_iPM)   (protein-weighted usage)",
     r_ipm, p_ipm, "tRNA adaptiveness vs translational demand"),
    (axes[1], "log10_fraw", "log₁₀(frac_raw)   (unweighted codon usage)",
     r_raw, p_raw, "tRNA adaptiveness vs raw codon usage"),
]:
    ax.scatter(codon_df["log10_w"], codon_df[ycol],
               s=38, alpha=0.75, color="#4C8BB5",
               edgecolors="white", linewidths=0.5)
    m, b = np.polyfit(codon_df["log10_w"], codon_df[ycol], 1)
    xs = np.linspace(codon_df["log10_w"].min(),
                     codon_df["log10_w"].max(), 50)
    ax.plot(xs, m * xs + b, color="#c0392b", linewidth=1.2, label="OLS fit")
    texts = [ax.text(row["log10_w"], row[ycol],
                     f"{row['codon']}({row['aa']})",
                     fontsize=6, color="#333333")
             for _, row in codon_df.iterrows()]
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.4),
                expand_text=(1.15, 1.3), expand_points=(1.2, 1.4),
                force_text=(0.5, 0.7))
    ax.set_xlabel("log₁₀(w_rel)   (global tAI weight)", fontsize=10)
    ax.set_ylabel(ylab, fontsize=10)
    ax.set_title(f"{title}\nPearson r (log-log) = {rr:+.3f}   "
                 f"p = {pp:.1e}   n = {len(codon_df)}", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="lower right")

fig.tight_layout()
fig.savefig(f"{OUT_DIR}/codon_w_vs_usage.pdf", dpi=300)
plt.close(fig)
print("Saved: codon_w_vs_usage.pdf")

# =============================================================================
# 12. Final report — printed so it can be redirected to a .txt file for
#     downstream writing / figure captions. Run as:
#         python Translation_Residual_L2_elongation.py > Translation_Residual.txt
# =============================================================================

report = f"""
================================================================================
Level 2 — Translation Elongation Efficiency in JCVI-Syn1.0
Report generated from Translation_Residual_L2_elongation.py
================================================================================

GOAL
----
Quantify how much of the log10(iPM) ~ log10(mRNA TPM) residual (the r ≈ 0.6
transcriptome-proteome correlation in Syn1) is explained by codon-level
elongation efficiency, without touching the Level 1 initiation model.

DATA
----
  - Transcriptomics + proteomics table: syn1_genes_transcriptomics_proteomics.csv
  - Protein-coding loci:                {len(faa_locus_tags)} (from syn1_proteins.faa)
  - CDSs extracted and frame-valid:     {len(cds_seqs)}
  - tRNAs in omics table:               {len(trna_sorted)}
  - Genome:                             JCVI-Syn1.0 (CP002027.1), 1,078,809 bp, circular
  - Genetic code:                       Mycoplasma table 4 (UGA = Trp)

PRELIMINARY — tRNA distribution and expression
-----------------------------------------------
  Illumina vs PacBio tRNA TPM correlation (n = 30):
    Pearson (log10) r = 0.618   p = 2.75e-04
    Spearman        r = 0.578   p = 8.13e-04
  Measured tRNA abundance spans ~2 orders of magnitude, justifying use of
  measured tRNA TPM rather than gene copy number in tAI.

  4 multi-tRNA operons identified:
    OP_00381  (+)  n=9   Arg,Pro,Ala,Met,Met,Ser,Met,Asp,Phe
    OP_00277  (-)  n=5   Leu,Lys,Gln,Tyr,Thr
    OP_00363  (-)  n=4   Thr,Val,Glu,Asn
    OP_00341  (-)  n=2   Trp,Sec (Sup-UCA reading UGA)

CAI — codon adaptation index
----------------------------
  Reference set: top 20% of protein-coding genes by iPM_mean
                 ({n_ref} genes out of {len(prot)} with iPM data)
  Correlations (n = {len(ok)}):
    CAI vs log10(avg_sense_TPM):   Pearson r = {r_pt:+.3f}  p = {p_pt:.2e}
    CAI vs log10(iPM_mean):        Pearson r = {r_pi:+.3f}  p = {p_pi:.2e}
                                   Spearman r = {r_si:+.3f}  p = {p_si:.2e}
    CAI vs proteome residual:      Pearson r = {r_res_p:+.3f}  p = {p_res_p:.2e}
  Mean CAI in reference set: {ok.loc[ok['in_reference'], 'CAI'].mean():.3f}
  Mean CAI outside:          {ok.loc[~ok['in_reference'], 'CAI'].mean():.3f}

tAI — decoding matrix (Maier et al. 2011, Mpn s-values + sG:U = 0.41)
---------------------------------------------------------------------
  s-values used:
      Watson-Crick                                 s = 0.00
      G34 : U3                                     s = 0.41   (dos Reis default)
      U34 : G3                                     s = 0.68
      I34 : U3   I34 : C3   I34 : A3               s = 0.00, 0.28, 0.99
      L34 : A3   (lysidine, Ile2-CAU -> AUA)       s = 0.89
      Four-way wobble  U34 : U3 / U34 : C3         s = 0.70, 0.95
  Four-way wobble applied to U34 tRNAs of single-decoder 4-codon boxes:
      Ala-UGC, Pro-UGG, Val-UAC   (paper)
      Gly-UCC, Leu-UAG, Ser-UGA, Thr-UGU  (Syn1 extension, same logic)
  Inosine I34 applied only to Arg-ACG (the one bacterial A-to-I case).

  Decoder coverage of the 62 sense codons: 61 reachable, 1 orphan.
  Orphan codon: CGG (Arg) — no cognate tRNA in Syn1.
      Arg-ACG (I34) reads CGU/CGC/CGA but I:G is disallowed.
      Arg-UCU reads the AGR box, never CGN at all.
  Usage of CGG in Syn1 CDSs:  24 / 307,942 codons (0.008%), the rarest
  Arg codon by ~40x. Syn1 strongly prefers AGA (2.07%) for Arg, so the
  CGG orphan has negligible functional impact.
  All Arg codons (Syn1 coded proteome):
      CGU 0.33%   CGC 0.05%   CGA 0.05%   CGG 0.01%
      AGA 2.07%   AGG 0.06%
  For per-gene tAI, CGG is assigned w = geomean of defined w =
  {geo_fallback:.3f} (global) / {geo_fallback_aa:.3f} (per-aa), per dos Reis 2004.

  tAI codon weight distribution (w_rel, global normalization):
      Met/Trp forced to 1.00 (single-codon families under Myco code).
      Thr-ACU reaches 1.00 because the cognate tRNA dominates Illumina TPM.
      Cys, His, Gly codons drop below 0.01 — their cognate tRNAs are the
      least-expressed in the Illumina data, which crushes them under
      global w = W/max(W) normalization.

CORRELATIONS (n = {len(ok2)})
-----------------------------
                            vs log10(TPM)  vs log10(iPM)  vs residual
  CAI                         {r_pt:+.3f}         {r_pi:+.3f}         {r_res_p:+.3f}
  tAI (global W/maxW)         {rt_p:+.3f}         {ri_p:+.3f}         {rr_p:+.3f}
  tAI (per-aa family)         {r_taa_tpm:+.3f}         {r_taa_ipm:+.3f}         {r_taa_res:+.3f}

  Cross-metric:
    CAI  vs tAI (global):  r = {rct:+.3f}
    CAI  vs tAI (per-aa):  r = {r_taa_cai:+.3f}
    tAI  vs tAI (per-aa):  r = {r_taa_tai:+.3f}

REGRESSION ΔR² (dependent variable = log10(iPM_mean))
-----------------------------------------------------
  Baseline model:  log10(iPM) ~ log10(TPM)
      R² = {r2_b:.4f}   (Pearson r = {np.sqrt(r2_b):.3f}, n = {len(ok2)})

  Single-metric additions:
      + CAI                                   R² = {r2_c:.4f}   ΔR² = {r2_c - r2_b:+.4f}
      + tAI (global)                          R² = {r2_t:.4f}   ΔR² = {r2_t - r2_b:+.4f}
      + tAI (per-aa)                          R² = {r2_taa:.4f}   ΔR² = {r2_taa - r2_b3:+.4f}

  Joint additions:
      + CAI + tAI (global)                    R² = {r2_ct:.4f}   ΔR² = {r2_ct - r2_b:+.4f}
      + CAI + tAI (per-aa)                    R² = {r2_cta:.4f}   ΔR² = {r2_cta - r2_b3:+.4f}

  Best joint-model coefficients (CAI + tAI_global):
      intercept = {bct[0]:+.3f}
      β_logTPM  = {bct[1]:+.3f}
      β_CAI     = {bct[2]:+.3f}
      β_tAI     = {bct[3]:+.3f}

  ΔR² of tAI(global) beyond CAI      = {r2_ct - r2_c:+.4f}
  ΔR² of tAI(per-aa) beyond CAI      = {r2_cta - r2_c:+.4f}

INTERPRETATION
--------------
  1. CAI carries the dominant codon-level signal (+0.079 ΔR², ~23% of
     baseline R²). Proteome-based reference selection (top 20% by iPM)
     finds exactly the codons that predict protein abundance; mRNA level
     is uninformative for CAI (r ~ 0), so the ΔR² is nearly orthogonal
     to the baseline transcriptome term.

  2. Measured-tRNA-TPM tAI adds only ~+0.012 beyond CAI, regardless of
     normalization. Global-w tAI even has a negative correlation with
     mRNA TPM (r = {rt_p:+.3f}) — an artifact of a few super-abundant tRNAs
     (Thr-ACU, Met) dominating the denominator. Per-amino-acid
     normalization flips the TPM correlation to the expected sign
     (r = {r_taa_tpm:+.3f}) and roughly doubles standalone ΔR², but the
     beyond-CAI contribution stays near +0.012.

  3. CAI and tAI are essentially uncorrelated (r = {rct:+.3f} global,
     r = {r_taa_cai:+.3f} per-aa). In most bacteria CAI and tAI correlate
     strongly; the decoupling in Syn1 suggests that Illumina tRNA TPM
     is a noisy proxy for functional decoding capacity — aminoacylation
     status, anticodon modifications, and ribosome-local tRNA delivery
     are not captured by steady-state RNA-seq.

  4. For the Level 1 / Level 2 / Level 3 decomposition, CAI should be
     kept as the Level 2 elongation term. tAI can be retained as a
     sensitivity check but does not earn an independent slot in the
     explanatory model.

OUTPUTS
-------
  residual_analysis/
    trna_genome_annotation.csv
    trna_genome_map.pdf
    trna_TPM_illumina_vs_pacbio.csv
    trna_illumina_vs_pacbio_TPM.pdf
    gene_CAI.csv
    gene_CAI_omics_merged.csv
    CAI_vs_omics.pdf
    syn1_trnascan.out
    syn1_trnascan.ss
    tAI_tRNA_table.csv
    tAI_codon_weights.csv       (+ count_raw, frac_raw, frac_iPM)
    gene_tAI.csv                (both tAI and tAI_aa)
    gene_tAI_omics_merged.csv
    tAI_vs_omics.pdf
    tAI_aa_vs_omics.pdf

================================================================================
END OF LEVEL 2 REPORT
================================================================================
"""
print(report)
# with open(f"{OUT_DIR}/Translation_Residual_L2_report.txt", "w") as fh:
#     fh.write(report)
# print(f"Saved: {OUT_DIR}/Translation_Residual_L2_report.txt")