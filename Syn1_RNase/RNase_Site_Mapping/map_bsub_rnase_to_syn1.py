#!/usr/bin/env python
"""
map_bsub_rnase_to_syn1.py
=========================

Transfer experimentally-mapped *B. subtilis* RNase III and RNase Y cleavage
sites (Taggart, Charbonnier, ... Gene-Wei Li, "A high-resolution view of RNA
endonuclease cleavage in Bacillus subtilis", Mol Cell 2025) onto JCVI-Syn1 via
reciprocal-best-hit (RBH) protein homology, then predict the corresponding
Syn1 cleavage geometry:

  * RNase III -- folds the homologous Syn1 gene with ViennaRNA and scans for
    paired cut sites carrying the canonical RNase III 2-nt 3' overhang
    (RNase III is a dsRNA-specific endonuclease, so cuts must sit in a helix).

  * RNase Y -- projects each B. subtilis site by its transcript fraction into
    the Syn1 homolog, then scores the *downstream* secondary structure
    (RNase Y cleaves single-stranded RNA just 5' of a stable downstream stem).

This is the "RNase mapping" half of Syn1_RNase/End_Annotation/Peaks_Annotation.ipynb
(cells 45-76), split into a standalone, re-runnable pipeline.

Inputs (all read by absolute path):
  - Taggart SI tables (RNase III/Y cleavage sites): End_Annotation/map_RNase_Bsubtilis/
        Table S3 = RNase III high-confidence (73, has 'mRNA index')
        Table S4 = RNase Y  high-confidence (669)
        Table S5 = RNase III loose          (174)   <- used for transfer
        Table S6 = RNase Y  loose           (1477)  <- used for transfer
  - B. subtilis 168 GenBank  : End_Annotation/map_RNase_Bsubtilis/B_subtilis_168_NC_000964.3.gb
  - Syn1 GenBank (CP002027.1): Genomes_Input/syn1.gb
  - Syn1 genome FASTA        : Genomes_Input/syn1_genome.fasta

Outputs -> ./output/ (homology/, rnaseIII/, rnaseY/).

Run in the RNAseq conda env (ViennaRNA `RNA` + biopython); BLAST+ on PATH:
    /home/enguang/anaconda3/envs/RNAseq/bin/python map_bsub_rnase_to_syn1.py
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

# ----------------------------------------------------------------------------
# Paths / config
# ----------------------------------------------------------------------------
PROJECT = Path("/data/enguang/Transcriptomics/Minimal_Cells_Transcriptomics_Proteomics")

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"          # self-contained: B. subtilis inputs copied from End_Annotation/

TAGGART_XLSX = INPUTS / "TaggartSupplementalTables_241231.xlsx"
BSUB_GB = INPUTS / "B_subtilis_168_NC_000964.3.gb"

SYN1_GB = PROJECT / "Genomes_Input" / "syn1.gb"
SYN1_FASTA = PROJECT / "Genomes_Input" / "syn1_genome.fasta"

OUT = HERE / "output"
OUT_HOM = OUT / "homology"
OUT_R3 = OUT / "rnaseIII"
OUT_RY = OUT / "rnaseY"
for d in (OUT_R3, OUT_RY):          # OUT_HOM/blast created only when building locally
    d.mkdir(parents=True, exist_ok=True)

# Taggart sheet roles (loose sets feed the homology transfer)
SHEET_R3_LOOSE = "Table S5"
SHEET_RY_LOOSE = "Table S6"

# RBH homology filters -- kept IDENTICAL to Genomes_Input/Homology_Build.py (the
# canonical repo homology builder).  Same reciprocal-BLASTP + RBH method, same
# pident/coverage cutoffs, PLUS evalue<=1e-5 and bitscore>=50.  The extra cutoffs
# only drop ~5 twilight-zone RBH pairs (bitscore<50); atpA/0792 and every RNase
# III/Y candidate gene are unaffected.  Verified: this yields the same 349 strict
# RBH pairs as Genomes_Input/homology_syn1_bsub/bsub_syn1_rbh_homology_table.tsv.
MIN_PIDENT = 25.0
MIN_QCOV = 50.0
MIN_SCOV = 50.0
MIN_EVALUE = 1e-5
MIN_BITSCORE = 50.0
BLAST_EVALUE = "1e-5"
BLAST_MAX_TARGET = "5"

# Homology source.  Prefer the canonical pre-built RBH table from
# Genomes_Input/homology_syn1_bsub/ (produced by Homology_Build.py -- the repo's
# single homology source of truth) so we don't re-run BLAST or maintain a second
# homology implementation.  Falls back to building locally (identical method +
# filters above) if the canonical files are absent.
CANONICAL_HOM_DIR = PROJECT / "Genomes_Input" / "homology_syn1_bsub"
USE_PREBUILT_HOMOLOGY = True

# RNase III ViennaRNA cleavage scan (cell 61)
R3_FLANK = 0          # fold the whole homologous gene
R3_TEMPERATURE_C = 25
R3_NO_LP = True
R3_MIN_STEMRUN = 4
R3_MIN_DIST = 20
R3_MAX_DIST = 250
R3_W = 6              # neighborhood half-window for cross-pair evaluation
R3_DIRECT_OFFSET = 6
R3_OVERHANG = 2
R3_OVERHANG_WINDOW = 0
R3_TOP_K_PER_REGION = 5
R3_STORE_STRUCT = True   # keep the dot-bracket fold (needed to draw the structure later)

# Homology-anchored RNase III cut finder (Stage 4b): project each B. subtilis cut
# into Syn1 by transcript fraction, then read out the local fold at that position.
R3_ANCHOR_UP = 60            # nt upstream of the projected cut to include in the local fold
R3_ANCHOR_DOWN = 60          # nt downstream
R3_ANCHOR_PARTNER_WIN = 5    # search +/- this many nt for the cut's base-paired partner

# RNase Y projection + downstream-structure scoring (cells 74, 76)
RY_PROJECTION_HALF_WINDOW = 30
RY_KEEP_ONLY_RBH = True
RY_MULTI_HOMOLOG_POLICY = "best_bitscore"   # or "all"
RY_FOLD_UP_LEN = 20
RY_FOLD_DOWN_LEN = 70
RY_DOWNSTREAM_MIN = 5
RY_DOWNSTREAM_MAX = 50
RY_PRIMARY_WIN = (5, 44)     # focused downstream 40-mer, midpoint +24.5 (paper-matched)
RY_CENTER_WIN = (-20, 19)    # cleavage-centered 40-mer
RY_MIN_STEM_BP = 4
RY_MAX_BULGE_RUN = 3
RY_MIN_HAIRPIN_LOOP = 3
RY_MAX_HAIRPIN_LOOP = 15
RY_CUT_LOCAL_WIN = (-10, 9)
RY_MIN_DOWNSTREAM_AU = 0.50
RY_PRIMARY_MFE_STRONG = -7.0
RY_PRIMARY_DELTA_MFE_STRONG = -1.0


# ----------------------------------------------------------------------------
# Generic helpers
# ----------------------------------------------------------------------------
COMP = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def to_rna(seq_dna: str) -> str:
    return seq_dna.replace("T", "U").replace("t", "u")


def normalize_strand(x) -> str:
    s = str(x).strip()
    if s in {"+", "plus", "Plus", "1", "forward", "Forward"}:
        return "+"
    if s in {"-", "minus", "Minus", "-1", "reverse", "Reverse"}:
        return "-"
    if s in {".", "", "nan", "None"}:
        return "."
    raise ValueError(f"Unrecognized strand value: {x!r}")


def safe_int(x):
    if pd.isna(x):
        return np.nan
    try:
        return int(float(x))
    except Exception:
        return np.nan


def safe_float(x):
    if pd.isna(x):
        return np.nan
    try:
        return float(x)
    except Exception:
        return np.nan


def clean_str(x) -> str:
    return "" if pd.isna(x) else str(x).strip()


def load_fasta_as_dict(path: Path) -> dict[str, str]:
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(str(path), "fasta")}


# ----------------------------------------------------------------------------
# Stage 1 -- B. subtilis / Syn1 cleavage-site annotation (cells 51, 52)
# ----------------------------------------------------------------------------
def parse_genbank_genes(gb_path: Path) -> pd.DataFrame:
    """One row per `gene` feature, products/gene-names backfilled from CDS/RNA."""
    rows = []
    for record in SeqIO.parse(str(gb_path), "genbank"):
        seqid = record.id
        product_map: dict[str, str] = {}
        gene_name_map: dict[str, str] = {}
        for feat in record.features:
            q = feat.qualifiers
            lt = q.get("locus_tag", [""])[0]
            if not lt:
                continue
            if q.get("gene", [""])[0] and lt not in gene_name_map:
                gene_name_map[lt] = q["gene"][0]
            if q.get("product", [""])[0] and lt not in product_map:
                product_map[lt] = q["product"][0]
        for feat in record.features:
            if feat.type != "gene":
                continue
            start = int(feat.location.start) + 1        # 1-based inclusive
            end = int(feat.location.end)                # inclusive
            sv = feat.location.strand
            strand = "+" if sv == 1 else "-" if sv == -1 else "."
            q = feat.qualifiers
            gene = q.get("gene", [""])[0]
            lt = q.get("locus_tag", [""])[0]
            old_lt = q.get("old_locus_tag", [""])[0]
            if not gene and lt in gene_name_map:
                gene = gene_name_map[lt]
            product = product_map.get(lt, "")
            parts = [x for x in [gene, lt, old_lt] if x]
            name = "/".join(parts) if parts else (lt or "")
            rows.append({
                "seqid": seqid, "start": start, "end": end, "strand": strand,
                "gene": gene, "locus_tag": lt, "old_locus_tag": old_lt,
                "product": product, "name": name, "gene_len": end - start + 1,
            })
    return pd.DataFrame(rows)


def _genes_to_string(df: pd.DataFrame, field: str = "name") -> str:
    vals = [str(x) for x in df[field].fillna("").tolist() if str(x) and str(x) != "."]
    seen, uniq = set(), []
    for v in vals:
        if v not in seen:
            uniq.append(v); seen.add(v)
    return ",".join(uniq)


def _nearest_flanking(genes_df, pos, strand, field):
    sub = genes_df[genes_df["strand"] == strand] if strand in {"+", "-"} else genes_df
    left = sub[sub["end"] < pos]
    right = sub[sub["start"] > pos]
    l = str(left.sort_values("end").iloc[-1][field]) if not left.empty else ""
    r = str(right.sort_values("start").iloc[0][field]) if not right.empty else ""
    return l, r


def _position_within_gene(pos, start, end, strand) -> dict:
    glen = end - start + 1
    if glen <= 0:
        return dict(offset_from_5prime=np.nan, offset_from_3prime=np.nan, fraction_from_5prime=np.nan)
    if strand == "+":
        o5, o3 = pos - start, end - pos
    elif strand == "-":
        o5, o3 = end - pos, pos - start
    else:
        o5 = o3 = np.nan
    return dict(offset_from_5prime=o5, offset_from_3prime=o3,
                fraction_from_5prime=(o5 / glen if pd.notna(o5) else np.nan))


def annotate_and_summarize(enzyme: str, sheet: str, site_prefix: str,
                           genes_use: pd.DataFrame,
                           out_sites: Path, out_genes: Path):
    """Annotate Taggart cleavage sites to B. subtilis genes (sense-intragenic
    prioritized) and collapse to a per-gene summary.  Shared by RNase III/Y."""
    sites = pd.read_excel(TAGGART_XLSX, sheet_name=sheet, skiprows=1)
    required = ["Position of identified 5' end", "Strand", "Local sequence",
                "Endonuclease sensitivity score"]
    missing = [c for c in required if c not in sites.columns]
    if missing:
        raise ValueError(f"{enzyme} {sheet} missing columns: {missing}")

    sites = sites.copy()
    sites["Position of identified 5' end"] = pd.to_numeric(sites["Position of identified 5' end"], errors="coerce")
    sites["Endonuclease sensitivity score"] = pd.to_numeric(sites["Endonuclease sensitivity score"], errors="coerce")
    sites = sites.dropna(subset=["Position of identified 5' end", "Strand"]).reset_index(drop=True)
    sites["Position of identified 5' end"] = sites["Position of identified 5' end"].astype(int)
    sites["Strand"] = sites["Strand"].map(normalize_strand)
    sites["site_id"] = [f"{site_prefix}{i}" for i in range(len(sites))]

    annot_rows = []
    for _, row in sites.iterrows():
        pos = int(row["Position of identified 5' end"])
        strand = row["Strand"]
        overlap = genes_use[(genes_use["start"] <= pos) & (genes_use["end"] >= pos)]
        same = overlap[overlap["strand"] == strand]
        opp = overlap[overlap["strand"] != strand]
        if not same.empty:
            p = same.sort_values(["start", "end"]).iloc[0]; context = "sense_intragenic"
        elif not opp.empty:
            p = opp.sort_values(["start", "end"]).iloc[0]; context = "antisense_intragenic"
        else:
            p = None; context = "intergenic"
        if p is not None:
            pinfo = _position_within_gene(pos, int(p["start"]), int(p["end"]), p["strand"])
            prec = dict(primary_gene=p["gene"], primary_locus_tag=p["locus_tag"],
                        primary_old_locus_tag=p["old_locus_tag"], primary_product=p["product"],
                        primary_name=p["name"], primary_gene_start=int(p["start"]),
                        primary_gene_end=int(p["end"]), primary_gene_strand=p["strand"],
                        primary_gene_length=int(p["gene_len"]))
        else:
            pinfo = dict(offset_from_5prime=np.nan, offset_from_3prime=np.nan, fraction_from_5prime=np.nan)
            prec = dict(primary_gene="", primary_locus_tag="", primary_old_locus_tag="",
                        primary_product="", primary_name="", primary_gene_start=np.nan,
                        primary_gene_end=np.nan, primary_gene_strand="", primary_gene_length=np.nan)
        lg, rg = _nearest_flanking(genes_use, pos, strand, "name")
        lt, rt = _nearest_flanking(genes_use, pos, strand, "locus_tag")
        annot_rows.append({
            "site_id": row["site_id"], "pos": pos, "strand": strand,
            "local_sequence": row["Local sequence"], "sensitivity_score": row["Endonuclease sensitivity score"],
            "context": context, **prec,
            "offset_from_primary_gene_5prime": pinfo["offset_from_5prime"],
            "offset_from_primary_gene_3prime": pinfo["offset_from_3prime"],
            "fraction_from_primary_gene_5prime": pinfo["fraction_from_5prime"],
            "same_strand_overlapping_genes": _genes_to_string(same, "name"),
            "same_strand_overlapping_locus_tags": _genes_to_string(same, "locus_tag"),
            "opposite_strand_overlapping_genes": _genes_to_string(opp, "name"),
            "left_flanking_gene_same_strand": lg, "right_flanking_gene_same_strand": rg,
            "left_flanking_locus_tag_same_strand": lt, "right_flanking_locus_tag_same_strand": rt,
        })
    annot = pd.DataFrame(annot_rows).sort_values(["strand", "pos"]).reset_index(drop=True)
    annot.to_csv(out_sites, sep="\t", index=False)

    sense = annot[annot["context"] == "sense_intragenic"].copy()
    gene_rows = []
    for (pname, ploc, pgene, strand), g in sense.groupby(
            ["primary_name", "primary_locus_tag", "primary_gene", "strand"], dropna=False):
        g = g.sort_values("pos")
        positions = g["pos"].astype(int).tolist()
        sc = pd.to_numeric(g["sensitivity_score"], errors="coerce")
        frac = pd.to_numeric(g["fraction_from_primary_gene_5prime"], errors="coerce")
        gene_rows.append({
            "primary_name": pname, "primary_gene": pgene, "primary_locus_tag": ploc, "strand": strand,
            "n_candidate_sites": len(g),
            "site_ids": ",".join(map(str, g["site_id"].tolist())),
            "site_positions": ",".join(map(str, positions)),
            "min_position": min(positions), "max_position": max(positions),
            "mean_sensitivity_score": float(sc.mean()) if sc.notna().any() else np.nan,
            "median_sensitivity_score": float(sc.median()) if sc.notna().any() else np.nan,
            "max_sensitivity_score": float(sc.max()) if sc.notna().any() else np.nan,
            "mean_fraction_from_5prime": float(frac.mean()) if frac.notna().any() else np.nan,
            "min_fraction_from_5prime": float(frac.min()) if frac.notna().any() else np.nan,
            "max_fraction_from_5prime": float(frac.max()) if frac.notna().any() else np.nan,
            "local_sequences": "|".join(g["local_sequence"].astype(str).tolist()),
        })
    genes = pd.DataFrame(gene_rows)
    if not genes.empty:
        genes = genes.sort_values(["n_candidate_sites", "max_sensitivity_score"],
                                  ascending=[False, False]).reset_index(drop=True)
    genes.to_csv(out_genes, sep="\t", index=False)
    print(f"[{enzyme}] sites={len(annot)}  sense-intragenic={len(sense)}  genes={len(genes)}")
    return annot, genes


# ----------------------------------------------------------------------------
# Stage 2 -- reciprocal-best-hit protein homology (cells 55, 56)
# ----------------------------------------------------------------------------
def parse_genbank_proteins(gb_path: Path):
    fasta_records, ann_rows = [], []
    id_counts: dict[str, int] = {}          # disambiguate duplicate primary_ids (matches Homology_Build.py)
    for record in SeqIO.parse(str(gb_path), "genbank"):
        seqid = record.id
        for feat in record.features:
            if feat.type != "CDS":
                continue
            q = feat.qualifiers
            translation = q.get("translation", [""])[0]
            if not translation:
                continue
            lt = q.get("locus_tag", [""])[0]
            gene = q.get("gene", [""])[0]
            old_lt = q.get("old_locus_tag", [""])[0]
            protein_id = q.get("protein_id", [""])[0]
            product = q.get("product", [""])[0]
            s1 = int(feat.location.start) + 1
            e1 = int(feat.location.end)
            strand = "+" if feat.location.strand == 1 else "-" if feat.location.strand == -1 else "."
            base_id = lt or protein_id or f"{seqid}_{s1}_{e1}"
            n_seen = id_counts.get(base_id, 0)
            primary_id = base_id if n_seen == 0 else f"{base_id}_{n_seen + 1}"
            id_counts[base_id] = n_seen + 1
            header = f"{primary_id}|{gene}|{seqid}"
            fasta_records.append((header, translation))
            ann_rows.append(dict(primary_id=primary_id, locus_tag=lt, gene=gene,
                                 old_locus_tag=old_lt, protein_id=protein_id, product=product,
                                 seqid=seqid, start_1b=s1, end_1b=e1, strand=strand,
                                 protein_len=len(translation), fasta_header=header))
    return fasta_records, pd.DataFrame(ann_rows)


def write_faa(records, out_faa: Path):
    with out_faa.open("w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def build_homology_table() -> pd.DataFrame:
    """Parse proteins, run reciprocal BLASTP, return strict RBH homology table."""
    (OUT_HOM / "blast").mkdir(parents=True, exist_ok=True)
    bsub_records, bsub_ann = parse_genbank_proteins(BSUB_GB)
    syn1_records, syn1_ann = parse_genbank_proteins(SYN1_GB)
    bsub_faa = OUT_HOM / "bsub_proteins.faa"
    syn1_faa = OUT_HOM / "syn1_proteins.faa"
    write_faa(bsub_records, bsub_faa)
    write_faa(syn1_records, syn1_faa)
    bsub_ann.to_csv(OUT_HOM / "bsub_proteins_annotation.tsv", sep="\t", index=False)
    syn1_ann.to_csv(OUT_HOM / "syn1_proteins_annotation.tsv", sep="\t", index=False)
    print(f"[homology] B.subtilis proteins={len(bsub_ann)}  Syn1 proteins={len(syn1_ann)}")

    # reciprocal BLASTP
    blast_dir = OUT_HOM / "blast"
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    fmt = "6 " + " ".join(cols)
    fwd_tsv = OUT_HOM / "bsub_vs_syn1.tsv"
    rev_tsv = OUT_HOM / "syn1_vs_bsub.tsv"

    def run(cmd):
        print("  $", " ".join(str(c) for c in cmd))
        subprocess.run(cmd, check=True)

    run(["makeblastdb", "-in", str(bsub_faa), "-dbtype", "prot", "-out", str(blast_dir / "bsub_prot_db")])
    run(["makeblastdb", "-in", str(syn1_faa), "-dbtype", "prot", "-out", str(blast_dir / "syn1_prot_db")])
    with fwd_tsv.open("w") as fh:
        subprocess.run(["blastp", "-query", str(bsub_faa), "-db", str(blast_dir / "syn1_prot_db"),
                        "-evalue", BLAST_EVALUE, "-max_target_seqs", BLAST_MAX_TARGET,
                        "-outfmt", fmt], check=True, stdout=fh)
    with rev_tsv.open("w") as fh:
        subprocess.run(["blastp", "-query", str(syn1_faa), "-db", str(blast_dir / "bsub_prot_db"),
                        "-evalue", BLAST_EVALUE, "-max_target_seqs", BLAST_MAX_TARGET,
                        "-outfmt", fmt], check=True, stdout=fh)

    bsub_len = dict(zip(bsub_ann["fasta_header"], bsub_ann["protein_len"]))
    syn1_len = dict(zip(syn1_ann["fasta_header"], syn1_ann["protein_len"]))
    fwd = pd.read_csv(fwd_tsv, sep="\t", names=cols)
    rev = pd.read_csv(rev_tsv, sep="\t", names=cols)
    for df in (fwd, rev):
        for c in cols[2:]:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    fwd["q_len"] = fwd["qseqid"].map(bsub_len); fwd["s_len"] = fwd["sseqid"].map(syn1_len)
    fwd["qcov"] = 100.0 * fwd["length"] / fwd["q_len"]; fwd["scov"] = 100.0 * fwd["length"] / fwd["s_len"]
    rev["q_len"] = rev["qseqid"].map(syn1_len); rev["s_len"] = rev["sseqid"].map(bsub_len)
    rev["qcov"] = 100.0 * rev["length"] / rev["q_len"]; rev["scov"] = 100.0 * rev["length"] / rev["s_len"]

    # deterministic best-hit pick (stable sort + sseqid tiebreak), matches Homology_Build.py
    sk = ["qseqid", "bitscore", "evalue", "pident", "sseqid"]
    sa = [True, False, True, False, True]
    fwd_best = fwd.sort_values(sk, ascending=sa, kind="mergesort").groupby("qseqid", as_index=False).first()
    rev_best = rev.sort_values(sk, ascending=sa, kind="mergesort").groupby("qseqid", as_index=False).first()
    rev_lookup = dict(zip(rev_best["qseqid"], rev_best["sseqid"]))

    rbh_rows = []
    for r in fwd_best.itertuples(index=False):
        rbh_rows.append(dict(bsub_id=r.qseqid, syn1_id=r.sseqid,
                             is_rbh=(rev_lookup.get(r.sseqid) == r.qseqid),
                             pident=r.pident, align_len=r.length, bsub_len=r.q_len, syn1_len=r.s_len,
                             bsub_cov_pct=r.qcov, syn1_cov_pct=r.scov, evalue=r.evalue, bitscore=r.bitscore))
    rbh = pd.DataFrame(rbh_rows)

    bsub_meta = bsub_ann.rename(columns={"fasta_header": "bsub_id", "primary_id": "bsub_primary_id",
                                         "locus_tag": "bsub_locus_tag", "gene": "bsub_gene",
                                         "product": "bsub_product", "seqid": "bsub_seqid",
                                         "start_1b": "bsub_start_1b", "end_1b": "bsub_end_1b",
                                         "strand": "bsub_strand"})
    syn1_meta = syn1_ann.rename(columns={"fasta_header": "syn1_id", "primary_id": "syn1_primary_id",
                                         "locus_tag": "syn1_locus_tag", "gene": "syn1_gene",
                                         "product": "syn1_product", "seqid": "syn1_seqid",
                                         "start_1b": "syn1_start_1b", "end_1b": "syn1_end_1b",
                                         "strand": "syn1_strand"})
    rbh = rbh.merge(bsub_meta[["bsub_id", "bsub_primary_id", "bsub_locus_tag", "bsub_gene",
                               "bsub_product", "bsub_seqid", "bsub_start_1b", "bsub_end_1b",
                               "bsub_strand"]], on="bsub_id", how="left")
    rbh = rbh.merge(syn1_meta[["syn1_id", "syn1_primary_id", "syn1_locus_tag", "syn1_gene",
                               "syn1_product", "syn1_seqid", "syn1_start_1b", "syn1_end_1b",
                               "syn1_strand"]], on="syn1_id", how="left")
    rbh.to_csv(OUT_HOM / "bsub_syn1_besthits_table.tsv", sep="\t", index=False)

    rbh_strict = rbh[(rbh["is_rbh"]) & (rbh["pident"] >= MIN_PIDENT)
                     & (rbh["bsub_cov_pct"] >= MIN_QCOV) & (rbh["syn1_cov_pct"] >= MIN_SCOV)
                     & (rbh["evalue"] <= MIN_EVALUE) & (rbh["bitscore"] >= MIN_BITSCORE)].copy()
    rbh_strict = rbh_strict.sort_values(["pident", "bitscore", "bsub_cov_pct", "syn1_cov_pct"],
                                        ascending=[False, False, False, False]).reset_index(drop=True)
    rbh_strict.to_csv(OUT_HOM / "bsub_syn1_rbh_homology_table.tsv", sep="\t", index=False)
    print(f"[homology] best-hit pairs={len(rbh)}  strict RBH pairs={len(rbh_strict)}")
    return rbh_strict


def get_homology_table():
    """Return (strict RBH homology table, Syn1 protein annotation).

    Prefer the canonical pre-built table from Genomes_Input/homology_syn1_bsub/
    (Homology_Build.py output -- same RBH method + filters as build_homology_table
    below).  Fall back to building it locally if the canonical files are missing,
    so the pipeline still runs standalone."""
    rbh_path = CANONICAL_HOM_DIR / "bsub_syn1_rbh_homology_table.tsv"
    ann_path = CANONICAL_HOM_DIR / "syn1_proteins_annotation.tsv"
    if USE_PREBUILT_HOMOLOGY and rbh_path.exists() and ann_path.exists():
        rbh = pd.read_csv(rbh_path, sep="\t")
        syn1_ann = pd.read_csv(ann_path, sep="\t")
        print(f"[homology] using canonical pre-built table: {rbh_path} ({len(rbh)} strict RBH pairs)")
        return rbh, syn1_ann
    print("[homology] canonical table absent -> building locally (Homology_Build.py method)")
    rbh = build_homology_table()
    syn1_ann = pd.read_csv(OUT_HOM / "syn1_proteins_annotation.tsv", sep="\t")
    return rbh, syn1_ann


# ----------------------------------------------------------------------------
# Stage 3 -- transfer RNase III candidate genes to Syn1 (cells 58, 59)
# ----------------------------------------------------------------------------
RBH_COLS = ["bsub_id", "bsub_primary_id", "bsub_locus_tag", "bsub_gene", "bsub_product",
            "syn1_id", "syn1_primary_id", "syn1_locus_tag", "syn1_gene", "syn1_product",
            "pident", "align_len", "bsub_cov_pct", "syn1_cov_pct", "evalue", "bitscore"]


def _merge_with_rescue(left: pd.DataFrame, rbh: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in RBH_COLS if c in rbh.columns]
    m1 = left.merge(rbh[["bsub_locus_tag_key"] + cols].drop_duplicates("bsub_locus_tag_key"),
                    on="bsub_locus_tag_key", how="left")
    unmatched = m1["syn1_locus_tag"].isna() if "syn1_locus_tag" in m1.columns else pd.Series([True] * len(m1))
    un = m1[unmatched].drop(columns=[c for c in cols if c in m1.columns], errors="ignore")
    m2 = un.merge(rbh[["bsub_gene_key"] + cols].drop_duplicates("bsub_gene_key"),
                  on="bsub_gene_key", how="left")
    return pd.concat([m1[~unmatched], m2], ignore_index=True)


def transfer_rnaseIII(sites: pd.DataFrame, genes: pd.DataFrame, rbh: pd.DataFrame,
                      syn1_ann: pd.DataFrame):
    sites = sites.copy(); genes = genes.copy(); rbh = rbh.copy()
    sites["bsub_locus_tag_key"] = sites.get("primary_locus_tag", "").map(clean_str)
    sites["bsub_gene_key"] = sites.get("primary_gene", "").map(clean_str)
    genes["bsub_locus_tag_key"] = genes.get("primary_locus_tag", "").map(clean_str)
    genes["bsub_gene_key"] = genes.get("primary_name", "").map(lambda x: clean_str(x).split("/")[0])
    rbh["bsub_locus_tag_key"] = rbh["bsub_locus_tag"].map(clean_str)
    rbh["bsub_gene_key"] = rbh["bsub_gene"].map(clean_str)

    sites_final = _merge_with_rescue(sites, rbh)
    genes_final = _merge_with_rescue(genes, rbh)
    sites_final["has_syn1_homolog"] = sites_final["syn1_locus_tag"].notna()
    genes_final["has_syn1_homolog"] = genes_final["syn1_locus_tag"].notna()

    sites_final.to_csv(OUT_R3 / "rnaseIII_hits_transferred_to_syn1.tsv", sep="\t", index=False)
    genes_final.to_csv(OUT_R3 / "rnaseIII_genes_transferred_to_syn1.tsv", sep="\t", index=False)

    # add Syn1 coordinates (cell 59)
    ann = syn1_ann.rename(columns={"locus_tag": "syn1_locus_tag", "gene": "syn1_gene_annot",
                                   "product": "syn1_product_annot", "seqid": "syn1_seqid",
                                   "start_1b": "syn1_start_1b", "end_1b": "syn1_end_1b",
                                   "strand": "syn1_strand", "primary_id": "syn1_primary_id",
                                   "protein_id": "syn1_protein_id"})
    ann = ann.sort_values(["syn1_locus_tag", "syn1_start_1b"]).drop_duplicates("syn1_locus_tag")
    keep = ["syn1_locus_tag", "syn1_seqid", "syn1_start_1b", "syn1_end_1b", "syn1_strand",
            "syn1_primary_id", "syn1_protein_id", "syn1_gene_annot", "syn1_product_annot"]
    merged = genes_final.merge(ann[keep], on="syn1_locus_tag", how="left")
    merged.to_csv(OUT_R3 / "rnaseIII_genes_transferred_to_syn1_with_coords.tsv", sep="\t", index=False)

    n_hom = int(genes_final["has_syn1_homolog"].sum())
    n_coord = int(merged["syn1_start_1b"].notna().sum())
    print(f"[RNase III] candidate genes={len(genes)}  with Syn1 homolog={n_hom}  with coords={n_coord}")
    return merged


# ----------------------------------------------------------------------------
# Stage 4 -- predict RNase III cleavage geometry in Syn1 (cell 61, ViennaRNA)
# ----------------------------------------------------------------------------
def _norm_seqid(x: str) -> str:
    x = str(x)
    if "|" in x:
        parts = [p for p in x.split("|") if p]
        for p in parts:
            if p.startswith("CP") and "." in p:
                return p
        return max(parts, key=len) if parts else x
    return x


def _dotbracket_to_pairs(struct: str) -> list[int]:
    p = [0] * (len(struct) + 1)
    st = []
    for i, ch in enumerate(struct, start=1):
        if ch == "(":
            st.append(i)
        elif ch == ")" and st:
            j = st.pop(); p[i] = j; p[j] = i
    return p


def _local_stemrun(struct, ptab, idx) -> int:
    n = len(struct)
    if idx < 1 or idx > n or ptab[idx] == 0:
        return 0
    arm = struct[idx - 1]
    if arm not in ("(", ")"):
        return 0
    run = 1
    i = idx - 1
    while i >= 1 and struct[i - 1] == arm and ptab[i] != 0:
        run += 1; i -= 1
    i = idx + 1
    while i <= n and struct[i - 1] == arm and ptab[i] != 0:
        run += 1; i += 1
    return run


def _pairing_metrics(ptab, rel1, rel2, n, w, direct_offset):
    W1 = range(max(1, rel1 - w), min(n, rel1 + w) + 1)
    W2 = range(max(1, rel2 - w), min(n, rel2 + w) + 1)
    W1set, W2set = set(W1), set(W2)
    direct = (ptab[rel1] != 0) and (abs(ptab[rel1] - rel2) <= direct_offset)
    c12 = sum(1 for i in W1 if ptab[i] != 0 and ptab[i] in W2set)
    c21 = sum(1 for j in W2 if ptab[j] != 0 and ptab[j] in W1set)
    frac = 0.5 * (c12 / len(W1set) + c21 / len(W2set))
    return direct, c12, c21, frac


def _overhang_delta(ptab, i1, i2, overhang=2, window=0):
    n = len(ptab) - 1
    paired = lambda i: ptab[i] if 1 <= i <= n else 0
    best_delta, best = None, None
    for a in range(-window, window + 1):
        for b in range(-window, window + 1):
            x1, x2 = i1 + a, i2 + b
            y1, y2 = paired(x1), paired(x2)
            if y1 == 0 or y2 == 0:
                continue
            S = (y1 - y2) + (x1 - x2)
            delta = min(abs(S - 2 * overhang), abs(S + 2 * overhang))
            if best_delta is None or delta < best_delta:
                best_delta, best = delta, (x1, y1, x2, y2, a, b)
    return best_delta, best


def predict_rnaseIII_cleavage(cand: pd.DataFrame):
    import RNA
    genome = load_fasta_as_dict(SYN1_FASTA)
    md = RNA.md(); md.temperature = R3_TEMPERATURE_C
    if R3_NO_LP:
        md.noLP = 1

    cand = cand.dropna(subset=["syn1_locus_tag"]).copy()
    cand["syn1_seqid"] = cand["syn1_seqid"].map(_norm_seqid)
    cand["syn1_start_1b"] = cand["syn1_start_1b"].astype(int)
    cand["syn1_end_1b"] = cand["syn1_end_1b"].astype(int)
    cand["syn1_strand"] = cand["syn1_strand"].astype(str)

    optional_meta = ["syn1_gene", "syn1_locus_tag", "syn1_product",
                     "bsub_gene", "bsub_locus_tag", "primary_name", "primary_locus_tag"]
    all_rows = []
    for ridx, row in cand.reset_index(drop=True).iterrows():
        seqid = _norm_seqid(row["syn1_seqid"])
        g0, g1, strand = int(row["syn1_start_1b"]), int(row["syn1_end_1b"]), row["syn1_strand"]
        seq = genome[seqid]; G = len(seq)
        gs = max(1, g0 - R3_FLANK); ge = min(G, g1 + R3_FLANK)
        frag = seq[gs - 1:ge]
        orient = "minus" if str(strand).startswith("-") else "plus"
        rna = to_rna(frag if orient == "plus" else revcomp(frag))

        fc = RNA.fold_compound(rna, md)
        struct, mfe = fc.mfe()
        ptab = _dotbracket_to_pairs(struct)
        n = len(rna)
        g_from_rel = (lambda rel: gs + rel - 1) if orient == "plus" else (lambda rel: ge - rel + 1)
        paired_positions = [i for i in range(1, n + 1) if ptab[i] != 0]

        region_rows = []
        for i in paired_positions:
            for j in paired_positions:
                if j <= i:
                    continue
                dist = j - i
                if dist < R3_MIN_DIST or dist > R3_MAX_DIST:
                    continue
                sr_i = _local_stemrun(struct, ptab, i)
                sr_j = _local_stemrun(struct, ptab, j)
                if sr_i < R3_MIN_STEMRUN or sr_j < R3_MIN_STEMRUN:
                    continue
                direct, c12, c21, frac = _pairing_metrics(ptab, i, j, n, R3_W, R3_DIRECT_OFFSET)
                if not (direct or c12 >= 5 or frac >= 0.25):
                    continue
                od, best = _overhang_delta(ptab, i, j, R3_OVERHANG, R3_OVERHANG_WINDOW)
                if best is None:
                    continue
                x1, y1, x2, y2, a, b = best
                score = ((10 if od == 0 else max(0, 8 - od)) + c12 + int(round(frac * 10))
                         + sr_i + sr_j + (5 if direct else 0))
                out = dict(candidate_region_id=ridx, syn1_seqid=seqid, syn1_strand=strand,
                           fold_orient=orient, fold_g_start_1b=gs, fold_g_end_1b=ge, fold_len=n,
                           mfe_kcal_mol=float(mfe), rel_cut1=i, rel_cut2=j,
                           genomic_cut1=g_from_rel(i), genomic_cut2=g_from_rel(j), rel_distance=dist,
                           partner_cut1=int(ptab[i]), partner_cut2=int(ptab[j]),
                           cut1_stemrun=sr_i, cut2_stemrun=sr_j, direct_hit_pair=bool(direct),
                           cross_pairs_12=int(c12), cross_pairs_21=int(c21), cross_pair_frac=float(frac),
                           overhang=R3_OVERHANG, overhang_window=R3_OVERHANG_WINDOW,
                           overhang_best_delta=int(od), best_shift_a=int(a), best_shift_b=int(b),
                           best_rel_cut1=int(x1), best_rel_cut2=int(x2), best_partner1=int(y1),
                           best_partner2=int(y2), best_genomic_cut1=int(g_from_rel(x1)),
                           best_genomic_cut2=int(g_from_rel(x2)), composite_score=int(score))
                for c in optional_meta:
                    if c in row.index:
                        out[c] = row[c]
                if R3_STORE_STRUCT:
                    out["dotbracket"] = struct
                region_rows.append(out)
        if not region_rows:
            continue
        rdf = pd.DataFrame(region_rows).sort_values(
            ["composite_score", "overhang_best_delta", "cross_pair_frac", "cross_pairs_12",
             "cut1_stemrun", "cut2_stemrun"], ascending=[False, True, False, False, False, False]
        ).reset_index(drop=True)
        rdf["rank_within_region"] = range(1, len(rdf) + 1)
        all_rows.append(rdf.head(R3_TOP_K_PER_REGION))

    out_tsv = OUT_R3 / "rnaseIII_syn1_predicted_cleavage_pairs.tsv"
    if all_rows:
        outdf = pd.concat(all_rows, ignore_index=True)
        outdf.to_csv(out_tsv, sep="\t", index=False)
        print(f"[RNase III] regions scanned={len(cand)}  predicted cleavage pairs={len(outdf)}")
        return outdf
    print("[RNase III] no cleavage pairs passed structural filters")
    return pd.DataFrame()


# ----------------------------------------------------------------------------
# Stage 4b -- HOMOLOGY-ANCHORED RNase III cleavage sites
# ----------------------------------------------------------------------------
# predict_rnaseIII_cleavage() above folds the whole gene and picks the strongest
# stem *de novo* -- it ignores where the B. subtilis cut actually was.  This stage
# instead ANCHORS every call to a real B. subtilis site: project that site into
# the Syn1 homolog by transcript fraction (the same projector the RNase Y stage
# uses), fold a local window, and read out the secondary structure at the
# projected position.  RNase III is a dsRNA endonuclease, so a bona-fide site
# should sit in (or right at the base of) a helix.  Output gives the *exact*
# projected Syn1 cut coordinate per B. subtilis site -- this is what the panel-f
# atpA split line is anchored to.
def anchor_rnaseIII_cleavage(sites_annot: pd.DataFrame, genes_coords: pd.DataFrame):
    import RNA
    genome = load_fasta_as_dict(SYN1_FASTA)
    md = RNA.md(); md.temperature = R3_TEMPERATURE_C
    if R3_NO_LP:
        md.noLP = 1

    gc = genes_coords.dropna(subset=["syn1_locus_tag", "syn1_start_1b"]).drop_duplicates("primary_locus_tag")
    syn1_of = {clean_str(r["primary_locus_tag"]): r for _, r in gc.iterrows()}

    sub = sites_annot[sites_annot["context"] == "sense_intragenic"].copy()
    rows = []
    for _, s in sub.iterrows():
        bloc = clean_str(s["primary_locus_tag"])
        if bloc not in syn1_of:
            continue
        g = syn1_of[bloc]
        bstart, bend, bstr = int(s["primary_gene_start"]), int(s["primary_gene_end"]), s["primary_gene_strand"]
        frac = _tx_fraction(int(s["pos"]), bstart, bend, bstr)
        if pd.isna(frac):
            continue
        sstart, send, sstr = int(g["syn1_start_1b"]), int(g["syn1_end_1b"]), str(g["syn1_strand"])
        seqid = _norm_seqid(str(g["syn1_seqid"]))
        cut1 = _project_fraction(sstart, send, sstr, frac)
        if pd.isna(cut1) or seqid not in genome:
            continue
        wseq, cut0, fgs, fge = _fold_window(genome[seqid], int(cut1), sstr, R3_ANCHOR_UP, R3_ANCHOR_DOWN)
        if len(wseq) < 20 or pd.isna(cut0):
            continue
        struct, mfe = RNA.fold_compound(to_rna(wseq), md).mfe()
        ptab = _dotbracket_to_pairs(struct)
        rel = int(cut0) + 1                         # 1-based projected-cut index in the fold
        in_helix = ptab[rel] != 0
        eval_rel = rel
        if not in_helix:                            # nearest paired base within +/- window
            for d in range(1, R3_ANCHOR_PARTNER_WIN + 1):
                hit = next((c for c in (rel - d, rel + d) if 1 <= c <= len(struct) and ptab[c] != 0), None)
                if hit:
                    eval_rel = hit
                    break
        partner_rel = ptab[eval_rel]
        stemrun = _local_stemrun(struct, ptab, eval_rel) if partner_rel else 0
        g_from_rel = (lambda r: fgs + r - 1) if sstr == "+" else (lambda r: fge - r + 1)
        rows.append(dict(
            bsub_locus_tag=bloc, bsub_gene=clean_str(s["primary_gene"]), bsub_site_id=s["site_id"],
            bsub_site_pos=int(s["pos"]), bsub_fraction_from_5prime=round(float(frac), 4),
            sensitivity_score=safe_float(s["sensitivity_score"]),
            syn1_locus_tag=g["syn1_locus_tag"], syn1_gene=clean_str(g.get("syn1_gene_annot", "")),
            syn1_product=clean_str(g.get("syn1_product_annot", "")), syn1_seqid=seqid, syn1_strand=sstr,
            syn1_projected_cut_1b=int(cut1), cut_in_helix=bool(in_helix),
            cut_or_near_paired=bool(partner_rel != 0),
            paired_partner_genomic=int(g_from_rel(partner_rel)) if partner_rel else np.nan,
            cut_stemrun=int(stemrun), local_fold_start_1b=fgs, local_fold_end_1b=fge,
            local_mfe_kcal_mol=round(float(mfe), 2), pident=safe_float(g.get("pident", np.nan)),
            fold_struct=struct, fold_seq=to_rna(wseq), cut_idx0_in_fold=int(cut0),
        ))
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["syn1_locus_tag", "syn1_projected_cut_1b"]).reset_index(drop=True)
    out.to_csv(OUT_R3 / "rnaseIII_syn1_anchored_cleavage_sites.tsv", sep="\t", index=False)
    n_helix = int(out["cut_in_helix"].sum()) if not out.empty else 0
    n_gene = out["syn1_locus_tag"].nunique() if not out.empty else 0
    print(f"[RNase III] homology-anchored cuts={len(out)} in {n_gene} genes "
          f"({n_helix} land directly in a helix)")
    return out


# ----------------------------------------------------------------------------
# Stage 5 -- project RNase Y sites into Syn1 homologs (cell 74)
# ----------------------------------------------------------------------------
def _split_csv(x):
    return [t.strip() for t in str(x).split(",") if t.strip()] if pd.notna(x) and str(x).strip() else []


def _split_pipe(x):
    return [t.strip() for t in str(x).split("|") if t.strip()] if pd.notna(x) and str(x).strip() else []


def _tx_fraction(pos, start1, end1, strand) -> float:
    glen = end1 - start1 + 1
    if glen <= 1 or not (start1 <= pos <= end1):
        return np.nan
    if strand == "+":
        return (pos - start1) / (glen - 1)
    if strand == "-":
        return (end1 - pos) / (glen - 1)
    return np.nan


def _project_fraction(start1, end1, strand, frac) -> int:
    glen = end1 - start1 + 1
    if glen <= 0 or pd.isna(frac):
        return np.nan
    frac = min(max(float(frac), 0.0), 1.0)
    if strand == "+":
        return int(start1 + round(frac * max(1, glen - 1)))
    if strand == "-":
        return int(end1 - round(frac * max(1, glen - 1)))
    return np.nan


def project_rnaseY(gene_summary: pd.DataFrame, rbh: pd.DataFrame) -> pd.DataFrame:
    hom = rbh.copy()
    for c in ["bsub_start_1b", "bsub_end_1b", "syn1_start_1b", "syn1_end_1b",
              "align_len", "bsub_len", "syn1_len"]:
        hom[c] = hom[c].apply(safe_int)
    for c in ["pident", "bsub_cov_pct", "syn1_cov_pct", "evalue", "bitscore"]:
        hom[c] = hom[c].apply(safe_float)
    hom["bsub_strand"] = hom["bsub_strand"].map(normalize_strand)
    hom["syn1_strand"] = hom["syn1_strand"].map(normalize_strand)
    for c in ["bsub_id", "syn1_id", "bsub_primary_id", "bsub_locus_tag", "bsub_gene", "bsub_product",
              "bsub_seqid", "syn1_primary_id", "syn1_locus_tag", "syn1_gene", "syn1_product", "syn1_seqid"]:
        hom[c] = hom[c].map(clean_str)
    hom["is_rbh_bool"] = hom["is_rbh"].astype(str).str.strip().str.lower().isin({"1", "true", "t", "yes", "y"})
    if RY_KEEP_ONLY_RBH:
        hom = hom[hom["is_rbh_bool"]]
    hom = hom[(hom["pident"].fillna(-np.inf) >= MIN_PIDENT)
              & (hom["bsub_cov_pct"].fillna(-np.inf) >= MIN_QCOV)
              & (hom["syn1_cov_pct"].fillna(-np.inf) >= MIN_SCOV)].copy()

    # expand gene summary -> one row per site
    exp = []
    for _, row in gene_summary.iterrows():
        ids = _split_csv(row["site_ids"]); pos = _split_csv(row["site_positions"])
        loc = _split_pipe(row["local_sequences"])
        n = max(len(ids), len(pos), len(loc))
        if n == 0:
            continue
        ids += [""] * (n - len(ids)); pos += [""] * (n - len(pos)); loc += [""] * (n - len(loc))
        for sid, spos, lseq in zip(ids, pos, loc):
            p = safe_int(spos)
            if pd.isna(p):
                continue
            exp.append(dict(primary_name=clean_str(row["primary_name"]), primary_gene=clean_str(row["primary_gene"]),
                            primary_locus_tag=clean_str(row["primary_locus_tag"]),
                            bsub_site_strand=normalize_strand(row["strand"]),
                            bsub_site_id=str(sid).strip(), bsub_site_pos=int(p),
                            local_sequence=str(lseq).strip(),
                            mean_sensitivity_score=safe_float(row.get("mean_sensitivity_score", np.nan)),
                            median_sensitivity_score=safe_float(row.get("median_sensitivity_score", np.nan)),
                            max_sensitivity_score=safe_float(row.get("max_sensitivity_score", np.nan)),
                            n_candidate_sites_in_gene=safe_int(row.get("n_candidate_sites", np.nan))))
    sites_exp = pd.DataFrame(exp)
    print(f"[RNase Y] expanded sites={len(sites_exp)}  filtered homology rows={len(hom)}")

    # join expanded sites to homology (locus_tag -> gene -> primary_id -> id)
    parts = []
    for left_key, right_key in [("primary_locus_tag", "bsub_locus_tag"),
                                ("primary_gene", "bsub_gene"),
                                ("primary_name", "bsub_primary_id"),
                                ("primary_name", "bsub_id")]:
        a = sites_exp.copy(); a["join_key"] = a[left_key].map(clean_str)
        b = hom.copy(); b["join_key"] = b[right_key].map(clean_str)
        m = a.merge(b, on="join_key", how="inner", suffixes=("", "_hom"))
        m["join_method"] = f"{left_key} -> {right_key}"
        parts.append(m)
    matched = pd.concat(parts, ignore_index=True)
    dedup = ["bsub_site_id", "bsub_site_pos", "primary_locus_tag", "primary_gene", "primary_name",
             "syn1_id", "syn1_locus_tag", "syn1_gene", "syn1_start_1b", "syn1_end_1b", "syn1_strand"]
    matched = matched.drop_duplicates(subset=dedup)
    if RY_MULTI_HOMOLOG_POLICY == "best_bitscore":
        matched = (matched.sort_values(["bsub_site_id", "bsub_site_pos", "bitscore", "pident",
                                        "bsub_cov_pct", "syn1_cov_pct"],
                                       ascending=[True, True, False, False, False, False])
                          .drop_duplicates(["bsub_site_id", "bsub_site_pos"], keep="first"))

    matched["site_inside_bsub_gene"] = ((matched["bsub_site_pos"] >= matched["bsub_start_1b"])
                                        & (matched["bsub_site_pos"] <= matched["bsub_end_1b"]))
    proj = matched[matched["site_inside_bsub_gene"]].copy()
    proj["bsub_fraction_from_5prime"] = proj.apply(
        lambda r: _tx_fraction(int(r["bsub_site_pos"]), int(r["bsub_start_1b"]),
                               int(r["bsub_end_1b"]), r["bsub_strand"]), axis=1)
    proj = proj[proj["bsub_fraction_from_5prime"].notna()].copy()
    proj["syn1_projected_pos_1b"] = proj.apply(
        lambda r: _project_fraction(int(r["syn1_start_1b"]), int(r["syn1_end_1b"]),
                                    r["syn1_strand"], float(r["bsub_fraction_from_5prime"])), axis=1)
    proj["syn1_proj_win_start_1b"] = (proj["syn1_projected_pos_1b"] - RY_PROJECTION_HALF_WINDOW).clip(lower=1).astype(int)
    proj["syn1_proj_win_end_1b"] = (proj["syn1_projected_pos_1b"] + RY_PROJECTION_HALF_WINDOW).astype(int)
    proj = proj.reset_index(drop=True)
    proj["syn1_site_id"] = [f"Syn1_RY_{i+1}" for i in range(len(proj))]

    front = ["syn1_site_id", "bsub_site_id", "bsub_site_pos", "bsub_site_strand", "local_sequence",
             "primary_name", "primary_gene", "primary_locus_tag",
             "bsub_id", "bsub_primary_id", "bsub_locus_tag", "bsub_gene", "bsub_product", "bsub_seqid",
             "bsub_start_1b", "bsub_end_1b", "bsub_strand",
             "syn1_id", "syn1_primary_id", "syn1_locus_tag", "syn1_gene", "syn1_product", "syn1_seqid",
             "syn1_start_1b", "syn1_end_1b", "syn1_strand",
             "is_rbh", "pident", "align_len", "bsub_len", "syn1_len", "bsub_cov_pct", "syn1_cov_pct",
             "evalue", "bitscore", "join_method", "bsub_fraction_from_5prime", "syn1_projected_pos_1b",
             "syn1_proj_win_start_1b", "syn1_proj_win_end_1b", "n_candidate_sites_in_gene",
             "mean_sensitivity_score", "median_sensitivity_score", "max_sensitivity_score"]
    final = proj[[c for c in front if c in proj.columns]].copy()
    final = final.sort_values(["syn1_seqid", "syn1_strand", "syn1_projected_pos_1b", "bitscore"],
                              ascending=[True, True, True, False]).reset_index(drop=True)
    final.to_csv(OUT_RY / "rnaseY_syn1_projected_sites.tsv", sep="\t", index=False)
    print(f"[RNase Y] projected Syn1 sites={len(final)}  in {final['syn1_id'].nunique()} genes")
    return final


# ----------------------------------------------------------------------------
# Stage 6 -- RNase Y downstream secondary-structure support (cell 76)
# ----------------------------------------------------------------------------
def _subseq_1b(seq, s1, e1, strand="+"):
    n = len(seq); s1 = max(1, int(s1)); e1 = min(n, int(e1))
    if e1 < s1:
        return ""
    s = seq[s1 - 1:e1]
    return revcomp(s) if strand == "-" else s


def _fold_window(seq, cut1, strand, up_len, down_len):
    gs, ge = cut1 - up_len, cut1 + down_len
    sub = _subseq_1b(seq, gs, ge, strand)
    if strand == "+":
        cut_idx0 = cut1 - max(1, gs)
    elif strand == "-":
        cut_idx0 = min(len(seq), ge) - cut1
    else:
        return "", np.nan, np.nan, np.nan
    return sub, int(cut_idx0), max(1, gs), min(len(seq), ge)


def _rel_window(seq, cut_idx0, rel_start, rel_end):
    if pd.isna(cut_idx0):
        return "", np.nan, np.nan
    s0 = max(0, int(cut_idx0) + int(rel_start))
    e0 = min(len(seq) - 1, int(cut_idx0) + int(rel_end))
    if e0 < s0:
        return "", np.nan, np.nan
    return seq[s0:e0 + 1], s0, e0


def _au(seq):
    if not seq:
        return np.nan
    s = seq.upper().replace("T", "U")
    return (s.count("A") + s.count("U")) / len(s)


def _build_pair_map(struct):
    stack, pair = [], [-1] * len(struct)
    for i, ch in enumerate(struct):
        if ch == "(":
            stack.append(i)
        elif ch == ")" and stack:
            j = stack.pop(); pair[i] = j; pair[j] = i
    return pair


def _best_downstream_stem(struct, cut0, dmin, dmax):
    empty = dict(has_downstream_stem=False, best_stem_bp=0, best_stem_start_rel=np.nan,
                 best_stem_end_rel=np.nan, best_loop_len=np.nan, best_segment_struct="",
                 _i=None, _j=None, _ii=None, _jj=None)
    n = len(struct)
    if pd.isna(cut0):
        return empty
    s0 = max(0, int(cut0) + dmin); s1 = min(n - 1, int(cut0) + dmax)
    if s1 < s0:
        return empty
    pair = _build_pair_map(struct)
    best = dict(empty)
    for i in range(s0, s1 + 1):
        if struct[i] != "(":
            continue
        j = pair[i]
        if j == -1 or j <= i:
            continue
        bp, ii, jj = 1, i, j
        while True:
            ni, nj = ii + 1, jj - 1
            if ni >= nj:
                break
            if struct[ni] == "(" and pair[ni] == nj:
                bp += 1; ii, jj = ni, nj; continue
            found = False
            for di in range(0, RY_MAX_BULGE_RUN + 1):
                for dj in range(0, RY_MAX_BULGE_RUN + 1):
                    if di == 0 and dj == 0:
                        continue
                    ti, tj = ii + 1 + di, jj - 1 - dj
                    if ti < tj and 0 <= ti < n and 0 <= tj < n and struct[ti] == "(" and pair[ti] == tj:
                        bp += 1; ii, jj = ti, tj; found = True; break
                if found:
                    break
            if not found:
                break
        loop_len = (jj - ii - 1) if jj > ii else np.nan
        if bp >= RY_MIN_STEM_BP and (pd.isna(loop_len) or (RY_MIN_HAIRPIN_LOOP <= loop_len <= RY_MAX_HAIRPIN_LOOP)):
            if bp > best["best_stem_bp"]:
                best.update(has_downstream_stem=True, best_stem_bp=int(bp),
                            best_stem_start_rel=int(i - cut0), best_stem_end_rel=int(j - cut0),
                            best_loop_len=int(loop_len) if pd.notna(loop_len) else np.nan,
                            best_segment_struct=struct[i:j + 1], _i=i, _j=j, _ii=ii, _jj=jj)
    return best


def _classify_rnaseY(primary_mfe, delta_mfe, has_stem, best_bp, downstream_au):
    primary_ok = pd.notna(primary_mfe) and primary_mfe <= RY_PRIMARY_MFE_STRONG
    delta_ok = pd.notna(delta_mfe) and delta_mfe <= RY_PRIMARY_DELTA_MFE_STRONG
    stem_ok = bool(has_stem) and int(best_bp) >= RY_MIN_STEM_BP
    au_ok = pd.notna(downstream_au) and downstream_au >= RY_MIN_DOWNSTREAM_AU
    if primary_ok and delta_ok and stem_ok and au_ok:
        return "supported"
    if primary_ok and delta_ok:
        return "mfe_supported"
    if stem_ok:
        return "stem_only"
    return "not_supported"


def score_rnaseY_structure(proj: pd.DataFrame):
    import RNA
    genome = load_fasta_as_dict(SYN1_FASTA)

    def fold(seq):
        if not seq:
            return "", np.nan
        struct, mfe = RNA.fold_compound(seq.replace("T", "U")).mfe()
        return struct, float(mfe)

    proj = proj.copy()
    proj["syn1_projected_pos_1b"] = pd.to_numeric(proj["syn1_projected_pos_1b"], errors="coerce")
    proj = proj.dropna(subset=["syn1_projected_pos_1b", "syn1_seqid", "syn1_strand"])
    proj["syn1_projected_pos_1b"] = proj["syn1_projected_pos_1b"].astype(int)
    proj["syn1_strand"] = proj["syn1_strand"].map(normalize_strand)

    rows = []
    for _, r in proj.iterrows():
        seqid = _norm_seqid(str(r["syn1_seqid"]).strip())
        cut1 = int(r["syn1_projected_pos_1b"]); strand = r["syn1_strand"]
        if seqid not in genome:
            continue
        gseq = genome[seqid]
        broad, cut0, fgs, fge = _fold_window(gseq, cut1, strand, RY_FOLD_UP_LEN, RY_FOLD_DOWN_LEN)
        if len(broad) < 20:
            continue
        bstruct, bmfe = fold(broad)
        pseq, _, _ = _rel_window(broad, cut0, *RY_PRIMARY_WIN); pstruct, pmfe = fold(pseq); pau = _au(pseq)
        cseq, _, _ = _rel_window(broad, cut0, *RY_CENTER_WIN); cstruct, cmfe = fold(cseq); cau = _au(cseq)
        dmfe = float(pmfe - cmfe) if pd.notna(pmfe) and pd.notna(cmfe) else np.nan
        cut_local, _, _ = _rel_window(broad, cut0, *RY_CUT_LOCAL_WIN); cut_local_au = _au(cut_local)
        down_seq, _, _ = _rel_window(broad, cut0, RY_DOWNSTREAM_MIN, RY_DOWNSTREAM_MAX); down_au = _au(down_seq)
        stem = _best_downstream_stem(bstruct, cut0, RY_DOWNSTREAM_MIN, RY_DOWNSTREAM_MAX)
        if stem.get("has_downstream_stem") and stem["_i"] is not None:
            i, j, ii, jj = stem["_i"], stem["_j"], stem["_ii"], stem["_jj"]
            left_arm, loop_seq, right_arm = broad[i:ii + 1], broad[ii + 1:jj], broad[jj:j + 1]
        else:
            left_arm = loop_seq = right_arm = ""
        cls = _classify_rnaseY(pmfe, dmfe, stem["has_downstream_stem"], stem["best_stem_bp"], down_au)
        rank = 0.0
        if pd.notna(pmfe):
            rank += -pmfe
        if pd.notna(dmfe):
            rank += -dmfe
        rank += 1.5 * int(bool(stem["has_downstream_stem"])) + float(stem["best_stem_bp"])
        rank += 1.0 * int(pd.notna(down_au) and down_au >= RY_MIN_DOWNSTREAM_AU)

        out = dict(r)
        out.update(fold_region_start_1b=fgs, fold_region_end_1b=fge, fold_seq=broad,
                   cut_idx0_in_fold_seq=cut0, fold_struct=bstruct, fold_mfe_kcal_mol=bmfe,
                   primary_window_seq=pseq, primary_window_struct=pstruct,
                   primary_window_mfe_kcal_mol=pmfe, primary_window_au_content=pau,
                   center_window_seq=cseq, center_window_struct=cstruct,
                   center_window_mfe_kcal_mol=cmfe, center_window_au_content=cau,
                   delta_mfe_primary_minus_center=dmfe, cut_local_seq=cut_local,
                   cut_local_au_content=cut_local_au, downstream_seq=down_seq,
                   downstream_au_content=down_au, has_downstream_stem=stem["has_downstream_stem"],
                   best_stem_bp=stem["best_stem_bp"], best_stem_start_rel=stem["best_stem_start_rel"],
                   best_stem_end_rel=stem["best_stem_end_rel"], best_loop_len=stem["best_loop_len"],
                   best_segment_struct=stem["best_segment_struct"], best_left_arm=left_arm,
                   best_loop_seq=loop_seq, best_right_arm=right_arm, downstream_stem_class=cls,
                   rnaseY_structure_rank_score=rank)
        rows.append(out)

    stem_df = pd.DataFrame(rows)
    if stem_df.empty:
        print("[RNase Y] no sites scored")
        return stem_df
    order = {"supported": 0, "mfe_supported": 1, "stem_only": 2, "not_supported": 3}
    stem_df["class_rank"] = stem_df["downstream_stem_class"].map(order).fillna(99)
    stem_df = stem_df.sort_values(["class_rank", "rnaseY_structure_rank_score",
                                   "primary_window_mfe_kcal_mol", "best_stem_bp",
                                   "syn1_seqid", "syn1_projected_pos_1b"],
                                  ascending=[True, False, True, False, True, True]).reset_index(drop=True)
    supported = stem_df[stem_df["downstream_stem_class"].isin(["supported", "mfe_supported"])].copy()
    stem_df = stem_df.drop(columns=["class_rank"])
    supported = supported.drop(columns=["class_rank"])
    stem_df.to_csv(OUT_RY / "rnaseY_syn1_projected_sites_with_stem.tsv", sep="\t", index=False)
    supported.to_csv(OUT_RY / "rnaseY_syn1_projected_sites_supported_by_downstream_stem.tsv", sep="\t", index=False)
    print(f"[RNase Y] scored sites={len(stem_df)}  supported/mfe_supported={len(supported)}")
    return stem_df


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    print("=" * 72)
    print("B. subtilis RNase III / RNase Y  ->  Syn1  cleavage-site mapping")
    print("=" * 72)

    # Stage 1: annotate cleavage sites to B. subtilis genes
    bsub_genes = parse_genbank_genes(BSUB_GB)
    print(f"[input] B. subtilis gene features parsed: {len(bsub_genes)}")
    _, r3_genes = annotate_and_summarize(
        "RNase III", SHEET_R3_LOOSE, "R3_", bsub_genes,
        OUT_R3 / "rnaseIII_loose_sites_annotated.tsv", OUT_R3 / "rnaseIII_loose_gene_summary.tsv")
    _, ry_genes = annotate_and_summarize(
        "RNase Y", SHEET_RY_LOOSE, "RY_", bsub_genes,
        OUT_RY / "rnaseY_loose_sites_annotated.tsv", OUT_RY / "rnaseY_loose_gene_summary.tsv")

    # Stage 2: reciprocal-best-hit homology (canonical pre-built table if available)
    rbh, syn1_ann = get_homology_table()
    r3_sites = pd.read_csv(OUT_R3 / "rnaseIII_loose_sites_annotated.tsv", sep="\t")

    # Stage 3+4: RNase III transfer + cleavage prediction (de-novo + homology-anchored)
    r3_coords = transfer_rnaseIII(r3_sites, r3_genes, rbh, syn1_ann)
    predict_rnaseIII_cleavage(r3_coords)
    anchor_rnaseIII_cleavage(r3_sites, r3_coords)

    # Stage 5+6: RNase Y projection + downstream-structure scoring
    ry_proj = project_rnaseY(ry_genes, rbh)
    score_rnaseY_structure(ry_proj)

    print("\nDone. Outputs in", OUT)


if __name__ == "__main__":
    main()
