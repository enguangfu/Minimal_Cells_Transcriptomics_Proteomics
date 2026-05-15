from __future__ import annotations

import logging
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

# -----------------------
# User config
# -----------------------

MOTHER_FOLDER = Path("./")
HOST_ORGANISM = "syn1"
GUEST_ORGANISM = "mpn"

HOST_GB = Path("./syn1.gb")
GUEST_GB = Path("./Mpn.gb")

OUT_DIR = MOTHER_FOLDER / f"homology_{HOST_ORGANISM}_{GUEST_ORGANISM}"

GUEST_FAA = OUT_DIR / f"{GUEST_ORGANISM}_proteins.faa"
HOST_FAA = OUT_DIR / f"{HOST_ORGANISM}_proteins.faa"

GUEST_ANN = OUT_DIR / f"{GUEST_ORGANISM}_proteins_annotation.tsv"
HOST_ANN = OUT_DIR / f"{HOST_ORGANISM}_proteins_annotation.tsv"

# -----------------------
# Homology search config
# -----------------------

GUEST_VS_HOST = OUT_DIR / f"{GUEST_ORGANISM}_{HOST_ORGANISM}.tsv"
HOST_VS_GUEST = OUT_DIR / f"{HOST_ORGANISM}_{GUEST_ORGANISM}.tsv"

OUT_RBH = OUT_DIR / f"{GUEST_ORGANISM}_{HOST_ORGANISM}_rbh_homology_table.tsv"
OUT_ALL = OUT_DIR / f"{GUEST_ORGANISM}_{HOST_ORGANISM}_besthits_table.tsv"

# BLAST search parameters (must match the columns in BLAST_COLS below)
BLAST_DB_DIR = OUT_DIR / "blast_db"
BLAST_EVALUE = 1e-5
BLAST_MAX_TARGET_SEQS = 5
BLAST_THREADS = 4
FORCE_BLAST = False  # re-run BLAST even if output TSVs already exist

# Filtering thresholds for final RBH table
MIN_PIDENT = 25.0
MIN_QCOV = 50.0
MIN_SCOV = 50.0
MIN_EVALUE = 1e-5
MIN_BITSCORE = 50.0

BLAST_PARAMS_TAG = (
    f"blastp;evalue={BLAST_EVALUE};max_target_seqs={BLAST_MAX_TARGET_SEQS};"
    f"min_pident={MIN_PIDENT};min_qcov={MIN_QCOV};min_scov={MIN_SCOV};"
    f"max_evalue={MIN_EVALUE};min_bitscore={MIN_BITSCORE}"
)

REPORT_PATH = OUT_DIR / f"{GUEST_ORGANISM}_{HOST_ORGANISM}_homology_report.txt"

# Columns produced by `blastp -outfmt 6` (default)
BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]
NUMERIC_COLS = [
    "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]


# -----------------------
# Parser
# -----------------------
def parse_genbank_proteins(gb_path: Path):
    """Extract CDS translations + annotation rows from a GenBank file.

    Pseudogenes (CDS without a /translation qualifier) are skipped and logged.
    Duplicate primary_ids are disambiguated with a numeric suffix.
    """
    fasta_records = []
    ann_rows = []
    skipped = 0
    id_counts: dict[str, int] = {}

    for record in SeqIO.parse(str(gb_path), "genbank"):
        seqid = record.id

        for feat in record.features:
            if feat.type != "CDS":
                continue

            q = feat.qualifiers
            translation = q.get("translation", [""])[0]
            if not translation:
                skipped += 1
                continue

            locus_tag = q.get("locus_tag", [""])[0]
            gene = q.get("gene", [""])[0]
            old_locus_tag = q.get("old_locus_tag", [""])[0]
            protein_id = q.get("protein_id", [""])[0]
            product = q.get("product", [""])[0]

            # NB: for compound locations (joins / programmed frameshifts) start/end
            # span the whole joined range. Adequate for annotation, not for slicing.
            start_1b = int(feat.location.start) + 1
            end_1b = int(feat.location.end)
            strand = "+" if feat.location.strand == 1 else "-" if feat.location.strand == -1 else "."

            base_id = locus_tag or protein_id or f"{seqid}_{start_1b}_{end_1b}"
            n_seen = id_counts.get(base_id, 0)
            primary_id = base_id if n_seen == 0 else f"{base_id}_{n_seen + 1}"
            id_counts[base_id] = n_seen + 1

            header = f"{primary_id}|{gene}|{seqid}"
            fasta_records.append((header, translation))

            ann_rows.append({
                "primary_id": primary_id,
                "locus_tag": locus_tag,
                "gene": gene,
                "old_locus_tag": old_locus_tag,
                "protein_id": protein_id,
                "product": product,
                "seqid": seqid,
                "start_1b": start_1b,
                "end_1b": end_1b,
                "strand": strand,
                "protein_len": len(translation),
                "fasta_header": header,
            })

    if skipped:
        logging.info("%s: skipped %d CDS without /translation (pseudogenes)",
                     gb_path.name, skipped)

    ann_df = pd.DataFrame(ann_rows)
    return fasta_records, ann_df


def write_fasta(records, out_faa: Path) -> None:
    with out_faa.open("w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


# -----------------------
# BLAST runner
# -----------------------
def _require_blast_binaries() -> None:
    for tool in ("makeblastdb", "blastp"):
        if shutil.which(tool) is None:
            raise RuntimeError(f"Required executable not found in PATH: {tool}")


def make_blast_db(faa: Path, db_prefix: Path) -> None:
    """Build a protein BLAST database if it's missing or stale."""
    phr = db_prefix.with_suffix(db_prefix.suffix + ".phr")
    if phr.exists() and phr.stat().st_mtime >= faa.stat().st_mtime:
        logging.info("BLAST DB up to date: %s", db_prefix)
        return

    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["makeblastdb", "-in", str(faa), "-dbtype", "prot", "-out", str(db_prefix)]
    logging.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_blastp(query_faa: Path, db_prefix: Path, out_tsv: Path) -> None:
    """Run blastp -outfmt 6 with the columns expected by the downstream parser."""
    if out_tsv.exists() and not FORCE_BLAST:
        logging.info("BLAST output exists, skipping: %s (set FORCE_BLAST=True to redo)",
                     out_tsv)
        return

    outfmt = "6 " + " ".join(BLAST_COLS)
    cmd = [
        "blastp",
        "-query", str(query_faa),
        "-db", str(db_prefix),
        "-evalue", str(BLAST_EVALUE),
        "-max_target_seqs", str(BLAST_MAX_TARGET_SEQS),
        "-num_threads", str(BLAST_THREADS),
        "-outfmt", outfmt,
        "-out", str(out_tsv),
    ]
    logging.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


# -----------------------
# BLAST loader / RBH
# -----------------------
def load_blast_tsv(path: Path, q_len_map: dict, s_len_map: dict) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", names=BLAST_COLS)
    for c in NUMERIC_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Replace NaN evalue with +inf so ascending sort pushes missing values last.
    df["evalue"] = df["evalue"].fillna(np.inf)

    df["q_len"] = df["qseqid"].map(q_len_map)
    df["s_len"] = df["sseqid"].map(s_len_map)
    df["qcov"] = 100.0 * df["length"] / df["q_len"]
    df["scov"] = 100.0 * df["length"] / df["s_len"]
    return df


def best_hit_per_query(df: pd.DataFrame) -> pd.DataFrame:
    """Highest bitscore wins; evalue then pident break ties deterministically."""
    ordered = df.sort_values(
        ["qseqid", "bitscore", "evalue", "pident", "sseqid"],
        ascending=[True, False, True, False, True],
        kind="mergesort",
    )
    return ordered.groupby("qseqid", as_index=False).first()


def write_report(
    report_path: Path,
    guest_ann: pd.DataFrame,
    host_ann: pd.DataFrame,
    rbh_df: pd.DataFrame,
    strict: pd.DataFrame,
) -> None:
    """Write a human-readable summary of the homology-build run."""
    guest_tag, host_tag = GUEST_ORGANISM, HOST_ORGANISM
    n_guest = len(guest_ann)
    n_host = len(host_ann)
    n_best = len(rbh_df)
    n_rbh_all = int(rbh_df["is_rbh"].sum())
    n_strict = len(strict)

    pid = strict["pident"]
    bits = strict["bitscore"]
    gcov = strict[f"{guest_tag}_cov_pct"]
    hcov = strict[f"{host_tag}_cov_pct"]

    def fmt_stats(s: pd.Series) -> str:
        if s.empty:
            return "n/a"
        return (f"median={s.median():.1f}  mean={s.mean():.1f}  "
                f"min={s.min():.1f}  max={s.max():.1f}")

    # Identity distribution bins — useful for spotting twilight-zone mass
    if not pid.empty:
        bins = [0, 25, 30, 40, 60, 80, 100.01]
        labels = ["<25", "25-30", "30-40", "40-60", "60-80", ">=80"]
        id_hist = pd.cut(pid, bins=bins, labels=labels, right=False).value_counts().reindex(labels, fill_value=0)
        id_hist_str = "\n".join(f"    {lab:>6}%  {int(n)}" for lab, n in id_hist.items())
    else:
        id_hist_str = "    (empty)"

    guest_rbh_ids = set(strict[f"{guest_tag}_id"])
    host_rbh_ids = set(strict[f"{host_tag}_id"])
    guest_cov_frac = 100.0 * len(guest_rbh_ids) / n_guest if n_guest else 0.0
    host_cov_frac = 100.0 * len(host_rbh_ids) / n_host if n_host else 0.0

    lines = [
        "=" * 70,
        f" Homology build report: {guest_tag} vs {host_tag}",
        "=" * 70,
        f"Run (UTC):          {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        f"Inputs:             {GUEST_GB.name}, {HOST_GB.name}",
        f"Output directory:   {OUT_DIR}",
        "",
        "Parameters",
        "----------",
        f"  blastp evalue:        {BLAST_EVALUE}",
        f"  max_target_seqs:      {BLAST_MAX_TARGET_SEQS}",
        f"  threads:              {BLAST_THREADS}",
        f"  strict min %identity: {MIN_PIDENT}",
        f"  strict min qcov/scov: {MIN_QCOV} / {MIN_SCOV}",
        f"  strict max evalue:    {MIN_EVALUE}",
        f"  strict min bitscore:  {MIN_BITSCORE}",
        "",
        "Protein counts",
        "--------------",
        f"  {guest_tag} proteins (guest): {n_guest}",
        f"  {host_tag} proteins (host):  {n_host}",
        "",
        "Homology results",
        "----------------",
        f"  Best-hit pairs (guest -> host):   {n_best}",
        f"  Reciprocal best hits (any qual): {n_rbh_all}",
        f"  Strict RBH (pass all filters):   {n_strict}",
        f"  Guest proteins with strict RBH:  {len(guest_rbh_ids)} / {n_guest} ({guest_cov_frac:.1f}%)",
        f"  Host  proteins with strict RBH:  {len(host_rbh_ids)} / {n_host} ({host_cov_frac:.1f}%)",
        "",
        "Strict-RBH quality (distribution across kept pairs)",
        "---------------------------------------------------",
        f"  %identity:       {fmt_stats(pid)}",
        f"  bitscore:        {fmt_stats(bits)}",
        f"  {guest_tag} coverage %: {fmt_stats(gcov)}",
        f"  {host_tag} coverage %:  {fmt_stats(hcov)}",
        "",
        "  %identity histogram (strict RBH):",
        id_hist_str,
        "",
        "Outputs",
        "-------",
        f"  {GUEST_FAA.name}",
        f"  {HOST_FAA.name}",
        f"  {GUEST_ANN.name}",
        f"  {HOST_ANN.name}",
        f"  {GUEST_VS_HOST.name}  (blastp guest -> host)",
        f"  {HOST_VS_GUEST.name}  (blastp host -> guest)",
        f"  {OUT_ALL.name}        (all best-hit pairs + RBH flag)",
        f"  {OUT_RBH.name}        (strict RBH subset)",
        "",
        "Note: RBH with BLASTP is reliable for clear orthologs but can miss",
        "twilight-zone homologs (<30% identity) and domain-fusion cases.",
        "Consider eggNOG-mapper or HMMER/Pfam to recover those.",
        "=" * 70,
        "",
    ]

    text = "\n".join(lines)
    report_path.write_text(text)
    # Also echo to stdout so it shows up in job logs
    print(text)


def rename_annotation(ann: pd.DataFrame, tag: str) -> pd.DataFrame:
    cols = ["fasta_header", "primary_id", "locus_tag", "gene", "product",
            "seqid", "start_1b", "end_1b", "strand"]
    mapping = {
        "fasta_header": f"{tag}_id",
        "primary_id": f"{tag}_primary_id",
        "locus_tag": f"{tag}_locus_tag",
        "gene": f"{tag}_gene",
        "product": f"{tag}_product",
        "seqid": f"{tag}_seqid",
        "start_1b": f"{tag}_start_1b",
        "end_1b": f"{tag}_end_1b",
        "strand": f"{tag}_strand",
    }
    return ann[cols].rename(columns=mapping)


def build_rbh_table(
    fwd_best: pd.DataFrame,
    rev_best: pd.DataFrame,
    guest_tag: str,
    host_tag: str,
) -> pd.DataFrame:
    rev_lookup = dict(zip(rev_best["qseqid"], rev_best["sseqid"]))

    rows = []
    for row in fwd_best.itertuples(index=False):
        g_id = row.qseqid
        h_id = row.sseqid
        rows.append({
            f"{guest_tag}_id": g_id,
            f"{host_tag}_id": h_id,
            "is_rbh": rev_lookup.get(h_id) == g_id,
            "pident": row.pident,
            "align_len": row.length,
            f"{guest_tag}_len": row.q_len,
            f"{host_tag}_len": row.s_len,
            f"{guest_tag}_cov_pct": row.qcov,
            f"{host_tag}_cov_pct": row.scov,
            "evalue": row.evalue,
            "bitscore": row.bitscore,
        })
    return pd.DataFrame(rows)


# -----------------------
# Main
# -----------------------
def main() -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    for gb in (HOST_GB, GUEST_GB):
        if not gb.exists():
            logging.error("GenBank input missing: %s", gb)
            return 1

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ---- Extract proteins ----
    guest_records, guest_ann = parse_genbank_proteins(GUEST_GB)
    host_records, host_ann = parse_genbank_proteins(HOST_GB)

    write_fasta(guest_records, GUEST_FAA)
    write_fasta(host_records, HOST_FAA)
    guest_ann.to_csv(GUEST_ANN, sep="\t", index=False)
    host_ann.to_csv(HOST_ANN, sep="\t", index=False)

    logging.info("%s proteins: %d", GUEST_ORGANISM, len(guest_ann))
    logging.info("%s proteins: %d", HOST_ORGANISM, len(host_ann))
    for p in (GUEST_FAA, HOST_FAA, GUEST_ANN, HOST_ANN):
        logging.info("Wrote: %s", p)

    # ---- Reciprocal BLASTP ----
    _require_blast_binaries()
    guest_db = BLAST_DB_DIR / f"{GUEST_ORGANISM}_prot_db"
    host_db = BLAST_DB_DIR / f"{HOST_ORGANISM}_prot_db"
    make_blast_db(GUEST_FAA, guest_db)
    make_blast_db(HOST_FAA, host_db)
    run_blastp(GUEST_FAA, host_db, GUEST_VS_HOST)
    run_blastp(HOST_FAA, guest_db, HOST_VS_GUEST)

    guest_len = dict(zip(guest_ann["fasta_header"], guest_ann["protein_len"]))
    host_len = dict(zip(host_ann["fasta_header"], host_ann["protein_len"]))

    fwd = load_blast_tsv(GUEST_VS_HOST, guest_len, host_len)
    rev = load_blast_tsv(HOST_VS_GUEST, host_len, guest_len)

    fwd_best = best_hit_per_query(fwd)
    rev_best = best_hit_per_query(rev)

    rbh_df = build_rbh_table(fwd_best, rev_best, GUEST_ORGANISM, HOST_ORGANISM)

    guest_meta = rename_annotation(guest_ann, GUEST_ORGANISM)
    host_meta = rename_annotation(host_ann, HOST_ORGANISM)

    rbh_df = rbh_df.merge(guest_meta, on=f"{GUEST_ORGANISM}_id", how="left")
    rbh_df = rbh_df.merge(host_meta, on=f"{HOST_ORGANISM}_id", how="left")

    # Provenance: parameters + run timestamp for reproducibility
    rbh_df["source"] = f"blastp;{BLAST_PARAMS_TAG};run_utc={datetime.now(timezone.utc).isoformat(timespec='seconds')}"

    rbh_df.to_csv(OUT_ALL, sep="\t", index=False)

    strict = rbh_df[
        rbh_df["is_rbh"]
        & (rbh_df["pident"] >= MIN_PIDENT)
        & (rbh_df[f"{GUEST_ORGANISM}_cov_pct"] >= MIN_QCOV)
        & (rbh_df[f"{HOST_ORGANISM}_cov_pct"] >= MIN_SCOV)
        & (rbh_df["evalue"] <= MIN_EVALUE)
        & (rbh_df["bitscore"] >= MIN_BITSCORE)
    ].copy()

    strict = strict.sort_values(
        ["pident", "bitscore", f"{GUEST_ORGANISM}_cov_pct", f"{HOST_ORGANISM}_cov_pct"],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)

    strict.to_csv(OUT_RBH, sep="\t", index=False)

    logging.info("All best-hit pairs: %d", len(rbh_df))
    logging.info("Strict RBH homolog pairs: %d", len(strict))
    logging.info("Wrote: %s", OUT_ALL)
    logging.info("Wrote: %s", OUT_RBH)

    write_report(REPORT_PATH, guest_ann, host_ann, rbh_df, strict)
    logging.info("Wrote: %s", REPORT_PATH)
    return 0


if __name__ == "__main__":
    sys.exit(main())
