# Genome Reduction Pipeline — JCVI-Syn1.0 → JCVI-Syn3A

Three-step pipeline that maps the genome-reduction events between Syn1 and
Syn3A and produces an interactive circular visualization.

```
01_align.sh        nucmer + dnadiff alignment
02_analyze.py      analysis & multi-sheet Excel + narrative report
03_visualize.py    interactive HTML (Plotly polar plot)
```

---

## Inputs

All paths relative to the project root (one level above this folder):

| File | Purpose |
|---|---|
| `Genomes_Input/syn1_genome.fasta` | Syn1 chromosome (CP002027.1, 1,078,809 bp) |
| `Genomes_Input/syn3A_genome.fasta` | Syn3A chromosome (CP016816.2, 543,379 bp) |
| `Genomes_Input/syn1.genes.gff3` | Syn1 gene annotation |
| `Genomes_Input/syn3a_genome.gff3` | Syn3A gene annotation |

---

## Environment

The pipeline uses the `RNAseq` conda env recipe. Each script resolves its
binaries from `$PATH`, so just activate the env before running:

```bash
conda activate RNAseq
```

Required binaries (resolved via `which`): `nucmer`, `delta-filter`,
`show-coords`, `show-snps`, `dnadiff`, `samtools`, `bedtools`.
Required Python packages: `pandas`, `openpyxl`, `plotly`.

---

## Run

```bash
bash 01_align.sh        # ~1 min
python 02_analyze.py    # ~5 s
python 03_visualize.py  # ~2 s
```

After step 3, open `aln/analysis/genome_reduction_circle.html` in a browser.

---

## Outputs

The three deliverables live in `aln/analysis/` and are the only files most
readers ever need to open:

| File | What it is |
|---|---|
| `genome_reduction_summary.xlsx` | **Canonical analysis table.** 4 sheets — see below. |
| `genome_reduction_summary.txt`  | Narrative report (Parts A–D: deletions, insertions, relocations, biology). |
| `genome_reduction_circle.html`  | Interactive Plotly polar plot — outer ring = Syn1, inner ring = Syn3A. Hover any segment for gene content. |

### Excel sheet map

| Sheet | One row per … | Key columns |
|---|---|---|
| `events`           | change event | `block_index_syn1`, `S1`, `E1`, `S2`, `E2`, `LEN1`, `LEN2`, `PCT_IDY`, `Change Case`, `Syn1_genes`, `Syn3A_genes` |
| `short_insertions` | small (<1 kb) qdiff insertion | overlapping gene + nearest gene on each strand |
| `dnadiff_summary`  | headline metric | totals, identity, SNPs, indels, relocations, etc. |
| `legend`           | term / column / case | definitions |

Each `events` row is one of four `Change Case` values:

| Change Case | Meaning |
|---|---|
| `retained_ordered`   | Syn1 block kept in Syn3A at the position predicted by linear order. |
| `retained_relocated` | Syn1 block kept but moved (LIS outlier). Currently 1 block: `lap` (MMSYN1_0154). |
| `deleted`            | Syn1 region not represented in Syn3A. |
| `inserted`           | Syn3A region with no homolog in Syn1. The 1.1 kb entry carries `met14p` (JCVISYN3A_0931). |

`Syn1_genes` / `Syn3A_genes` are comma-separated, deduplicated locus-tag
suffixes. Format: `0001/dnaA_1` when there's a real gene name, just `0005`
when there isn't, `.` when intergenic.

---

## Layout reference

```
Genome_Reduction/
├── 01_align.sh              # nucmer + dnadiff
├── 02_analyze.py            # builds the Excel + .txt
├── 03_visualize.py          # builds the HTML
├── README.md                # this file
└── aln/
    ├── run.log              # 01_align.sh stdout/stderr
    ├── raw/                 # nucmer + dnadiff intermediates (rarely opened)
    │   ├── syn1_vs_syn3A.{delta,1delta}
    │   ├── syn1_vs_syn3A.coords.tsv     # alignment block table (96 rows)
    │   ├── syn1_vs_syn3A.snps.tsv       # per-base SNP/indel list
    │   ├── syn1_deleted_regions.bed     # raw deletion intervals
    │   ├── syn1_deleted_genes.tsv       # raw deletion ↔ gene join
    │   └── dnadiff_out.{report,qdiff,rdiff,1coords,mcoords,...}
    └── analysis/            # human-facing deliverables (3 files)
        ├── genome_reduction_summary.xlsx
        ├── genome_reduction_summary.txt
        └── genome_reduction_circle.html
```

If you only ever look at one thing, look at the Excel.
If you want the picture, open the HTML.
If you want the prose summary, open the .txt.

---

## Tunable parameters

In `01_align.sh`:

| Var | Default | What it controls |
|---|---|---|
| `THREADS` | 8 | nucmer / samtools parallelism |
| `MIN_DEL_BP` | 50 | drop deletions smaller than this |
| `MERGE_GAP` | 50 | close alignment gaps ≤ this before complement |
| `delta-filter -i 95 -l 100` | inline | identity / length cutoffs for retained blocks |

In `02_analyze.py`:

| Const | Default | What it controls |
|---|---|---|
| `MIN_INS_BP` | 10 | drop syn3A-only intervals smaller than this from the events sheet |
| `SHORT_INS_MAX` | 1000 | upper size for the `short_insertions` sheet |

---

## Re-running cleanly

`01_align.sh` writes only into `aln/raw/`. `02_analyze.py` and `03_visualize.py`
write only into `aln/analysis/`. To rebuild from scratch:

```bash
rm -rf aln && bash 01_align.sh && python 02_analyze.py && python 03_visualize.py
```
