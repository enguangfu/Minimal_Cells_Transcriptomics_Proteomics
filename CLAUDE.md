# Project: Minimal Cells Transcriptomics & Proteomics

Transcriptomics (PacBio, ONT, Illumina) and Proteomics Analysis of Synthetic Bacteria JCVI-syn1, and JCVI-syn3A

All paths are relative to the project root:
`/data/enguang/Transcriptomics/Minimal_Cells_Transcriptomics_Proteomics/`

---

## Directory Structure

The tree is grouped by organism / platform / role.

```
.
├── Genomes_Input/                       reference FASTAs + GFFs (both organisms)
│
├── Syn1_Transcriptomics/
│   ├── PacBio/
│   │   ├── PacBio_Raw/                  raw HiFi reads from NCBI
│   │   └── PacBio_Processing/           FLNC recovery → minimap2 → sorted BAM + depth bedGraphs
│   ├── Illumina/
│   │   ├── Illumina_Raw/                raw paired-end FASTQs
│   │   └── Illumina_Processing/         FastQC + bowtie2 + strand-split bedGraphs
│   ├── Isoforms_PacBio/                 PacBio isoform clusters (Cluster_Isoform.py)
│   └── Gene_TPM/                        per-gene TPM tables (syn1 PacBio + Illumina)
│
├── Syn3A_Transcriptomics/
│   ├── ONT/
│   │   ├── ONT_Raw/                     raw direct-RNA FASTQs
│   │   └── ONT_Processing/              minimap2 -ax map-ont → sorted BAM + depth bedGraphs
│   ├── Illumina/
│   │   ├── Illumina_Raw/                paired-end FASTQs (renamed from SRR pair)
│   │   └── Illumina_Processing/         FastQC + bowtie2 + strand-split bedGraphs
│   ├── Isoform_Cluster/                 ONT isoform clusters (Cluster_Isoform_Syn3A.py)
│   └── Gene_TPM/                        syn3A TPM (Illumina + ONT) + Palsson comparison
│
├── Syn1_Operon/                         syn1 operon segmentation, annotation, visualization
├── Syn3A_Operon/                        syn3A operon segmentation, annotation, visualization
│
├── Syn1_Syn3A_Proteomics/                proteomics tables for both organisms
├── Syn1_Corr_RNA_Proteins/               syn1 RNA × protein correlation analysis
├── Syn1_Novel_ORF/                      novel-ORF discovery from syn1 PacBio isoforms
│
└── Genome_Reduction/                    syn1 → syn3A comparison; the central downstream layer
```

### Genomes_Input

Reference genomes for both organisms. The Mycoplasma genetic code applies to both (UGA = Trp, not stop); 16S rRNA 3' tail (anti-SD): `ACCUCCUUU`.
- **Syn1 FASTA:** `syn1_genome.fasta` (CP002027.1, 1,078,809 bp, circular)
- **Syn1 GFF3:** `syn1.genes.gff3` — canonical gene annotation
- **Syn3A FASTA:** `syn3A_genome.fasta` (CP016816.2, 543,379 bp, circular)
- **Syn3A GFF3:** `syn3a_genome.gff3` — 493 `gene` + 3 `pseudogene` records (pseudogenes 0051, 0546, 0602 use feature type `pseudogene`; parse both).
- **Syn3A GenBank:** `syn3a.gb` — used to mark genes absent in Syn3A compared to Syn1.

### Syn1_Transcriptomics

- `PacBio/PacBio_Raw/` — raw PacBio HiFi reads from NCBI.
- `PacBio/PacBio_Processing/` — FLNC recovery (reorientation, primer removal, polyA trimming) → `minimap2 -ax map-hifi` → quality control.
  - **PacBio BAM:** `syn1.PacBio.FLNC.sorted.HQ.bam` (lives next to the script).
  - **PacBio depth bedGraphs:** `depth_bedgraph/syn1.PacBio.FLNC.HQ.{plus,minus,total}.bedGraph`.
- `Illumina/Illumina_Raw/` — paired-end FASTQs from SRA.
- `Illumina/Illumina_Processing/` — FastQC + MultiQC, bowtie2 paired-end alignment (dUTP / fr-firststrand → R2 = transcript strand), per-strand bedGraphs (`bam/`, `depth_bedgraph/`).
- `Isoforms_PacBio/` — `Cluster_Isoform.py` clusters PacBio isoforms.
  - **Isoform clusters:** `isoform_clusters_annotated.tsv` with columns `isoform_id, chrom, strand, start0, end0, pos5p0, pos3p0, n_reads, ...`. `pos5p0` = 5' end (0-based); use for 5' UTR extraction.
- `Gene_TPM/` — `Gene_Transcriptomics.py` computes per-gene sense/antisense TPM for syn1 PacBio + Illumina.
  - **Output:** `syn1_Illumina_PacBio_TPM_profiles.csv`.

### Syn3A_Transcriptomics

- `ONT/ONT_Raw/` — direct-RNA FASTQs.
- `ONT/ONT_Processing/` — `minimap2 -ax map-ont` (NOT splice; bacteria are intron-less) → sorted BAM + per-strand bedGraphs.
  - **ONT BAM:** `syn3A.ONT.rep1.sorted.bam`.
  - **ONT depth bedGraphs:** `depth_bedgraph/syn3A.ONT.rep1.{plus,minus,total}.bedGraph`.
- `Illumina/Illumina_Raw/` — paired-end FASTQs (`syn3A_rep1_{1,2}.fastq`).
  - Note: SRR19432056 = R1 mate, SRR19432057 = R2 mate of ONE paired-end run (the depositor uploaded the two mates as separate SRA accessions); `00_retrive_fastq.sh` downloads both and renames.
  - Library: Kapa Stranded RNA-Seq + Ribo-Zero Gram-Negative (dUTP / fr-firststrand).
- `Illumina/Illumina_Processing/` — FastQC, bowtie2-build (one-shot), bowtie2 paired-end alignment (dUTP), per-strand bedGraphs in `depth_bedgraph/`.
- `Isoform_Cluster/` — output of `Cluster_Isoform_Syn3A.py` (ONT isoform clusters), shared with the segmentation pipeline.
  - **Output:** `isoform_clusters_annotated.tsv`.
- `Gene_TPM/` — `Syn3A_TPM.py` computes per-gene TPM from the syn3A Illumina + ONT bedGraphs.
  - **Output:** `syn3A_TPM_Illumina_ONT.tsv` + two correlation PDFs (Illumina vs ONT; our Illumina vs Palsson).
  - `Processed_TPM_Palsson/` contains the reported TPMs from the iScience paper for cross-validation (Pearson r ≈ 0.998 vs our computation).

### Syn1_Operon

Operon segmentation + annotation + visualization for syn1, driven by PacBio isoforms.
- **Operons:** `operons.candidate_blocks.tsv` with columns `operon_id, chrom, strand, start0, end0, n_isoforms, n_reads_total, member_ids, tss, tts, sense_gene_names, ...`. Operon IDs use `OP_*` prefix.
- Visualization: `Operon_Visualization.py` produces `operon_plots/OP_*.pdf` and `_wdepth.pdf` variants.

### Syn3A_Operon

Same pipeline, syn3A flavour:
- `Operon_Segmentation_Syn3A.py` produces `operons.candidate_blocks.tsv` (operon IDs use `OP3A_*` prefix).
- `Operon_Annotation_Syn3A.py` builds ORF-coverage and multiplicity reports.
- `Operon_Visualization_syn3A.py` produces `operon_plots/`.

### Syn1_Syn3A_Proteomics

Proteomics tables for both organisms (used by `Syn1_Corr_RNA_Proteins/` and `Genome_Reduction/`).
- **Syn1 (absolute copy numbers):** `syn1_proteomics_localization_2026.csv` — use only when absolute quantities are explicitly needed.
- **Syn3A (relative iPM + absolute copy numbers):** `syn3a_proteomics_summary_2026.csv`.

### Syn1_Corr_RNA_Proteins

Syn1 RNA × protein correlation.
- `Transcription_Translation.py` joins syn1 PacBio/Illumina TPMs with proteomics.
- `Translation_Residual_L2_elongation.py` explains residuals.
- **Combined table:** `syn1_genes_transcriptomics_proteomics.csv` — relative iPM values only.

### Syn1_Novel_ORF

Abnormal-transcription / novel-ORF discovery from PacBio isoforms.

### Genome_Reduction

Compare how the syn1 → syn3A reduction reshapes transcription and translation. Three-step alignment + analysis pipeline plus four follow-on scripts:

- `01_align.sh` — nucmer + dnadiff → `aln/raw/`
- `02_analyze.py` — builds `aln/analysis/genome_reduction_summary.{xlsx,txt}` (canonical event table: deletions, insertions, relocations)
- `03_visualize.py` — interactive Plotly circular map
- `05_deletion_operon.py` — overlaps the 95 deletions with `Syn1_Operon/operons.candidate_blocks.tsv` at single-bp resolution. Two complementary classifications per operon:
  - **Truncation pattern (operon-span level):** `overlap_class` and `per_hit_classes` ∈ {fully_deleted, 5'_truncation_gene, 5'_truncation_UTR, 3'_truncation_gene, 3'_truncation_UTR, intra_truncated} or `multi:…`
  - **Deletion pattern (gene level):** `gene_deletion_pattern` ∈ {intact, all_deleted, leading_deleted, lagging_deleted, intra_deleted, fully_deleted}; per-gene retained/partial/fully columns. Output: `aln/analysis/operon_deletion_classification.tsv` + per-category visualization PDFs (with optional `_wdepth.pdf` versions interleaving `OP_*_wdepth.pdf` panels).
- `06_delete_gene.py` — for each retained syn1 gene, finds same-strand transcription-direction (upstream/downstream) and strand-agnostic (cw/ccw) neighbors in both syn1 and syn3A, computes `unaltered_cw_bps` / `unaltered_ccw_bps` in the syn1 frame (circular wrap), and adds `operon_change` ∈ {leading_gene_deleted, promoter_deleted, ""}. Output: `aln/analysis/retained_gene_context.tsv` + boxplot of TPM_fold_change by operon_change.
- `07_operon_change.py` — tests three biological questions against syn3A ONT reads:
  - **Q1/Q2** (preservation of syn1 operons in syn3A): for every consecutive retained pair, count spanning (strict) and bridging (loose) ONT reads. Operon-level verdict ∈ {preserved_strict, preserved_loose, split}.
  - **Q3** (new gene proximities): for every newly-adjacent same-strand syn3A pair (≥1 syn1 gene deleted between the ancestors), count spanning / bridging reads as co-transcription evidence.
  - Output: `operon_change/Q1Q2_pair_preservation.tsv`, `Q1Q2_operon_preservation.tsv`, `Q3_new_pair_candidates.tsv`, narrative `operon_change_summary.txt`, plus side-by-side syn1-vs-syn3A comparison PDFs in `operon_change/comparison_plots/{split_intact, leading_deleted, intra_deleted, lagging_deleted}/`.
- `Operon_Comparison_Syn1_Syn3A.py` — the side-by-side comparison plotter used by 07 (also importable standalone).
- `Compare_RNA_Protein.py` — joins TPM + proteomics tables across the two organisms.

#### Key outputs (cross-referenced often)

- `Genome_Reduction/aln/raw/syn1_deleted_regions.bed` — 95 syn1 deletion intervals (≥ 50 bp), the authoritative deletion list.
- `Genome_Reduction/aln/analysis/genome_reduction_summary.xlsx` — `events` sheet is the canonical per-block table (filter `Change Case == deleted` for the lost-gene set).
- `Genome_Reduction/aln/analysis/operon_deletion_classification.tsv` — operon × deletion crosstable (truncation + gene-deletion patterns).
- `Genome_Reduction/aln/analysis/retained_gene_context.tsv` — retained-gene neighbor and unaltered-bp table with the `operon_change` tag.
- `Genome_Reduction/syn1_vs_syn3a_RNA_protein.csv` — paired expression table: `TPM_mean_syn1, TPM_mean_syn3A_ONT, TPM_fold_change_ONT, TPM_mean_syn3a_Illumina, TPM_fold_change_Illumina, iPM_mean_syn1, iPM_mean_syn3a, iPM_fold_change, rna_type, sense_covering_ops`. Join key is `locus_syn1 = MMSYN1_NNNN`.
- `Genome_Reduction/operon_change/Q1Q2_*.tsv` and `Q3_new_pair_candidates.tsv` — ONT-spanning-evidence tests for operon preservation and new gene proximities.


---

## Syn1 ↔ Syn3A correspondence

- Locus tags: syn1 uses `MMSYN1_NNNN`, syn3A uses `JCVISYN3A_NNNN`. **Numeric suffix is preserved** across the two annotations, so `MMSYN1_0025 ↔ JCVISYN3A_0025` is reliable (verified, only ~4 syn1 small-feature drops and 3 syn3A additions break the 1:1).
- Syn3A is essentially a Syn1 subset: ~536 kb removed in 95 discrete deletions; retained sequence at 99.90% identity (36 SNPs, 12 indels); 0 inversions/translocations; 1 block relocation (gene `lap`, MMSYN1_0154 → syn3A ~311.7 kb); 1 truly novel CDS (`JCVISYN3A_0931`, met14p). Two extra rpmG L33 paralogs newly annotated (`JCVISYN3A_0930`, `0932`).

---


## Key refs

Creation of a Bacterial Cell Controlled by a Chemically Synthesized Genome, 2010, Science

Design and synthesis of a minimal bacterial genome, 2016, Science

Essential metabolism for a minimal cell, 2019, eLife




---

## Notes
### Coordinate conventions

- All internal coordinates are **0-based half-open** (`start0, end0`).
- GFFs are 1-based inclusive on disk → subtract 1 from `start` when loading (see `load_gff_genes` in `Genome_Reduction/06_delete_gene.py` for the canonical pattern).
- Genomes are circular — neighbor lookup and intergenic-distance calculations should wrap around the chromosome end.
