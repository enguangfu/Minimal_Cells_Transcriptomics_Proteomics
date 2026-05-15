# Project: Minimal Cells Transcriptomics & Proteomics

Transcriptomics (PacBio, ONT, Illumina) and Proteomics Analysis of Synthetic Bacteria JCVI-syn1, and JCVI-syn3A

All paths are relative to the project root:
`/data/enguang/Transcriptomics/Minimal_Cells_Transcriptomics_Proteomics/`

---

## Directory Structure in Order 

### PacBio_Raw

Retrieve PacBio fastq files of Syn1 from NCBI

### PacBio_Processing

- Recover FLNC reads from raw PacBio HiFi cDNA reads (reorientation, removal of primers, trimming of polyAs)
- Map to syn1 ref genome with *minimap2*
- Quality control
- Output: BAM file of isoforms mapped to syn1 genome

### Isoform_Cluster

Cluster the PacBio isoforms
- Output: clustered PacBio isoforms of syn1 transcriptome

### Operon_Annotation_Visualization

Segement, annotate, and visualize syn1's operons with PacBio isoforms
- Output: operons of syn1

### Transcriptomics_Quantification

Quantification of TPMs of different genes in syn1 (PacBio, Illumina) and syn3A (ONT)

### Proteomics_Quantification

Quantifications of iPMs of different proteins in syn1 and syn3A

### ONT_Raw

Retrieve ONT fastq files of syn3A from NCBI

### ONT_Processing

Map ONT isoforms to syn3A ref genome

### Syn3A_Illumina

Map Short Illumina reads to syn3A genome

### Novel_ORF_Discovery

Find abnormal transcription and novel ORFs in syn1 using PacBio isoforms 

### Genome_Reduction

Compare how genome reduction affect the transcription and translation in syn3A

#### Pipeline

Three-step alignment + analysis, plus two follow-on scripts:

- `01_align.sh` — nucmer + dnadiff → `aln/raw/`
- `02_analyze.py` — builds `aln/analysis/genome_reduction_summary.{xlsx,txt}` (canonical event table: deletions, insertions, relocations)
- `03_visualize.py` — interactive Plotly circular map
- `05_deletion_operon.py` — overlaps the 95 deletions with `operons.candidate_blocks.tsv` at single-bp resolution. Two complementary classifications per operon:
  - **Truncation pattern (operon-span level):** `overlap_class` and `per_hit_classes` ∈ {fully_deleted, 5'_truncation_gene, 5'_truncation_UTR, 3'_truncation_gene, 3'_truncation_UTR, intra_truncated} or `multi:…`
  - **Deletion pattern (gene level):** `gene_deletion_pattern` ∈ {intact, all_deleted, leading_deleted, lagging_deleted, intra_deleted, fully_deleted}; per-gene retained/partial/fully columns. Output: `aln/analysis/operon_deletion_classification.tsv` + per-category visualization PDFs (with optional `_wdepth.pdf` versions interleaving `OP_*_wdepth.pdf` panels).
- `06_delete_gene.py` — for each retained syn1 gene, finds same-strand transcription-direction (upstream/downstream) and strand-agnostic (cw/ccw) neighbors in both syn1 and syn3A, computes `unaltered_cw_bps` / `unaltered_ccw_bps` in the syn1 frame (circular wrap), and adds `operon_change` ∈ {leading_gene_deleted, promoter_deleted, ""}. Output: `aln/analysis/retained_gene_context.tsv` + boxplot of TPM_fold_change by operon_change.

#### Key outputs (cross-referenced often)

- `Genome_Reduction/aln/raw/syn1_deleted_regions.bed` — 95 syn1 deletion intervals (≥ 50 bp), the authoritative deletion list.
- `Genome_Reduction/aln/analysis/genome_reduction_summary.xlsx` — `events` sheet is the canonical per-block table (filter `Change Case == deleted` for the lost-gene set).
- `Genome_Reduction/aln/analysis/operon_deletion_classification.tsv` — operon × deletion crosstable (truncation + gene-deletion patterns).
- `Genome_Reduction/aln/analysis/retained_gene_context.tsv` — retained-gene neighbor and unaltered-bp table with the `operon_change` tag.
- `Genome_Reduction/syn1_vs_syn3a_RNA_protein.csv` — paired expression table: `TPM_mean_syn1, TPM_mean_syn3A_ONT, TPM_fold_change_ONT, TPM_mean_syn3a_Illumina, TPM_fold_change_Illumina, iPM_mean_syn1, iPM_mean_syn3a, iPM_fold_change, rna_type, sense_covering_ops`. Join key is `locus_syn1 = MMSYN1_NNNN`.

### Operon_Syn3A

Analyze operonal structure in syn3A using ONT: cluster ONT isoforms then do the segmentations.


---

## Key Data Files

### Genome Inputs
- **Syn1 FASTA:** `Genomes_Input/syn1_genome.fasta` (CP002027.1, 1,078,809 bp, circular)
- **Syn1 GFF3:** `Genomes_Input/syn1.genes.gff3` — canonical gene annotation
- **Syn3A GenBank:** `Genomes_Input/syn3a.gb` — used to mark genes absent in Syn3A
- **Syn3A FASTA:** `Genomes_Input/syn3A_genome.fasta` (CP016816.2, 543,379 bp, circular)                                                                                                                                                 
- **Syn3A GFF3:** `Genomes_Input/syn3a_genome.gff3` — 493 `gene` + 3 `pseudogene` records. Pseudogenes (e.g. 0051, 0546, 0602) have feature type `pseudogene`, not `gene` — parse both.  
- **Genetic code:** Mycoplasma (UGA = Trp, not stop)
- **16S rRNA 3' tail (anti-SD):** `ACCUCCUUU`
 
### Transcriptomics of Syn1
- **PacBio BAM:** `./syn1.PacBio.FLNC.sorted.HQ.bam`
  - Mapped PacBio isoforms to syn1 after quality control 
- **Isoforms:** `isoform_annotation/isoform_clusters_annotated.tsv`
  - Columns: `isoform_id, chrom, strand, start0, end0, pos5p0, pos3p0, n_reads, ...`
  - `pos5p0` = 5' end (0-based); use for 5' UTR extraction
<!-- - **Isoform membership:** `isoform_annotation/isoform_membership.tsv`
  - Columns: `isoform_id, chrom, strand, read_id, pos5p0, pos3p0, start0, end0` -->
<!-- - **PacBio depth bedGraphs:** `PacBio_Processing/depth_bedgraph/`
  - `syn1.PacBio.FLNC.HQ.plus.bedGraph`
  - `syn1.PacBio.FLNC.HQ.minus.bedGraph`
  - Format: 4-col `chrom, start0, end0, depth` -->
- **Operons:** `Operon_Annotation_Visualization/operons.candidate_blocks.tsv`
  - Columns: `operon_id, chrom, strand, start0, end0, n_isoforms, n_reads_total, member_ids, tss, tts, sense_gene_names, ...`
- **Quantified Gene TPMs:** `Transcriptomoics_Quantification/syn1_Illumina_PacBio_TPM_profiles.csv`

### Transcriptomics of Syn3A
- **Direct ONT BAM:** `ONT_Processing/syn3A.ONT.rep1.sorted.bam`
- **Illumina from iScience paper:** `Syn3A_Illumia/bam/syn3A_rep1.sorted.bam`
- **Operons:** `Operon_Syn3A/operons.candidate_blocks.tsv`


### Proteomics of Syn1 and Syn3A
- **Absolute protein copy numbers:** `Proteomics_Quantification/syn1_proteomics_localization_2026.csv` — absolute copy numbers per protein for Syn1 (use only when absolute quantities are explicitly needed)
- **Transcriptomics + Proteomics table:** `Gene_Transcritpomics_Proteomics/syn1_genes_transcriptomics_proteomics.csv` — relative iPM values only
- **syn3A proteomics:**: `Proteomics_Quantification/syn3a_proteomics_summary_2026` - relative iPM and absolute copy numbers per protein for Syn3A.

### Syn1 ↔ Syn3A correspondence  

- Locus tags: syn1 uses `MMSYN1_NNNN`, syn3A uses `JCVISYN3A_NNNN`. **Numeric suffix is preserved** across the two annotations, so `MMSYN1_0025 ↔ JCVISYN3A_0025` is reliable (verified, only ~4 syn1 small-feature drops and 3 syn3A additions break the 1:1).                                                    
- Syn3A is essentially a Syn1 subset: ~536 kb removed in 95 discrete deletions; retained sequence at 99.90% identity (36 SNPs, 12 indels); 0 inversions/translocations; 1 block relocation (gene `lap`, MMSYN1_0154 → syn3A ~311.7 kb); 1 truly novel CDS (`JCVISYN3A_0931`, met14p). Two extra rpmG L33 paralogs newly annotated (`JCVISYN3A_0930`, `0932`).                                                                     
- **Mycoplasma genetic code applies to both.**   

---


## Key refs

Creation of a Bacterial Cell Controlled by a Chemically Synthesized Genome, 2010, Science

Design and synthesis of a minimal bacterial genome, 2016, Science

Essential metabolism for a minimal cell, 2019, eLife


<!-- ## Analysis Scripts
- `Gene_Transcritpomics_Proteomics/Translation_Residual.py` — three-level analysis to explain r=0.6 transcriptome–proteome correlation
- `Operon_Annotation_Visualization/Operon_Visualization.py` — operon visualization pipeline -->


---

## Notes
### Coordinate conventions

- All internal coordinates are **0-based half-open** (`start0, end0`).
- GFFs are 1-based inclusive on disk → subtract 1 from `start` when loading (see `load_gff_genes` in `Genome_Reduction/06_delete_gene.py` for the canonical pattern).
- Genomes are circular — neighbor lookup and intergenic-distance calculations should wrap around the chromosome end.
