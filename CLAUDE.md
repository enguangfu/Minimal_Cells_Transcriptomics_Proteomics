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
├── Syn1_RNase/                          R2: RNA-processing / ribonuclease analysis + B.subtilis→syn1 RNase-site mapping
│
├── Genome_Reduction/                    syn1 → syn3A comparison; the central downstream layer
│
└── Manuscript/                          LaTeX manuscript (sections/, figures/, SI/) + references.bib
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
  - Use wild-type Syn3A log phase in Sandberg *et al.*, iScience paper.
  - Note: SRR19432056 = R1 mate, SRR19432057 = R2 mate of ONE paired-end run (the depositor uploaded the two mates as separate SRA accessions); `00_retrive_fastq.sh` downloads both and renames.
  - Library: Kapa Stranded RNA-Seq + Ribo-Zero Gram-Negative (dUTP / fr-firststrand).
  - Read length of fixed 51 nts.
- `Illumina/Illumina_Processing/` — FastQC, bowtie2-build (one-shot), bowtie2 paired-end alignment (dUTP), per-strand bedGraphs in `depth_bedgraph/`.
- `Isoform_Cluster/` — output of `Cluster_Isoform_Syn3A.py` (ONT isoform clusters), shared with the segmentation pipeline.
  - **Output:** `isoform_clusters_annotated.tsv`.
- `Gene_TPM/` — `Syn3A_TPM.py` computes per-gene TPM from the syn3A Illumina + ONT bedGraphs.
  - **Output:** `syn3A_TPM_Illumina_ONT.tsv` + two correlation PDFs (Illumina vs ONT; our Illumina vs Palsson).
  - `Processed_TPM_Palsson/` contains the reported TPMs from the iScience paper for cross-validation (Pearson r ≈ 0.998 vs our computation).

### Syn1_Operon

Operon segmentation + annotation + visualization for syn1, driven by PacBio isoforms.
- `Operon_Segmentation.py` — containment clustering + read-evidence co-transcription merge → **`operons.candidate_blocks.tsv`** (canonical **459-operon** map; `operon_id, chrom, strand, start0, end0, n_isoforms, n_reads_total, member_ids, tss, tts, sense_gene_names, segmentation_type, ...`; `OP_*` IDs).
- `Operon_Annotation.py` — analysis-only: canonical-operon promoter (−10/−35, via `promoter_motif.py`) and terminator (TransTermHP) signatures + R1 figures; persists `annotation/canonical/{promoter_minus10_classification,terminator_tts_classification,operon_utr}.tsv`.
- `build_operon_xlsx.py` — pure assembly → **`operon.xlsx`** (SI Supplementary Data S1): joins the operon map + the promoter/terminator signature tsvs + `protein_complexes.xlsx` (curated complex→loci); promoter/terminator filled for the 127 canonical operons only.
- `Operon_Visualization.py` — `operon_plots/OP_*.pdf` (+ `_wdepth.pdf`).

### Syn3A_Operon

Same pipeline, syn3A flavour:
- `Operon_Segmentation_Syn3A.py` produces `operons.candidate_blocks.tsv` (operon IDs use `OP3A_*` prefix).
- `Operon_Annotation_Syn3A.py` builds ORF-coverage and multiplicity reports.
- `Operon_Visualization_syn3A.py` produces `operon_plots/`.

### Syn1_Syn3A_Proteomics

Proteomics tables for both organisms (used by `Syn1_Corr_RNA_Proteins/` and `Genome_Reduction/`).
- **Syn1 (absolute copy numbers):** `syn1_proteomics_localization_2026.csv` — use only when absolute quantities are explicitly needed.
- **Syn3A (relative iPM + absolute copy numbers in 2026):** `syn3a_proteomics_summary_2026.csv`.
- **syn3A (tertiary function annotation + proteomics number in 2019 and 2026):** `syn3A_proteome_annotated.xlsx` sheet Syn3A_Proteome
- **Annotation report:** `report_annotation_stats_syn3A.py` → builds `syn3A_proteome_annotated.xlsx` (derived, reordered + Protein Sequence + Exp. Ptn. Cnt 2019/2026) and `syn3A_tertiary_function_composition.html` (self-contained interactive: clickable composition bars → filterable/sortable protein table with sticky first 3 cols + CSV/Excel download). Pairs validate against `Syn3A_annotation/function_hierachy.tsv` (controlled vocab).

### Syn1_Corr_RNA_Proteins

Syn1 RNA × protein correlation.
- `Transcription_Translation.py` joins syn1 PacBio/Illumina TPMs with proteomics.
- `Translation_Residual_L2_elongation.py` explains residuals.
- **Combined table:** `syn1_genes_transcriptomics_proteomics.csv` — relative iPM values only.

### Syn1_Novel_ORF

Abnormal-transcription / novel-ORF discovery from PacBio isoforms.

### Syn1_RNase

R2 (RNA-processing / ribonuclease) analysis: endpoint-context erosion + ribonuclease inventory + B.subtilis→syn1 cleavage-site mapping.
- `RNA_Processing.py` — endpoint-context erosion from PacBio isoforms (per-isoform 5′/3′ intragenic-vs-intergenic labelling → 4 categories; ORF start-without-stop pass). Output → `RNase/isoform_endpoint_context.tsv`, `RNA_Processing.txt`.
- `R2_figure_panels.py` + `R2_legend_strip.py` → `R2_panels/` (Fig 2 born-at-size panels).
- `fold_3prime_terminator.py` — 3′-terminus 2° structure (ViennaRNA) for the 0178 intrinsic terminator (panel d).
- `RNase_Site_Mapping/` — self-contained B.subtilis→syn1 RNase-site transfer (inputs in `inputs/`; run in the **RNAseq** conda env — needs ViennaRNA `RNA`, biopython, BLAST+):
  - `map_bsub_rnase_to_syn1.py` — reciprocal-best-hit BLASTP homology (reuses `Genomes_Input/homology_syn1_bsub/`) transfers Taggart RNase III/Y sites; whole-gene ViennaRNA fold + homology-anchored facing-duplex test (18 RNase III genes; 0/5 paired confirm a conserved duplex). → `output/rnaseIII/rnaseIII_syn1_predicted_cleavage_pairs.tsv`.
  - `render_rnaseIII_structure.py` — local 2° structure at the homology-mapped cuts → `output/rnaseIII/stems/`.

### Genome_Reduction

Compare how the syn1 → syn3A reduction reshapes transcription. The pipeline is
numbered in dependency order; each downstream script writes to its own folder
**directly under `Genome_Reduction/`** (only `01–03` write to `aln/`). The
central concept is the **deletion junction** (05): every deletion is reframed as
a junction between the two retained operon fragments it joins, and the
co-expression scripts (06, 07) validate that structure against syn3A reads.

- `01_align.sh` → `02_analyze.py` → `03_visualize.py` — nucmer/dnadiff align syn3A→syn1, build the canonical event table (`aln/analysis/genome_reduction_summary.{xlsx,txt}`; `events` sheet, filter `Change Case == deleted`), and a circular Plotly map.
- `04_deletion_overlaid_operon.py` — overlap the 95 deletions with the 459 operons at single-bp resolution → per-operon truncation (`overlap_class`) and gene-deletion (`gene_deletion_pattern`) classes.
- `05_deletion_junction.py` — **the junction taxonomy**: `strand_relationship` ∈ {tandem, convergent, divergent, intra_operon} and (tandem only) `junction_type` ∈ {fusion, decapitation, readthrough_extension, clean_excision}, from facing-regulator loss; consistency-checked vs 04.
- `06_single_operon_coexpression.py` / `07_operon_pair_coexpression.py` — read-based co-transcription tests (ONT spanning/bridging + Illumina gap depth): operon-internal pairs (06; pristine controls + intra_operon) and cross-junction operon pairs (07; stratified by `junction_type`, fusion = predict yes, clean_excision = negative control).
- `08_delete_gene.py` — per retained gene: neighbors + unaltered-bp context + **`gene_impact_class`** (transcriptional impact by promoter-source change, integrating 04–07; precedence `promoter_lost > promoter_disconnected > new_promoter_fusion > readthrough_exposed > promoter_proximity_changed > context_only > unaffected`).
- `09_Compare_RNA_Protein.py` — builds the paired **mean-normalized** RNA(TPM)+protein(iPM) change table and owns the RNA/TPM figures + outliers; computes PTR (relIPM/relTPM proxy) and per-curated-category (Primary/Secondary/Tertiary) shares + story plots. **TPM policy:** Illumina for both organisms (ONT only QC scatters + non-coding syn3A fallback).
- `10_Compare_Ptn.py` — protein(iPM) counterpart (reads 09's CSV, no recompute): iPM figures, PTR-by-Secondary, proteome story plots, and `macromolecule_complex_abundance.tsv` (limiting-subunit RNAP + degradosome estimates). Run after 09.
- `coexpression_common.py` — shared 06/07 primitives (loaders, ONT spanning/bridging counters, `test_pair`); `Operon_Comparison_Syn1_Syn3A.py` — importable syn1↔syn3A comparison plotters.

#### Key outputs (cross-referenced often)

- `Genome_Reduction/aln/raw/syn1_deleted_regions.bed` — 95 syn1 deletion intervals (≥ 50 bp), the authoritative deletion list.
- `Genome_Reduction/aln/analysis/genome_reduction_summary.xlsx` — `events` sheet is the canonical per-block table (filter `Change Case == deleted` for the lost-gene set).
- `Genome_Reduction/deletion_overlaid_operon/operon_deletion_classification.tsv` — operon × deletion crosstable (truncation + gene-deletion patterns); read by 05/06/07/08.
- `Genome_Reduction/deletion_junction/deletion_junctions.tsv` — one row per deletion: operon_L/R, strand_relationship, junction_type, facing-regulator loss.
- `Genome_Reduction/delete_gene/retained_gene_context.tsv` — retained-gene neighbor / unaltered-bp table with the `gene_impact_class` column (the gene-level reduction-impact classification).
- `Genome_Reduction/Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv` — paired expression table (mean-normalized; TPM = Illumina): `relTPM_syn1, rank_syn1_TPM, relTPM_syn3a, rank_syn3a_TPM, TPM_fold_change, TPM_abs_change, relIPM_syn1, relIPM_syn3a, iPM_fold_change, iPM_abs_change, PTR_syn1, PTR_syn3a, PTR_fold_change, rna_type, sense_covering_ops`. Join key `locus_syn1 = MMSYN1_NNNN`.
- `Genome_Reduction/{single_operon_coexpression,operon_pair_coexpression}/` — ONT read-evidence for operon-internal and cross-junction co-transcription.

> **Run note:** 06/07 import `Operon_Comparison_Syn1_Syn3A.py`; run them with `Genome_Reduction/` as the working directory.

---

## Manuscript

LaTeX manuscript (`Manuscript/main.tex`, Nature Microbiology); per-section sources in `sections/{results,methods}/`. `MANUSCRIPT.md` (project root) is the drafting guide. SI assembly under `Manuscript/SI/`:
- `build_S4_qc.py` — zips the four library RNA-QC PDFs + a README → `Supplementary_Data_S4_QC.zip`.
- `gather_SI_excel.sh` — copies the three data workbooks (operon / syn1_omics / genome_reduction `.xlsx`) → `Supplementary_Data_S{1,2,3}_*.xlsx`.


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
- GFFs are 1-based inclusive on disk → subtract 1 from `start` when loading (see `load_gff_genes` in `Genome_Reduction/08_delete_gene.py` or `coexpression_common.py` for the canonical pattern).
- Genomes are circular — neighbor lookup and intergenic-distance calculations should wrap around the chromosome end.

### Deletion-corrected cross-organism comparisons

When comparing rel-units between syn1 and syn3A (per-category mRNA-pool shares,
complex abundance, etc.), renormalize the syn1 side to the **retained-gene pool**
(loci kept in syn3A). syn1's full-pool mean is diluted by the ~420 deleted genes,
which inflates retained-gene `rel*` values and biases the syn3A/syn1 ratio. Applied
inside `09_Compare_RNA_Protein.py`'s `_function_category_tpm_analysis` and both
story plots (composition + dumbbell) in 09 and 10. The main `TPM_FC_vs_absChange`
and `_impact_class_boxplot` still use the original full-pool normalization.

### Output conventions

Read OUTPUT.md.