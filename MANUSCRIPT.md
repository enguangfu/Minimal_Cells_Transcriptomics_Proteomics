# MANUSCRIPT.md — drafting guide for Methods & Results

This file is the **single source of truth** the LLM uses to draft Methods and
Results paragraphs. Each subsection here mirrors a `.tex` file under
`Manuscript/sections/{results,methods}/`. Fill in the bullets; the LLM will
expand them into prose, quoting numbers verbatim from the files you list.

---

## 0. Style guide (read before drafting any subsection)

**Target journal:** Nature Microbiology  
**Length budget:** Results ≈ 3000; Methods ≈ unlimited (unlimited); Intro + Discussion ≈ 500  
**Tense / voice:**  
- Methods: past tense, passive acceptable ("Reads were aligned with …").
- Results: past tense for what we did, present for what the data show ("syn3A retains 99.90% identity …").  

**Person:** "we" not allowed in Results and Methods  
**Abbreviations:** define on first use in each of Abstract / Main / Methods.  
**Figure references:** Fig. 2b   
**Citations:** numeric, `\cite{key}`; bibfile `Manuscript/references.bib`.  
**Organism names:** "JCVI-syn1.0" first use then "syn1"; Same for syn3A  
**Number formatting:** sizes/coords — bp for exact coordinates (`1,078,809 bp`), kb/Mb for sizes (`536 kb`, `1.08 Mb`); read/feature counts — number + nonbreaking space + `k`/`M` (`2.6~M`, `267~k`), exact counts use `{,}` separators (`82{,}000`); p-values `P = 0.003`.  
**Software versions:** name + exact version used in the scripts, e.g. minimap2 v2.30, bowtie2 v2.5.5, samtools v1.22.1, FastQC, MultiQC (citation keys TBD).  
**Units / symbols:** micro prefix as math `$\mu$` (`$\mu$l`, `$\mu$g`, `$\mu$L`), not the raw `µ` glyph; temperatures as `\textdegree C` (not `°C`/`ºC`/`℃`); `®`/`™` as `\textregistered{}`/`\texttrademark`.  
**Species / gene / locus:** species names via macros `\Mmy` / `\mmy` (not raw `\emph{...}`); gene names italic via `\gene{name}{num}` (italic name / plain locus); locus tags (`MMSYN1_NNNN`, `JCVISYN3A_NNNN`) in **plain text** via `\locus` / `\locusA`.  
**Things to NEVER do:** _<e.g. don't write "novel" without qualification, don't claim causation from correlation, …>_  
- em-dash (`---`): use comma or parenthese instead. (en-dash `--` for numeric ranges, e.g. `0.88--3.07~kb`, is correct and retained)
- bullet points
- hedging-as-polteness such as "I think", "I feel"  

**Line Breaking in Latex:** One sentence per line for easier tracking  

### Exemplar paragraph (tone we want to match)
> _<paste one paragraph from a paper whose voice you want to emulate, or one of your own already-good paragraphs>_

---

## How to fill in each subsection

```
**Claim:** one sentence — the takeaway a reader should remember.
**Logic:** why this analysis answers the claim (1–3 sentences of scientific reasoning).
**Analysis:** scripts + key parameters (paths relative to repo root).
**Outputs:** tables/figures with paths.
**Numbers to cite:** the 2–5 values that MUST appear verbatim (n, %, p, fold-change…).
**Figure panels:** which panels of which manuscript figure this maps to.
**Conclusion:** what we conclude, and what we explicitly do *not* claim.
**Caveats / hedges:** limitations the prose must acknowledge.
**Notes for LLM:** anything special (e.g. "cite Sandberg 2023 here", "do not call ONT quantitative").
```

The **Numbers to cite** line is the single most important field — if listed, they get quoted; if absent, the LLM may pick the wrong column from the output table.

---

# RESULTS

## R1 — Gene Co-transcriptions as Operons in Syn1 from PacBio Long-read RNA Sequencing
**Tex file:** `Manuscript/sections/results/operons.tex`
**Figure:** `Manuscript/figures/operon.pdf`


How transcription and further RNA processing can complexity the transcriptome even of the reduced organism

FLNC reads clustered to isoforms: 2.6M by 10 bps to 205k clustered; Use 10 for visual of 50 for operon backbone;  distribution of sharpness 5' and 3' ends

Find the longest RNA isoform for representative of operons; then merge and rescue to get total 480 operons

Statistics on operon size: operon can cover anti-sense genes; longest is rPtn operon; operon length distribution;

Promoter and terminator signatures: -10 box significant; -35 box not; terminator predicted as hairpins + polyU

\textit{Protein complexes in operons}

tRNAs in operons

Put one operon here: DCW operon?

**Claim:**
**Logic:**
**Analysis:**
- Scripts: `Syn1_Operon/…`, `Syn1_Transcriptomics/Isoforms_PacBio/Cluster_Isoform.py`
- Key params:
**Outputs:**
- `Syn1_Operon/operons.candidate_blocks.tsv`
- `Syn1_Operon/operon_plots/`
**Numbers to cite:**
**Figure panels:**
**Conclusion:**
**Caveats / hedges:**
**Notes for LLM:**

---

## R2 — Pervasive and biased RNA processing further complexifies the transcriptome
**Tex file:** `Manuscript/sections/results/motifs_RNase.tex`
**Figure:** `Manuscript/figures/rnase.pdf`

**Claim:**
**Logic:**
**Analysis:**
- Scripts:
- Key params:
**Outputs:**
**Numbers to cite:**
**Figure panels:**
**Conclusion:**
**Caveats / hedges:**
**Notes for LLM:**

---

## R3 — Correlation between transcriptomics and proteomics in the reduced organism
**Tex file:** `Manuscript/sections/results/corr_RNA_ptn.tex`
**Figure:** _<which?>_

**Claim:**
**Logic:**
**Analysis:**
- Scripts: `Syn1_Corr_RNA_Proteins/Transcription_Translation.py`, `Translation_Residual_L2_elongation.py`
- Key params:
**Outputs:**
- `Syn1_Corr_RNA_Proteins/syn1_genes_transcriptomics_proteomics.csv`
- `Syn1_Corr_RNA_Proteins/residual_analysis/…`
**Numbers to cite:**
**Figure panels:**
**Conclusion:**
**Caveats / hedges:**
**Notes for LLM:**

---

## R4 — Novel transcription and translation activities
**Tex file:** `Manuscript/sections/results/novel.tex`
**Figure:** `Manuscript/figures/novel.pdf` (+ `novel-si.pdf`)

**Claim:**
**Logic:**
**Analysis:**
- Scripts: `Syn1_Novel_ORF/…`
- Key params:
**Outputs:**
**Numbers to cite:**
**Figure panels:**
**Conclusion:**
**Caveats / hedges:**
**Notes for LLM:**

---

## R5 — Genome reduction to the minimal cell, JCVI-syn3A
**Tex file:** `Manuscript/sections/results/genome_reduction.tex`
**Figure:** `Manuscript/figures/reduction.pdf`

This is a large, multi-panel result — break into sub-stories so each maps to a
single panel / paragraph.

### R5.1 — Reduction landscape (deletions/insertions/relocations)
**Claim:**
**Logic:**
**Analysis:** `Genome_Reduction/01_align.sh`, `02_analyze.py`, `03_visualize.py`
**Outputs:** `aln/raw/syn1_deleted_regions.bed`, `aln/analysis/genome_reduction_summary.xlsx`
**Numbers to cite:** _(e.g. 95 deletions; ~536 kb removed; 99.90% identity; 36 SNPs; 12 indels)_
**Figure panels:**
**Conclusion:**
**Caveats / hedges:**

### R5.2 — Deletion × operon overlay (truncation + gene-deletion patterns)
**Claim:**
**Logic:**
**Analysis:** `Genome_Reduction/04_deletion_overlaid_operon.py`
**Outputs:** `deletion_overlaid_operon/operon_deletion_classification.tsv`
**Numbers to cite:**
**Figure panels:**
**Conclusion:**

### R5.3 — Deletion-junction taxonomy
**Claim:**
**Logic:**
**Analysis:** `Genome_Reduction/05_deletion_junction.py`
**Outputs:** `deletion_junction/deletion_junctions.tsv`, `deletion_junction_summary.txt`
**Numbers to cite:** _(counts per `strand_relationship` × `junction_type`)_
**Figure panels:**
**Conclusion:**

### R5.4 — Co-expression validation (single-operon + cross-junction pairs)
**Claim:**
**Logic:**
**Analysis:** `06_single_operon_coexpression.py`, `07_operon_pair_coexpression.py`, `coexpression_common.py`
**Outputs:** `single_operon_coexpression/`, `operon_pair_coexpression/`
**Numbers to cite:**
**Figure panels:**
**Conclusion:**

### R5.5 — Per-gene impact classification
**Claim:**
**Logic:**
**Analysis:** `Genome_Reduction/08_delete_gene.py`
**Outputs:** `delete_gene/retained_gene_context.tsv` (`gene_impact_class` column)
**Numbers to cite:** _(class counts; promoter_lost vs others)_
**Figure panels:**
**Conclusion:**

### R5.6 — RNA-level remodeling (TPM)
**Claim:**
**Logic:**
**Analysis:** `Genome_Reduction/09_Compare_RNA_Protein.py`
**Outputs:**
- `Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv`
- `Compare_RNA_Protein/TPM_FC_vs_absChange.pdf`, `TPM_FC_by_impact_class.pdf`
- `Compare_RNA_Protein/mRNA_pool_composition_by_secondary.pdf`, `tertiary_share_change_dumbbell.pdf`
- `deleted_gene_occupancy.txt`
**Numbers to cite:** _(deleted-gene transcriptome share; promoter_lost median FC; r-protein share shift)_
**Figure panels:**
**Conclusion:**
**Caveats:** TPM = Illumina for quantification; ONT only for non-coding fallback / QC. Mean-normalized; deletion-corrected for cross-organism shares.

### R5.7 — Protein-level remodeling (iPM) and PTR
**Claim:**
**Logic:**
**Analysis:** `Genome_Reduction/10_Compare_Ptn.py`
**Outputs:**
- `Compare_RNA_Protein/iPM_*` plots, `iPM_pool_composition_by_secondary.pdf`, `iPM_tertiary_share_change_dumbbell.pdf`
- `Compare_RNA_Protein/PTR_TPMfc_vs_iPMfc.pdf`, `PTR_by_category_boxplot.pdf`
- `Compare_RNA_Protein/macromolecule_complex_abundance.tsv`
- `protein_{upgrade,downgrade}/`
**Numbers to cite:** _(RNAP ~21–35% down; degradosome ~36–68% up; n PTR outliers per Secondary category)_
**Figure panels:**
**Conclusion:**
**Caveats:** PTR is a steady-state proxy, NOT Ribo-seq TE. r-proteins excluded from PTR (digestion bias).

---

# METHODS

## M1 — Illumina short-read sequencing of syn1.0 transcriptome  
**Tex file:** `Manuscript/sections/methods/illumina_syn1.tex`  
**Analysis:** `Syn1_Transcriptomics/Illumina/Illumina_Processing/`  
**Key params to mention:** FastQC + MultiQC; bowtie2 paired-end; dUTP / fr-firststrand (R2 = transcript strand); per-strand bedGraph.  
**Inputs:** SRA accessions in `Syn1_Transcriptomics/Illumina/Illumina_Raw/00_retrive_fastq.sh`  
**Notes for LLM:** Subsubsection **Illumina MiSeq Read Processing, Mapping to the Genome** and **Quantification of TPM from Sequencing Depth** need written. 

---

## M2 — PacBio long-read sequencing of syn1.0 transcriptome  
**Tex file:** `Manuscript/sections/methods/pacbio_syn1.tex`  
**Analysis:** `Syn1_Transcriptomics/PacBio/PacBio_Processing/`  
**Key params:** FLNC recovery (reorientation, primer removal, polyA trimming); `minimap2 -ax map-hifi`; HQ filtering.  
**Outputs:** `syn1.PacBio.FLNC.sorted.HQ.bam`, `depth_bedgraph/…`  
**Inputs:** SRA accessions in `Syn1_Transcriptomics/PacBio/PacBio_Raw/00_retrieve_fastq.sh`  
**Notes for LLM:** This section mostly finished and polished.  

---

## M3 — Operon identification from PacBio long-read transcriptomics in syn1.0
**Tex file:** `Manuscript/sections/methods/operon_analysis.tex`  
**Analysis:** `Syn1_Transcriptomics/Isoforms_PacBio/Cluster_Isoform.py`, `Syn1_Operon/…`  
**Key params:** clustering thresholds, min reads, TSS/TTS calling rule.  
**Outputs:** `isoform_clusters_annotated.tsv`, `operons.candidate_blocks.tsv`  
**Notes for LLM:** Subsubsection **Locate TSS, TTS, and RNA Cleavage Sites from PacBio RNASeq** need changed and written.  

---

## M4 — RNA processing and ribonucleases  
**Tex file:** `Manuscript/sections/methods/RNA_processing.tex`  
**Analysis:** _<scripts>_  
**Inputs:** `Genomes_Input/Motif_Identifications.xlsx`  
**Notes for LLM:** The whole section needs my further analysis.  

---

## M5 — Proteomics of syn1 and syn3A  
**Tex file:** `Manuscript/sections/methods/proteomics_syn1_syn3A.tex`  
**Analysis / source:** `Syn1_Syn3A_Proteomics/` — `syn1_proteomics_localization_2026.csv`, `syn3a_proteomics_summary_2026.csv`, `syn3A_proteome_annotated.xlsx`.  
**Key params:** absolute copy numbers vs iPM; 2019 vs 2026 measurements; tertiary function annotation built by `report_annotation_stats_syn3A.py`.  
**Notes for LLM:** distinguish what we measured vs reused from prior work; cite original datasets; subsubsection **Relative Protein Quantification** and **Localization of Proteome** need written.

---

## M6 — Correlation between transcriptome and proteome
**Tex file:** `Manuscript/sections/methods/corr_transcriptome_proteome.tex`  
**Analysis:** `Syn1_Corr_RNA_Proteins/Transcription_Translation.py`, `Translation_Residual_L2_elongation.py`  
**Key params:** TPM source (PacBio vs Illumina); relative iPM; residual model covariates (CAI, tAI, TIR, etc.).  
**Outputs:** `syn1_genes_transcriptomics_proteomics.csv`, `residual_analysis/`  
**Notes for LLM:** Entire session needs written.  

---

## M7 — Novel open-reading frames from long-read transcripts
**Tex file:** `Manuscript/sections/methods/novel_orf.tex`  
**Analysis:** `Syn1_Novel_ORF/…`   
**Key params:** Mycoplasma genetic code (UGA = Trp); anti-SD `ACCUCCUUU`; ORF length/SD-distance thresholds.  
**Notes for LLM:** Entire session needs written.  

---

## M8 — Oxford Nanopore (ONT) and Illumina sequencing of syn3A transcriptome
**Tex file:** `Manuscript/sections/methods/ont_illumina_syn3A.tex`  
**Analysis:** `Syn3A_Transcriptomics/ONT/ONT_Processing/` (ONT) + `Syn3A_Transcriptomics/Illumina/Illumina_Processing/` (Illumina)  
**Key params:** ONT direct-RNA, `minimap2 -ax map-ont` (NOT splice, bacteria are intron-less), per-strand depth; Illumina syn3A paired-end bowtie2 (dUTP / fr-firststrand), per-strand bedGraph.  
**Inputs:** ONT raw `Syn3A_Transcriptomics/ONT/ONT_Raw/`; Illumina syn3A SRA accessions (SRR19432056/57 mate pair) via `Syn3A_Transcriptomics/Illumina/Illumina_Raw/00_retrive_fastq.sh`.  
**Notes for LLM:** file now covers BOTH ONT and Illumina-syn3A mapping (two subsubsections); rRNA operons at ~55,460 and ~343,267 bp drive the multi-mapping fraction.  

---

## M9 — Genome reduction
**Tex file:** `Manuscript/sections/methods/genome_reduction.tex`  
**Analysis pipeline (in order):** `Genome_Reduction/01_align.sh` → `02_analyze.py` → `03_visualize.py` → `04_deletion_overlaid_operon.py` → `05_deletion_junction.py` → `06_single_operon_coexpression.py` → `07_operon_pair_coexpression.py` → `08_delete_gene.py` → `09_Compare_RNA_Protein.py` → `10_Compare_Ptn.py`.

**Key params / definitions to spell out in Methods:**
- Coordinate convention: 0-based half-open; circular wrap.
- Locus-tag correspondence: `MMSYN1_NNNN ↔ JCVISYN3A_NNNN` (numeric suffix preserved).
- Junction taxonomy: `strand_relationship` {tandem, convergent, divergent, intra_operon}; `junction_type` (tandem only) {fusion, decapitation, readthrough_extension, clean_excision}.
- `gene_impact_class` precedence: promoter_lost > promoter_disconnected > new_promoter_fusion > readthrough_exposed > promoter_proximity_changed > context_only > unaffected.
- Mean-normalization (`rel*`) and **retained-pool deletion correction** for cross-organism shares.
- TPM platform policy: Illumina for both organisms (coding); ONT for syn3A non-coding only.
- PTR definition: relIPM/relTPM; `PTR_fold_change = iPM_FC / TPM_FC`; explicitly a steady-state proxy, not Ribo-seq TE.
- Co-expression test (06/07): ONT spanning/bridging reads + Illumina gap depth; thresholds in `coexpression_common.py`.

**Notes for LLM:** This is the longest Methods subsection — draft it in the same order as the pipeline.

---