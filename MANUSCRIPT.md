# MANUSCRIPT.md — drafting guide for Methods & Results

This file is the **single source of truth** the LLM uses to draft Methods and
Results paragraphs. Each subsection here mirrors a `.tex` file under
`Manuscript/sections/{results,methods}/`. Fill in the bullets; the LLM will
expand them into prose, quoting numbers verbatim from the files you list.

TO DO Reminder:

- Fill out the logics in Results and generate the Results.
- Regenerate and format the figures.
- Reorganize Operon_Visual Jupyter Notebook before publication.

---

## 0. Style guide (read before drafting any subsection)

**Target journal:** Nature Microbiology  
**Length budget:** Results ≈ 3000; Methods ≈ unlimited (unlimited); Intro + Discussion ≈ 500  
**Tense / voice:**  
- Methods: past tense, passive acceptable ("Reads were aligned with …").
- Results: past tense for what we did, present for what the data show ("syn3A retains 99.90% identity …").  

**Person:** "we" not allowed in Results and Methods  
**Abbreviations:** define on first use in each of Abstract / Main / Methods. Do NOT define abbr in figure legends.
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

## Overview of RESULTS

Six sections and six multi-panel figures; each section will be of 500 words.

Chain of logics for each section; use one or multiple paragraphs for each logic.

---

## R1 — Gene Co-transcriptions as Operons in Syn1 from PacBio Long-read RNA Sequencing

**Tex file:** `Manuscript/sections/results/operons.tex`  

### One-sentence Summary
**480 operons were identified using PacBio long-read RNAseq, with transcription signatures located.**

### Figure
**Figure:** `Manuscript/figures/operon.pdf`

- Panel a: Two-gene operon co-transcription and the following RNA processing.
- Panel b: Number of sense genes per operon in syn1.
- Panel c: Lengths per operon in syn1.
- Panel d: Transcription promoter and terminator.
- Panel e: Macromolecular complex ... operonal structure in syn1.

### Chain of Logics

#### L1.1: Co-transcription and further RNA processing complexify the transcriptome even of the reduced bacterium.

- **Logic:** Co-transcription with multiple transcription start and terminator sites generates multiple transcription units; RNA processing from endo- or exo-ribonucleases can digest the transcripts to even more **RNA isoforms**.
- **Analysis:** None
- **Outputs:** None
- **Numbers to cite:**  None
- **Figure panels:** a  
- **Conclusion:** None
- **Caveats:** None
- **Notes for LLM:** This is a descriptive part.

#### L1.2: 2.6 M full-length PacBio RNA seq clustered into 267k isoform clusters.

- **Logic:** PacBio raw cDNA reads were processed and quality controlled to output 2.6 M RNA reads. Clustering was applied to suppress the noise.
- **Analysis:** 
  - Processing: `Syn1_Transcriptomics/PacBio/PacBio_Processing/PacBio_Processing.py`
  - Clustering: `Syn1_Transcriptomics/PacBio/Isoforms_PacBio/Cluster_Isoform.py`
- **Outputs:** 
  - BAM file: `Syn1_Transcriptomics/PacBio/PacBio_Processing/syn1.PacBio.FLNC.sorted.HQ.bam`
  - Isoform clusters: `Syn1_Transcriptomics/PacBio/Isoforms_PacBio/isoform_clusters_annotated.xlsx`
- **Numbers to cite:**  cluster eps = 10 bp, low MAD
- **Figure panels:** None
- **Conclusion:** 267k isoform clusters with sharp ends serve as solid foundation for operon calling.
- **Caveats:** None
- **Notes for LLM:** More details presented in Methods M2.

#### L1.3: 480 operons were mapped by full-length PacBio RNA seq.

- **Logic:** Unique longest isoform clusters as evidence of gene co-transcription were constructed by containment to cover 316 initial operons. Overlap between operons was solved. Uncovered genes were rescued. Finally, 480 operons for 911 genes in syn1. The statistics on the size and length of operons were reported.
- **Analysis:** 
  - Operon segmentation: `Syn1_Operon/Operon_Segmentation.ipynb`
  - Operon annotation: `Syn1_Operon/Operon_Annotation.ipynb`
- **Outputs:** 
  - Operons: `Syn1_Operon/operons.candidate_blocks.tsv`
- **Numbers to cite:**  MIN_READS threshold = 50. mean and median of length, sense gene count, anti-sense gene count; largest operon of ribosomal proteins
- **Figure panels:** b,c
- **Conclusion:** None
- **Caveats:** None
- **Notes for LLM:** More details are in Methods M3

#### L1.4: Transcription signatures located for operons.

- **Logic:** Exact transcription start site (TSS) and termination sites (TTS) were located for canonical operons whose boundaries were intergenic.
- **Analysis:** 
  - Operon segmentation: `Syn1_Operon/Operon_Annotation.ipynb`
- **Outputs:** 
  - Transcription signatures: `Syn1_Operon/annotation/canonical/`
- **Numbers to cite:** None
- **Figure panels:** d
- **Conclusion:** Signatures are both consistent with previous knowledge: -10 box of TSS has TANAAT, -35 just AT rich; TTS as intrinsic terminators as hairpin + polyU
- **Caveats:** The TSS, TTS sites are only for canonical operons; we might need to refine for all cases.
- **Notes for LLM:** Corresponding Method needs to be finished.

#### L1.5: One instance of polycistronic operons

- **Logic:** The choice not decided yet: could be rPtn operons, or other complexes
- **Analysis:** None
- **Outputs:** 
- **Numbers to cite:** None
- **Figure panels:** e
- **Conclusion:** 
- **Caveats:** 
- **Notes for LLM:** Use for instance to catch eyes of readers

---

## R2 — Pervasive and biased RNA processing further complexifies the transcriptome
**Tex file:** `Manuscript/sections/results/RNase.tex`

### One-sentence Summary
**Pervasive and biased RNA processing caused truncated RNA isoforms with more 3' erosion.**

### Figure
**Figure:** `Manuscript/figures/rnase.pdf`

- Panel a: RNA isoform distribution for gene 0154/lap with more 5' erosion.
- Panel b: RNA isoform distribution for gene 0178 with more 3' erosion.
- Panel c: RNA isoform truncation categories.
- Panel d: Biased RNA Processing schematics.
- Panel e: RNA isoform distributions for ATP synthase operon.

### Chain of Logics

#### L2.1: Distinct RNA isoforms distributions found for operons.

- **Logic:** Truncated isoforms compared to the full transcription units exist for operons because of the RNA processing. Distinct patterns of truncations can be found, using genes 0154 and 0178 as examples.
- **Analysis:** None
- **Outputs:** None
- **Numbers to cite:**  None
- **Figure panels:** a,b
- **Conclusion:** None
- **Caveats:** None
- **Notes for LLM:** None

#### L2.2: Transcriptome-wide, significantly more 3' erosion found.

- **Logic:** Overlaid RNA isoforms to gene ORFs to find significantly more 3' erosions.
- **Analysis:** `Syn1_Operon/RNA_Processing.py`
- **Outputs:** 
  - `Syn1_Operon/RNA_Processing.txt`
  - plots in `Syn1_Operon/RNase/`
- **Numbers to cite:**  None
- **Figure panels:** c
- **Conclusion:** None
- **Caveats:** None
- **Notes for LLM:** Method **RNA Processing Analysis** needs to be finished.

#### L2.3: 3' Exo ribonucleases can contribute to the biased RNA digestion. (Putative)

- **Logic:** RNA digestion can be initiated not only by endo-ribonuclease, but also 3'~5' exo-ribonucleases.
- **Analysis:** None
- **Outputs:** None
- **Numbers to cite:**  None
- **Figure panels:** d
- **Conclusion:** None
- **Caveats:** None
- **Notes for LLM:** This analysis needs to be implemented by checking the order of RNA secondary structure at the 3' end.

#### L2.4: ATP synthase operon is co-expressed in one-go but cut at $\alpha$ subunit.

- **Logic:** Macromolecular complexes' gene co-expression can be altered by RNA processing. ATP synthase's RNA isoform distribution has a clear pattern of isolation at $\alpha$ subunits, which was identified as endo RNase III cleavage site.
- **Analysis:** 
- **Outputs:** 
  - plot: `Syn1_Operon/ATP_Synthase_wdepth.pdf`
- **Numbers to cite:**  None
- **Figure panels:** e
- **Conclusion:** The RNase complexifies the subunit synthesis of complexes.
- **Caveats:** None
- **Notes for LLM:** Needs to check for other complexes; visualize the RNA secondary structure at $\alpha$

---

## R3 — High Correlation between transcriptomics and proteomics in the reduced organism.

**Tex file:** `Manuscript/sections/results/corr_RNA_ptn.tex`

### One-sentence Summary
**High correlation found between transcriptome and proteome.**

### Figure
**Figure:** `Manuscript/figures/correlation.pdf`

- Panel a: Log10 distribution of proteins copy numbers as cytoplasm, membrane, lipoprotein and secreted proteins.
- Panel b: Correlation between Illumina TPM and Mass-spec iPM.
- Panel c: PacBio TPM and Illumina TPM correlation of syn1.
- Panel d: Predicted Translation Initiation Rate (TIR) poorly correlates with residual.
- Panel e: Codon Adaptation Index (CAI) decently correlates with residual.
- Panel f: Correlation coefficient R improvements as only cytosolic proteins or entire proteome and with CAI or not.
- Panel g: Intrinsic protein half-lives distribution transferred from Mpn.
- Panel h: Protein half-lives poorly correlate with the residual.

### Chain of Logics

#### L3.1: Distribution of protein copy numbers in syn1

- **Logic:** Starting from the iBAQ after the mapping using SpectroNaut software, iPM was calculated, which was then converted to absolute protein copy numbers with protein dry mass per cell; Cytoplasmic protein has a median copy number of 47, while only 10 for membrane proteins.
- **Analysis:** 
  - Copy number quantification: `Syn1_Syn3A_Proteomics/Protein_Quantification_Localization.ipynb`
- **Outputs:** 
  - syn1 proteomics: `Syn1_Syn3A_Proteomics/syn1_proteomics_localization_2026.csv`
- **Numbers to cite:**  Median copy number of protein at locations
- **Figure panels:** a
- **Conclusion:** Poor coverage of membrane proteins because of the protease digestion; 
- **Caveats:** None
- **Notes for LLM:** Method **Relative Protein Quantification** and **Localization of Proteome** needs to be finished; **Absolute Intracellular Protein Quantification** needs to be polished.

#### L3.2: Using Illumina TPM as standard of transcriptome quantification

- **Logic:** PacBio and Illumina TPMs of syn1 were correlated to get r of 0.62; no significant TPM and length bias was found.
- **Analysis:** `Syn1_Transcriptomics/Gene_TPM/Gene_Transcriptomics.py`
- **Outputs:** plots in `Syn1_Transcriptomics/Gene_TPM/`
- **Numbers to cite:**  r of 0.62
- **Figure panels:** c
- **Conclusion:** As convention, Illumina TPMs were used to do correlation.
- **Caveats:** None
- **Notes for LLM:** Method **Illumina MiSeq Read Processing, Mapping to the Genome** under M1 and two **Quantification of TPM from Sequencing Depth** under M1 and M2 need finished.

#### L3.3: Decent correlation found between transcriptome and proteome for syn1

- **Logic:** Pearson r of 0.7 found between two omics for cytosolic proteins; lower r for all since poor coverage of membrane proteins.
- **Analysis:** `Syn1_Corr_RNA_Proteins/Transcription_Translation.py`
- **Outputs:** 
  - Same name Txt file
- **Numbers to cite:**  Pearson r's
- **Figure panels:** b
- **Conclusion:** Decent correlation.
- **Caveats:** None
- **Notes for LLM:** None

#### L3.4: Predicted TIR had low correlation with residuals between two omics

- **Logic:** Translation initiation rates predicted for all gene ORFs by OSTIR were evaluated to validate how much of the residuals can be explained by translation initiation.
- **Analysis:** `Syn1_Corr_RNA_Proteins/Translation_Residual_L1_initiation.py`
- **Outputs:** 
  - Same name Txt file
  - plots under `Syn1_Corr_RNA_Proteins/residual_analysis`
- **Numbers to cite:**  Pearson r's
- **Figure panels:** d
- **Conclusion:** No improvement on correlations, could be because the prediction itself was inaccurate.
- **Caveats:** None
- **Notes for LLM:** Method **Translation Initiation Rate Prediction** needs to be finished.

#### L3.5: Translation elongation factor improved the correlations

- **Logic:** CAI as metric for translation elongation efficiency significantly improved the correlation.
- **Analysis:** `Syn1_Corr_RNA_Proteins/Translation_Residual_L2_elongation.py`
- **Outputs:** 
  - Same name Txt file
  - plots under `Syn1_Corr_RNA_Proteins/residual_analysis`
- **Numbers to cite:**  Pearson r's
- **Figure panels:** e,f
- **Conclusion:** Translation elongation affected the protein biosynthesis.
- **Caveats:** None
- **Notes for LLM:** Method **Codon Adaptation Index (CAI)** needs to be finished.

#### L3.6: Protein degradation had low correlation with residuals.

- **Logic:** Intrinsic protein degradation as a result of protease activities in syn1 was transferred from Mpn by finding the reciprocal homologs; Lowest half-life was 4.7 hours, way longer than the half-life of 1 hour for syn1.
- **Analysis:** `Syn1_Corr_RNA_Proteins/Translation_Residual_L3_degradation.py`
- **Outputs:** 
  - Same name Txt file
  - plots under `Syn1_Corr_RNA_Proteins/residual_analysis`
- **Numbers to cite:**  half-life distributions; Pearson r's
- **Figure panels:** g,h
- **Conclusion:** Intrinsic protein degradation is way slower than protein synthesis and doubling in syn1, thus having a minor effect on proteome abundances.
- **Caveats:** Only a subset of proteins found reciprocal homologs; the intrinsic half-lives were corrected by protease (Lon or FtsH) abundances.
- **Notes for LLM:** Method **Protein Degradation Rate Mapping from Mpn** needs to be finished.

---

## R4 — Novel transcription and translation activities of the synthetic bacterium, Syn1

**Tex file:** `Manuscript/sections/results/novel.tex`

### One-sentence Summary
**Long-read RNA seq reveals anti-sense and intergenic transcription.**

### Figure
**Figure:** `Manuscript/figures/novel.pdf`

- Panel a: Distribution of anti-sense percentage of all isoforms
- Panel b: Schematics of read-through, embedded, and spurious promoter to explain anti-sense coverage in operons.
- Panel c: Distribution of 5' and 3' untranslated regions (UTR) in operons
- Panel d: Truly intergenic transcription between 0154 and 0155.
- Panel e: New ORF encoding peptide of length 118 aas found after gene 0592.

### Chain of Logics

#### L4.1: Minor percentage of isoforms have anti-sense coverages that can be categorized into three cases.

- **Logic:** 2% of the isoforms have antisense transcription coverage. Two thirds came from the spurious promoters in the AT-rich genome; the transcription read-through can cause anti-sense transcription at the end or embedded inside the operons.
- **Analysis:** 
  - `Syn1_Novel_ORF/Abnormal_Transcripts.py`
- **Outputs:** 
  - Same name Txt file
  - plots under `Syn1_Novel_ORF`
- **Numbers to cite:**  2%; cases of anti-sense coverage
- **Figure panels:** a,b
- **Conclusion:** Full-length RNA isoforms reveal new cases of anti-sense transcription as read-throughs.
- **Caveats:** None
- **Notes for LLM:** Method **Novel Transcription and Translation** needs to be finished and polished.

#### L4.2: Distribution of 5' and 3' UTR lengths.

- **Logic:** UTR lengths were evaluated for all canonical operons with median of tens of nucleotides.
- **Analysis:** 
  - `Syn1_Operon/Operon_Annotation.ipynb`
- **Outputs:** 
  - Output inside Notebook
- **Numbers to cite:**  Median, and maximum
- **Figure panels:** c
- **Conclusion:** Median values were consistent with previous papers; high outliers are due to the anti-sense transcription or truncations.
- **Caveats:** None
- **Notes for LLM:** None

#### L4.3: One truly isolated intergenic transcription.

- **Logic:** One truly isolated intergenic transcription was found between genes 0154 and 0155.
- **Analysis:** Operon coverage analysis in `Operon_Segmentation.ipynb`
- **Outputs:** 
  - `Syn1_Operon/operons.candidate_blocks.tsv`
- **Numbers to cite:**  None
- **Figure panels:** d
- **Conclusion:** None
- **Caveats:** None
- **Notes for LLM:** None

#### L4.4: Novel peptide identified by enumerating all possible ORFs in isoforms having high abnormal fraction.

- **Logic:** Driven by curiosity, the possible translation of the abnormal fraction was checked by enumerating all ORFs using OSTIR.
- **Analysis:** 
  - `Syn1_Novel_ORF/Novel_translation.ipynb`
- **Outputs:** 
  - `Syn1_Operon/operons.candidate_blocks.tsv`
- **Numbers to cite:**  None
- **Figure panels:** e
- **Conclusion:** Two predicted ORFs were identified in Mass-spec proteome, and both were located near less annotated genes. Also, these two regions were deleted in syn3A.
- **Caveats:** None
- **Notes for LLM:** See Method **Novel Open-Reading Frames from Full-length RNA Isoforms** for details; the Analysis Jupyter Notebook can be cleaned up.

#### L4.5: Traits of transcription of synthetic elements.

- **Logic:** Yeast vector and watermarks were inserted into syn1's synthetic genome. Strikingly, one of the Yeast genes, 0918, was heavily transcribed in an anti-sense way but not translated at all. The watermarks were minimally transcribed as noise.
- **Analysis:** None
- **Outputs:** None
- **Numbers to cite:**  None
- **Figure panels:** None
- **Conclusion:** None
- **Notes for LLM:** A separate Python script can be created to finalize this analysis.

---

## R5 — Operonal structure changes to the minimal cell, JCVI-syn3A
**Tex file:** `Manuscript/sections/results/reduction_operons.tex`

### One-sentence Summary
**Halving the genome was a gene-order-preserving deletion campaign that excised whole operons, decapitated some retained operons by deleting their promoters, and fused only a small number of new cross-junction transcription units.**

### Figure
**Figure:** `Manuscript/figures/genome_reduction.pdf`

- Panel a: Schematics of genome reduction from syn1 to syn3A.
- Panel b: Fusion of new operons.
- Panel c: Box plot of gene expressions to highlight decapitated ones having lower values.
- Panel d: the HupA operon, whose true promoter, located inside gene 0349, was deleted.
- Panel e: Gene essentiality evaluation for those trace-expressed genes that are still essential.

### Chain of Logics

#### L5.1: Reduction from syn1 to syn3A was a gene-order-preserving deletion campaign.

- **Logic:** Aligning syn3A back onto syn1 (nucmer/dnadiff) reframes minimization as a set of discrete deletions; counting the cuts, measuring retained-sequence identity, and testing for inversions/translocations/relocations distinguishes whether reduction rewired the genome or simply removed pieces of it.
- **Analysis:** `Genome_Reduction/01_align.sh` -> `02_analyze.py` -> `03_visualize.py`
- **Outputs:**
  - `Genome_Reduction/aln/raw/syn1_deleted_regions.bed`
  - `Genome_Reduction/aln/analysis/genome_reduction_summary.{xlsx,txt}`
- **Numbers to cite:** 1,078,809 -> 543,379 bp; 95 deletions (>= 50 bp); 536,543 bp (~536 kb) removed; mean 5,647 bp; largest 71,578 bp; 99.90% identity; 36 SNPs; 12 indels; 0 inversions; 0 translocations; 1 relocation (lap / MMSYN1_0154, ~110 kb downstream); 6 insertions / 1,324 bp; 1 novel CDS JCVISYN3A_0931 (met14p).
- **Figure panels:** a
- **Conclusion:** Half the genome was excised in discrete cuts while retained DNA stayed essentially identical and gene order was preserved; expression differences at retained genes are therefore not attributable to sequence divergence.
- **Caveats:** dnadiff reports 118 raw reference-side "insertion" events; the filtered >= 50 bp BED set (95) is the authoritative deletion list.
- **Notes for LLM:** Method M9. Use this to frame the whole section as a structural/regulatory (not sequence-level) story.

#### L5.2: Deletions overlaid on syn1's 480 operons show whole-operon excision dominating over partial truncation.

- **Logic:** Intersecting the 95 deletions with the 480 syn1 operons at single-bp resolution classifies how each operon was hit, separating operons removed wholesale from those left partially truncated; the truncations are what create the junction effects in L5.3 and L5.4.
- **Analysis:** `Genome_Reduction/04_deletion_overlaid_operon.py`
- **Outputs:** `Genome_Reduction/deletion_overlaid_operon/operon_deletion_classification.tsv`
- **Numbers to cite:** span-level overlap_class (n=480): intact 180, fully_deleted 179, 3'_truncation_gene 50, 5'_truncation_gene 32, intra_truncated 18; gene-level gene_deletion_pattern: all_deleted 238, intact 190, leading_deleted 23, lagging_deleted 20, intra_deleted 9; 414 syn1 genes overlapped by a deletion.
- **Figure panels:** a
- **Conclusion:** Reduction preferentially removed entire operons; the minority of partial truncations (5' vs 3') sets up the junction taxonomy.
- **Caveats:** the two axes (span-level truncation vs gene-level deletion) differ by design; 180 vs 190 "intact" reflects operons whose genes are all kept but whose UTR/flank was nicked.
- **Notes for LLM:** None.

#### L5.3: Same-strand deletion junctions can fuse new transcription units, but true fusion is rare.

- **Logic:** Each deletion is recast as a junction between the nearest retained operons on either side; relative orientation (tandem/convergent/divergent) and facing-regulator loss decide whether a new co-transcribed unit can form, and ONT spanning/bridging reads test whether the new cross-junction gene pair is actually co-transcribed.
- **Analysis:**
  - junction taxonomy: `Genome_Reduction/05_deletion_junction.py`
  - read validation: `06_single_operon_coexpression.py`, `07_operon_pair_coexpression.py`, `coexpression_common.py`
- **Outputs:**
  - `Genome_Reduction/deletion_junction/deletion_junctions.tsv`, `deletion_junction_summary.txt`
  - `Genome_Reduction/operon_pair_coexpression/`, `single_operon_coexpression/`
- **Numbers to cite:** 95 junctions: tandem 55, convergent 19, divergent 15, intra_operon 6; tandem junction_type: fusion 3, decapitation 10, readthrough_extension 8, clean_excision 34; cross-junction co-transcription (loose): fusion 67% (n=3) vs clean_excision 15% (negative control); pristine single-operon baseline preserved_loose 65% (51 testable); fusion exemplar DEL_014 OP_00045 -> OP_00053 (MMSYN1_0082 -> MMSYN1_0094), n_span=2, n_bridge=37.
- **Figure panels:** b
- **Conclusion:** Operon fusion is real but rare (3 events); the dominant junction outcome is clean excision of whole operon(s) between intact neighbors.
- **Caveats:** ONT depth is low, so most positive calls are loose-bridge rather than strict-spanning; convergent/divergent junctions are opposite-strand and not expected to co-transcribe.
- **Notes for LLM:** This logic, not "all new operons fused," supports the reworded one-sentence summary.

#### L5.4: Decapitated operons that lost their own promoter are the one class that robustly drops in expression; HupA is the showcase.

- **Logic:** Classifying every retained gene by promoter-source change isolates operons whose own promoter was deleted (promoter_lost / decapitation); their syn3A TPM is compared against the other impact classes to test whether promoter loss, not sequence change, predicts lower expression.
- **Analysis:**
  - per-gene impact: `Genome_Reduction/08_delete_gene.py`
  - expression: `09_Compare_RNA_Protein.py`
- **Outputs:**
  - `Genome_Reduction/delete_gene/retained_gene_context.tsv` (`gene_impact_class` column)
  - `Genome_Reduction/Compare_RNA_Protein/TPM_FC_by_impact_class.pdf`
- **Numbers to cite:** gene_impact_class: promoter_lost 35, promoter_disconnected 8, new_promoter_fusion 3, readthrough_exposed 14, promoter_proximity_changed 17, context_only 43, unaffected 377; HupA (MMSYN1_0350) relTPM 6.68 -> 0.13 (FC 0.020), relIPM 6.48 -> 0.092 (FC 0.014); other decapitated drops rpmE/L31 FC 0.010, rpsU/S21 FC 0.044.
- **Figure panels:** c, d
- **Conclusion:** Promoter-source loss drives the largest expression decreases; promoter_lost is the only impact class robustly down in TPM.
- **Caveats:** the class is assigned at operon level; 8 junctions lose only UTR (genes intact); the 05-vs-04 consistency check flags 2 flank operons as all_deleted.
- **Notes for LLM:** HupA's true promoter sits inside deleted gene MMSYN1_0349 (panel d).

#### L5.5: A few trace-expressed retained genes remain essential.

- **Logic:** Crossing gene essentiality against syn3A expression surfaces genes that are essential yet barely transcribed, i.e. retained through minimization despite minimal expression.
- **Analysis:** TBD (not produced by the 01-10 pipeline).
- **Outputs:** TBD
- **Numbers to cite:** TBD
- **Figure panels:** e
- **Conclusion:** TBD
- **Caveats:** essentiality calls are inherited from the syn3A design literature, not measured here.
- **Notes for LLM:** ANALYSIS NOT YET DONE. Keep the logic; an essentiality x trace-expression script must be written and the essentiality source supplied before this can be drafted.

---

## R6 — Transcriptome and Proteome Changes to minimal cell, Syn3A
**Tex file:** `Manuscript/sections/results/reduction_omics.tex`

### One-sentence Summary
**More transcription on ribosomal protein operons suppresses the expression of enzymatic proteins in central metabolism.**

### Figure
**Figure:** `Manuscript/figures/reduction.pdf`

- Panel a: mRNA pool compositions in syn1 and syn3A as secondary protein functions.
- Panel b: Significant mRNA pool share changes from syn1 to syn3A as tertiary functions.
- Panel c: Transcription and translation changes of RNAP, degradosome and enzymes in central metabolism.
- Panel d: Flux comparison of ATP/GTP generation between syn1 and syn3A.

### Chain of Logics

#### L6.1: The ~418 deleted loci carried about a fifth of syn1's coding expression, freeing pool capacity.

- **Logic:** Quantifying the share of syn1's transcriptome and proteome contributed by loci absent from syn3A measures how much expression budget minimization freed for reallocation, and which RNA classes were lost.
- **Analysis:** `Genome_Reduction/09_Compare_RNA_Protein.py`
- **Outputs:**
  - `Compare_RNA_Protein/deleted_gene_occupancy.txt`
  - `Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv`
- **Numbers to cite:** 418 deleted loci (911 -> 496); by RNA type mRNA 382, pseudo 33, ncRNA 2, tRNA 1; deleted share = 21.78% of the syn1 mRNA pool (14.22% of all-RNA TPM), 22.25% of the iPM proteome; top deleted by TPM lacZ, pdhA/pdhB, ald; unclear-function proteins occupy only ~3%.
- **Figure panels:** a
- **Conclusion:** Minimization removed ~1/5 of the coding transcriptome and proteome, concentrated in dispensable metabolism, leaving pool capacity that syn3A redistributes.
- **Caveats:** shares are raw syn1 TPM/iPM; cross-organism comparisons in L6.2-L6.4 are mean-normalized and deletion-corrected to the retained-gene pool.
- **Notes for LLM:** "Non-essential" in the heading is shorthand for syn3A-deleted loci.

#### L6.2: The retained mRNA pool reallocates toward the translation machinery.

- **Logic:** After renormalizing to retained genes, per-function mRNA-pool shares are compared between organisms to test whether the freed capacity is taken up by ribosome/translation operons rather than spread evenly.
- **Analysis:** `Genome_Reduction/09_Compare_RNA_Protein.py`
- **Outputs:**
  - `Compare_RNA_Protein/TPM_change_by_{secondary,tertiary}.tsv`
  - `Compare_RNA_Protein/mRNA_pool_composition_by_secondary.pdf`, `tertiary_share_change_dumbbell.pdf`
- **Numbers to cite:** largest r-protein TPM absolute gains rpsK +6.04, rplO +5.04, infA +4.97, rplX +4.83, rplN +4.33; rpoA +2.69; per-category median FC + Mann-Whitney p from TPM_change_by_*.tsv.
- **Figure panels:** a, b
- **Conclusion:** The transcriptome shifts toward translation/ribosome biogenesis at the expense of central metabolism.
- **Caveats:** the shift is a pool-level reallocation, not uniform; several decapitated r-proteins (rpmF/L32, rpmE/L31, rpsU/S21) crash (links to L5.4).
- **Notes for LLM:** None.

#### L6.3: RNA polymerase is downregulated while the degradosome is upregulated in syn3A.

- **Logic:** Estimating each machine's assembled abundance from its limiting (lowest-stoichiometry) subunit compares transcription capacity against RNA-turnover capacity; opposite movement is coherent with syn3A's longer cell cycle.
- **Analysis:** `Genome_Reduction/10_Compare_Ptn.py`
- **Outputs:**
  - `Compare_RNA_Protein/macromolecule_complex_abundance.tsv`
  - `Compare_RNA_Protein/PTR_TPMfc_vs_iPMfc.pdf`, `PTR_by_category_boxplot.pdf`
- **Numbers to cite:** RNAP MIN(rpoA/2, rpoC, rpoB) TPM FC 0.65 (~35% down), iPM FC 0.79 (~21% down); degradosome MIN(rny, rnjA, yhaM+rnr) TPM FC 1.68 (~68% up), iPM FC 1.36 (~36% up); syn3A cell cycle 105 vs 60 min.
- **Figure panels:** c
- **Conclusion:** Lower transcription capacity plus higher RNA turnover is coherent with slower growth.
- **Caveats:** limiting-subunit estimate; PTR is a steady-state proxy, not Ribo-seq TE; r-proteins excluded from PTR (digestion bias).
- **Notes for LLM:** None.

#### L6.4: Central-metabolism enzymes are coordinately downgraded, predicting lower ATP/GTP generation.

- **Logic:** Glycolytic and energy-generating enzymes are tracked at both RNA and protein level; coordinated downgrade across the pathway predicts reduced ATP/GTP flux, to be confirmed against a metabolic-flux comparison.
- **Analysis:**
  - RNA/protein evidence: `09_Compare_RNA_Protein.py` / `10_Compare_Ptn.py`
  - flux comparison: TBD (not in the 01-10 pipeline)
- **Outputs:**
  - `Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv` (enzyme FCs)
  - flux output: TBD
- **Numbers to cite (omics, TPM_FC / iPM_FC, both down):** gapDH 0.44/0.63, eno 0.37/0.38, pgk 0.35/0.44, pdhC 0.45/0.34, ackA 0.19/0.44, pta 0.27/0.26, ldh 0.66/0.70, pyk -/0.34; flux numbers TBD.
- **Figure panels:** d
- **Conclusion:** Central-carbon and acetate-pathway enzymes drop in concert, predicting suppressed ATP/GTP output (flux quantification pending).
- **Caveats:** the flux claim is currently inferred from enzyme abundance only.
- **Notes for LLM:** FLUX ANALYSIS (panel d) NOT YET DONE. Keep the logic; the ATP/GTP flux comparison needs a metabolic model and must be finished before drafting the flux claim. Model/source to be supplied.

---

# METHODS

## M1 — Illumina short-read sequencing of syn1.0 transcriptome  
**Tex file:** `Manuscript/sections/methods/illumina_syn1.tex`  
**Analysis:** 
- Illumina mapping: `Syn1_Transcriptomics/Illumina/Illumina_Processing/01_quality_control.sh`, `02_alignment_seqdepth.sh`  
- Gene TPM: `Syn1_Transcriptomics/Gene_TPM/Gene_Transcriptomics.py`
**Key params to mention:** FastQC + MultiQC; bowtie2 paired-end; dUTP / fr-firststrand (R2 = transcript strand); per-strand bedGraph.  
**Inputs:** SRA accessions in `Syn1_Transcriptomics/Illumina/Illumina_Raw/00_retrive_fastq.sh` 
**Notes for LLM:** Subsubsection **Illumina MiSeq Read Processing, Mapping to the Genome** and **Quantification of TPM from Sequencing Depth** need to be written. 

---

## M2 — PacBio long-read sequencing of syn1.0 transcriptome  
**Tex file:** `Manuscript/sections/methods/pacbio_syn1.tex`  
**Analysis:** `Syn1_Transcriptomics/PacBio/PacBio_Processing/`  
**Key params:** FLNC recovery (reorientation, primer removal, polyA trimming); `minimap2 -ax map-hifi`; HQ filtering.  
**Outputs:** `syn1.PacBio.FLNC.sorted.HQ.bam`, `depth_bedgraph/…`  
**Inputs:** SRA accessions in `Syn1_Transcriptomics/PacBio/PacBio_Raw/00_retrieve_fastq.sh`  
**Notes for LLM:** This section is mostly finished and polished.  

---

## M3 — Operon identification from PacBio long-read transcriptomics in syn1.0
**Tex file:** `Manuscript/sections/methods/operon_analysis.tex`  
**Analysis:** `Syn1_Transcriptomics/Isoforms_PacBio/Cluster_Isoform.py`, `Syn1_Operon/…`  
**Key params:** clustering thresholds, min reads, TSS/TTS calling rule.  
**Outputs:** `isoform_clusters_annotated.tsv`, `operons.candidate_blocks.tsv`  
**Notes for LLM:** Subsubsection **Locate TSS, TTS, and RNA Cleavage Sites from PacBio RNASeq** needs to be changed and written.  

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
**Notes for LLM:** distinguish what we measured vs reused from prior work; cite original datasets; subsubsection **Relative Protein Quantification** and **Localization of Proteome** need to be written.

---

## M6 — Correlation between transcriptome and proteome
**Tex file:** `Manuscript/sections/methods/corr_transcriptome_proteome.tex`  
**Analysis:** `Syn1_Corr_RNA_Proteins/Transcription_Translation.py`, `Translation_Residual_L2_elongation.py`  
**Key params:** TPM source (PacBio vs Illumina); relative iPM; residual model covariates (CAI, tAI, TIR, etc.).  
**Outputs:** `syn1_genes_transcriptomics_proteomics.csv`, `residual_analysis/`  
**Notes for LLM:** Entire section needs to be written.  

---

## M7 — Novel Transcription and Translation
**Tex file:** `Manuscript/sections/methods/novel_orf.tex`  
**Analysis:** `Syn1_Novel_ORF/…`   
**Key params:** Mycoplasma genetic code (UGA = Trp); anti-SD `ACCUCCUUU`; ORF length/SD-distance thresholds.  
**Notes for LLM:** Entire section needs to be written.  

---

## M8 — Oxford Nanopore (ONT) and Illumina sequencing of syn3A transcriptome
**Tex file:** `Manuscript/sections/methods/ont_illumina_syn3A.tex`  
**Analysis:** `Syn3A_Transcriptomics/ONT/ONT_Processing/` (ONT) + `Syn3A_Transcriptomics/Illumina/Illumina_Processing/` (Illumina)  
**Key params:** ONT direct-RNA, `minimap2 -ax map-ont` (NOT splice, bacteria are intron-less), per-strand depth; Illumina syn3A paired-end bowtie2 (dUTP / fr-firststrand), per-strand bedGraph.  
**Inputs:** ONT raw `Syn3A_Transcriptomics/ONT/ONT_Raw/`; Illumina syn3A SRA accessions (SRR19432056/57 mate pair) via `Syn3A_Transcriptomics/Illumina/Illumina_Raw/00_retrive_fastq.sh`.  
**Notes for LLM:** This section is mostly polished. The file now covers BOTH ONT and Illumina-syn3A mapping (two subsubsections); rRNA operons at ~55,460 and ~343,267 bp drive the multi-mapping fraction.  

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

**Notes for LLM:** This is the longest Methods subsection — draft it in the same order as the pipeline; More details were recorded in the CLAUDE.md file.

---