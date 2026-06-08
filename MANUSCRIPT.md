# MANUSCRIPT.md — drafting guide for Methods & Results

This file is the **single source of truth** the LLM uses to draft Methods and
Results paragraphs. Each subsection here mirrors a `.tex` file under
`Manuscript/sections/{results,methods}/`. Fill in the bullets; the LLM will
expand them into prose, quoting numbers verbatim from the files you list.

TO DO Reminder:

**A. Pending analyses (these block their Results text + figure panels):**
- [~] R2 / M4 — RNA-processing/RNase analysis. DONE: L2.2 erosion (panel d) + L2.4 ATP-synthase α cleavage (panel f); R2 prose (`RNase.tex`) + M4 Methods (`RNA_processing.tex`) drafted. **L2.3 SCOPED 2026-06-08** to ship with: (1) 2° structures of 3 examples — lap 5', 0178 3', atpA/α RNase III site (user-built, ViennaRNA); (2) panel e proposed-hypothesis schematic (proteomics capacity + kinetic argument); (3) the B.subtilis→Syn1 homology / α RNase III site visualization (user). These three are mostly user-built → then wrap R2. The other-membrane-complex scan is also deferred (optional).
  - **DEFERRED (backlog) — genome-wide 3'-end 2°-structure test (L2.3 future work):** test the 3'→5'-exo signature transcriptome-wide. **Design is SETTLED — do NOT use the confounded version:** do NOT compare intragenic vs intergenic 3' ends (terminators are trivially structured → no evidence). Use within-class controls instead: (i) **terminal accessibility** — mean base-pairing prob of the last ~5 nt (ViennaRNA partition function) at intragenic 3' ends vs **dinucleotide-shuffled** windows (composition-matched null) → predict the 3' termini are more single-stranded (a loadable 3'-OH); (ii) **3'-vs-5' mirror-asymmetry meta-profile** — pairing probability aligned at the endpoint for intragenic 3' vs 5' ends → predict a structure peak just UPSTREAM (5') of 3' ends (the YhaM-stall stem) mirrored DOWNSTREAM (3') of 5' ends, both with accessible termini. N = 6,302 intragenic 3' / 9,992 5' unique positions (from `Syn1_RNase/RNase/isoform_endpoint_context.tsv`); ViennaRNA 2.6.4 in the RNAseq env; inputs ready.
- [ ] R5 L5.5 (panel e) — **DEFERRED** essentiality × trace-expression script; supply the syn3A essentiality source. 
- [ ] R6 L6.5 (no panel) — ATP/GTP flux comparison; supply the metabolic model/source. (panel d/e now = rPtn operon + tRNA junction, L6.3.)
- [x] R4 L4.5 — synthetic-element transcription script (yeast gene 0918 antisense; watermarks).
- [ ] R1 L1.5 (panel e) — decide which polycistronic operon to showcase (r-protein operon vs another complex).
- [x] R1 L1.4 (panel e,f,g) - RESOLVED 2026-06-08: the "canonical operons with no TransTermHP terminator" were a matching bug in `Operon_Annotation.py:find_terminators_near_tts` — the commented strict rule tested the terminator `end0` on BOTH strands, so '-'-strand terminators (whose 3' boundary is `start0`) were dropped (e.g. OP_00099/0178 missed its conf-100 TERM 82). Fixed to a strand-correct 3'-boundary window (`_term_3p`), re-ran: 97/127 canonical operons have a TTS terminator (was 98 via a midpoint hack; OP_00061 correctly dropped, OP_00099 now mapped). operons.tex L25 updated (97/127, 98% within 10 nt); term_* figures regenerated. Biology: a 3' terminator hairpin does NOT preclude 3' erosion (RNase R reads through structure) → L2.1/L2.3 reframe from "unstructured 3' end" to "limiting 3'→5' read-through capacity."

**B. Methods still to write:**
- [~] M4 — "RNA processing and ribonucleases" DRAFTED (ribonuclease inventory + RNA-processing endpoint analysis + RNase III/Y homology table; versions pinned). Pending: the B.subtilis→Syn1 RNase-site-mapping + 3'-end 2°-structure subsubsection (commented TODO, after tomorrow's analysis).

**C. Results prose (draft from logics, §0 style — one section at a time):**
- [x] Generate all six Results sections. R5 (`reduction_operons.tex`) drafted (L5.1–L5.4, 459-operon numbers); R6 (`reduction_omics.tex`) drafted (L6.1–L6.5, 459-era numbers; L6.3 = the 11 kb rPtn operon + tRNA-junction paragraph, panels d/e). Only blocked items remain: R5-L5.5 (panel e essentiality) and R6-L6.5 flux (no panel) — both await analyses A.
- [x] R2 (`RNase.tex`) DRAFTED 2026-06-07 (L2.1–L2.4 prose, panels a–f, ribonuclease table; L2.3 exo mechanism kept as a hedged prediction pending tomorrow's 2°-structure analysis). Compiles clean (durand_rnases_2018, redko_minimal_2013, janssen_tmrna_2012). **ALL SIX RESULTS SECTIONS + R2 NOW DRAFTED.**

**D. Figures:**
- [ ] Regenerate + format the six multi-panel figures in Illustrator.

**E. Numbers to re-verify (flagged during Methods drafting):**
- [x] M7 — recounted on revised (post Apr-22) clusters: 837 abnormal isoforms -> 29,443 candidate ORFs -> top-100 -> 48 unique / 47 proteotypic; 2 MS peptides reproduced (NOVEL_PEP_002; NOVEL_PEP_043 = old 030).
- [x] M8 — RE-VERIFIED (Syn3A_TPM.py rerun): the r=0.998 is **our depth-based Illumina TPM vs Palsson/Sandberg reported Illumina TPM** (n=458; Spearman 0.998), NOT Illumina vs ONT. Illumina vs ONT agreement is only Pearson r=0.570 / Spearman 0.558 (n=496). The earlier note conflated the two. syn3A Illumina alignment % still not locally re-verifiable (no logs).

**F. Style / front-and-back matter:**
- [x] §0 — fill the "Exemplar paragraph" field and the "Things to NEVER do" list.
- [ ] Draft Abstract, Introduction, and Discussion (~500 words each; not yet in this file).
- [x] Reorganize the Operon_Visual Jupyter notebook before publication. (Progress: `Operon_Annotation.py` cleaned to analysis-only and confirmed to run end-to-end; the per-operon `plot_one_operon` driver was split out into a new `Operon_Visualization.ipynb`.)
- [ ] Final LaTeX pass (resolve overfull \hbox lines / typesetting).

**H. Operon segmentation — RESOLVED (new segmentation finalized 2026-06-04):**
- [x] Reproducibility: the current `Operon_Segmentation.py`, run in conda env `RNAseq` (has pysam), reproduces the canonical `operons.candidate_blocks.tsv` (**459 operons**) byte-for-byte. The earlier 480-vs-483 gap is gone — the finalized map is 459 operons. Stage-by-stage counts are now persisted to `Syn1_Operon/Operon_Segmentation.txt` (mirrors `Operon_Annotation.txt`).
- [x] Dedup bug fixed: `dedup_operon_gene_lists` runs as the final pass before `to_csv`, so the tsv's `sense_gene_count` already equals the unique-loci count (max 21, 0 mismatches across all 459 operons). No consumption-time recount needed.
- [x] Co-transcription merge: the Step-5a overlap / coverage-hole merge is now read-evidence based (≥50 strand-specific PacBio bridging reads across the junction AND ≥0.5 gap/flank depth continuity), so the 11 kb rPtn supercluster (0652–0672) stays intact (61 candidate junctions → 43 pass → 37 merged operons). Running the script from `Syn1_Operon/` is safe — it regenerates the identical canonical map.
- [ ] Could use promoter and terminator predictions to judge the merging of truncated operons; then we can expand the statistics of transcription signatures to all operons; tell if internal promoter and terminator.
- [ ] Also checking the operon gene coverage.

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
**Referring to files:** For now, use files names in the Git repo; replaced with SI file names later.

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

**SI file:** operons.xlsx (boundaries, promoters, terminators, gene coverage, protein complex annotation)

### One-sentence Summary
**459 operons were identified using PacBio long-read RNAseq, with transcription signatures located.**

### Figure
**Figure:** `Manuscript/figures/operon.pdf`

- Panel a: Two-gene operon co-transcription and subsequent RNA processing (schematic).
- Panel b: Sense genes per operon.
- Panel c: Operon length distribution.
- Panel d: Promoter and terminator on a genome axis: -35/-10 logos at canonical TSSs and the dnaA-dnaN (MMSYN1_0001-0002) terminator.
- Panel e: Terminator loop-length distribution.
- Panel f: Terminator stem-length distribution.
- Panel g: Terminator 3' poly-U tail logo.
- Panel h: Operonal structure of a macromolecular-complex cluster (L1.5 showcase).

### Chain of Logics

#### L1.1: Co-transcription and further RNA processing complexify the transcriptome even of the reduced bacterium.

- **Logic:** Co-transcription with multiple transcription start and terminator sites generates multiple transcription units; RNA processing from endo- or exo-ribonucleases can digest the transcripts to even more **RNA isoforms**.
- **Analysis:** None
- **Outputs:** None
- **Numbers to cite:**  None
- **Figure panels:** a
- **Conclusion:** None
- **Caveats:** None
- **Notes for LLM:** This is a descriptive part. DONE — L1.1 drafted in `operons.tex`.

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
- **Notes for LLM:** More details presented in Methods M2. DONE — L1.2 drafted in `operons.tex`.

#### L1.3: 459 operons were mapped by full-length PacBio RNA seq.

- **Logic:** Unique longest isoform clusters as evidence of gene co-transcription were constructed by containment to cover 313 initial operons. Overlapping operons and clustering coverage holes were resolved by a read-level co-transcription test (PacBio bridging + depth continuity), collapsing the set to 275 isoform-derived operons. Uncovered genes were rescued. Finally, 459 operons for 911 genes in syn1. The statistics on the size and length of operons were reported.
- **Analysis:** 
  - Operon segmentation: `Syn1_Operon/Operon_Segmentation.ipynb`
  - Operon annotation: `Syn1_Operon/Operon_Annotation.ipynb`
- **Outputs:** 
  - Operons: `Syn1_Operon/operons.candidate_blocks.tsv`
- **Numbers to cite:**  
  - MIN_READS threshold = 50;
  - Sense genes per operon: median 1, mean 2.0 (911/459)
  - On average, 911/459 genes per operon
  - Operon length: median 1,650 bp, mean 2,049 bp
  - 69 operons (15.0%) enclose ≥1 antisense gene
  - Largest operon of ribosomal proteins and translation machineries: MMSYN1_0652–MMSYN1_0672 (21 genes, ~11 kb; OP_00341, 10,954 bp)
  - Nine operons cover only anti-sense genes or intergenic regions (0 sense genes)
- **Figure panels:** b,c
- **Conclusion:** None
- **Caveats:** RESOLVED — the segmentation now applies `dedup_operon_gene_lists` so `sense_gene_count` equals the unique-loci count (max 21, 0 mismatches across all 459 operons); panel b / mean are no longer inflated.
- **Notes for LLM:** More details are in Methods M3. DONE — L1.3 drafted in `operons.tex` (segmentation-algorithm paragraph + size-statistics paragraph).

#### L1.4: Transcription signatures located for operons.

- **Logic:** Exact transcription start site (TSS) and termination sites (TTS) were located for canonical operons whose boundaries were intergenic.
- **Analysis:** 
  - Annotation (analysis-only; runs end-to-end via `MPLBACKEND=Agg python Operon_Annotation.py`): `Syn1_Operon/Operon_Annotation.py`
- **Outputs:** 
  - `Syn1_Operon/Operon_Annotation.txt` (terminator statistics section: stem/loop/poly-U tract lengths, stem G+C, TTS→poly-U distance)
  - Promoter logos: `Syn1_Operon/annotation/canonical/` (`promoter_logo_minus35.pdf`, `promoter_logo_minus10.pdf`, `promoter_logo_tss.pdf`)
  - Terminator hairpins: `tts_hairpins/TERM_*.pdf` (clean 1×1 in, logomaker classic colors) + `tts_hairpins/all_hairpins.pdf` (labeled QC grid)
  - Terminator stat strips (each 7/3 × 7/9 in): `term_stem_length.pdf`, `term_loop_length.pdf`, `term_tail3_logo.pdf`, `term_stem_gc.pdf`, `term_tts_polyU_distance.pdf`
- **Numbers to cite:** 
  - 127 canonical operons (isoform_operon with TSS+TTS both intergenic)
  - -10 box: `TANAAT` hexamer 87/127 (69%); extended `TNNTANAAT` 52/127 (41%)
  - -35 box: no fixed hexamer, broadly AT-rich
  - terminators near TTS: 98/127 (77%)
  - stem median 8 bp (4–19); loop median 4 nt (3–7); poly-U tract median 5 nt
  - stem G+C: median 43%, mean 44% vs genome 24% (G+C-enriched, not absolute GC-rich)
  - TTS to 3' poly-U end: median 3 nt beyond, IQR 1–4, 97% within ±10 nt
- **Figure panels:** d (promoter + terminator logos), e (terminator stat strips: stem length, loop length, poly-U tail logo, stem G+C, TTS→poly-U distance)
- **Conclusion:** Signatures are both consistent with previous knowledge: -10 box of TSS has TANAAT, -35 just AT rich; TTS as intrinsic terminators (G+C-enriched hairpin + poly-U) whose 3' poly-U end coincides with the mapped TTS.
- **Caveats:** The TSS, TTS sites are only for canonical operons; we might need to refine for all cases. Stem G+C 40% is enriched relative to the AT-rich genome but not "GC-rich" in absolute terms — phrase as G+C-enriched.
- **Notes for LLM:** Corresponding Method finished (M3). DONE — L1.4 drafted in `operons.tex` (promoter sentence + terminator sentences with stem G+C and TTS→poly-U distance). Panel d: -35 / -10 promoter logos regenerated to OUTPUT.md spec (1×1 in, x-tick labels removed, enlarged ATCG). Panel e candidate = the five terminator stat strips (user is assembling); tex keeps Fig. d refs until panel letters are fixed.

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

- Panel a: Legends of erosions of RNA isoforms and of RNA processing molecules (shared top strip: 4 erosion categories + 5 RNase/ribosome/tmRNA icons). (7, 7/10) — `Syn1_RNase/R2_panels/R2a_legend_strip.pdf` (icons also as standalone `R2_icon_*.pdf`).
- Panel b: RNA isoform distribution for gene 0178 with more 3' erosion. (7/2, 7/4) — `Syn1_RNase/R2_panels/R2b_0178_3p_erosion.pdf`
- Panel c: RNA isoform distribution for gene 0154/lap with more 5' erosion. (7/2, 7/4) — `Syn1_RNase/R2_panels/R2c_lap_5p_erosion.pdf`
- Panel d: 3' secondary structure for gene 0178 (7/4, 7/4) — `Syn1_RNase/R2_panels/R2d_0178_3prime_structure.pdf` (from `fold_3prime_terminator.py`, mirrored from `terminator_3p/`)
- Panel e: RNA isoform truncation categories. (7/4, 7/4) — `Syn1_RNase/R2_panels/R2e_truncation_categories.pdf`
- Panel f: Biased RNA Processing schematics: endo and exo from 3'. (7/2, 7/4) — Illustrator (no matplotlib file)
- Panel g: RNA isoform distributions for ATP synthase operon — isoforms split into two regions at atpA/α (0792) where RNase III cuts; coloured by 5'-block (a,c,b,δ; teal) vs 3'-block (γ,β,ε; orange), gene arrows tinted by block, depth steps down at the α cut. (7, 7/3) — `Syn1_RNase/R2_panels/R2g_atp_synthase.pdf` (matplotlib: gene arrows + isoforms + depth only; the F1/F0 scheme, SD strengths, subunit labels, "RNase III on α" scissors are added in Illustrator).

### Chain of Logics

#### L2.1: Distinct RNA isoforms distributions found for operons.

- **Logic:** Truncated isoforms compared to the full transcription units exist for operons because of the RNA processing. Distinct patterns of truncations can be found, using genes 0154 and 0178 as examples. 0178 has structured 3' end as shown in d, but RNase R can digest through the dsRNA structured region.
- **Analysis:** None
- **Outputs:** None
- **Numbers to cite:**  None
- **Figure panels:** b,c,d
- **Conclusion:** None
- **Caveats:** None
- **Notes for LLM:** None

#### L2.2: Transcriptome-wide, significantly more 3' erosion found.

- **Logic:** Overlaid RNA isoforms to gene ORFs to find significantly more 3' erosions.
- **Analysis:** `Syn1_RNase/RNA_Processing.py` (endpoint-context: each isoform 5'/3' end labelled intragenic vs intergenic; + ORF start-without-stop)
- **Outputs:** 
  - `Syn1_RNase/RNA_Processing.txt`
  - plots + tables in `Syn1_RNase/RNase/`
- **Numbers to cite:** (current clustering, n_reads>=10, 20,885 isoforms) full-length/unprocessed 8.3% isoforms / 27.0% reads; 3'-intragenic-only 30.3% iso / 47.5% reads (the dominant non-canonical class); both-eroded 41.9% iso; 5'-intragenic-only 19.5% iso; 5'/3' intragenic ratio 0.40 by reads (3' erosion dominates); 42.2% of contained ORFs have a start codon but no stop codon (read-weighted 38.9%).
- **Figure panels:** e
- **Conclusion:** Most full-length isoforms are processed; endpoints fall inside ORFs far more often at the 3' end, evidence of biased 3'-directed erosion.
- **Caveats:** intragenic-endpoint logic counts a 3' end inside a gene body as processing (cannot be a terminator); does not distinguish endo- from exo-nucleolytic origin.
- **Notes for LLM:** Make this logic short and to the point.

#### L2.3 (HYPOTHESIS): The 3'-bias reflects an asymmetric exonucleolytic *clearance* bottleneck — limited 3'→5' read-through (scarce RNase R, stalling YhaM) and ribosome trapping on non-stop ends — not biased endonucleolytic cutting.

- **Logic (the chain):** Endonucleolytic cleavage (RNase Y/III) is *symmetric* — each cut yields one upstream fragment (intact 5' + a NEW intragenic 3' end) and one downstream fragment (NEW intragenic 5' end + intact 3' terminator) — so endo cutting alone cannot create a 3'/5' asymmetry; equal numbers of each are born. The observed bias must therefore arise from differential *clearance* of these fragments. (i) Downstream / 5'-eroded fragments are cleared rapidly by the abundant 5'→3' machinery (RNase J1 + J2) → low steady state → rarely captured. (ii) Upstream / 3'-eroded fragments can only be *fully* erased by RNase R — the lone 3'→5' exo that reads THROUGH structure (it needs an unstructured ss 3' overhang to load, then unwinds via intrinsic helicase activity); the abundant YhaM is a "generator not finisher" — it trims the ss 3' tail but STALLS at the first stem base, manufacturing 3' ends rather than removing them. So 3'-eroded intermediates are both over-produced (YhaM stalls) and under-cleared (RNase R scarce) → they accumulate and dominate the steady-state long-read pool. **Two reinforcing mechanisms:** (A) **lower 3'→5' read-through capacity** (RNase R 36 vs RNase J1+J2 ~234 copies; YhaM 117 stalls at structure); (B) **ribosome trapping** at non-stop 3' ends (the 42% start-but-no-stop ORFs from L2.2) physically blocks exo entry until tmRNA–SmpB rescue, which is limiting (SmpB 14).
- **Analysis — SCOPED 2026-06-08:** the genome-wide structure test is DEFERRED (full design parked in TODO §A); R2 ships with three worked examples + the kinetic argument. Three strands of support:
  1. **2° structure of the three examples** (user-built, ViennaRNA 2.6.4 in RNAseq env): the 5' erosion region of the lap operon (panel b), the 3' erosion region of the 0178 operon (panel c), and the atpA/α RNase III cleavage site (panel f) — qualitative illustration that the eroded 3' ends sit at accessible / stall-competent structures.
  2. **Ribonuclease capacity asymmetry (proteomics, in hand, Table S1):** 5'→3' (RNase J1+J2 ~234) vs 3'→5' read-through (RNase R 36); YhaM 117 (stalls), SmpB 14 → panel e proposed-hypothesis schematic.
  3. **Non-stop / ribosome-trapping link (in hand):** the L2.2 42%-no-stop ORFs × the limiting tmRNA–SmpB capacity.
- **Outputs:** Table S1 ribonuclease abundances (in hand); the three example 2° structures (user-built); genome-wide ViennaRNA test DEFERRED.
- **Numbers to cite:** 5'→3' RNase J2 142 + J1 92 = ~234 copies vs 3'→5' read-through RNase R **36**; YhaM 117 (Mn2+-dependent, stalls at stem base); SmpB 14; tmRNA 2.5% of non-rRNA. RNase R loads on a ss 3' overhang ≥7 nt (optimal ≥10). [structure-enrichment numbers PENDING]
- **Figure panels:** f (biased-processing schematic).
- **Conclusion:** The 3'-erosion bias is best explained by a 3'→5' clearance bottleneck (scarce read-through RNase R; abundant but structure-stalling YhaM) compounded by ribosome trapping on non-stop products under limiting trans-translation — not by biased endonucleolytic cutting.
- **Caveats:** copy number is a capacity *proxy*, not measured flux (RNase R is processive, so 36 copies aren't negligible); repeated endo cuts can substitute for RNase R in clearing structured RNA; the ribosome-trapping arm applies only to translated (CDS) 3' ends; long-read capture biases.
- **Notes for LLM:** Literature backing (web-searched, verified): **RNase R** requires an unstructured ss 3' overhang to LOAD but then degrades THROUGH structure via intrinsic helicase activity (Chu et al. 2017 PubMed 29036353; "How RNase R Degrades Structured RNA" JBC 2016) — so its "unstructured 3'" requirement is at INITIATION only; **YhaM** is Mn2+-dependent, prefers a ss 3' docking site and STALLS at the base of secondary structure (Oussenko et al. 2002, J Bacteriol 184:6250); RNase R is also a principal nuclease of non-stop mRNA decay with trans-translation (need a citation before asserting in text). **Table S1 RNase R row should be nuanced** to "needs ss 3' overhang to load, then reads through structure" (not "only unstructured"). The tex L2.3 paragraph can now state the capacity-asymmetry + ribosome-trapping HYPOTHESIS with the Table S1 numbers; the structure evidence is the three example folds (qualitative), with the genome-wide test framed as future work. Candidate bib keys: durand_rnases_2018, durand_three_2012, redko_minimal_2013, janssen_tmrna_2012.

#### L2.4: ATP synthase operon is co-expressed in one-go but cut at $\alpha$ subunit.

- **Logic:** Macromolecular complexes' gene co-expression can be altered by RNA processing. ATP synthase's RNA isoform distribution has a clear pattern of isolation at $\alpha$ subunits, which was identified as endo RNase III cleavage site. Comment on the other membrane complexes.
- **Analysis:** `Syn1_RNase/R2_figure_panels.py` (panel_f).
- **Outputs:** 
  - `Syn1_RNase/R2_panels/R2f_atp_synthase.pdf`
- **Numbers to cite:** the atp operon (minus strand) is segmented into two overlapping operons that meet AT atpA/$\alpha$ (MMSYN1_0792): the 5'-block OP_00395 (0797–atpH/0793 + 5' of $\alpha$; 12 member isoforms, top isoform 1,161 reads) and the 3'-block OP_00394 (3' of $\alpha$ + atpG/atpD/atpC = 0791–0789; 6 members); the RNase III cut at $\alpha$ (~933.78 kb) drops the minus-strand depth from ~9k to ~half across the junction.
- **Figure panels:** g
- **Conclusion:** The RNase complexifies the subunit synthesis of complexes.
- **Caveats:** None
- **Notes for LLM:** Panel-f matplotlib (gene arrows + 2-region-coloured isoforms + depth) DONE 2026-06-07; the F1/F0 scheme, SD strengths, subunit labels and "RNase III on $\alpha$" scissors are Illustrator. Still TODO: check OTHER membrane complexes for the same pattern; visualize the RNA secondary structure at $\alpha$.

---

## R3 — High Correlation between transcriptome and proteome in the reduced organism.

**Tex file:** `Manuscript/sections/results/corr_RNA_ptn.tex`

**SI file:** `Syn1_Corr_RNA_Proteins/syn1_omics.xlsx` (all 911 genes; columns: locusTag, gene_name, rna_type, gene_product, protein_localization, TPM_illumina, TPM_PacBio, iPM_mean, protein_copy_number, TIR, CAI, protein_halflife_h; gene name/product taken from the syn3A proteome where an ortholog exists). Built by `Syn1_Corr_RNA_Proteins/R3_figure_panels.py`.

### One-sentence Summary
**High correlation found between transcriptome and proteome.**

### Figure
**Figure:** `Manuscript/figures/correlation.pdf`

- Panel a: Per-protein copy-number distribution by localization (cytoplasmic, lipoprotein, membrane, extracellular).
- Panel b: Illumina sense TPM vs proteome iPM (log10).
- Panel c: PacBio vs Illumina sense TPM (log10).
- Panel d: Predicted TIR vs proteome residual.
- Panel e: CAI vs proteome residual.
- Panel f: Model Pearson R for whole proteome and cytoplasmic proteins, with/without CAI.
- Panel g: Intrinsic protein half-lives transferred from Mpn.
- Panel h: Protein half-life vs proteome residual.

### Chain of Logics

#### L3.1: Distribution of protein copy numbers in syn1

- **Logic:** Starting from the iBAQ after the mapping using SpectroNaut software, iPM was calculated, which was then converted to absolute protein copy numbers with protein dry mass per cell; Cytoplasmic protein has a median copy number of 47, while only 10 for membrane proteins.
- **Analysis:** 
  - Copy number quantification: `Syn1_Syn3A_Proteomics/Protein_Quantification_Localization.ipynb`
- **Outputs:** 
  - syn1 proteomics: `Syn1_Syn3A_Proteomics/syn1_proteomics_localization_2026.csv`
- **Numbers to cite:**  721/828 annotated proteins detected (87%); median copies/cell: cytoplasmic 47 (n=516), lipoprotein 21 (n=68), membrane 10 (n=126), extracellular 3 (n=11); most abundant = EF-Tu (tuf, MMSYN1_0151) ~7,200 copies
- **Figure panels:** a
- **Conclusion:** Poor coverage of membrane proteins because of the protease digestion; 
- **Caveats:** None
- **Notes for LLM:** Method M5 done (Relative/Absolute Protein Quantification + Localization). DONE — Results prose drafted in `corr_RNA_ptn.tex` (L3.1 copy-number paragraph).

#### L3.2: Using Illumina TPM as standard of transcriptome quantification

- **Logic:** PacBio and Illumina TPMs of syn1 were correlated to get r of 0.62; no significant TPM and length bias was found.
- **Analysis:** `Syn1_Transcriptomics/Gene_TPM/Gene_Transcriptomics.py`
- **Outputs:** plots in `Syn1_Transcriptomics/Gene_TPM/`
- **Numbers to cite:**  PacBio vs Illumina sense TPM Pearson r=0.62 (log10, n=884 at TPM≥0.5; figure threshold low_threshold=0.5); no abundance/length bias
- **Figure panels:** c
- **Conclusion:** As convention, Illumina TPMs were used to do correlation.
- **Caveats:** Do NOT say "Illumina TPM used since correlate with iPM better."
- **Notes for LLM:** Methods M1/M2/M8 done. DONE — Results prose drafted in `corr_RNA_ptn.tex` (L3.2). r=0.62 verified = the panel-c figure value (TPM≥0.5); all-feature r=0.61, mRNA-only r=0.68.

#### L3.3: Decent correlation found between transcriptome and proteome for syn1

- **Logic:** Pearson r of 0.7 found between two omics for cytosolic proteins; lower r for all since poor coverage of membrane proteins.
- **Analysis:** `Syn1_Corr_RNA_Proteins/Transcription_Translation.py`
- **Outputs:** 
  - Same name Txt file
- **Numbers to cite:**  all proteins Pearson r=0.61 (R²=0.38, n=717; Spearman 0.67); cytosolic-only r=0.70 (R²=0.49, n=512; Spearman 0.75)
- **Figure panels:** b
- **Conclusion:** Decent correlation.
- **Caveats:** None
- **Notes for LLM:** DONE — Results prose drafted in `corr_RNA_ptn.tex` (L3.3). The old tex heading's "0.67 for cytosolic" was STALE — cytosolic Pearson is 0.70 (0.6971); corrected in the prose.

#### L3.4: Predicted TIR had low correlation with residuals between two omics

- **Logic:** Translation initiation rates predicted for all gene ORFs by OSTIR were evaluated to validate how much of the residuals can be explained by translation initiation.
- **Analysis:** `Syn1_Corr_RNA_Proteins/Translation_Residual_L1_initiation.py`
- **Outputs:** 
  - Same name Txt file
  - plots under `Syn1_Corr_RNA_Proteins/residual_analysis`
- **Numbers to cite:**  residual vs log10(TIR) Pearson r=0.17; ΔR² from TIR +0.020 (2% of variance); baseline on n=566 UTR-genes r=0.59 / R²=0.35
- **Figure panels:** d
- **Conclusion:** No improvement on correlations, could be because the prediction itself was inaccurate.
- **Caveats:** None
- **Notes for LLM:** Method M6 (TIR) done. DONE — Results prose drafted in `corr_RNA_ptn.tex` (L3.4); cite OSTIR `roots_ostir_2021`.

#### L3.5: Translation elongation factor improved the correlations

- **Logic:** CAI as metric for translation elongation efficiency significantly improved the correlation.
- **Analysis:** `Syn1_Corr_RNA_Proteins/Translation_Residual_L2_elongation.py`
- **Outputs:** 
  - Same name Txt file
  - plots under `Syn1_Corr_RNA_Proteins/residual_analysis`
- **Numbers to cite:**  CAI vs residual Pearson r=0.36; ΔR² from CAI +0.080 (+21% of baseline); cytosolic ΔR² +0.052; CAI reference set = top 20% by iPM = 144/721 genes
- **Figure panels:** e,f
- **Conclusion:** Translation elongation affected the protein biosynthesis.
- **Caveats:** None
- **Notes for LLM:** Method M6 (CAI) done. DONE — Results prose drafted in `corr_RNA_ptn.tex` (L3.5); cite CAI `sharp_codon_1987`.

#### L3.6: Protein degradation had low correlation with residuals.

- **Logic:** Intrinsic protein degradation as a result of protease activities in syn1 was transferred from Mpn by finding the reciprocal homologs; Lowest half-life was 4.7 hours, way longer than the half-life of 1 hour for syn1.
- **Analysis:** 
  - `Syn1_Corr_RNA_Proteins/Translation_Residual_L3_degradation.py`
  - Homology build: `Genomes_Input/Homology_Build.py`
- **Outputs:** 
  - Same name Txt file
  - plots under `Syn1_Corr_RNA_Proteins/residual_analysis/`
  - Homology under `Genomes_Input/homology_syn1_mpn/`
- **Numbers to cite:**  245 Mpn-mapped half-lives; Syn1-corrected median 32 h, shortest 4.7 h vs ~1 h doubling; residual vs log10(t½) Pearson r=−0.12 (P=0.06, n=242); ΔR² from half-life ≤0.01
- **Figure panels:** g,h
- **Conclusion:** Intrinsic protein degradation is way slower than protein synthesis and doubling in syn1, thus having a minor effect on proteome abundances.
- **Caveats:** Only a subset of proteins found reciprocal homologs; the intrinsic half-lives were corrected by protease (Lon or FtsH) abundances.
- **Notes for LLM:** Method M6 (degradation) done. DONE — Results prose drafted in `corr_RNA_ptn.tex` (L3.6); cite Mpn half-lives `burgos_protein_2020`.

---

## R4 — Novel transcription and translation activities in the synthetic bacterium, Syn1

**Tex file:** `Manuscript/sections/results/novel.tex`

### One-sentence Summary
**Long-read RNA seq reveals anti-sense and intergenic transcription.**

### References:
- Syn1, Science, 2010
- antisense RNAs as noise, Science Advances, 2016

### Figure
**Figure:** `Manuscript/figures/novel.pdf`

<!-- - Panel a: Distribution of anti-sense percentage of all isoforms. (7/4,7/4) -->
- Panel a: Schematics of read-through, embedded, and spurious promoter to explain anti-sense coverage in operons. (7/2,7/4)
- Panel b: Isoform length and read distribution for all three cases of anti-sense transcription. (7/4, 7/4)
- Panel c: Isoform read distribution for all three cases of anti-sense transcription. (7/4, 7/4)
- Panel d: Unexpected anti-sense transcription inside inserted Yeast vector. (7/2,7/4)
- Panel e: Intergenic coverage of all RNA isoforms (7/4,7/4)
- Panel f: Distribution of 5' and 3' untranslated regions (UTR) in operons. (7/4, 7/4)
- Panel g: Truly intergenic transcription between 0154 and 0155. (7/2, 7/4)
- Panel h: New ORF encoding peptide of length 118 aas found after gene 0592. (7/2, 7/4)

### Chain of Logics

#### L4.1: Minor percentage of isoforms have anti-sense coverages that can be categorized into three cases.

- **Logic:** 1.4% of the isoforms have antisense transcription coverage. About two thirds came from the spurious promoters in the AT-rich genome; the transcription read-through can cause anti-sense transcription at the end or embedded inside the operons.
- **Analysis:** 
  - `Syn1_Novel_ORF/Novel_translation.ipynb`
  - `Syn1_Novel_ORF/Abnormal_Transcripts.py`
- **Outputs:** 
  - Same name Txt file
- **Numbers to cite:** 1.4% antisense (267 isoforms -> 89 clusters); spurious 59 (66%), read-through 30 (34%, incl. 4 embedded)
- **Figure panels:** a
- **Conclusion:** Full-length RNA isoforms reveal new cases of anti-sense transcription as read-throughs.
- **Caveats:** None
- **Notes for LLM:** DONE — `novel.tex` drafted (L4.1–L4.6) against the final 8-panel figure (`R4_dist_panels.py` b/c/e/f + `R4_track_panels.py` a/d/g/h). Cites gibson_creation_2010, llorens-rico_bacterial_2016, roots_ostir_2021.

#### L4.2: Spurious promoter has highest isoform read support out of three cases, but still minor compared to sense-transcription.

- **Logic:** Spurious promoter has most isoform read support; embedded case tend to be in longer length; all three cases were minor in read compared to sense transciption, except one case of 0918, yeast gene.
- **Analysis:** 
  - `Syn1_Novel_ORF/Novel_translation.ipynb`
- **Outputs:** 
  - Same name Txt file
  - plots under `Syn1_Novel_ORF`
- **Numbers to cite:** per-case counts 59 spurious / 26 read-through / 4 embedded; spurious has the highest read support (max ~16.6k reads), embedded the longest isoforms (median ~2.2 kb) (Abnormal_Transcripts.txt)
- **Figure panels:** b,c
- **Conclusion:** All three antisense classes are minor versus sense transcription; read-through has the highest typical (median) support, embedded the longest isoforms, and one spurious-promoter locus (his3/0918) is the sole high-abundance exception.
- **Caveats:** None
- **Notes for LLM:** Panel a,b, and c are in the same row. DONE.

#### L4.3: Unexpected transcription of yeast vector gene 0918.

- **Logic:** Yeast vector elements were carried into syn1's synthetic genome during assembly. Strikingly, the yeast selection marker his3/0918 was heavily transcribed antisense (depth >30k) but not translated. The antisense isoforms initiate from a spurious promoter just upstream, which is itself deleted in syn3A; the his3 body sits at the deletion boundary, so the non-coding transcription is removed by minimization.
- **Analysis:** `Syn1_Novel_ORF/R4_track_panels.py` (panel d): isoform table + PacBio plus-strand depth + deletion overlay from `Genome_Reduction/aln/raw/syn1_deleted_regions.bed`. -10 box at the antisense TSS scored by `R4_track_panels.py:quantify_novel_promoters()` (same algorithm as canonical operons, via `Syn1_Operon/promoter_motif.py`).
- **Outputs:** `R4_panels/panel_d_his3_antisense.pdf`; `R4_panels/novel_promoter_minus10.txt`
- **Numbers to cite:** his3/0918 antisense depth >30,000; 28 antisense isoforms (top 16.6k reads); the driving spurious-promoter region deleted in syn3A; the antisense TSS (pos5p0 27522) carries a perfect -10 hexamer TAAAAT (TANAAT consensus, 0 mismatch; core_6mer tier) with an AT-rich -35 (CTTTGAA), confirming a genuine spurious promoter. Watermarks W1-W4 (located by exact sequence): length-weighted mean PacBio depth ~283/+ , ~360/- vs genome-wide average ~2133/2051 (6-8x lower); covered ORFs all hypothetical/watermark calls (plus real 0590) -> minimally transcribed noise.
- **Figure panels:** d
- **Conclusion:** The yeast marker his3 is heavily transcribed antisense yet untranslated, and its driving spurious promoter (which carries a canonical TAAAAT -10 box, i.e. a real but mislocated sigma-factor promoter) is deleted in syn3A; the four watermarks are only minimally transcribed (noise).
- **Notes for LLM:** DONE. his3 itself is RETAINED at the deletion boundary -- only its upstream promoter region is deleted, so do not claim the gene is deleted. Watermark expression quantified in `Abnormal_Transcripts.py` -> `watermark_expression.txt` (located by sequence, not in panel d); one sentence added to L4.3 in novel.tex.

#### L4.4: RNA isoforms has much more intergenic coverage than anti-sense.

- **Logic:** UTR and regions between genes contributed to the intergenic coverage. lengths were evaluated for all canonical operons with median of tens of nucleotides.
- **Analysis:** 
  - `Syn1_Operon/Operon_Annotation.ipynb`
- **Outputs:** 
  - Output inside Notebook
- **Numbers to cite:**  Median, and maximum
- **Figure panels:** e,f
- **Conclusion:** Median values were consistent with previous papers; high outliers are due to the anti-sense transcription or truncations.
- **Caveats:** None
- **Notes for LLM:**  Panel d,e and f in the same row.

#### L4.5: One truly isolated intergenic transcription.

- **Logic:** One truly isolated intergenic transcription was found between genes 0154 and 0155.
- **Analysis:** Operon coverage analysis in `Operon_Segmentation.ipynb`; -10 box at the intergenic TSS scored by `R4_track_panels.py:quantify_novel_promoters()` (via `Syn1_Operon/promoter_motif.py`).
- **Outputs:** 
  - `Syn1_Operon/operons.candidate_blocks.tsv`
  - `R4_panels/novel_promoter_minus10.txt`
- **Numbers to cite:**  intergenic TSS (pos5p0 199379) has NO canonical -10 box (best TANAAT match TAAAAA, 1 mismatch; no_minus10 tier; -35 TATTGTA) — contrast with the his3/0918 spurious promoter.
- **Figure panels:** g
- **Conclusion:** A single genuinely intergenic transcript exists (between lap/0154 and pseudo/0155, ~980 reads, 199,379-200,123); it lies in a syn3A-deleted region, while lap occupies the retained gap and was relocated. Unlike the his3 spurious promoter, its TSS lacks a recognizable -10 box (1 mismatch from TANAAT), consistent with pervasive rather than promoter-driven transcription.
- **Caveats:** None
- **Notes for LLM:** Panel g and h in the same row. DONE.


#### L4.6: Two novel peptides identified by enumerating all possible ORFs in isoforms having high abnormal fraction.

- **Logic:** The possible translation of the abnormal RNA isoforms was checked by enumerating all ORFs using OSTIR on the abnormal RNA isoforms. 
- **Analysis:** 
  - `Syn1_Novel_ORF/Novel_translation.ipynb`
- **Outputs:** 
  - `Syn1_Operon/operons.candidate_blocks.tsv`
- **Numbers to cite:** 837 abnormal isoforms -> ~29,000 candidate ORFs -> top 100 -> 48 unique / 47 proteotypic; 2 MS-confirmed (NOVEL_PEP_002, 118 aa intergenic near 0592; NOVEL_PEP_043 = old 030, 225 aa, 54-aa N-term extension of 0768), both deleted in syn3A
- **Figure panels:** h
- **Conclusion:** Two predicted ORFs were identified in Mass-spec proteome, and both were located near less annotated genes. Also, these two regions were deleted in syn3A.
- **Caveats:** Only top 100 ORFs were selected to do the new proteomics search, thus we cannot assure if all ORFs were translated or not (leave this question to the reviewers); the new canonical cluster isoforms gave new ORF candidiates, which were highly similar to the old ones that searched against raw proteoimcs.
- **Notes for LLM:** See Method **Novel Open-Reading Frames from Full-length RNA Isoforms** for details; the Analysis Jupyter Notebook can be cleaned up.

---

## R5 — Operonal structure changes to the minimal cell, JCVI-syn3A
**Tex file:** `Manuscript/sections/results/reduction_operons.tex`

### One-sentence Summary
**Halving the genome was a gene-order-preserving deletion campaign that excised whole operons, decapitated some retained operons by deleting their promoters, and fused only a small number of new cross-junction transcription units.**

### Figure
**Figure:** `Manuscript/figures/genome_reduction.pdf`

- Panel a: Schematics of genome reduction from syn1 to syn3A. (7/2, 7/2)
- Panel b: 0083 and rpsT/0082 co-expressed in syn1. (7/2, 7/4)
- Panel c: 0094 and 0082 co-expressed in syn3A. (7/2, 7/4)
- Panel d: Box plot of gene expressions to highlight decapitated ones having lower values. (7/4, 7/4)
- Panel e: the HupA operon, whose true promoter, located inside gene 0349, was deleted. (21/4, 7/4)
- Panel f: Gene essentiality evaluation for those trace-expressed genes that are still essential.

Extended Figure:

- rpsO situation same as rpsT/0082

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

#### L5.2: Deletions overlaid on syn1's 459 operons show whole-operon excision dominating over partial truncation.

- **Logic:** Intersecting the 95 deletions with the 459 syn1 operons at single-bp resolution classifies how each operon was hit, separating operons removed wholesale from those left partially truncated; the truncations are what create the junction effects in L5.3 and L5.4.
- **Analysis:** `Genome_Reduction/04_deletion_overlaid_operon.py`
- **Outputs:** `Genome_Reduction/deletion_overlaid_operon/operon_deletion_classification.tsv`
- **Numbers to cite:** span-level overlap_class (n=459): fully_deleted 181, intact 162, 3'_truncation_gene 47, 5'_truncation_gene 28, intra_truncated 17 (plus 9 UTR-only and 15 multi-hit); gene-level gene_deletion_pattern: all_deleted 235 (51.2%), intact 172 (37.5%), leading_deleted 21 (4.6%), lagging_deleted 20 (4.4%), intra_deleted 11 (2.4%); 414 syn1 genes overlapped by a deletion.
- **Figure panels:** a
- **Conclusion:** Reduction preferentially removed entire operons; the minority of partial truncations (5' vs 3') sets up the junction taxonomy.
- **Caveats:** the two axes (span-level truncation vs gene-level deletion) differ by design; 162 vs 172 "intact" reflects operons whose genes are all kept but whose UTR/flank was nicked.
- **Notes for LLM:** None.

#### L5.3: Same-strand deletion junctions can fuse new transcription units, but true fusion is rare.

- **Logic:** Each deletion is recast as a junction between the nearest retained operons on either side; relative orientation (tandem/convergent/divergent) and facing-regulator loss decide whether a new co-transcribed unit can form, and ONT spanning/bridging reads test whether the new cross-junction gene pair is actually co-transcribed.
- **Analysis:**
  - junction taxonomy: `Genome_Reduction/05_deletion_junction.py`
  - read validation: `06_single_operon_coexpression.py`, `07_operon_pair_coexpression.py`, `coexpression_common.py`
- **Outputs:**
  - `Genome_Reduction/deletion_junction/deletion_junctions.tsv`, `deletion_junction_summary.txt`
  - `Genome_Reduction/operon_pair_coexpression/`, `single_operon_coexpression/`
- **Numbers to cite:** 95 junctions: tandem 53, convergent 19, divergent 15, intra_operon 8; tandem junction_type: fusion 3, decapitation 9, readthrough_extension 11, clean_excision 30; cross-junction co-transcription (loose): fusion 67% (2/3) vs clean_excision 10% (3/30, negative control); pristine single-operon baseline preserved_loose 60% (45 testable, 111 pairs); fusion exemplar DEL_014 OP_00043 -> OP_00050 (MMSYN1_0094 -> MMSYN1_0082 = rpsT/S20), n_span=2, n_bridge=37, however the TPM FC was still low as 0.074 for rpsT/0082 since the fused promoter of 0094 is weak (in syn1, 0082 was co-transcribed with 0083 instead); a second r-protein rpsO/S15 (MMSYN1_0294) followed the same route (lost its own promoter, gained a weak fused one) and likewise collapsed to TPM FC 0.036. Both rpsT/0082 and rpsO/0294 are gene_impact_class new_promoter_fusion (from 08).
- **Figure panels:** b,c in the same column
- **Conclusion:** Operon fusion is real but rare (3 events); the dominant junction outcome is clean excision of whole operon(s) between intact neighbors.
- **Caveats:** ONT depth is low, so most positive calls are loose-bridge rather than strict-spanning; convergent/divergent junctions are opposite-strand and not expected to co-transcribe.
- **Notes for LLM:** This logic, not "all new operons fused," supports the reworded one-sentence summary. The weak-fused-promoter expression story (rpsT/0082, rpsO/0294 both crash despite structural fusion) lives HERE in L5.3, not L5.4, so the two paragraphs do not overlap; L5.4 covers only pure decapitation (hupA).

#### L5.4: Decapitated operons that lost their own promoter are the one class that robustly drops in expression; HupA is the showcase.

- **Logic:** Classifying every retained gene by promoter-source change isolates operons whose own promoter was deleted (promoter_lost / decapitation); their syn3A TPM is compared against the other impact classes to test whether promoter loss, not sequence change, predicts lower expression.
- **Analysis:**
  - per-gene impact: `Genome_Reduction/08_delete_gene.py`
  - expression: `09_Compare_RNA_Protein.py`
- **Outputs:**
  - `Genome_Reduction/delete_gene/retained_gene_context.tsv` (`gene_impact_class` column)
  - `Genome_Reduction/Compare_RNA_Protein/TPM_FC_by_impact_class.pdf`
- **Numbers to cite:** gene_impact_class (retained genes): promoter_lost 42, promoter_disconnected 6, new_promoter_fusion 3, readthrough_exposed 24, promoter_proximity_changed 17, context_only 45, unaffected 360; promoter_lost is the only class robustly down in TPM (median FC 0.44, Mann-Whitney p=2.7e-4 vs unaffected median 0.76); HupA (MMSYN1_0350) relTPM 6.68 -> 0.13 (FC 0.020), relIPM 6.48 -> 0.092 (FC 0.014); HupA operon -10 box = perfect TANAAT (TATAAT), extended TNNTANAAT match, strong_9mer tier, -10 window 441019-441024 inside deleted DEL_050 (440092-441059) [from R5_panels/R5_panel_stats.txt via promoter_motif.scan_minus10]. NOTE: the weak-fused-promoter r-proteins rpsT/0082 and rpsO/0294 belong to L5.3 (new_promoter_fusion), NOT here. rpmE/L31 (0137) and rpsU/S21 (0482) are gene_impact_class unaffected (operon structure intact) with only mild TPM dips FC 0.198 and 0.413, so they are NOT decapitation cases either.
- **Figure panels:** d,e in the same row
- **Conclusion:** Promoter-source loss drives the largest expression decreases; promoter_lost is the only impact class robustly down in TPM.
- **Caveats:** the class is assigned at operon level; 8 junctions lose only UTR (genes intact); the 05-vs-04 consistency check flags 2 flank operons as all_deleted.
- **Notes for LLM:** HupA's true promoter sits inside deleted gene MMSYN1_0349 (panel e); HupA (promoter_lost) is pure decapitation, no replacement promoter, and is the ONLY decapitation gene featured in L5.4. The weak-fused-promoter r-proteins (rpsT/0082, rpsO/0294) are covered in L5.3, not here. rpmE/L31 and rpsU/S21 were NOT affected by any deletion (operon structure intact), so do not cite them as decapitated.

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
**Figure:** `Manuscript/figures/reduction_omics.pdf`

- Panel a: mRNA pool compositions in syn1 and syn3A as secondary protein functions.
- Panel b: Significant mRNA pool share changes from syn1 to syn3A as tertiary functions.
- Panel c: Transcription and translation changes of RNAP, degradosome and enzymes in central metabolism.
- Panel d: The giant ~11 kb ribosomal-protein operon OP_00341 (MMSYN1_0652–0672), one polycistron, no internal terminator; gene track + RNA isoforms + depth; supplies 12% of the syn1 → 34% of the syn3A coding mRNA pool. (14/3, 7/3). [L6.3. Replaces the blocked ATP/GTP flux panel; the flux prediction stays text-only in L6.5.]
- Panel e: The syn3A tRNA operon (Thr/Val/Glu/Asn, MMSYN1_0678–0681) relocated by the flanking deletions to ~770 bp upstream of the rPtn operon; Illumina + ONT depth across the silent inter-operon gap show no read-through (the two operons stay independent). Broken minus-strand axis. (7/3, 7/3). [L6.3]

### Chain of Logics

#### L6.1: The ~418 deleted loci carried about a fifth of syn1's coding expression, freeing pool capacity.

- **Logic:** Quantifying the share of syn1's transcriptome and proteome contributed by loci absent from syn3A measures how much expression budget minimization freed for reallocation, and which RNA classes were lost.
- **Analysis:** `Genome_Reduction/09_Compare_RNA_Protein.py`
- **Outputs:**
  - `Compare_RNA_Protein/deleted_gene_occupancy.txt`
  - `Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv`
- **Numbers to cite:** 418 deleted loci (911 -> 496); by RNA type mRNA 382, pseudo 33, ncRNA 2, tRNA 1; deleted share = 21.78% of the syn1 mRNA pool, 22.25% of the iPM proteome; top deleted by TPM lacZ, pdhA/pdhB, ald; unclear-function proteins occupy only ~3%.
- **Figure panels:** a
- **Conclusion:** Minimization removed ~1/5 of the coding transcriptome and proteome, concentrated in dispensable metabolism, leaving pool capacity that syn3A redistributes.
- **Caveats:** shares are raw syn1 TPM/iPM; cross-organism comparisons in L6.2-L6.4 are mean-normalized and deletion-corrected to the retained-gene pool.
- **Notes for LLM:** "Non-essential" in the heading is shorthand for syn3A-deleted loci; Do NOT mention the deleted share of all-RNA TPM since rRNA were depleted during RNA sample preparation.

#### L6.2: The retained mRNA pool reallocates toward the translation machinery.

- **Logic:** After renormalizing to retained genes, per-function mRNA-pool shares are compared between organisms to test whether the freed capacity is taken up by ribosome/translation operons rather than spread evenly.
- **Analysis:** `Genome_Reduction/09_Compare_RNA_Protein.py`
- **Outputs:**
  - `Compare_RNA_Protein/TPM_change_by_{secondary,tertiary}.tsv`
  - `Compare_RNA_Protein/mRNA_pool_composition_by_secondary.pdf`, `tertiary_share_change_dumbbell.pdf`
- **Numbers to cite:** largest r-protein TPM absolute gains rpsK +6.04, rplO +5.04, infA +4.97, rplX +4.83, rplN +4.33; rpoA +2.69; per-category median FC + Mann-Whitney p from TPM_change_by_*.tsv.
- **Figure panels:** a, b
- **Conclusion:** The transcriptome shifts toward translation/ribosome biogenesis at the expense of central metabolism.
- **Caveats:** the shift is a pool-level reallocation, not uniform. Two routes make some r-proteins buck the up-trend: (i) the fusion-affected rpsT/S20 (0082) and rpsO/S15 (0294) crash because their deleted promoter was replaced by a weak fused one (gene_impact_class new_promoter_fusion; links to L5.3/L5.4); (ii) three structurally-INTACT r-proteins rpmF/L32 (0526), rpmE/L31 (0137), rpsU/S21 (0482) also drop sharply in transcript (TPM FC 0.126 / 0.010 / 0.044) yet are gene_impact_class **unaffected — NOT decapitated** (their operon structure is untouched), and their protein is in fact maintained or up (iPM FC 36.2 / 1.37 / 1.98). [CORRECTION 2026-06-05: the earlier "decapitated rpmF/rpmE/rpsU" claim was wrong — none are decapitation cases; rpsT/rpsO are the only deletion-hit r-proteins, both new_promoter_fusion.]
- **Notes for LLM:** The 11~kb rPtn operon now has its own paragraph (L6.3) — do NOT fold it back in here.

#### L6.3: A single 11 kb ribosomal-protein operon (OP_00341) triples its mRNA-pool share and dominates syn3A transcription, expressed from its own retained promoter despite a newly adjacent tRNA operon.

- **Logic:** The reallocation toward translation (L6.2) is concentrated in one polycistron — OP_00341 (MMSYN1_0652–0672) — so its mRNA-pool occupancy, transcript structure, and new syn3A genomic neighbourhood are examined to see how the dominant translation unit is expressed and whether the deletion that relocated an upstream tRNA operon couples the two.
- **Analysis:** `Genome_Reduction/09_Compare_RNA_Protein.py` (pool shares); `Genome_Reduction/R6_panel_e_trna_rptn.py` + `coexpression_common.py` (the tRNA-junction co-expression test, the 06/07 method).
- **Outputs:**
  - `Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv`
  - `R6_panels/R6_stats.txt`
- **Numbers to cite:** OP_00341 = 21 genes, ~11 kb (10,954 bp), minus strand, one polycistron with NO internal terminator; coding mRNA-pool share 12.1% (syn1) → 34.0% (syn3A), share FC 2.80, per-gene relTPM FC 1.48; full-length ~11 kb reads rare (1–2, PacBio read-length limit); depth = 5' polarity gradient (~90k at 5') with a sharp internal step at tx~2100 (endonucleolytic cut, likely RNase Y/degradosome which is up in syn3A, L6.4). New upstream neighbour after DEL_074 (5,509 bp; dhaK/0673–0676) + DEL_075 (912 bp; 0677): co-directional 4-tRNA operon MMSYN1_0678–0681 = Thr/Val/Glu/Asn; TSS(806176)→nearest deletion 179 bp (promoter intact); TSS→tRNA-3' 7,193 bp (syn1) → 772 bp (syn3A). Co-expression test rpsJ/0672 ↔ tRNA cluster: ONT 0/3084 spanning reads; Illumina true inter-operon middle (419784–420350) mean depth 27 = 1.2% of flanking → SPLIT (not co-transcribed).
- **Figure panels:** d, e
- **Conclusion:** The 11 kb r-protein operon, expressed from its own retained promoter as one endonucleolytically processed transcript, carries about a third of the syn3A coding mRNA pool; the deletion that parked a tRNA operon within ~770 bp upstream changed its neighbour but not its regulation, and the two operons stay transcriptionally independent — so the upregulation is the intact promoter plus the L6.2 pool reallocation, NOT tRNA read-through.
- **Caveats:** the ~11 kb full-length isoform is undersampled by PacBio read-length, so single-unit structure is inferred from continuous depth + no internal terminator, not from many full-span reads; the internal step is a depth/3'-end signature consistent with RNase Y cleavage, not a mapped cut site.
- **Notes for LLM:** This is the separate rPtn paragraph promised in L6.2. The tRNA-upstream result is a NEGATIVE co-expression finding — state plainly the operons stay independent (no read-through). Panels d (operon structure) and e (tRNA junction) both belong here.

#### L6.4: RNA polymerase is downregulated while the degradosome is upregulated in syn3A.

- **Logic:** Estimating gene-expression machinery RNAP and degradosome's assembled abundance from its limiting (lowest-stoichiometry) subunit compares transcription capacity against RNA-turnover capacity; opposite movement is coherent with syn3A's longer cell cycle.
- **Analysis:** `Genome_Reduction/10_Compare_Ptn.py`
- **Outputs:**
  - `Compare_RNA_Protein/macromolecule_complex_abundance.tsv`
  - `Compare_RNA_Protein/PTR_TPMfc_vs_iPMfc.pdf`, `PTR_by_category_boxplot.pdf`
- **Numbers to cite:** RNAP MIN(rpoA/2, rpoC, rpoB) TPM FC 0.65 (~35% down), iPM FC 0.79 (~21% down); degradosome MIN(rny, rnjA, yhaM+rnr) TPM FC 1.68 (~68% up), iPM FC 1.36 (~36% up); syn3A cell cycle 105 vs 60 min.
- **Figure panels:** c
- **Conclusion:** Lower transcription capacity plus higher RNA turnover is coherent with slower growth.
- **Caveats:** limiting-subunit estimate; PTR is a steady-state proxy, not Ribo-seq TE; r-proteins excluded from PTR (digestion bias).
- **Notes for LLM:** Degradosome is RNase Y, the endo-ribonucleases based, not RNase-J based; No need to mention the caveates at the end of this logic.

#### L6.5: Central-metabolism enzymes are coordinately downgraded, predicting lower ATP/GTP generation.

- **Logic:** Glycolytic and energy-generating enzymes are tracked at both RNA and protein level; coordinated downgrade across the pathway predicts reduced ATP/GTP flux, to be confirmed against a metabolic-flux comparison.
- **Analysis:**
  - RNA/protein evidence: `09_Compare_RNA_Protein.py` / `10_Compare_Ptn.py`
  - flux comparison: TBD (not in the 01-10 pipeline)
- **Outputs:**
  - `Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv` (enzyme FCs)
  - flux output: TBD
- **Numbers to cite (omics, TPM_FC / iPM_FC, both down):** gapDH 0.44/0.63, eno 0.37/0.38, pgk 0.35/0.44, pdhC 0.45/0.34, ackA 0.19/0.44, pta 0.27/0.26, ldh 0.66/0.70, pyk -/0.34; flux numbers TBD.
- **Figure panels:** (none — flux is text-only; the metabolic-flux panel was dropped and panel d reassigned to the giant rPtn operon, see L6.2)
- **Conclusion:** Central-carbon and acetate-pathway enzymes drop in concert, predicting suppressed ATP/GTP output (flux quantification pending).
- **Caveats:** the flux claim is currently inferred from enzyme abundance only.
- **Notes for LLM:** FLUX ANALYSIS NOT YET DONE (needs a metabolic model) — prose stops at the prediction with NO panel. Panel d is now the giant ribosomal-protein operon OP_00341 (belongs to L6.2). Update gapDH→GapA in the numbers above if reused.

---

# METHODS

## M1 — Illumina short-read sequencing of syn1.0 transcriptome  
**Tex file:** `Manuscript/sections/methods/illumina_syn1.tex`  
**Analysis:** 
- Illumina mapping: `Syn1_Transcriptomics/Illumina/Illumina_Processing/01_quality_control.sh`, `02_alignment_seqdepth.sh`  
- Gene TPM: `Syn1_Transcriptomics/Gene_TPM/Gene_Transcriptomics.py`
**Key params to mention:** 3 SRA datasets (SRR35996296/297 = technical reps of one RNA sample, SRR35996298 = second biological sample), 2x251 nt MiSeq, Kapa Hyper Stranded mRNA kit (dUTP / fr-firststrand, R2 = transcript strand); FastQC + MultiQC, no trimming; bowtie2 v2.5.5 default paired-end (99.49-99.56% overall alignment); samtools v1.22.1; per-strand bedGraph; depth-based TPM = length-normalized mean depth / (sense+antisense total) x 1e6, two-step replicate averaging (tech reps r=0.98; bio samples r=0.92-0.94).  
**Inputs:** SRA accessions in `Syn1_Transcriptomics/Illumina/Illumina_Raw/00_retrive_fastq.sh` 
**Outputs:** `Illumina_Processing/depth_bedgraph/SRR3599629{6,7,8}.{plus,minus}.bedGraph`; `Gene_TPM/syn1_illumina_TPM_profiles.tsv`, `syn1_Illumina_PacBio_TPM_profiles.csv`  
**Notes for LLM:** DONE — both subsubsections written (Illumina MiSeq read processing/mapping; TPM from sequencing depth). Stale commented NextFlow/BWA-MEM block removed (actual pipeline is bowtie2). RNA prep, QC, and library-prep subsubsections were already drafted. 

---

## M2 — PacBio long-read sequencing of syn1.0 transcriptome  
**Tex file:** `Manuscript/sections/methods/pacbio_syn1.tex`  
**Analysis:** `Syn1_Transcriptomics/PacBio/PacBio_Processing/`  
**Key params:** 3 technical reps (SRR36012641/642/643) pooled to 2.95 M HiFi reads; custom FLNC recovery (reorientation via H1/BCRC, primer trim, polyA trim) -> 2.62 M; `minimap2 -ax map-hifi --secondary=no` v2.30 (CP002027.1); pysam per-read HQ filter (MAPQ>=20, aln-frac>=0.7, clip<=0.3, |qlen-refspan|<=100, concatemer flag) -> 99.3% retained (2.60 M HQ); samtools v1.22.1 per-strand bedGraph; clustering (`Cluster_Isoform.py`, complete-linkage Chebyshev eps=10 bp): 621 k tuples -> 267 k clusters; depth-based TPM as in M1 (single pooled library, no replicate averaging).  
**Outputs:** `PacBio_Processing/syn1.PacBio.FLNC.sorted.HQ.bam`, `depth_bedgraph/syn1.PacBio.FLNC.HQ.{plus,minus,total}.bedGraph`; `Isoforms_PacBio/isoform_clusters_annotated.tsv`; `Gene_TPM/syn1_pacbio_TPM_profiles.tsv`  
**Inputs:** SRA accessions in `Syn1_Transcriptomics/PacBio/PacBio_Raw/00_retrieve_fastq.sh`  
**Notes for LLM:** DONE — polished. Fixed wrong/inconsistent numbers (retention 99.6->99.3%, min MAPQ 36->20, 612->621 k tuples, n>=10 discard 990 k/38%->874 k/34%, n>=50 1.4->1.5 M/+29->+22%), typos (Quanlity, missing space), added minimap2+samtools citations, and wrote the empty TPM subsubsection. Left untouched: sub-percent before/after-filter rounding (aligned-fraction 0.998/0.994, soft-clip 0.2%/0.6%, read length 3.07/3.93 kb) — negligible, flag if you want them aligned to the after-filter log.  

---

## M3 — Operon identification from PacBio long-read transcriptomics in syn1.0
**Tex file:** `Manuscript/sections/methods/operon_analysis.tex`  
**Analysis:** `Syn1_Transcriptomics/Isoforms_PacBio/Cluster_Isoform.py`, `Syn1_Operon/…`  
**Key params:** clustering thresholds, min reads, TSS/TTS calling rule.  
**Outputs:** `isoform_clusters_annotated.tsv`, `operons.candidate_blocks.tsv`  
**Notes for LLM:** Subsubsection **Locate Transcription Promoter and Terminator Signatures** finished.  

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
**Key params:** Spectronaut iBAQ -> iPM (iPM_i = 1e6 * iBAQ_i / sum iBAQ_j) per rep, mean across 3 reps; absolute copy number = (iPM/1e6) x total proteins/cell (syn1 ~127 k from dry mass 12.8 fg x 58.2% protein / avg MW); localization via DeepTMHMM (TMRs) + SignalP 6, priority signal-peptide > TMR > cytoplasmic; 2019 vs 2026 measurements; tertiary function annotation built by `report_annotation_stats_syn3A.py`.  
**Numbers to cite:** syn1 detected 721/828 (87.1%); median copy number cytoplasmic 47, lipoprotein 21, membrane 10, extracellular 3 (n = 516/68/126/11).  
**Notes for LLM:** DONE — wrote **Relative Protein Quantification** (iBAQ->iPM) and **Localization of Proteome** (DeepTMHMM+SignalP, cytoplasmic/membrane/lipoprotein/extracellular), polished **Absolute Intracellular Protein Quantification** (fixed detected 735->721, "weght"->"weight", removed empty "(see )" ref). Added DeepTMHMM (`hallgren_deeptmhmm_2022`) + SignalP (`teufel_signalp_2022`) citations. Sample-prep / LC-MS-MS subsubsections were already drafted. Verified syn3A 2026 detection = 446/455 (the old 449 was wrong); per author, the syn3A 2019-vs-2026 comparison paragraph (measured-vs-reused) was deleted from the tex (not discussed for now).

---

## M6 — Correlation between transcriptome and proteome
**Tex file:** `Manuscript/sections/methods/corr_transcriptome_proteome.tex`  
**Analysis:** `Transcription_Translation.py` (base correlation), `Translation_Residual_L1_initiation.py` (TIR/OSTIR), `Translation_Residual_L2_elongation.py` (CAI), `Translation_Residual_L3_degradation.py` + `Genomes_Input/Homology_Build.py` (degradation).  
**Key params:** Illumina sense TPM vs iPM, log10 OLS, residual = log10(iPM) − fit; TIR via OSTIR (anti-SD ACCUCCUUU, 30-nt windows, read-weighted); CAI (Sharp & Li, ref set = top 20% by iPM); Mpn half-lives (Burgos 2020) via reciprocal-best-hit blastp, protease-abundance correction (Lon 0.84 / FtsH 2.08, Mpn from Maier 2011).  
**Outputs:** `syn1_genes_transcriptomics_proteomics.csv`, `residual_analysis/`  
**Numbers to cite:** correlation all r=0.61 (R²=0.38, n=717) / cytoplasmic r=0.70 (R²=0.49, n=512); ΔR²: TIR +0.020 (2%), CAI +0.080 (+21%), degradation ≤0.01; shortest Mpn-mapped half-life 4.7 h (median 32 h) vs ~1 h doubling.  
**Notes for LLM:** DONE — intro (base correlation + residual definition) + 3 subsubsections written (TIR, CAI, degradation). Per author, **tAI left out** of Methods (sensitivity check only). Citations added: OSTIR `roots_ostir_2021`, CAI `sharp_codon_1987`, Burgos `burgos_protein_2020`, Maier `maier_quantification_2011`. CAI is the one layer that improved the fit; TIR and degradation explained little.  

---

## M7 — Novel Transcription and Translation
**Tex file:** `Manuscript/sections/methods/novel_orf.tex`  
**Analysis:** `Syn1_Novel_ORF/Abnormal_Transcripts.py` (antisense classes), `Syn1_Novel_ORF/Novel_translation.ipynb` (novel ORFs).  
**Key params:** antisense labeling base-by-base vs gene model, 3 classes (spurious-promoter / read-through / embedded); OSTIR all-start-codon scan (anti-SD ACCUCCUUU, genetic code 4 / UGA=Trp, ORFs >=15 aa); synthesis-rate rank (reads x TIR), top 100, in-silico trypsin (1 missed cleavage, 7-25 aa), Spectronaut re-search of augmented DB.  
**Numbers to cite:** 1.4% antisense (267 isoforms -> 89 clusters); spurious-promoter 59 (66%), read-through 30 (34%, incl. 4 embedded); 837 abnormal isoforms -> ~29,000 candidate ORFs -> top 100 -> 48 unique -> 47 with proteotypic peptides; 2 ORFs confirmed by MS (near MMSYN1_0592 [revised NOVEL_PEP_002] and the 54-aa N-term extension of MMSYN1_0768 [revised NOVEL_PEP_043, = NOVEL_PEP_030 in the old-cluster MS search Excel]), both deleted in syn3A. NOTE: the Spectronaut MS search ran on the old-cluster candidate DB; both confirmed peptides are reproduced in the revised top-100.  
**Notes for LLM:** DONE (polished + renumbered on revised clusters). Added antisense subsubsection (L4.1, 3 classes) + polished novel-ORF subsubsection. anti-SD ACCTCCTTT->ACCUCCUUU. RECOUNTED on the revised (post Apr-22) clusters via Novel_translation.py: 837 abnormal isoforms -> 29,443 candidate ORFs -> top-100 -> 48 unique / 47 proteotypic (no longer "~26,000 / 44"). The "2 confirmed in MS" is the R4-L4.6 result (kept out of Methods); both peptides reproduced in the revised top-100 (NOVEL_PEP_002; NOVEL_PEP_043 = old 030). NOTE: the Spectronaut search itself ran on the old-cluster DB.  

---

## M8 — Oxford Nanopore (ONT) and Illumina sequencing of syn3A transcriptome
**Tex file:** `Manuscript/sections/methods/ont_illumina_syn3A.tex`  
**Analysis:** `Syn3A_Transcriptomics/ONT/ONT_Processing/` (ONT) + `Syn3A_Transcriptomics/Illumina/Illumina_Processing/` (Illumina)  
**Key params:** ONT direct-RNA, `minimap2 -ax map-ont` (NOT splice, bacteria are intron-less), per-strand depth; Illumina syn3A paired-end bowtie2 (dUTP / fr-firststrand), per-strand bedGraph.  
**Inputs:** ONT raw `Syn3A_Transcriptomics/ONT/ONT_Raw/`; Illumina syn3A SRA accessions (SRR19432056/57 mate pair) via `Syn3A_Transcriptomics/Illumina/Illumina_Raw/00_retrive_fastq.sh`.  
**Notes for LLM:** DONE (polished). Seven subsubsections: ONT RNA sample prep, RNA QC, ONT library prep + direct sequencing, raw-read processing/base-calling, ONT→syn3A mapping, Illumina→syn3A mapping, and "Quantification of TPM from Sequencing Depth" (added; mirrors M1/M2). TPM subsubsection: depth-based TPM (same definition as M1) computed for both Illumina + ONT per-strand bedGraphs over 496 loci (493 gene + 3 pseudogene, both feature types parsed); single replicate each → no averaging; output `syn3A_TPM_Illumina_ONT.tsv`; Illumina = quantitative track, ONT = orthogonal check. RE-VERIFIED (Syn3A_TPM.py rerun): our depth-based Illumina TPM validated against Palsson/Sandberg-reported Illumina TPM at Pearson r=0.998 / Spearman 0.998 on log10 (n=458); the Illumina-vs-ONT cross-platform agreement is only Pearson r=0.570 / Spearman 0.558 (n=496). Do NOT report 0.998 as the Illumina/ONT correlation. (script: `Syn3A_Transcriptomics/Gene_TPM/Syn3A_TPM.py`; 496 loci re-counted from the GFF3.) ONT mapping numbers re-verified against `syn3A.ONT.rep1.sorted.bam.qc_report.txt`: 734.08k reads → 559k (76.2%) primary mapped, mean 383 nt (49--2,858), Q31.2, 98.2% ≥Q20; 175k unmapped (23.9%, mean Q20.7, 42.0% <Q15, 61.0% <300nt); 20.2k secondary (3.6%) at the two rRNA operons (~55,460 & ~343,267 bp); 7.3k supplementary. minimap2 v2.30 `-ax map-ont -p 0.99 --MD` (NOT splice; old splice-preset numbers kept commented). Illumina syn3A: bowtie2 v2.5.5 default paired-end, 98.88% overall (83.55% concordant-once, 10.43% multi = rRNA operons, 6.02% discordant/unpaired), dUTP/fr-firststrand strand split — alignment % not locally re-verifiable (no logs/ in OneDrive copy; trusted from prior run). Polish pass: units → mL/$\mu$L, centrifugal force → `$\times$g` (matches M1), removed a double space. Compiles clean (24 pp).  

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

**Notes for LLM:** DONE — polished. This is the longest Methods subsection, drafted in the same order as the pipeline; More details were recorded in the CLAUDE.md file.

---