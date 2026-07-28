# MANUSCRIPT.md — drafting guide for Methods & Results

This file is the **single source of truth** the LLM uses to draft Methods and
Results paragraphs. Each subsection here mirrors a `.tex` file under
`Manuscript/sections/{results,methods}/`. Fill in the bullets; the LLM will
expand them into prose, quoting numbers verbatim from the files you list.

Now Extended/SI displays:
- FigS1: RNA seq comparison
- TableS1: statistics on RNA seqs
- TableS2: Rnases
- FigS2: corr between mRNA and protein
- TableS3: parameters on calculating absolute protein copy numbers
- FigS3: intergenic transcription between 0884/0885 and novel peptide for 0768

## Manuscript Revision

Now i have discussed the entire manuscript in details with PI; following is the final TODO list before sending to all co-authors. I will ask you to revise some analysis and redraw panels; the tex files now will be edited on overleaf.

### ⬜ Open / remaining

Questions from Troy:

beatify the operon drawing notebook; learn from the David 2022 format

Things for the experimentalists:
Quality report of ONT: second Syn1 run

Try to visualize the counts of isoform arrows.

Discussion:
**Write as flowing prose that FOLLOWS the three Intro questions (introduction.tex:14-17), NOT a literal numbered Q&A. Target ~500 words (style budget). Do NOT re-report that transcript is the strongest determinant of protein abundance — that is R3's result; cut the repetition (current draft says it twice, discussion.tex:5 + :11).**

Opening (1 sentence): full-length sequencing gave an operon-resolved transcriptome of syn1 and matched transcriptome/proteome in both cells; set up the three questions below.

Q1 — How is gene expression organized across the synthetic genome?
- PacBio resolved **459 operons** in syn1, each carrying matched promoter/terminator signatures; the transcriptome is richer than the reduced gene count implies.
- Primary transcripts are pervasively reshaped by RNA processing, with markedly **more 3' than 5' erosion**.
- RETAIN (item 1 — the interpretation, not just the observation): the 3'-bias is a **3'->5' exonucleolytic CLEARANCE bottleneck** (scarce RNase R that reads through structure; YhaM stalls at stems; ribosome trapping on non-stop ends), not biased endonucleolytic cutting — a hypothesis that direct cleavage-site mapping will test (discussion.tex:9). (You can mention the interpretation, but it is just a hypothesis so keep it short)

Q2 — What did minimization do to the transcription AND translation program of the retained genes? (ALL observations live here now)
- Removed synthetic/mis-annotation baggage first: antisense + intergenic transcription traced to inherited *M. mycoides* mis-annotation and synthetic-construction artifacts, and those regions were among the first excised in syn3A.
- Structurally conservative (gene order preserved, whole operons excised) yet functionally transformative: a few retained operons were **decapitated** when their own promoter was deleted — clearest = **HupA** (promoter inside gpsA/0349, deleted -> HupA mRNA collapses 7.6x -> 0.13x avg depth).
- **MOVED HERE from Q3 (reallocation observations):** the retained transcriptome shifted toward the translation machinery —
   - the single **11 kb ribosomal-protein operon tripled its mRNA-pool share**;
   - **RNA polymerase and central-carbon / central-metabolic enzymes declined**;
   - the **RNA degradosome rose**.
- RETAIN (item 2 — two mechanistic hypotheses for the rProtein-operon rise): (a) the **tRNA operon** that minimization brought adjacent may raise local RNA-polymerase availability without read-through; (b) **HupA loss** may have shifted genome-wide supercoiling toward strong promoters (discussion.tex:21).

Q3 — How does minimization contribute to the phenotype differences? (INTERPRET the Q2 changes as causes)
1. **Longer cell cycle** (Syn3A ~2 h vs Syn1 ~1 h; pelletier_genetic_2021): lower transcription + lower energy output (RNAP + central metabolism down) coupled to faster RNA turnover (degradosome up) is coherent with slower growth. *Tentative add-on you flagged "will check":* under laboratory adaptive evolution Syn3A raises rProteins further and grows faster (ribosome content sets growth rate; sandberg_adaptive_2023), and prunes synthetic baggage (tetM/0913 silenced, echoing the his3/0918 over-transcription already lost) — reuse the drafted-but-commented paragraph (discussion.tex:23-25) ONLY if supportable.
2. **Denser cytosol** (gilbert_generating_2021): the giant rProtein operon UP together with several other rProtein genes DOWN imbalances rProtein production -> disrupts ribosome assembly -> surplus subunits crowd the cytosol.
3. **HupA** (as noted in Q2): its collapse accounts for Syn3A losing the persistent chromosome contacts Syn1 maintains (gilbert_dynamics_2023, gilbert_generating_2021).

Limitations:
- Direct ONT gave lower coverage and shorter reads than PacBio cDNA Iso-Seq, plausibly because pervasive 3' RNA processing degrades these transcripts — direct-RNA of degradation-prone messages would benefit from a dedicated RNA-protection protocol.
- **No PacBio library for Syn3A**: syn3A operon-level changes are inferred from short ONT/Illumina reads and await PacBio confirmation (Methods).
- Expression is framed as mRNA-pool share and excludes rRNA.
- Ribosomal- and membrane-protein abundances are the least reliable; **ribosome profiling** could sharpen them and test whether the elevated rProtein transcripts yield proportionally more ribosomes.

Implications:
- RETAIN (item 3 — cross-organism framing): reallocating resources between transcription and translation echoes the cellular-economy balancing seen when other genomes are streamlined (e.g. *B. subtilis*, michalik_bacillus_2021) and when synthetic cells are built from scratch — a general principle, not a syn3A quirk.
- Delete this since repeat with item 3: Emphasize that **allocation of transcription/translation resources** is a key consideration for anyone streamlining or designing minimal genomes.
- This map of gene co-expression + transcript quantification establishes the RNA-level **foundation for a 4D whole-cell model (4DWCM)** of the minimal cell (thornburg_bringing_2026).



### ✅ Finished

General directive (applied to all figures): Always denote syn1 or syn3A in the operon plot; normalize depth to avg depth in each organism.

- Fig 4 SI panel (Fig S3, intergenic 0884/0885) [DONE 2026-07-07]: `Syn1_Novel_ORF/FigS3_intergenic_0884_0885.py` -> `R4_panels/FigS3_intergenic_0884_0885.pdf` (-> `figures/si-intergenic.pdf`). **Combined Syn1-vs-Syn3A**, absolute syn1 coords: gene track (0884+ / 0885- / 0886-) + 6 mirror depth tracks (Syn1 PacBio / ONT1 / ONT2 / Illumina-avg, then Syn3A ONT / Illumina mapped through the retained block onto the syn1 axis; deleted block = blank on syn3A tracks). + up (blue) / - down (orange), each ÷ its own both-strand genome-mean. Deletion band shaded red ("deleted in Syn3A"), the ~180 bp intergenic gap shaded yellow. **INTERPRETATION (author, corrected from my "decapitation"):** a discrete ~180 bp + transcript sits in the 0884/0885 gap NOT overlapping either gene (0884 body is bare; 0884 = "conserved hypothetical protein"); a candidate UNANNOTATED transcription unit. Abundant in Illumina (21.8x) + ONT1 (21.6x) but under-sampled by size-selected PacBio (0.15x, short transcript). syn3A: 0884 (hypothetical) deleted, but the intergenic unit RETAINED though ~8x lower (Illumina 2.8x, ONT 1.5x). Wired into `novel.tex` (sifigure Fig S3; reframed reference sentence). Compiles clean (50 pp). NOTE Illumina norm uses replicate-WEIGHTED-avg mean (~292, matches R5/R6 panels), not Table S1 pooled 766.

- Fig 5 panel d + Fig 6 panel b [DONE 2026-07-07]:
  - Fig 5 panel d (hupA, `Genome_Reduction/R5_figure_panels.py::panel_d`): resized to (7, 7/3), added a **Syn1 ONT isoform track** beneath the PacBio one (new `load_syn1_ont_plus_isoforms` reads `syn1.ONT.merged.sorted.bam` via samtools; + strand reads over hupA clustered by 5'/3' ends). ONT 5' ends pile up at the TSS inside 0349: 3,313 reads over hupA, 1,540 start at/upstream of the deletion junction. Confirms the promoter lies inside gpsA/0349 (deleted in syn3A; hupA depth 7.6x -> 0.13x). **5' quantification (PacBio vs ONT):** PacBio pins the TSS sharply -- 79% of 2,130 reads at a single 5' base (441031, = annotated TSS, anchors the -10 box 441019-441024); ONT 5' ends fall inside 0349 too but shifted ~15-30 nt downstream (modal +14, median +29 nt), the expected 5' digestion of direct-RNA. Results + caption say PacBio fixes the TSS sharply / ONT corroborates just downstream. `_junction_panel` extended with an optional `ont_iso` track (panel c unaffected). Results + caption in `reduction_operons.tex` updated (PacBio+ONT place the TSS inside 0349).
  - Fig 6 panel b (`R6_figure_panels.py::panel_b`): prepended "mRNA" to both axis labels -> "mRNA TPM fold change" / "mRNA TPM absolute change".

- Protein copy-number: syn3A Methods + Table S3 + Results diff [DONE 2026-07-07]:
  - Checked the dry masses: NOT identical. syn1 = 12.8 fg dry mass x 58.2% protein = 7.48 fg protein/cell -> ~127,800 proteins/cell; syn3A = 10.2 fg x 54.727% = 5.56 fg -> ~100,300 proteins/cell. EF-Tu(0151): syn1 7,173 vs syn3A 5,089 copies. Constants in `Syn1_Syn3A_Proteomics/Protein_Quantification_Localization.ipynb` cells 16 (syn1) + 35 (syn3A).
  - Methods (`proteomics_syn1_syn3A.tex`): added syn3A relative-detection sentence (446/455) + syn3A absolute-quantification paragraph. **Table S3** (`tbl:copynum-params`, all parameters incl. cell volume um^3 + protein density um^-3, mimics OLD `dry_mass_param.tex` via booktabs+minipage since threeparttable isn't loaded) now lives in `results/corr_RNA_ptn.tex` right after Fig S2 (moved out of Methods per author). Numbering unchanged: S1 seq-summary, S2 rnase, S3 copynum-params (R3 is after R2 in doc order; verified in .aux, table on p.14 right after Fig S2 p.13).
  - Results (`corr_RNA_ptn.tex`): added the copy-number difference sentence. **Direction fix**: syn3A median copy number (~66) is nearly DOUBLE syn1 (~31), because a comparable protein pool is spread over fewer genes (446 vs 721) even though syn1 has 1.27x more total protein molecules. (The TODO had the direction reversed.)
  - Also reconciled the now-stale Fig 3 / Fig S2 captions + in-text panel refs to the rebuilt figures. **Panel order (match in Illustrator):** main correlation.pdf a/d = copy-number dist, b/e = mRNA-iPM corr, c/f = half-life (a-c syn1, d-f syn3A); SI si-correlation a-d = syn1 TIR/CAI/half-life-residual/R, e-h = syn3A same. Table S3 also lists cell volume (um^3) + protein density (um^-3) per organism. Full manuscript compiles clean (49 pp, 0 undefined).

- **Fig 3 COMPLETE [DONE 2026-07-07].** Panel order finalized (a/d copy-number, b/e corr, c/f half-life). Legend polish: copy-number title compacted to "(median, n=unique)", corr+half-life legends -> 5 pt, localization labels shortened to Cyto/Lipo/Mem/Extra (defined in the caption). Proteome residual annotated in Illustrator (author). Remaining copy-number-assumption question moved to Deferred.
- Fig 3 rebuilt with the syn1(blue)/syn3A(red) organism convention [DONE 2026-07-07]:
  - Main figure = 6 panels (3 metrics x 2 organisms): mRNA-vs-iPM corr (7/3 sq; FILLED syn1 / OPEN syn3A circles), copy-number distribution by localization (7/3 x 7/6), half-life distribution (7/3 x 7/6; blue vs red bars). Every panel tagged with a bold blue "Syn1.0" / red "Syn3A" left-title (Fig 5/6 convention); localization palette SHARED so the two rows compare like-for-like. Corr xlabel -> "mRNA Illumina TPM (log10)".
  - SI figure (si-correlation) = 8 panels (4 metrics x 2 organisms): residual vs TIR / CAI / half-life + R-improvement; syn1 dots blue, syn3A dots red, black fit lines.
  - Dropped from R3: PacBio-vs-Illumina corr + length/abundance-bias panels (now live in R0 si-rnaseq). Also recoloured R0 Fig S1's gray scatter to blue(syn1)/red(syn3A) per author.
  - Scripts: `Syn1_Corr_RNA_Proteins/R3_figure_panels.py`, `Syn3A_Corr_RNA_Proteins/Corr_RNA_Protein_Syn3A.py`, `RNAseq_Comparison/FigS1_panels.py` (all rerun, same numbers). Combine correlation.pdf / si-correlation.pdf / si-rnaseq.pdf manually in Illustrator.

- Intro / R0 / Fig S1 / syn1 ONT [DONE 2026-07-06]:
  - Intro now names all three RNA-seq platforms (PacBio, ONT, Illumina); the detailed technique review moved to Methods `rnaseq_comparison.tex` (from OLD intro.tex L27-30; 11 technique refs ported into references.bib). "Add a few details in Intro" dropped (no longer needed).
  - New opening Results section `rna_seq.tex` (R0) reporting read length + per-strand coverage for all libraries, with Table S1 (sequencing summary).
  - Fig S1: `RNAseq_Comparison/{compute_platform_TPM.py,FigS1_panels.py}` -> 4 born-at-size panels in `FigS1_panels/` (a = syn1 4x4 TPM SPLOM; b/c = PacBio-vs-Illumina length/abundance bias, ONT dropped per author; d = syn3A ONT-vs-Illumina). Combine si-rnaseq.pdf manually in Illustrator; `sifigure` env fixed so this=Fig S1, R3 si-correlation=Fig S2. (This absorbs the old Fig 3 item "move Illumina/PacBio TPM corr panels to a new Fig S1".)
  - syn1 ONT redo: new `Syn1_Transcriptomics/ONT/` pipeline mirroring the Syn3A flow. Orientation verified (rep1 as-is 88.5% map; rep2 native 3'->5' reversed 78.1% map). Two runs differ a lot (per-gene sense TPM r=0.82, different rRNA kits) -> DO NOT average; merged BAM kept only as a browser track. Methods: `ont_syn1.tex` + combined `ont_syn1_syn3A.tex` + `illumina_syn3a.tex`.
  - Fig 2 panel g: SD strength added as digits.
  - Fig 6 panel b: up/down-regulated labels added; panel d: text expanded to 5.3x avg depth.

- R2: 
  - X Remove panel d since secondary prediction and terminator signature not consistent
  - X Add a panel as f to show the spatial organization of ATP synthase
  - X Delete the atpA mRNA secondary structure prediction
  - X DONE 2026-06-30 — panel g remade (R2_figure_panels.py -> R2_panels/R2g_atp_synthase.pdf): per-subunit gene colours + locusNum/subunit two-line labels, F0/F1 isoform colours, sqrt-of-reads line thickness, depth normalised to x-mean, RNase III cut removed, single 5'-sorted isoform stack; final size 7/3 (not 7/4).
  - X Done Compress panel g to size (7, 7/4); make depth track longer; denote thickness to isoform counts; color each gene/subunit; draw RNA isoforms in two colors, one for trans-membrane F0 and one for peripheral membrane F1; No internal stop in syn1 given no predicted terminator;
  - X Done RGB color for subunit: 0797 subunit black, a (61,132,181), c (174, 174, 174), b (202,112,199), delta (234, 52, 38), alpha (219, 120, 66), gamma (244, 193, 66), beta (158, 214, 126), epsilon (76, 124, 49) 
  - X DoneRGB for trans-membrane part (74, 124, 179) for peripheral part (0, 146, 69)
  - X Done remove legends of 5' and 3' block
  - X Done No internal promoter predicted from TransTermHP within this region; depth drop from delta to alpha can from RNAP deattach or RNase III digestion; drop from gamma to beta can be from
  - X Remove corresponding Methods section of finding dsRNA stem on atpA gene

- R3: The correlation analysis is sort of standard and we want to suppress this discussion even more: [DONE 2026-07-01]
  - X a. Put several panels to Supplemental figures; main correlation.pdf trimmed to 5 syn1 panels (a copynum, b TPMvsiPM, c PacBiovsIllumina, d R±CAI, e half-life dist); the residual/bias panels moved to new Figure S1 (`si-correlation.pdf`, panels a-l), inserted right after Fig 3 via the `sifigure` env in macros.tex (numbers as S1, restores main counter; verified compile).
  - X b. repeat the analysis on syn3A, and report the values; new `Syn3A_Corr_RNA_Proteins/Corr_RNA_Protein_Syn3A.py` (run in RNAseq env, needs ostir) -> `R3_panels_syn3A/panel_{f..l}.pdf` + `syn3A_genes_transcriptomics_proteomics.csv` + `R3_syn3A.txt`. **syn3A values:** TPM-iPM r=0.63 (all, n=446) / 0.65 (cyto, n=352); CAI dR2=+0.127 (0.39->0.52), r=0.46; TIR dR2=+0.015, r=0.15; R±CAI all 0.63->0.72, cyto 0.65->0.73; half-life median 5.7h (shortest 0.9h, 6.5% below the 105-min doubling). Method diffs: TIR gene-level OSTIR from scratch (no isoforms); CAI recomputed w/ syn3A top-20%-iPM ref set; half-life reused by locus suffix + re-scaled by syn3A Lon(518)/FtsH(260) & V=0.0335 fL (factors 0.16/1.78). syn1 bias panels (S1 a,b) added to `Syn1_Corr_RNA_Proteins/R3_figure_panels.py`. Methods: `corr_transcriptome_proteome.tex` "Repeating the analysis in Syn3A" subsection.
  - X c. discuss a little bit why R is not 1 as Abner emailed me; sentence drafted in `corr_RNA_ptn.tex` (finite measurement reproducibility + genuine elongation-driven offsets set the ceiling) — currently COMMENTED OUT by author.
  - NOTE: the old "no length/abundance bias" claim was softened to "weak dependence (r=0.22 / -0.23, both |r|<0.25)" — the S1 a,b panels now SHOW the weak bias, so the text and figure agree.
  - NOTE: panel S1j marks the 105-min doubling line + legend "shortest 0.9 h, 6.5% below". Results now states the shorter syn3A half-lives track the **Lon 5.4x / FtsH 1.2x higher concentration** in syn3A vs syn1 (concentration = copies/volume; syn1 Lon 128/FtsH 292 @ V=0.0446 fL, syn3A Lon 518/FtsH 260 @ V=0.0335 fL; computed in Corr_RNA_Protein_Syn3A.py).
  - NOTE: both R3 figures now included at **native (born-at-size) size** (`\includegraphics[scale=1]`), not `width=\linewidth`, to preserve the 7 pt panel fonts. To host the 7 in (177.8 mm) figures, `main.tex` geometry changed to `left=15mm, right=15mm` on A4 => **text block = 180 mm** (Nature's max figure width; top/bottom kept 1 in). Both figures now fit natively, centred, 0 overfull. (Nature provides no LaTeX template; A4 or Letter both fine for submission. Side effect: body text lines are now 180 mm wide, i.e. wider than before.)

- R4: [DONE 2026-06-30 — R4_track_panels.py; b/c untouched]
  - X panel a: matched the "Isoform span (5'->3')" label to the b/c xlabel font (real set_xlabel, 7pt) and lowered the case rows (ylim -0.2..2.9) so spurious/read-through/embedded align row-for-row with the b/c ridge baselines; isoform-span arrow moved into the bottom margin. (Earlier "move isoform/xlabels to top" idea dropped.)
  - X panel d: dotted (dashed + lightened) every gene antisense to the shown isoforms (his3/0918 AND 0917) via the shared draw_gene_track (also applies to g/h). REBUILT as a dedicated panel_d_his3(): 4 tracks (genes | Syn1 + isoforms | Syn1 + depth | Syn3A Illumina + depth, RED) on a deletion-junction RELATIVE x-axis = syn1_pos - 27638 (deletion end; his3/0918 positive, deleted upstream negative); junction syn1 27638 <-> syn3A 18715; xlabel "Relative genome position (bp)". Story: Syn1 antisense over-transcription starts at the spurious promoter (rel ~-116, inside the deleted region) and runs across his3; in Syn3A it collapses to background.
  - X depth normalization (Syn1 tracks d+g AND the Syn3A track): x TOTAL (plus+minus) genome-mean coverage, to match novel.tex's "genome-wide average" basis (Syn1 total 4184, Syn3A total 2355). his3 antisense = Syn1 ~8x vs Syn3A ~0.23x (== the 0.23x quoted in novel.tex L4.3); panel g now ~0.5x. NOTE for novel.tex: his3 Syn1 "depth exceeding 30,000" -> "~8x the genome mean".
  - X panels g, h xlabel -> "Syn1 Genome Position (kb)".
  
- R5: [figures DONE 2026-07-01 — layout now a/b/c/d; details in the R5 figure-list NOTE below]
  - X panel a: resize to 7/3, 7/3 (panel_a).
  - X Add the junction-reformation panel (14/3, 7/3): the DEL_014 fusion where one junction reforms the operons at both ends -> panel b.
  - X Merge panel b and c -> single panel b (panel_b): rpsT/0082 partner switch on the shared-0082 5' transcript axis (syn1 0083/0082 co-transcription + syn3A 0094/0082 fusion; the lone 2-read ONT isoform spanning 0094+0082 highlighted).
  - X panel e -> now panel d (panel_d, hupA): syn1 + syn3A depth tracks, renormalized to genome-mean coverage, Illumina both organisms, syn3A mapped through the retained blocks so deleted regions read as gaps.
  - X panel d (violin replacement): rather than a TPM-FC distribution, replaced with the two-decapitated-central-carbon-operons showcase -> panel c (pdh/acetate OP_00121 + PTS OP_00122; operon-spanning PacBio isoforms; log-y Illumina depth; highlights promoter_lost cases beyond hupA). The impact-class violin is kept as panel_impact -> SI.
  - X panels b/c/d x-axis unified to "Relative transcript position (nt)".
  - X prose reframed 2026-07-01: section retitled "Genome reduction decapitates the operons of key proteins"; leads with promoter-loss -> key-protein suppression (HupA + central-carbon in two dedicated paragraphs), fusion (rpsT/rpsO) framed as the rare, ineffective reciprocal. The operon-truncation / junction-taxonomy / gene-impact-class COUNTS were moved out of Results into Methods (`genome_reduction.tex`: 235/172/52 operon overlay; 53/19/15/8 + 30/11/9/3 junctions; 2/3 vs 3/30, baseline 27/45 co-transcription; 360/45/42/24/17/6/3 impact classes) + full tables in `\sdreduction`. rpsO TPM FC restored to 0.144 (protein unchanged/up, transcript-only).

- R6: [prose + panels DONE 2026-07-01/02]
  - X New story line [prose DONE 2026-07-01, see the R6 Figure NOTE below]: P1 removed half genes + 1/5 pool (panel a); P2 correlation still high but majority genes lower since translation up, esp. the 21-rPtn operon (new panel b); P3 RNAP/degradosome/central-carbon (panel c); P4 21-rPtn operon intact yet upstream now a tRNA operon (panels d, e).
  - X a: add title: mRNA pool share in Illustrator [Illustrator step — panel unchanged in code]
  - X new panel b: FC-vs-absChange (14/3, 7/2), black base dots, 51 rPtns green, correlation inset top-left (r=0.84, 66% below diagonal). `R6_figure_panels.py::panel_b`.
  - X panel c: flipped to landscape (7, 7/6) — vertical lollipop, entities on x, FC on log-y.
  - X panel c: legend "Transcript" -> "mRNA".
  - X naming DECIDED 2026-07-02: keep **relative TPM** (NOT "normalized TPM" — TPM is already normalized, so that name is ambiguous). Defined inline at first use in R6 P2 (reduction_omics.tex): "each gene's TPM relative to the average gene TPM, which counteracts the different gene numbers of the two cells"; relative iPM = its protein analogue. P1 normalization explanation shortened to a one-clause pointer (mean-normalized to the retained pool; Methods).
  - X recheck the normalization scheme (retained-pool mean-normalization) — DONE 2026-07-02; comparisons are mean-normalized to the retained pool, with the relative-TPM definition added inline (see the relative-TPM naming item above).
  - X d and e merged into one clean panel — DONE 2026-07-01: the rPtn operon's upstream-neighbour swap (dnaK -> tRNAs) with still no co-expression is now a single transcript-axis panel d.

## Deferred / future work

Intentionally deferred analyses; the manuscript already ships without them.

- **L1.4** Annotate promoter, terminator and RNA processing signatures for non-canonical operons; Mannual search for more examples of 3' and 5' digestions
- **L2.3 genome-wide 3'-end 2°-structure test.** Settled design (do NOT use the confounded intragenic-vs-intergenic comparison, since terminators are trivially structured): (i) terminal accessibility = mean base-pairing probability of the last ~5 nt (ViennaRNA partition function) at intragenic 3' ends vs dinucleotide-shuffled, composition-matched windows; (ii) a 3'-vs-5' mirror-asymmetry meta-profile of pairing probability aligned at the endpoint. N = 6,302 intragenic 3' / 9,992 5' unique positions (`Syn1_RNase/RNase/isoform_endpoint_context.tsv`); ViennaRNA 2.6.4 in the RNAseq env.
- **R5 L5.5 (panel e).** Essentiality × trace-expression: write the script and supply a syn3A essentiality source.
- **R6 L6.5 (no panel).** ATP/GTP flux comparison: needs a metabolic model.
- **R1 L1.5 (panel e).** Decide which polycistronic operon to showcase.
- **R4.** RNA-polymerase-conflict angle (Ju et al., Nat. Microbiol.).
- **R3 copy-number assumption consistency.** syn1 uses protein mass fraction 58.2% (Razin 1963) while syn3A uses 54.727% (Breuer 2019); syn3A gDW = 1.0161e-14 g has a commented-out alternative 5.4e-15 (~1.9x) in `Protein_Quantification_Localization.ipynb`. Optional: unify onto one basis. Dry masses are NOT identical, and Table S3 documents the values as used, so the manuscript ships without this.
- Make **deliverables** for following 4DWCM: omics in syn3A (Done); Genome Visualization of Operons (TSS, TTS), ORFs

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
**Referring to files:** cite the three Supplementary Data workbooks by macro (defined in `sections/macros.tex`), never by repo path: `\sdoperon` = Supplementary Data S1 (operon.xlsx), `\sdomics` = S2 (syn1_omics.xlsx), `\sdreduction` = S3 (genome_reduction.xlsx), `\sdqc` = S4 (the experimentalists' RNA-sample QC PDFs — Qubit + TapeStation — for the Syn1 PacBio/Illumina and Syn3A ONT/Illumina libraries; to be delivered as a single zip indexed by a README.txt). In-silico read-QC (FastQC/MultiQC) was DROPPED from Methods (those reports are not in the SI); the meaningful "no adapter/quality trimming" statement was kept. Repo working-filenames (`*.tsv/.bed/.bedGraph/.gff3`) and internal pipeline script names were stripped from Results/Methods; each section's "SI:" pointer is now bold "Supplementary Data S#." All three are described in `sections/SI.tex` (holds `\label{SI}`; the library QC plots fold in there too, so the old `\fig{qc-*}` refs were removed). `isoform_endpoint_context.tsv` is NOT an SI file (pointer removed). Build the actual `.xlsx` workbooks from the repo tables before submission.

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

**SI file:** `\sdoperon` = operon.xlsx (boundaries, member genes + coverage, TSS/TTS, promoter/terminator signatures, macromolecular-complex annotation) — cited in-text as "Supplementary Data S1". **Built by `Syn1_Operon/build_operon_xlsx.py`** (pure assembly, no recompute): joins `operons.candidate_blocks.tsv` + the promoter (`promoter_minus10_classification.tsv`) and terminator (`terminator_tts_classification.tsv`) signature tsvs + `protein_complexes.xlsx`. Promoter/terminator are filled for the **127 canonical operons only** (TSS+TTS intergenic); non-canonical rows blank. `Operon_Annotation.py` now also persists the per-operon terminator table (with stem/loop/poly-U geometry), symmetric with the promoter table — re-run it before rebuilding. Sheets: `Operons` (459 rows) + `Protein_complexes` (26).

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
- **Caveats:** the segmentation now applies `dedup_operon_gene_lists` so `sense_gene_count` equals the unique-loci count (max 21, 0 mismatches across all 459 operons); panel b / mean are no longer inflated.

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

#### L1.5: One instance of polycistronic operons

- **Logic:** The choice not decided yet: could be rPtn operons, or other complexes
- **Analysis:** None
- **Outputs:** 
- **Numbers to cite:** None
- **Figure panels:** e
- **Conclusion:** 
- **Caveats:** 

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
- Panel d: RNA isoform truncation categories. (7/4, 7/4) — `Syn1_RNase/R2_panels/R2e_truncation_categories.pdf`
- Panel e: Biased RNA Processing schematics: endo and exo from 3'. (7/2, 7/4) — Illustrator (no matplotlib file)
- Panel f: Subunit composition and spatial arrangement of ATP synthase - Illustrator
- Panel g: RNA isoform distributions for ATP synthase operon — isoforms split into two regions at atpA/α (0792) where RNase III cuts; coloured by 5'-block (a,c,b,δ; teal) vs 3'-block (γ,β,ε; orange), gene arrows tinted by block, depth steps down at the α cut. (7, 7/3) — `Syn1_RNase/R2_panels/R2g_atp_synthase.pdf` (matplotlib: gene arrows + isoforms + depth only; the F1/F0 scheme, SD strengths, subunit labels, "RNase III on α" scissors are added in Illustrator).

### Chain of Logics

#### L2.1: Distinct RNA isoforms distributions found for operons.

- **Logic:** Truncated isoforms compared to the full transcription units exist for operons because of the RNA processing. Distinct patterns of truncations can be found, using genes 0154 and 0178 as examples. 0178 has structured 3' end as shown in d, but RNase R can digest through the dsRNA structured region.
- **Analysis:** per-operon erosion-category composition from `Syn1_RNase/R2_panels/R2_panels.txt` (`R2_figure_panels.py`; operon-contained isoforms, n_reads≥2) — these are the figure-panel numbers, NOT the isoform_endpoint_context.tsv n≥10 set.
- **Outputs:** None
- **Numbers to cite:** **0178/OP_00099 (−, neopullulanase): 85 isoforms / 3,056 reads** — 4 unprocessed (2,001 reads, 65.5%), 49 iso 3′-eroded-only (29.5% reads), 27 iso 5′-eroded-only (4.6%); 5′/3′ eroded-read ratio 0.17. **lap/0154/OP_00078 (+, leucyl-aminopeptidase): 134 isoforms / 2,035 reads** — 10 unprocessed (1,011 reads, 49.7%), 92 iso 5′-eroded-only (41.8% reads), 29 iso 3′-eroded-only (8.2%); ratio 4.95. Opposite polarities.
- **Figure panels:** b,c,d
- **Conclusion:** None
- **Caveats:** None

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

#### L2.3 (HYPOTHESIS): The 3'-bias reflects an asymmetric exonucleolytic *clearance* bottleneck — limited 3'→5' read-through (scarce RNase R, stalling YhaM) and ribosome trapping on non-stop ends — not biased endonucleolytic cutting.

- **Logic (the chain):** Endonucleolytic cleavage (RNase Y/III) is *symmetric* — each cut yields one upstream fragment (intact 5' + a NEW intragenic 3' end) and one downstream fragment (NEW intragenic 5' end + intact 3' terminator) — so endo cutting alone cannot create a 3'/5' asymmetry; equal numbers of each are born. The observed bias must therefore arise from differential *clearance* of these fragments. (i) Downstream / 5'-eroded fragments are cleared rapidly by the abundant 5'→3' machinery (RNase J1 + J2) → low steady state → rarely captured. (ii) Upstream / 3'-eroded fragments can only be *fully* erased by RNase R — the lone 3'→5' exo that reads THROUGH structure (it needs an unstructured ss 3' overhang to load, then unwinds via intrinsic helicase activity); the abundant YhaM is a "generator not finisher" — it trims the ss 3' tail but STALLS at the first stem base, manufacturing 3' ends rather than removing them. So 3'-eroded intermediates are both over-produced (YhaM stalls) and under-cleared (RNase R scarce) → they accumulate and dominate the steady-state long-read pool. **Two reinforcing mechanisms:** (A) **lower 3'→5' read-through capacity** (RNase R 36 vs RNase J1+J2 ~234 copies; YhaM 117 stalls at structure); (B) **ribosome trapping** at non-stop 3' ends (the 42% start-but-no-stop ORFs from L2.2) physically blocks exo entry until tmRNA–SmpB rescue, which is limiting (SmpB 14).
- **Analysis:** the genome-wide structure test is deferred (design in Deferred / future work above); R2 ships with three worked examples + the kinetic argument. Three strands of support:
  1. **2° structure of the three examples** (user-built, ViennaRNA 2.6.4 in RNAseq env): the 5' erosion region of the lap operon (panel b), the 3' erosion region of the 0178 operon (panel c), and the atpA/α RNase III cleavage site (panel f) — qualitative illustration that the eroded 3' ends sit at accessible / stall-competent structures.
  2. **Ribonuclease capacity asymmetry (proteomics, in hand, Table S1):** 5'→3' (RNase J1+J2 ~234) vs 3'→5' read-through (RNase R 36); YhaM 117 (stalls), SmpB 14 → panel e proposed-hypothesis schematic.
  3. **Non-stop / ribosome-trapping link (in hand):** the L2.2 42%-no-stop ORFs × the limiting tmRNA–SmpB capacity.
- **Outputs:** Table S1 ribonuclease abundances (in hand); the three example 2° structures (user-built); genome-wide ViennaRNA test DEFERRED.
- **Numbers to cite:** 5'→3' RNase J2 142 + J1 92 = ~234 copies vs 3'→5' read-through RNase R **36**; YhaM 117 (Mn2+-dependent, stalls at stem base); SmpB 14; tmRNA 2.5% of non-rRNA. RNase R loads on a ss 3' overhang ≥7 nt (optimal ≥10). [structure-enrichment numbers PENDING]
- **Figure panels:** f (biased-processing schematic).
- **Conclusion:** The 3'-erosion bias is best explained by a 3'→5' clearance bottleneck (scarce read-through RNase R; abundant but structure-stalling YhaM) compounded by ribosome trapping on non-stop products under limiting trans-translation — not by biased endonucleolytic cutting.
- **Caveats:** copy number is a capacity *proxy*, not measured flux (RNase R is processive, so 36 copies aren't negligible); repeated endo cuts can substitute for RNase R in clearing structured RNA; the ribosome-trapping arm applies only to translated (CDS) 3' ends; long-read capture biases.

#### L2.4: ATP synthase operon is co-expressed in one-go but cut at $\alpha$ subunit.

- **Logic:** Macromolecular complexes' gene co-expression can be altered by RNA processing. ATP synthase's RNA isoform distribution has a clear pattern of isolation at $\alpha$ subunits, which was identified as endo RNase III cleavage site. Comment on the other membrane complexes.
- **Analysis:** `Syn1_RNase/R2_figure_panels.py` (panel_f).
- **Outputs:** 
  - `Syn1_RNase/R2_panels/R2g_atp_synthase.pdf`
  - `Syn1_RNase/RNase_Site_Mapping/output/rnaseIII/stems/R2_MMSYN1_0792_atpA_rnaseIII_stem.pdf` (local structure at the two homology-mapped α cuts; each at a stem, no duplex connector)
- **Numbers to cite:** the atp operon (minus strand) is segmented into two overlapping operons that meet AT atpA/$\alpha$ (MMSYN1_0792): the 5'-block OP_00395 (0797–atpH/0793 + 5' of $\alpha$; 12 member isoforms, top isoform 1,161 reads) and the 3'-block OP_00394 (3' of $\alpha$ + atpG/atpD/atpC = 0791–0789; 6 members); the minus-strand depth drops from ~9k to ~half across the junction. The cleavage is assigned to RNase III by homology to B. subtilis atpA (BSU_36830; 60.3% identity, reciprocal best hit). The two B. subtilis RNase III sites are transferred onto Syn1 by **homology mapping** (NOT a structurally confirmed duplex): the whole Syn1 gene is folded ONCE (RNAfold default), each B. subtilis cut is projected by transcript fraction, and the base-pairing is read from that single fold without re-folding; a duplex is "confirmed" only if the two cuts FACE across one helix (both cross-distances ≤ 4 nt — calibrated on B. subtilis atpA, whose two cuts pair ~105 nt apart yet face at cross-dist 3/3). **Mapped Syn1 cut sites = 932,767 / 932,881.** RESULT: **0/5 paired-site genes (incl. atpA) reproduce a clean cross-paired duplex** — Syn1 atpA's two cuts fold into SEPARATE local stem-loops (cross-dist 104/104), so the B. subtilis long stem is NOT structurally conserved. Each mapped cut still sits at a (local) stem, so the inset shows the local structure at the two mapped sites, not a staggered duplex. Source: `Syn1_RNase/RNase_Site_Mapping/output/rnaseIII/rnaseIII_syn1_predicted_cleavage_pairs.tsv`.
- **Figure panels:** g
- **Conclusion:** The RNase complexifies the subunit synthesis of complexes.
- **Caveats:** None

---

## R3 — High Correlation between transcriptome and proteome in the reduced organism.

**Tex file:** `Manuscript/sections/results/corr_RNA_ptn.tex`

**SI file (Supplementary Data S2, `\sdomics`):** `Syn1_Corr_RNA_Proteins/syn1_omics.xlsx` (all 911 genes; columns: locusTag, gene_name, rna_type, gene_product, protein_localization, TPM_illumina, TPM_PacBio, iPM_mean, protein_copy_number, TIR, CAI, protein_halflife_h; gene name/product taken from the syn3A proteome where an ortholog exists). Built by `Syn1_Corr_RNA_Proteins/R3_figure_panels.py`.

### One-sentence Summary
**High correlation found between transcriptome and proteome.**

### Figure

**Organism colour convention (Fig 5/6): Syn1.0 = blue `#3182bd`, Syn3A = red `#c0392b`.**
Every panel carries a bold left-title organism tag in that colour. The localization
palette (Okabe-Ito: cytoplasmic blue / lipoprotein green / membrane vermillion /
extracellular purple) is SHARED across both organisms so the two rows compare
like-for-like; organism is carried by the tag plus a second channel per panel
(corr = filled Syn1 / open Syn3A circles; half-life = blue vs red bars). Rebuilt
2026-07-07 by `Syn1_Corr_RNA_Proteins/R3_figure_panels.py` (syn1) +
`Syn3A_Corr_RNA_Proteins/Corr_RNA_Protein_Syn3A.py` (syn3A).

**Main figure** `Manuscript/figures/correlation.pdf` — 6 panels (3 metrics x 2 organisms):
- Syn1 mRNA TPM (log10, Illumina) vs proteome iPM (log10), FILLED circles by localization. (7/3, 7/3) — `R3_panels/panel_b_TPM_vs_iPM.pdf`
- Syn3A same, OPEN circles. (7/3, 7/3) — `R3_panels_syn3A/panel_g_TPM_vs_iPM.pdf`
- Syn1 copy-number distribution by localization. (7/3, 7/6) — `R3_panels/panel_a_copynumber_by_localization.pdf`
- Syn3A copy-number distribution by localization. (7/3, 7/6) — `R3_panels_syn3A/panel_f_copynumber_by_localization.pdf`
- Syn1 intrinsic half-life distribution (Mpn-transferred), blue bars. (7/3, 7/6) — `R3_panels/panel_g_halflife_distribution.pdf`
- Syn3A intrinsic half-life distribution (Mpn re-scaled), red bars. (7/3, 7/6) — `R3_panels_syn3A/panel_j_halflife_distribution.pdf`

**SI figure** `Manuscript/figures/si-correlation.pdf` — 8 panels (4 metrics x 2 organisms), each (7/4, 7/4); syn1 dots blue, syn3A dots red, black fit lines:
- Syn1 / Syn3A proteome residual vs TIR (log10) — `panel_d_TIR_vs_residual` / `panel_h_TIR_vs_residual`
- Syn1 / Syn3A proteome residual vs CAI — `panel_e_CAI_vs_residual` / `panel_i_CAI_vs_residual`
- Syn1 / Syn3A proteome residual vs half-life (log10) — `panel_h_halflife_vs_residual` / `panel_k_halflife_vs_residual`
- Syn1 / Syn3A Pearson R for whole proteome and cytoplasmic proteins, with/without CAI — `panel_f_R_improvement` / `panel_l_R_improvement`

The PacBio-vs-Illumina correlation + length/abundance-bias panels were dropped from R3 (now live in the R0 platform-comparison figure, `RNAseq_Comparison/`, si-rnaseq.pdf). The R0 Fig S1 gray scatter was also recoloured to the same blue(syn1)/red(syn3A) convention.

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

#### L3.2: Using Illumina TPM as standard of transcriptome quantification

- **Logic:** PacBio and Illumina TPMs of syn1 were correlated to get r of 0.62; no significant TPM and length bias was found.
- **Analysis:** `Syn1_Transcriptomics/Gene_TPM/Gene_Transcriptomics.py`
- **Outputs:** plots in `Syn1_Transcriptomics/Gene_TPM/`
- **Numbers to cite:**  PacBio vs Illumina sense TPM Pearson r=0.62 (log10, n=884 at TPM≥0.5; figure threshold low_threshold=0.5); no abundance/length bias
- **Figure panels:** c
- **Conclusion:** As convention, Illumina TPMs were used to do correlation.
- **Caveats:** Do NOT say "Illumina TPM used since correlate with iPM better."

#### L3.3: Decent correlation found between transcriptome and proteome for syn1

- **Logic:** Pearson r of 0.7 found between two omics for cytosolic proteins; lower r for all since poor coverage of membrane proteins.
- **Analysis:** `Syn1_Corr_RNA_Proteins/Transcription_Translation.py`
- **Outputs:** 
  - Same name Txt file
- **Numbers to cite:**  all proteins Pearson r=0.61 (R²=0.38, n=717; Spearman 0.67); cytosolic-only r=0.70 (R²=0.49, n=512; Spearman 0.75)
- **Figure panels:** b
- **Conclusion:** Decent correlation.
- **Caveats:** None

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
- **Numbers to cite:** 1.4% antisense (267 isoforms -> 89 clusters); spurious 59 (66%), read-through 30 (34%, incl. 4 embedded). **Cluster↔operon link:** 38/89 clusters transcribe an antisense gene enclosed by R1's 69 antisense-containing operons; R1's 9 sense-gene-less operons = 8 antisense-only (all match a spurious-promoter cluster, his3 the largest) + 1 purely intergenic (OP_00079 = the L4.5 isolated unit). Source: `Syn1_Novel_ORF/novel_tex_todos.py` → `novel_tex_todos.txt`.
- **Figure panels:** a
- **Conclusion:** Full-length RNA isoforms reveal new cases of anti-sense transcription as read-throughs.
- **Caveats:** None

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

#### L4.3: Unexpected transcription of yeast vector gene 0918.

- **Logic:** Yeast vector elements were carried into syn1's synthetic genome during assembly. Strikingly, the yeast selection marker his3/0918 was heavily transcribed antisense (depth >30k) but not translated. The antisense isoforms initiate from a spurious promoter just upstream. **his3/0918 is RETAINED in syn3A** (JCVISYN3A_0918; gene_impact_class context_only); its conspicuous antisense over-transcription collapses to background, the spurious upstream promoter evidently removed when an adjacent deletion truncated the operon's 3' end (the over-transcription, not the gene, is what minimization eliminates).
- **Analysis:** `Syn1_Novel_ORF/R4_track_panels.py` (panel d): isoform table + PacBio plus-strand depth + deletion overlay from `Genome_Reduction/aln/raw/syn1_deleted_regions.bed`. -10 box at the antisense TSS scored by `R4_track_panels.py:quantify_novel_promoters()` (same algorithm as canonical operons, via `Syn1_Operon/promoter_motif.py`).
- **Outputs:** `R4_panels/panel_d_his3_antisense.pdf`; `R4_panels/novel_promoter_minus10.txt`
- **Numbers to cite:** his3/0918 antisense depth >30,000; 28 antisense isoforms (top 16.6k reads); **syn3A Illumina at JCVISYN3A_0918 (18,716-19,378): antisense 0.23× / sense 0.28× genome mean — his3 RETAINED, the syn1 antisense standout gone** (`novel_tex_todos.txt`); the antisense TSS (pos5p0 27522) carries a perfect -10 hexamer TAAAAT (TANAAT consensus, 0 mismatch; core_6mer tier) with an AT-rich -35 (CTTTGAA), confirming a genuine spurious promoter. Watermarks W1-W4 (located by exact sequence): length-weighted mean PacBio depth ~283/+ , ~360/- vs genome-wide average ~2133/2051 (6-8x lower); covered ORFs all hypothetical/watermark calls (plus real 0590) -> minimally transcribed noise.
- **Figure panels:** d
- **Conclusion:** The yeast marker his3 is heavily transcribed antisense yet untranslated, and its driving spurious promoter (which carries a canonical TAAAAT -10 box, i.e. a real but mislocated sigma-factor promoter) is deleted in syn3A; the four watermarks are only minimally transcribed (noise).

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


#### L4.6: Two novel peptides identified by enumerating all possible ORFs in isoforms having high abnormal fraction.

- **Logic:** The possible translation of the abnormal RNA isoforms was checked by enumerating all ORFs using OSTIR on the abnormal RNA isoforms. 
- **Analysis:** 
  - `Syn1_Novel_ORF/Novel_translation.ipynb`
- **Outputs:** 
  - `Syn1_Operon/operons.candidate_blocks.tsv`
- **Numbers to cite:** 837 abnormal isoforms -> ~29,000 candidate ORFs -> top 100 -> 48 unique / 47 proteotypic; 2 MS-confirmed (NOVEL_PEP_002, 118 aa intergenic near 0592; NOVEL_PEP_043 = old 030, 225 aa, 54-aa N-term extension of 0768), both deleted in syn3A. **Context:** NOVEL_PEP_002 (728,399-728,756) flanked by the 0591 mmyCImod PSEUDOGENE + function-unknown 0592 cdsf; NOVEL_PEP_043 extends mmyCIVR/0768 within the mmyCIV restriction-modification cluster (`novel_tex_todos.py`)
- **Figure panels:** h
- **Conclusion:** Two predicted ORFs were identified in Mass-spec proteome, and both were located near less annotated genes. Also, these two regions were deleted in syn3A.
- **Caveats:** Only top 100 ORFs were selected to do the new proteomics search, thus we cannot assure if all ORFs were translated or not (leave this question to the reviewers); the new canonical cluster isoforms gave new ORF candidiates, which were highly similar to the old ones that searched against raw proteoimcs.

---

## R5 — Genome reduction decapitates the operons of key proteins
**Tex file:** `Manuscript/sections/results/reduction_operons.tex`
**Section title (current):** "Genome reduction decapitates the operons of key proteins" [retitled 2026-07-01; was "Operonal structure changes to the minimal cell, JCVI-syn3A"]

### One-sentence Summary
**Because the retained backbone is 99.90% identical, the deletions that most strongly suppress syn3A's retained proteins act regulatorily, by deleting an operon's own promoter (decapitation) — silencing functionally central genes such as the nucleoid protein HupA and the central-carbon enzymes — while the reciprocal event, fusing two operons into one weakly-driven unit, is rare and ineffective (rpsT, rpsO).**
<!-- Prose reframed 2026-07-01 to lead with promoter-loss -> key-protein suppression; the deletion/junction/impact-class taxonomy COUNTS were moved out of Results into Methods (genome_reduction.tex) + Supplementary Data S3. Old framing kept below for reference: "Halving the genome was a gene-order-preserving deletion campaign that excised whole operons, decapitated some retained operons by deleting their promoters, and fused only a small number of new cross-junction transcription units." -->

### Figure
**Figure:** `Manuscript/figures/genome_reduction.pdf`

- Panel a: Schematics of genome reduction from syn1 to syn3A. (7/3, 7/3)
- Panel b: 0083 and rpsT/0082 co-expressed in syn1 and 0094 and 0082 co-expressed in syn3A (14/3, 7/3)
- Panel c: Gene deletion disrupted expression of enzymes in central-carbon metabolism (7, 7/3)
- Panel d: the HupA operon, whose true promoter, located inside gene 0349, was deleted. (7, 7/4)
<!-- - Panel f: Gene essentiality evaluation for those trace-expressed genes that are still essential. -->

- NOTE [DONE 2026-07-01]: panels built by `Genome_Reduction/R5_figure_panels.py` (RNAseq env). a=panel_a; b=panel_b (rpsT/0082 partner switch — old syn1-b + syn3A-c merged onto the shared rpsT/0082 5' transcript axis; PacBio+ONT isoforms, the lone 2-read ONT isoform spanning 0094+0082 highlighted, Illumina depth both); c=panel_c (**two decapitated central-carbon operons**: pdh/acetate OP_00121 + PTS OP_00122; PacBio operon-spanning isoforms only, so each operon reads as a separate co-transcribed stack; **log-y Illumina depth** so the low PTS operon shows alongside the high pdh operon; per-region + per-gene avg depths exported to `R5_panels/R5_panel_stats.txt`); d=panel_d (hupA). The impact-class violin (former panel d) is kept as `panel_impact` for the SI. Depth = Illumina for both organisms (PacBio under-samples short genes, ONT is 3'-biased); syn3A depth mapped through retained blocks so deletions read as gaps. Shared 4-track junction plotter `_junction_panel(logy=)`. Stale panel PDFs (R5b_rpsT_operon_syn1, R5c_fusion_DEL014, R5d_pdh_operon, R5e_hupA_operon, R5bc_rpsT_fusion, R5d_TPM_FC_by_impact_class) are superseded.

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

#### L5.2: Deletions overlaid on syn1's 459 operons show whole-operon excision dominating over partial truncation.

- **Logic:** Intersecting the 95 deletions with the 459 syn1 operons at single-bp resolution classifies how each operon was hit, separating operons removed wholesale from those left partially truncated; the truncations are what create the junction effects in L5.3 and L5.4.
- **Analysis:** `Genome_Reduction/04_deletion_overlaid_operon.py`
- **Outputs:** `Genome_Reduction/deletion_overlaid_operon/operon_deletion_classification.tsv`
- **Numbers to cite:** span-level overlap_class (n=459): fully_deleted 181, intact 162, 3'_truncation_gene 47, 5'_truncation_gene 28, intra_truncated 17 (plus 9 UTR-only and 15 multi-hit); gene-level gene_deletion_pattern: all_deleted 235 (51.2%), intact 172 (37.5%), leading_deleted 21 (4.6%), lagging_deleted 20 (4.4%), intra_deleted 11 (2.4%); 414 syn1 genes overlapped by a deletion.
- **Figure panels:** a
- **Conclusion:** Reduction preferentially removed entire operons; the minority of partial truncations (5' vs 3') sets up the junction taxonomy.
- **Caveats:** the two axes (span-level truncation vs gene-level deletion) differ by design; 162 vs 172 "intact" reflects operons whose genes are all kept but whose UTR/flank was nicked.

#### L5.3: Same-strand deletion junctions can fuse new transcription units, but true fusion is rare.

- **Logic:** Each deletion is recast as a junction between the nearest retained operons on either side; relative orientation (tandem/convergent/divergent) and facing-regulator loss decide whether a new co-transcribed unit can form, and ONT spanning/bridging reads test whether the new cross-junction gene pair is actually co-transcribed.
- **Analysis:**
  - junction taxonomy: `Genome_Reduction/05_deletion_junction.py`
  - read validation: `06_single_operon_coexpression.py`, `07_operon_pair_coexpression.py`, `coexpression_common.py`
- **Outputs:**
  - `Genome_Reduction/deletion_junction/deletion_junctions.tsv`, `deletion_junction_summary.txt`
  - `Genome_Reduction/operon_pair_coexpression/`, `single_operon_coexpression/`
- **Numbers to cite:** 95 junctions: tandem 53, convergent 19, divergent 15, intra_operon 8; tandem junction_type: fusion 3, decapitation 9, readthrough_extension 11, clean_excision 30; cross-junction co-transcription (loose): fusion 67% (2/3) vs clean_excision 10% (3/30, negative control); pristine single-operon baseline preserved_loose 60% (45 testable, 111 pairs); fusion exemplar DEL_014 OP_00043 -> OP_00050 (MMSYN1_0094 -> MMSYN1_0082 = rpsT/S20), n_span=2, n_bridge=37, however the TPM FC was still low as 0.074 for rpsT/0082 since the fused promoter of 0094 is weak (in syn1, 0082 was co-transcribed with 0083 instead); a second r-protein rpsO/S15 (MMSYN1_0294) followed the same route (lost its own promoter, gained a weak fused one) and likewise dropped in transcript to **TPM FC 0.144** (relTPM 3.05 -> 0.44; its PROTEIN is unchanged/up, iPM FC 1.26 — so "collapse" is transcript-only). [CORRECTED 2026-07-01: earlier note said TPM FC 0.036, which was wrong; verified 0.144 from syn1_vs_syn3a_RNA_protein.tsv.] Both rpsT/0082 and rpsO/0294 are gene_impact_class new_promoter_fusion (from 08).
- **Figure panels:** b (rpsT/0082 partner switch; old syn1-b + syn3A-c merged onto the shared rpsT/0082 5' axis)
- **Conclusion:** Operon fusion is real but rare (3 events); the dominant junction outcome is clean excision of whole operon(s) between intact neighbors.
- **Caveats:** ONT depth is low, so most positive calls are loose-bridge rather than strict-spanning; convergent/divergent junctions are opposite-strand and not expected to co-transcribe.

#### L5.4: Decapitated operons that lost their own promoter are the one class that robustly drops in expression; HupA is the showcase.

- **Logic:** Classifying every retained gene by promoter-source change isolates operons whose own promoter was deleted (promoter_lost / decapitation); their syn3A TPM is compared against the other impact classes to test whether promoter loss, not sequence change, predicts lower expression.
- **Analysis:**
  - per-gene impact: `Genome_Reduction/08_delete_gene.py`
  - expression: `09_Compare_RNA_Protein.py`
- **Outputs:**
  - `Genome_Reduction/delete_gene/retained_gene_context.tsv` (`gene_impact_class` column)
  - `Genome_Reduction/Compare_RNA_Protein/TPM_FC_by_impact_class.pdf`
- **Numbers to cite:** gene_impact_class (retained genes): promoter_lost 42, promoter_disconnected 6, new_promoter_fusion 3, readthrough_exposed 24, promoter_proximity_changed 17, context_only 45, unaffected 360; promoter_lost is the only class robustly down in TPM (median FC 0.44, Mann-Whitney p=2.7e-4 vs unaffected median 0.76); HupA (MMSYN1_0350) relTPM 6.68 -> 0.13 (FC 0.020), relIPM 6.48 -> 0.092 (FC 0.014); HupA operon -10 box = perfect TANAAT (TATAAT), extended TNNTANAAT match, strong_9mer tier, -10 window 441019-441024 inside deleted DEL_050 (440092-441059) [from R5_panels/R5_panel_stats.txt via promoter_motif.scan_minus10]. NOTE: the weak-fused-promoter r-proteins rpsT/0082 and rpsO/0294 belong to L5.3 (new_promoter_fusion), NOT here. rpmE/L31 (0137) and rpsU/S21 (0482) are gene_impact_class unaffected (operon structure intact) with only mild TPM dips FC 0.198 and 0.413, so they are NOT decapitation cases either.
- **Central-carbon showcase (panel c):** two adjacent decapitated operons. OP_00121 (pdhC/0227-lpdA/0228-pta/0229-ackA/0230): DEL 288391-292905 removed promoter (TSS 291897) + PDH E1 subunits pdhA/0225 & pdhB/0226; region avg depth syn1 9.5x -> syn3A 3.0x (FC 0.32), 5'->3' per-gene gradient FC 0.43/0.35/0.27/0.18. OP_00122 (ptsI/0233-crr/0234-0235): DEL 298422-300803 removed promoter (TSS 300106) + 0231 & coaD/0232; region avg 1.4x -> 0.6x (FC 0.45). No PacBio isoform spans both operons (separate co-transcribed units). Full per-gene depths in `R5_panels/R5_panel_stats.txt`.
- **Figure panels:** c (two central-carbon operons), d (hupA). The impact-class violin (former panel d) moved to SI (`panel_impact`).
- **Conclusion:** Promoter-source loss drives the largest expression decreases; promoter_lost is the only impact class robustly down in TPM.
- **Caveats:** the class is assigned at operon level; 8 junctions lose only UTR (genes intact); the 05-vs-04 consistency check flags 2 flank operons as all_deleted.

#### L5.5: A few trace-expressed retained genes remain essential.

- **Logic:** Crossing gene essentiality against syn3A expression surfaces genes that are essential yet barely transcribed, i.e. retained through minimization despite minimal expression.
- **Analysis:** TBD (not produced by the 01-10 pipeline).
- **Outputs:** TBD
- **Numbers to cite:** TBD
- **Figure panels:** deferred — no panel in the current a-d layout (was the dropped essentiality panel).
- **Conclusion:** TBD
- **Caveats:** essentiality calls are inherited from the syn3A design literature, not measured here.

---

## R6 — Transcriptome and Proteome Changes to minimal cell, Syn3A
**Tex file:** `Manuscript/sections/results/reduction_omics.tex`

### One-sentence Summary
**More transcription on ribosomal protein operons suppresses the expression of enzymatic proteins in central metabolism.**

### Figure
**Figure:** `Manuscript/figures/reduction_omics.pdf`

- Panel a: mRNA pool compositions in syn1 and syn3A as secondary protein functions. [P1] (title "mRNA pool share" added in Illustrator)
- Panel b [NEW 2026-07-01, replaces the old tertiary-share dumbbell]: mRNA fold change (syn3A/syn1, x, log) vs absolute change (y, symlog) for the retained coding pool; base dots BLACK (alpha ramped by syn1 baseline), the 51 ribosomal proteins GREEN; INSET top-left = syn1-vs-syn3A relative-mRNA log-log correlation with y=x diagonal (r=0.84, 66% of genes below the diagonal = "majority shift down"). Size (14/3, 7/2). Built by `R6_figure_panels.py::panel_b`. [P2]
- Panel c: mRNA + protein fold change of RNAP, degradosome and central-carbon enzymes. FLIPPED 2026-07-01 from portrait (7/3, 7/2) to LANDSCAPE (7, 7/6): vertical lollipop, entities on x (rotated labels, family-coloured), fold change on log-y, ref line at 1; legend relabelled "Transcript" -> "mRNA". `R6_figure_panels.py::panel_c`. [P3]
- Panel d [MERGED 2026-07-01: old d (operon structure) + old e (tRNA junction) → one panel; old `R6_panel_e_trna_rptn.py`/`R6e_trna_rptn_syn3A.pdf` SUPERSEDED]: the 21-gene ~11 kb rPtn operon (rpsJ/0672→secY/0652, minus) + its swapped upstream neighbour, both cells on ONE transcript axis anchored on the shared operon 5′ end / TSS (806176; rel 0; rpsJ starts at ~+77 after the 5′ UTR; operon body positive, upstream negative). [origin changed rpsJ-5′ → operon-TSS 2026-07-01; title "21-gene ribosomal-protein operon" removed] 4 tracks top→down: syn1 genes (upstream dhaK/0673) | syn1 Illumina depth (×mean) | syn3A genes (upstream 4-tRNA operon 0678–0681, pulled to ~770 bp) | syn3A Illumina depth (×mean). Retained operon aligns 1:1; only the upstream neighbour changed; silent inter-operon gap shaded (no read-through). Every gene arrow labelled by gene name, alternating above/below the axis; the operon-body gray shade was removed (2026-07-01 per author). The promoter −10 box and the mRNA-pool-share text were computed and then TAKEN OFF the plot per author (2026-07-01) — both are EXPORT/prose only, not annotated on the figure. Size (7, 7/3). `R6_figure_panels.py::panel_d`. Exports to `R6_panels/R6de_rPtn_operon_depth.txt`: mean normalized operon-body depth (rel 0..10853) **syn1 5.36× vs syn3A 8.45× genome-mean (FC 1.58)**; coding mRNA-pool share **syn1 12.49% → syn3A 34.06% (share FC 2.73)**; retained promoter **−10 box TAGAAT** (canonical TANAAT, core_6mer, OP_00341 TSS 806176, retained in both — DEL_074 starts 179 bp upstream). [P4]

- NOTE [PROSE RESTRUCTURE 2026-07-01]: `reduction_omics.tex` reordered into 4 paragraph-blocks (P3 runs as 2 paras): **P1** removal + retained-pool composition leans to translation (Transl 53.1→64.9%, CCM 17.0→10.7%; panel a); **P2** hierarchy conserved (relTPM r=0.84 n=443; relIPM r=0.87 n=423) yet 66% of retained genes fell — freed capacity funnelled into ribosomal proteins (rpsK/rplO/rplX/rplN) + the one 21-gene operon (~1/3 of coding pool); selective exceptions rpsT/rpsO collapse + intact-but-down rpmF/rpmE/rpsU (panel b); **P3** RNAP↓ (0.65/0.79) vs degradosome↑ (1.68/1.36) coherent w/ 105-vs-60-min cell cycle, then the central-carbon enzyme cascade (glycolysis/pdh/acetate FCs) + ATP/GTP prediction (panel c); **P4** the 11 kb operon dominates from its own retained promoter, upstream neighbour swapped dhaK→tRNA operon, still no co-transcription (panel d, now the MERGED d+e). 66%-below-diagonal + r=0.84/0.87 now in `R6_panels/R6_stats.txt` ("[L6.1b / panel b]"); operon-body depth in `R6de_rPtn_operon_depth.txt`. Compiles clean. Deferred per author: fbaA/Cell-2022 glycolysis expansion. [d+e MERGED 2026-07-01 per author]

### Chain of Logics

#### L6.1: The ~418 deleted loci carried about a fifth of syn1's coding expression, freeing pool capacity.

- **Logic:** Quantifying the share of syn1's transcriptome and proteome contributed by loci absent from syn3A measures how much expression budget minimization freed for reallocation, and which RNA classes were lost.
- **Analysis:** `Genome_Reduction/09_Compare_RNA_Protein.py`
- **Outputs:**
  - `Compare_RNA_Protein/deleted_gene_occupancy.txt`
  - `Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv`
- **Numbers to cite:** 418 deleted loci (911 -> 496); by RNA type mRNA 382, pseudo 33, ncRNA 2, tRNA 1; deleted share = 21.78% of the syn1 mRNA pool, 22.25% of the iPM proteome; top deleted by TPM lacZ, pdhA/pdhB, ald; unclear-function proteins occupy only ~3%.
- **Figure panels:** a
- **Conclusion:** Minimization removed ~1/5 of the coding transcriptome and proteome, yet the deleted loci are nearly half (418/911) of the genes, so the removed genes were well below average in expression; the loss is concentrated in dispensable metabolism, leaving pool capacity that syn3A redistributes.
- **Caveats:** shares are raw syn1 TPM/iPM; cross-organism comparisons in L6.2-L6.4 are mean-normalized and deletion-corrected to the retained-gene pool.

#### L6.1b: The retained-gene expression hierarchy is strongly conserved across both layers, so the reallocation is concentrated in a few discrete movers.

- **Logic:** Before detailing the reallocation, the overall similarity of the retained-gene expression landscape between organisms is quantified, to show minimization preserved the ancestral program globally and that the changes below are concentrated exceptions (HupA down, the rPtn operon up), not a genome-wide drift.
- **Analysis:** `Genome_Reduction/09_Compare_RNA_Protein.py`
- **Outputs:**
  - `Compare_RNA_Protein/Compare_RNA_Protein.txt` (section "CROSS-ORGANISM CONSERVATION OF THE RETAINED-GENE EXPRESSION LANDSCAPE")
  - `Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv`
- **Numbers to cite:** transcriptome relTPM syn1 vs syn3A Pearson(log10) r=0.841 (Spearman 0.876), n=443 protein-coding (rna_type==mRNA); proteome relIPM syn1 vs syn3A Pearson(log10) r=0.867 (Spearman 0.890), n=423. Mean-normalized relative units, genes retained in both cells.
- **Figure panels:** none (text-only).
- **Conclusion:** The relative expression hierarchy of retained genes is largely conserved at both the RNA (r=0.84) and protein (r=0.87) layers; the reallocation in L6.2-L6.4 rides on this conserved backdrop as a small number of large, discrete shifts, chiefly HupA's collapse (L5.4) and the 11 kb r-protein operon's rise (L6.3).
- **Caveats:** the correlation is dominated by the bulk of unchanged genes; pool-composition shifts (L6.2) are driven by a minority of high-abundance loci, so high global conservation and concentrated reallocation are consistent.

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

#### L6.3: A single 11 kb ribosomal-protein operon (OP_00341) triples its mRNA-pool share and dominates syn3A transcription, expressed from its own retained promoter despite a newly adjacent tRNA operon.

- **Logic:** The reallocation toward translation (L6.2) is concentrated in one polycistron — OP_00341 (MMSYN1_0652–0672) — so its mRNA-pool occupancy, transcript structure, and new syn3A genomic neighbourhood are examined to see how the dominant translation unit is expressed and whether the deletion that relocated an upstream tRNA operon couples the two.
- **Analysis:** `Genome_Reduction/09_Compare_RNA_Protein.py` (pool shares); `Genome_Reduction/R6_panel_e_trna_rptn.py` + `coexpression_common.py` (the tRNA-junction co-expression test, the 06/07 method).
- **Outputs:**
  - `Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv`
  - `R6_panels/R6_stats.txt`
- **Numbers to cite:** OP_00341 = 21 genes, ~11 kb (10,954 bp), minus strand, one polycistron with NO internal terminator; coding mRNA-pool share 12.1% (syn1) → 34.0% (syn3A), share FC 2.80, per-gene relTPM FC 1.48; full-length ~11 kb reads rare (1–2, PacBio read-length limit); depth = 5' polarity gradient (~90k at 5') with a sharp internal step at tx~2100 (endonucleolytic cut, likely RNase Y/degradosome which is up in syn3A, L6.4). New upstream neighbour after DEL_074 (5,509 bp; dhaK/0673–0676) + DEL_075 (912 bp; 0677): co-directional 4-tRNA operon MMSYN1_0678–0681 = Thr/Val/Glu/Asn; TSS(806176)→nearest deletion 179 bp (promoter intact); **−10 box = TAGAAT** (canonical TANAAT, core_6mer tier, 0 mm; via promoter_motif.scan_minus10; retained in both cells); panel-d operon-body mean Illumina depth 5.36× (syn1) → 8.45× (syn3A) genome-mean, FC 1.58; TSS→tRNA-3' 7,193 bp (syn1) → 772 bp (syn3A). Co-expression test rpsJ/0672 ↔ tRNA cluster: ONT 0/3084 spanning reads; Illumina true inter-operon middle (419784–420350) mean depth 27 = 1.2% of flanking → SPLIT (not co-transcribed).
- **Figure panels:** d, e
- **Conclusion:** The 11 kb r-protein operon, expressed from its own retained promoter as one endonucleolytically processed transcript, carries about a third of the syn3A coding mRNA pool; the deletion that parked a tRNA operon within ~770 bp upstream changed its neighbour but not its regulation, and the two operons stay transcriptionally independent — so the upregulation is the intact promoter plus the L6.2 pool reallocation, NOT tRNA read-through.
- **Caveats:** the ~11 kb full-length isoform is undersampled by PacBio read-length, so single-unit structure is inferred from continuous depth + no internal terminator, not from many full-span reads; the internal step is a depth/3'-end signature consistent with RNase Y cleavage, not a mapped cut site.

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

---
# Introduction

**Tex file:** `Manuscript/sections/introduction.tex` — DRAFTED (180 words). "we" allowed (banned only in Results/Methods).
Funnel written as a PAIR with the Discussion: every gap raised here is answered there, and nothing is answered that was not first asked.
1. Broad: bacterial genes run in operons and RNA processing organizes the transcriptome above the single gene (Jacob-Monod lineage).
2. Organisms + gap: syn1 + syn3A NAMED up front as genetically modified and reduced versions of naturally occurring M. mycoides, defined at the genome-design + protein-essentiality level, their RNA-level organization (operons, processing, RNA→protein) uncharted; para 1 now ENDS with the three phenotype differences (slower division, denser cytosol, lost persistent chromosome contacts ← low HupA) the Discussion answers.
3. Three questions = the three Discussion clusters: (a) how is expression organized at the RNA level; (b) does transcription predict the proteome at minimal complexity; (c) what does minimization do to the expression program.
4. Approach (para 3, expanded 2026-06-28): long-read resolves co-transcription/isoforms/antisense, Illumina = quantification standard; PacBio (syn1) + ONT (syn3A) + matched proteomics → operon-resolved map for both organisms + the reduction comparison; foundation for whole-cell modeling.
Cites: `jacob_genetic_1961`, `gibson_creation_2010`, `hutchison_design_2016`, `breuer_essential_2019`, `pelletier_genetic_2021`, `gilbert_generating_2021`, `stark_rna_2019`, `byrne_realizing_2019`, `mattick_deciphering_2024`, `yan_smrt-cappable-seq_2018`, `grunberger_nanopore_2022`, `thornburg_fundamental_2022`, `thornburg_bringing_2026`.

# Discussion

**Tex file:** `Manuscript/sections/discussion.tex` — DRAFTED. Unheaded flowing prose, 6 short paragraphs; de-recapped into 3 cross-cutting ideas (NOT one bullet per R-section, the recap anti-pattern). NOTE 2026-06-28: Intro + Discussion now over the 500-word budget after the phenotype/long-read expansions (organisms framing + phenotype differences in Intro; phenotype-reasons + rProtein-imbalance in Discussion para 4); budget deferred for the first draft per author.

Complementarity map (Intro gap → Discussion answer):

| Intro raises | Discussion answers | from |
|---|---|---|
| transcriptome architecture uncharted | 459-operon map + pervasive 3′-biased processing | R1+R2 |
| does transcription predict the proteome | yes, dominant; residual partly = elongation; syn1 oddities are artifacts AND were deleted first | R3+R4 |
| what minimization does to expression | structurally conservative, functionally reallocates toward translation machinery; RNAP + central metabolism down; coherent w/ slow growth | R5+R6 |
| phenotype diffs (genome contacts, cell cycle, dense cytosol) | HupA decapitation → lost chromosome contacts; RNAP+metabolism↓ / degradosome↑ → slower division; rProtein imbalance → disrupted ribosome assembly → compact cytosol | R5+R6 |

6 paragraphs:
1. Take-home: transcription = primary layer shaping the proteome; processing + minimization = modifiers.
2. Architecture (R1+R2): 459 operons w/ matched promoters/terminators; pervasive 3′>5′ erosion; 3′→5′ exonucleolytic clearance-bottleneck HYPOTHESIS (awaits cleavage-site mapping).
3. Transcription→proteome + the loop (R3+R4): transcript level dominant, elongation explains part of residual; antisense/intergenic = mis-annotation + synthetic artifacts, and those regions were deleted first in the reduction.
4. Reduction reallocates AND explains the phenotype differences (R5+R6): para OPENS with the explicit phenotype-reasons statement, then — gene-order-preserving, whole-operon excision, promoter-loss decapitation (HupA showcase) → lost chromosome contacts; retained pool shifts to translation machinery (11 kb rProtein operon triples its share), RNAP + central metabolism down + degradosome up → coherent w/ longer cell cycle; rProtein imbalance (operon up + several other rProteins down) → disrupted ribosome assembly → compact cytosol. [rProtein-imbalance sentence REINSTATED & reframed 2026-06-28 as imbalance→ribosome-assembly→compact-cytosol (NOT the old aggregation framing), to answer the syn3A denser-cytosol phenotype; hedged could/plausible.]
5. Evolution of the minimal cell (Sandberg ALE, `sandberg_adaptive_2023`): our \syna~Illumina = the wild-type exponential ancestor = the evolutionary starting point; under selection for growth \syna~raises rProteins further still (ribosome content → protein-synthetic capacity → growth rate, the growth law); evolution also prunes synthetic baggage — the dispensable tetM/0913 resistance cassette (top tenth of the transcriptome, Illumina sense TPM 5,299, rank 49/496) is silenced under selection, echoing the his3/0918 over-transcription minimization had already silenced. Moger-Reischer (Nature 2023, fitness cost + recovery) deliberately NOT cited (no expression data). The ribosome-profiling caveat (do elevated rProtein transcripts → proportionally more ribosomes?) lives in the limitations paragraph, not here.
6. Limitations + outlook: homology sites await direct cleavage-site mapping (conserved enzyme ≠ conserved RNA structure); ONT truncation → PacBio confirmation; mRNA-pool-share excludes rRNA; rProtein/membrane abundances least reliable (ribosome profiling, which also tests whether elevated rProtein transcripts → ribosomes); closer = gene-coexpression into a 4D whole-cell model of the minimal cell.

Claims SOFTENED vs the raw outline (do not revert): "transcription alone determines covariance" → "dominant determinant"; "explain the longer cell cycle" → "coherent with". The rProtein-imbalance sentence (once CUT as the SPECULATIVE rProtein→aggregation line) was REINSTATED 2026-06-28 in a reframed, hedged form (imbalance → disrupted ribosome assembly → surplus subunits crowd cytosol → compact cytosol; "could"/"plausible"), now needed to answer the denser-cytosol phenotype raised in the new Intro.
Standalone "Relationship to previous literature" section dropped; literature woven into the interpretation (Heard: that return-to-the-Intro is what makes the two complementary).

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

---

## M2 — PacBio long-read sequencing of syn1.0 transcriptome  
**Tex file:** `Manuscript/sections/methods/pacbio_syn1.tex`  
**Analysis:** `Syn1_Transcriptomics/PacBio/PacBio_Processing/`  
**Key params:** 3 technical reps (SRR36012641/642/643) pooled to 2.95 M HiFi reads; custom FLNC recovery (reorientation via H1/BCRC, primer trim, polyA trim) -> 2.62 M; `minimap2 -ax map-hifi --secondary=no` v2.30 (CP002027.1); pysam per-read HQ filter (MAPQ>=20, aln-frac>=0.7, clip<=0.3, |qlen-refspan|<=100, concatemer flag) -> 99.3% retained (2.60 M HQ); samtools v1.22.1 per-strand bedGraph; clustering (`Cluster_Isoform.py`, complete-linkage Chebyshev eps=10 bp): 621 k tuples -> 267 k clusters; depth-based TPM as in M1 (single pooled library, no replicate averaging).  
**Outputs:** `PacBio_Processing/syn1.PacBio.FLNC.sorted.HQ.bam`, `depth_bedgraph/syn1.PacBio.FLNC.HQ.{plus,minus,total}.bedGraph`; `Isoforms_PacBio/isoform_clusters_annotated.tsv`; `Gene_TPM/syn1_pacbio_TPM_profiles.tsv`  
**Inputs:** SRA accessions in `Syn1_Transcriptomics/PacBio/PacBio_Raw/00_retrieve_fastq.sh`  

---

## M3 — Operon identification from PacBio long-read transcriptomics in syn1.0
**Tex file:** `Manuscript/sections/methods/operon_analysis.tex`  
**Analysis:** `Syn1_Transcriptomics/Isoforms_PacBio/Cluster_Isoform.py`, `Syn1_Operon/…`  
**Key params:** clustering thresholds, min reads, TSS/TTS calling rule.  
**Outputs:** `isoform_clusters_annotated.tsv`, `operons.candidate_blocks.tsv`  

---

## M4 — RNA processing and ribonucleases  
**Tex file:** `Manuscript/sections/methods/RNA_processing.tex`  
**Analysis:** _<scripts>_  
**Inputs:** `Genomes_Input/Motif_Identifications.xlsx`  

---

## M5 — Proteomics of syn1 and syn3A  
**Tex file:** `Manuscript/sections/methods/proteomics_syn1_syn3A.tex`  
**Analysis / source:** `Syn1_Syn3A_Proteomics/` — `syn1_proteomics_localization_2026.csv`, `syn3a_proteomics_summary_2026.csv`, `syn3A_proteome_annotated.xlsx`.  
**Key params:** Spectronaut iBAQ -> iPM (iPM_i = 1e6 * iBAQ_i / sum iBAQ_j) per rep, mean across 3 reps; absolute copy number = (iPM/1e6) x total proteins/cell (syn1 ~127 k from dry mass 12.8 fg x 58.2% protein / avg MW); localization via DeepTMHMM (TMRs) + SignalP 6, priority signal-peptide > TMR > cytoplasmic; 2019 vs 2026 measurements; tertiary function annotation built by `report_annotation_stats_syn3A.py`.  
**Numbers to cite:** syn1 detected 721/828 (87.1%); median copy number cytoplasmic 47, lipoprotein 21, membrane 10, extracellular 3 (n = 516/68/126/11).  

---

## M6 — Correlation between transcriptome and proteome
**Tex file:** `Manuscript/sections/methods/corr_transcriptome_proteome.tex`  
**Analysis:** `Transcription_Translation.py` (base correlation), `Translation_Residual_L1_initiation.py` (TIR/OSTIR), `Translation_Residual_L2_elongation.py` (CAI), `Translation_Residual_L3_degradation.py` + `Genomes_Input/Homology_Build.py` (degradation).  
**Key params:** Illumina sense TPM vs iPM, log10 OLS, residual = log10(iPM) − fit; TIR via OSTIR (anti-SD ACCUCCUUU, 30-nt windows, read-weighted); CAI (Sharp & Li, ref set = top 20% by iPM); Mpn half-lives (Burgos 2020) via reciprocal-best-hit blastp, protease-abundance correction (Lon 0.84 / FtsH 2.08, Mpn from Maier 2011).  
**Outputs:** `syn1_genes_transcriptomics_proteomics.csv`, `residual_analysis/`  
**Numbers to cite:** correlation all r=0.61 (R²=0.38, n=717) / cytoplasmic r=0.70 (R²=0.49, n=512); ΔR²: TIR +0.020 (2%), CAI +0.080 (+21%), degradation ≤0.01; shortest Mpn-mapped half-life 4.7 h (median 32 h) vs ~1 h doubling.  

---

## M7 — Novel Transcription and Translation
**Tex file:** `Manuscript/sections/methods/novel_orf.tex`  
**Analysis:** `Syn1_Novel_ORF/Abnormal_Transcripts.py` (antisense classes), `Syn1_Novel_ORF/Novel_translation.ipynb` (novel ORFs).  
**Key params:** antisense labeling base-by-base vs gene model, 3 classes (spurious-promoter / read-through / embedded); OSTIR all-start-codon scan (anti-SD ACCUCCUUU, genetic code 4 / UGA=Trp, ORFs >=15 aa); synthesis-rate rank (reads x TIR), top 100, in-silico trypsin (1 missed cleavage, 7-25 aa), Spectronaut re-search of augmented DB.  
**Numbers to cite:** 1.4% antisense (267 isoforms -> 89 clusters); spurious-promoter 59 (66%), read-through 30 (34%, incl. 4 embedded); 837 abnormal isoforms -> ~29,000 candidate ORFs -> top 100 -> 48 unique -> 47 with proteotypic peptides; 2 ORFs confirmed by MS (near MMSYN1_0592 [revised NOVEL_PEP_002] and the 54-aa N-term extension of MMSYN1_0768 [revised NOVEL_PEP_043, = NOVEL_PEP_030 in the old-cluster MS search Excel]), both deleted in syn3A. NOTE: the Spectronaut MS search ran on the old-cluster candidate DB; both confirmed peptides are reproduced in the revised top-100.  

---

## M8 — Oxford Nanopore (ONT) and Illumina sequencing of syn3A transcriptome
**Tex file:** `Manuscript/sections/methods/ont_illumina_syn3A.tex`  
**Analysis:** `Syn3A_Transcriptomics/ONT/ONT_Processing/` (ONT) + `Syn3A_Transcriptomics/Illumina/Illumina_Processing/` (Illumina)  
**Key params:** ONT direct-RNA, `minimap2 -ax map-ont` (NOT splice, bacteria are intron-less), per-strand depth; Illumina syn3A paired-end bowtie2 (dUTP / fr-firststrand), per-strand bedGraph.  
**Inputs:** ONT raw `Syn3A_Transcriptomics/ONT/ONT_Raw/`; Illumina syn3A SRA accessions (SRR19432056/57 mate pair) via `Syn3A_Transcriptomics/Illumina/Illumina_Raw/00_retrive_fastq.sh`.  

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


---