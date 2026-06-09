# Full-length transcriptomics and proteomics of synthetic minimal cells

Analysis code accompanying the manuscript:

> **Full-length transcriptomics and proteomics reveal how genome minimization reshapes gene expression in synthetic bacteria**

We combine full-length long-read RNA sequencing (PacBio for **JCVI-syn1.0**, Oxford Nanopore for **JCVI-syn3A**), short-read Illumina quantification, and matched proteomics to build a genome-wide, operon-resolved map of transcription, RNA processing, and translation for both synthetic cells, and to compare how halving the genome reshapes that map.

This repository holds the complete pipeline, from raw-read retrieval to the figures and Supplementary Data of the paper.

---

## Repository structure

```
.
├── Genomes_Input/            reference genomes (FASTA + GFF) for both organisms
│
├── Syn1_Transcriptomics/     syn1 PacBio + Illumina: read processing, isoforms, per-gene TPM
├── Syn3A_Transcriptomics/    syn3A ONT + Illumina: read processing, isoforms, per-gene TPM
│
├── Syn1_Operon/              syn1 operon segmentation, annotation, visualization (459 operons)
├── Syn3A_Operon/             syn3A operon segmentation + annotation
│
├── Syn1_Syn3A_Proteomics/    proteomics tables (relative iPM + absolute copy numbers) for both
├── Syn1_Corr_RNA_Proteins/   syn1 transcriptome × proteome correlation
├── Syn1_RNase/               RNA-processing / ribonuclease analysis (3' erosion, RNase-site mapping)
├── Syn1_Novel_ORF/           novel-ORF / antisense / intergenic transcription discovery
│
├── Genome_Reduction/         syn1 → syn3A comparison (deletions, operon remodeling, expression)
│
├── Manuscript/               LaTeX manuscript + figures + Supplementary Data assembly
└── env/                      conda environment specification
```

A detailed, file-level map of every script and its outputs lives in [`CLAUDE.md`](CLAUDE.md); output conventions are in [`OUTPUT.md`](OUTPUT.md).

---

## Environment

All analysis runs in a single conda environment (Python 3.10):

```bash
conda env create -f env/RNAseq.yml
conda activate RNAseq
```

It bundles the Python scientific stack (NumPy, pandas, SciPy, Matplotlib, Biopython, pysam) and the main bioinformatics tools: **minimap2** 2.30, **bowtie2** 2.5.5, **samtools** 1.23, **BLAST+** 2.17, and **ViennaRNA** 2.6.4. Two tools are installed separately — **SRA-Toolkit** (read download) and **TransTermHP** (terminator prediction) — see [`env/extra_softwares.txt`](env/extra_softwares.txt).

---

## Workflow

Data flows from raw reads (left/top) through per-organism processing and into the cross-organism comparison and the manuscript:

```mermaid
flowchart TD
    GEN["Genomes_Input<br/>reference FASTA + GFF"]

    subgraph SYN1["Syn1 transcriptome"]
        PB["PacBio HiFi"] --> PBP["PacBio_Processing"] --> ISO["Isoforms_PacBio"]
        ILL1["Illumina"] --> ILP1["Illumina_Processing"]
        PBP --> TPM1["Gene_TPM"]
        ILP1 --> TPM1
    end

    subgraph SYN3["Syn3A transcriptome"]
        ONT["ONT direct-RNA"] --> ONTP["ONT_Processing"]
        ILL3["Illumina"] --> ILP3["Illumina_Processing"]
        ONTP --> TPM3["Gene_TPM"]
        ILP3 --> TPM3
    end

    PROT["Syn1_Syn3A_Proteomics"]

    ISO --> OPN["Syn1_Operon<br/>459 operons"]
    ISO --> RNASE["Syn1_RNase"]
    ISO --> NOV["Syn1_Novel_ORF"]
    TPM1 --> CORR["Syn1_Corr_RNA_Proteins"]
    PROT --> CORR
    OPN --> GR["Genome_Reduction<br/>syn1 → syn3A"]
    PROT --> GR
    TPM1 --> GR
    TPM3 --> GR

    OPN --> MS["Manuscript<br/>figures + Supplementary Data S1–S4"]
    RNASE --> MS
    NOV --> MS
    CORR --> MS
    GR --> MS

    GEN -.-> SYN1
    GEN -.-> SYN3
    GEN -.-> OPN
    GEN -.-> GR
```

### Run order

Run the stages in this order; each folder's scripts read the outputs of the stages above it.

1. **Retrieve raw reads** — bash scripts in the `*_Raw/` folders download FASTQs from the NCBI SRA (`*_Transcriptomics/{PacBio,Illumina,ONT}/*_Raw/`).
2. **Process reads** — the `*_Processing/` folders map reads (minimap2 for long reads, bowtie2 for Illumina), sort/index with samtools, and emit per-strand depth bedGraphs.
3. **Cluster isoforms** — `Isoforms_PacBio/` (syn1) and `Isoform_Cluster/` (syn3A) collapse full-length reads into isoform clusters.
4. **Per-gene TPM** — `Gene_TPM/` computes sense/antisense TPM per gene from the depth tracks.
5. **Proteome** — `Syn1_Syn3A_Proteomics/` builds the per-protein relative (iPM) and absolute abundance tables.
6. **Operons** — `Syn1_Operon/` (and `Syn3A_Operon/`) segment and annotate operons from the isoforms, with promoter and terminator signatures.
7. **Per-organism analyses** — `Syn1_Corr_RNA_Proteins/` (RNA↔protein correlation), `Syn1_RNase/` (RNA processing + ribonuclease cleavage-site mapping), `Syn1_Novel_ORF/` (antisense / intergenic / novel ORFs).
8. **Genome reduction** — `Genome_Reduction/` runs scripts `01`→`10` in numeric order to recast the syn1→syn3A deletions as operon junctions and quantify the transcriptome/proteome reallocation. Run with `Genome_Reduction/` as the working directory.

> Most Python scripts carry their full method, parameters, and a result summary in a header docstring, and write a companion `.txt` log next to their outputs.

---

## Key outputs (Supplementary Data)

| | File | Built by |
|---|---|---|
| **S1** | `operon.xlsx` — per-operon table (boundaries, signals, complexes) | `Syn1_Operon/build_operon_xlsx.py` |
| **S2** | `syn1_omics.xlsx` — paired transcriptome + proteome for 911 syn1 genes | `Syn1_Corr_RNA_Proteins/` |
| **S3** | `genome_reduction.xlsx` — deletions, junctions, per-gene expression change | `Genome_Reduction/` |
| **S4** | `Supplementary_Data_S4_QC.zip` — RNA-sample QC reports | `Manuscript/SI/build_S4_qc.py` |

---

## Data availability

Raw sequencing data are deposited in the NCBI SRA under **BioProject [PRJNA1359397](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1359397)** (syn1 Illumina + PacBio, syn3A ONT). The syn3A Illumina data are from Sandberg *et al.* (2023). Reference genomes: **CP002027.1** (JCVI-syn1.0) and **CP016816.2** (JCVI-syn3A). Mass-spectrometry proteomics are on the MassIVE repository under **MSV000099558**.

---

## Citation

If you use this code or data, please cite the manuscript above. (Full citation details will be added upon publication.)
