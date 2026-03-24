# Transcriptomics and Proteomics of Minimal Cells

--- 

## File Architecture and Analysis Pipeline

File folder display here

- Run bash script under PacBio_Raw to retrieve fastq files from NCBI
- Run Jupyter Notebook under PacBio_Processing to process and map PacBio cDNA reads to Syn1 geonme
- Run bash script under Illumin_Raw to retrieve fastq files from NCBI
- Run bash scripts under Illumina_Processing to map Illumina to Syn1 genome
- Run Jupyter notebook under Proteomics to quantify Syn1 and Syn3a proteome
- Run Jupyter Notebook under Gene_Transcriptomics_Proteomics to correlate transcriptome and proteome

 
## Computation Environment

Install Conda environment with yml file under /env/

Install other packages in /env/extra_softwares.txt