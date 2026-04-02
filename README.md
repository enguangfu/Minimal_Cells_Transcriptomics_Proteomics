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

## Abner's Notes/Thoughts ##
### Overall ###
Changed hardcoded path variables to not need to be changed when run by other users.

### PacBio_Raw ###
Fixed yml file to add dependencies. The script will now pick up from where it left off if it terminated early. 

### PacBio_Processing ###
Fixed path issues in first chunk of code
Need to change the comments to something readable (not-run, etc...)
Edit bash 01 to run without path issues

### Illumina_Raw ###
Edited fasterq path

### Illumina_Processing ###
Works great!

### novel_orf_discovery ###
Edited some paths to make congruent with current directory structure
Some questions on methodology:
1. Why were the 5' and 3' binning thresholds set to different numbers (20bp and 30bp, respectively) for isoform definition?
2. What are the distinctions between "Strict" and "Exploratory" abnormal isoforms? Are these retained in the novel/abnormal ORF analysis?
3. When scoring the abnormal/novel isoforms, were any filters included to exclude lower scoring isoforms from the trypsan digestion step?
4. Can we generate some descriptive statistics on the start_spread_mad_bp and end_spread_mad_bp metrics for the isoforms? Maybe a histogram of both?
5. Would it be possible to generate the number of isoforms present before and after they were merged with one another based on the 15bp threshold?

### .gitignore ###
Added all results files to the gitignore to not crowd the repository


## Computation Environment

Install Conda environment with yml file under /env/

Install other packages in /env/extra_softwares.txt

**NOTE:** All other software in the extra_softwares.txt file are now included in the yml file. No need to include this extra file.