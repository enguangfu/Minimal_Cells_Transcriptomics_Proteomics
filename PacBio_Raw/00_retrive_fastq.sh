# 1) Define the list of accessions of three PacBio datasets of Syn1
SRR_ACCESSIONS=("SRR36012643" "SRR36012642" "SRR36012641")

FASTERQ_PATH="/home/enguang/software/sratoolkit.3.4.0-ubuntu64/bin/fasterq-dump"  # Adjust this path if fasterq-dump is not in your PATH

# 2) Loop through each accession
for ID in "${SRR_ACCESSIONS[@]}"; do
    echo "Processing $ID..."
    
    # Run fasterq-dump (using -e 16 for multi-threading speed)
    "$FASTERQ_PATH" --split-files --progress -e 16 "$ID"
    
    # Optional: Compress immediately to save disk space
    # gzip "${ID}"*.fastq
done

# 3) Rename the files
# Assuming your files are named SRR36012643.fastq, etc.
mv SRR36012643.fastq PacBio.hifi_reads.rep1.fastq
mv SRR36012642.fastq PacBio.hifi_reads.rep2.fastq
mv SRR36012641.fastq PacBio.hifi_reads.rep3.fastq

# 4) Compress them into fastq.gz
# -p 16: use 16 threads
# -k: keep original files (optional, remove if you want to save space)
pigz -p 16  PacBio.hifi_reads.rep1.fastq
pigz -p 16  PacBio.hifi_reads.rep2.fastq
pigz -p 16  PacBio.hifi_reads.rep3.fastq

# 5) Concatenate them into one merged file
cat PacBio.hifi_reads.rep1.fastq.gz \
    PacBio.hifi_reads.rep2.fastq.gz \
    PacBio.hifi_reads.rep3.fastq.gz > merged.hifi_reads.fastq.gz