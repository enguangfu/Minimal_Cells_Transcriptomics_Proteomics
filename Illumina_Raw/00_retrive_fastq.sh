# 1) Define the list of accessions of three Illumina datasets of Syn1
SRR_ACCESSIONS=("SRR35996298" "SRR35996297" "SRR35996296")

FASTERQ_PATH="/home/enguang/software/sratoolkit.3.4.0-ubuntu64/bin/fasterq-dump"  # Adjust this path if fasterq-dump is not in your PATH

# 2) Loop through each accession
for ID in "${SRR_ACCESSIONS[@]}"; do
    echo "Processing $ID..."
    
    # Run fasterq-dump (using -e 16 for multi-threading speed)
    "$FASTERQ_PATH" --split-files --progress -e 16 "$ID"
    
    # Optional: Compress immediately to save disk space
    # gzip "${ID}"*.fastq
done
