# 1) Define the list of accessions of three Illumina datasets of Syn1
SRR_ACCESSIONS=("SRR35996298" "SRR35996297" "SRR35996296")

FASTERQ_PATH=$(which fasterq-dump 2>/dev/null)
if [ -z "$FASTERQ_PATH" ]; then
    echo "Error: fasterq-dump not found in the active conda environment. Please activate the appropriate environment."
    exit 1
fi

# 2) Loop through each accession
for ID in "${SRR_ACCESSIONS[@]}"; do
    echo "Processing $ID..."
    
    # Run fasterq-dump (using -e 16 for multi-threading speed)
    "$FASTERQ_PATH" --split-files --progress -e 16 "$ID"
    
    # Optional: Compress immediately to save disk space
    # gzip "${ID}"*.fastq
done
