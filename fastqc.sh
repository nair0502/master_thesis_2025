#!/bin/bash
export PATH=$PATH:/dss/dsshome1/08/ge85jek2/miniconda3/envs/ncbi/bin

# Directory containing the trimmed FASTQ files
input_dir="vaishnavi/vaishnavi_masterthesis/trimmed"
# Output directory for FastQC reports
output_dir="vaishnavi/vaishnavi_masterthesis/fastqc_trimmed_results"

# Ensure output directory exists
mkdir -p "$output_dir"

# Process trimmed files with FastQC
echo "Starting FastQC analysis on trimmed files in $input_dir..."
echo "FastQC results will be saved in $output_dir."

# Loop through paired-end trimmed files
for fq1 in "${input_dir}"/*_trimmed_1.fastq; do
    fq2="${fq1/_trimmed_1.fastq/_trimmed_2.fastq}" # Find the corresponding _trimmed_2.fastq file
    base_name=$(basename "$fq1" "_trimmed_1.fastq")
    
    if [[ -f "$fq2" ]]; then
        # Process paired-end files
        echo "Processing paired-end files: $fq1 and $fq2"
        fastqc -o "$output_dir" "$fq1" "$fq2"
        echo "Processed: $fq1 and $fq2"
    else
        echo "Warning: Missing paired file for $fq1. Skipping..."
    fi
done

# Loop through single-end trimmed files
for fq in "${input_dir}"/*_trimmed.fastq; do
    # Skip paired-end files
    if [[ "$fq" == *_trimmed_1.fastq || "$fq" == *_trimmed_2.fastq ]]; then
        continue
    fi

    # Process single-end files
    base_name=$(basename "$fq" "_trimmed.fastq")
    echo "Processing single-end file: $fq"
    fastqc -o "$output_dir" "$fq"
    echo "Processed: $fq"
done

echo "FastQC analysis complete. Reports are saved in $output_dir."

