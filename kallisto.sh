#!/bin/bash

# Add Kallisto to the PATH
export PATH=$PATH:/miniconda3/envs/ncbi/bin/kallisto

# Input files and output directory
index_path="vaishnavi/vaishnavi_masterthesis/tomato_index.idx"
input_dir="vaishnavi/vaishnavi_masterthesis/trimmed"
output_dir="vaishnavi/vaishnavi_masterthesis/kallisto_output"
log_file="${output_dir}/kallisto_log.txt"

# Ensure output directory exists
mkdir -p "$output_dir"
touch "$log_file"

# Start processing
echo "Starting Kallisto quantification..." | tee -a "$log_file"

# Process paired-end reads
for read1 in "${input_dir}"/*_trimmed_1.fastq; do
    base=$(basename "$read1" "_trimmed_1.fastq")
    read2="${input_dir}/${base}_trimmed_2.fastq"
    sample_output_dir="${output_dir}/${base}"

    if [[ -f "$read1" && -f "$read2" ]]; then
        echo "Processing paired-end sample: $base" | tee -a "$log_file"
        mkdir -p "$sample_output_dir"
        
        # Run Kallisto for paired-end reads
        kallisto quant -i "$index_path" -o "$sample_output_dir" "$read1" "$read2" >> "$log_file" 2>&1
        echo "Finished: $base" | tee -a "$log_file"
    else
        echo "Warning: Missing paired file for sample $base. Skipping..." | tee -a "$log_file"
    fi
done

# Process single-end reads
for read in "${input_dir}"/*_trimmed.fastq; do
    # Skip paired-end files
    if [[ "$read" == *_trimmed_1.fastq || "$read" == *_trimmed_2.fastq ]]; then
        continue
    fi

    base=$(basename "$read" "_trimmed.fastq")
    sample_output_dir="${output_dir}/${base}"

    echo "Processing single-end sample: $base" | tee -a "$log_file"
    mkdir -p "$sample_output_dir"

    # Run Kallisto for single-end reads with default fragment length and SD
    kallisto quant -i "$index_path" -o "$sample_output_dir" --single -l 200 -s 20 "$read" >> "$log_file" 2>&1
    echo "Finished: $base" | tee -a "$log_file"
done

echo "Kallisto quantification complete. Log file: $log_file" | tee -a "$log_file"

