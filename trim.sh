#!/bin/bash
# Add fastp to PATH
export PATH=$PATH:/dss/dsshome1/08/ge85jek2/miniconda3/envs/ncbi/bin

# Directory containing the FASTQ files
input_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/vaishnavi/vaishnavi_masterthesis/fastq"
# Output directory for trimmed FASTQ files
output_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/vaishnavi/vaishnavi_masterthesis/trimmed"

# Ensure output directory exists
mkdir -p "$output_dir"

# Detect paired-end files
echo "Processing paired-end files..."
for fq1 in "${input_dir}"/*_1.fastq; do
    fq2="${fq1/_1.fastq/_2.fastq}"
    if [[ -f "$fq1" && -f "$fq2" ]]; then
        # Extract base name without _1.fastq
        base_name=$(basename "$fq1" "_1.fastq")
        out1="${output_dir}/${base_name}_trimmed_1.fastq"
        out2="${output_dir}/${base_name}_trimmed_2.fastq"

        # Run fastp for paired-end data
        fastp -i "$fq1" -I "$fq2" -o "$out1" -O "$out2"
        echo "Trimmed: ${base_name}_1.fastq and ${base_name}_2.fastq"
    else
        echo "Missing paired file for: $fq1, skipping..."
    fi
done

# Detect single-end files
echo "Processing single-end files..."
for fq in "${input_dir}"/*.fastq; do
    # Skip paired-end files
    if [[ "$fq" == *_1.fastq || "$fq" == *_2.fastq ]]; then
        continue
    fi

    # Extract base name without extension
    base_name=$(basename "$fq" .fastq)
    out="${output_dir}/${base_name}_trimmed.fastq"

    # Run fastp for single-end data
    fastp -i "$fq" -o "$out"
    echo "Trimmed: ${base_name}.fastq"
done

