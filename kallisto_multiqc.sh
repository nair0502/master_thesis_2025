#!/bin/bash
# Add MultiQC to PATH
export PATH=$PATH:/dss/dsshome1/08/ge85jek2/miniconda3/envs/multiqc/bin/multiqc

# Input directory containing Kallisto results
Kallisto_results_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/vaishnavi/vaishnavi_masterthesis/kallisto_output"

# Output directory for MultiQC report
multiqc_output_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/vaishnavi/vaishnavi_masterthesis/kallisto_multiqc_report"

# Ensure output directory exists
mkdir -p "$multiqc_output_dir"

# Run MultiQC on the Kallisto results directory
echo "Running MultiQC on Kallisto results..."
multiqc "$Kallisto_results_dir" -o "$multiqc_output_dir"

# Completion message
echo "MultiQC analysis complete."
echo "MultiQC report is saved in: $multiqc_output_dir/multiqc_report.html"

