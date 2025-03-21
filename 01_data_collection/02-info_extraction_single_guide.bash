#!/bin/bash
# 02-info_extraction_single_guide.bash

# -----------------------------------------------------------
# SLURM SBATCH Directives: Job submission options for the cluster.
# -----------------------------------------------------------
#SBATCH --job-name=UltraSeq_pipeline_single_guide
#SBATCH --mail-user=xhq@stanford.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50g 
#SBATCH --time=2-24:00:00
#SBATCH --account=mwinslow
#SBATCH --partition=batch

# -----------------------------------------------------------
# Environment Setup: Load configuration, modules, and Conda environment.
# -----------------------------------------------------------
source ../config.sh

module load adapterremoval/2.3.1
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate UltraSeq

# -----------------------------------------------------------
# Directories and Input Files Setup
# -----------------------------------------------------------
working_dir="$PROJECT_DIR/01_data_collection"
input_data_info_address="$working_dir/data/NGS_address"
# guide_ref="$working_dir/data/guide_reference.csv"
guide_ref="$working_dir/data/guide_reference-GS_single_guide.csv"
python_script_dir="$working_dir/main_code"

# Explanation for input files:
#   input_data_info_address:
#     A text file where each line is comma-separated with three fields:
#       1. The path to the FASTQ file for read 1.
#       2. The path to the FASTQ file for read 2.
#       3. The sample ID (used for naming output subdirectories).
#
#   guide_ref:
#     A CSV file serving as a reference for guide sequences with at least two files:
#        1. A column called Gene
#        1. A column called gRNA (which is effectively gRNA_complete)

# Define additional input/output directories.
Step1_address="$working_dir/data/Merging"
Step2_address="$working_dir/data/Bartender"
Step5_address="$working_dir/data/Processed_data"

# Create output directories if they don't exist.
mkdir -p "$Step1_address"
mkdir -p "$Step2_address"
mkdir -p "$Step5_address"

# -----------------------------------------------------------
# Main Processing Loop: Process each sample listed in the input data file.
# -----------------------------------------------------------
while read -r line; do
   # Parse the comma-separated fields.
   r1=$(echo "$line" | cut -d',' -f1)
   r2=$(echo "$line" | cut -d',' -f2)
   sampleID=$(echo "$line" | cut -d',' -f3)
   
   # --- Step 1: Sequence Trimming and Merging using AdapterRemoval ---
   # Create a subdirectory for the sample's merged reads.
   temp_folder1="$Step1_address/$sampleID"
   mkdir -p "$temp_folder1"
   
   # Run AdapterRemoval to trim adapters, collapse overlapping paired reads,
   # and compress the output.
   AdapterRemoval --file1 "$r1" --file2 "$r2" \
      --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
      --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
      --basename "$temp_folder1/Merged" --collapse --gzip
   echo "For sample ${sampleID}, sequence merging is finished."

   # --- Step 2: Generate Bartender Input ---
   # Create a subdirectory for bartender outputs for this sample.
   temp_folder2="$Step2_address/$sampleID"
   mkdir -p "$temp_folder2"
   
   # Run the Python script to extract sgRNA and clonal barcode sequences from the merged FASTQ file.
   python3 "$python_script_dir/single_guide_parsing.py" --a "$temp_folder1/Merged.collapsed.gz" --o "$temp_folder2"

   # --- Step 3: Clustering of sgRNA ---
   # Cluster sgRNA sequences using bartender with specified parameters.
   bartender_single_com -z -1 -d 2 -l 5 -f "$temp_folder2/gRNA.bartender" -o "$temp_folder2/gRNA"
   
   # Create a subdirectory for clonal barcode clustering outputs.
   temp_folder3="$Step2_address/$sampleID/Clonal_barcode"
   mkdir -p "$temp_folder3"
   
   # Run the Python script to aggregate sgRNA data with clonal barcode information.
   python3 "$python_script_dir/single_guide_aggregate_sgRNA.py" --a1 "$temp_folder2/gRNA_barcode.csv" \
      --a2 "$temp_folder2/gRNA_cluster.csv" --a3 "$guide_ref" \
      --a5 "$temp_folder2/gRNA.bartender" --a6 "$temp_folder2/clonalbarcode.bartender" \
      --o "$temp_folder2/"

   # --- Step 4: Clustering of Clonal Barcodes ---
   # For each line in the Bartender input address file (for barcode associated with same sgRNA), perform clustering.
   while read -r line2; do 
      # Remove the ".bartender" extension to create the output filename.
      new_name=${line2/.bartender/}
      bartender_single_com -z -1 -d 1 -l 5 -f "$line2" -o "$new_name"
   done < "$temp_folder2/Bartender_input_address"

   # --- Step 5: Combine Data for the Sample ---
   # Create a subdirectory to store the processed data for the sample.
   temp_folder4="$Step5_address/$sampleID"
   mkdir -p "$temp_folder4"
   
   # Run the Python script to combine clonal barcode data for this sample.
   python3 "$python_script_dir/single_guide_aggregate_barcode.py" --a "$temp_folder2" --o "$temp_folder4/"
   
   # Clean up: Remove the Clonal_barcode folder to free up disk space.
   rm -r "$temp_folder3"
done < "$input_data_info_address"

# --- Final Step: Aggregate Data Across All Samples ---
# Combine the processed data from all samples into a single final output.
python3 "$python_script_dir/single_guide_aggregate_sample.py" --o "$Step5_address/"

# Report job statistics using SLURM's sacct command.
sacct --format=JobID,JobName,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j "$SLURM_JOBID"

# --- Clean SLURM Output Log ---
# This command removes the log sections between specified patterns from the SLURM output file
# Note: The SLURM output file is assumed to be named as "slurm-<jobID>.out"
sed '/Loading barcodes from the file/,/Running bartender/d' slurm-${SLURM_JOBID}.out | \
sed '/Trimming paired end reads/,/Processed a total ../d' > slurm-${SLURM_JOBID}_clean.out
