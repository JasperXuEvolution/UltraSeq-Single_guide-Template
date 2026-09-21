#!/bin/bash
# 02a-info_extraction_per_sample.bash
#
# PARALLEL per-sample extraction as a SLURM JOB ARRAY.
# Each array task processes ONE sample (one line of data/NGS_address), running the
# full AdapterRemoval -> parsing -> Bartender clustering -> per-sample aggregation.
# The final cross-sample merge is done separately by 02b-aggregate_samples.bash.
#
# Do not sbatch this directly by hand unless you set --array yourself; use
# submit_02_parallel.sh, which sizes the array and chains the aggregation job.

# -----------------------------------------------------------
# SLURM directives (per-task resources; one task = one sample)
# -----------------------------------------------------------
#SBATCH --job-name=02a-extract_per_sample
#SBATCH --mail-user=your.email@example.com
#SBATCH --mail-type=END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50g
#SBATCH --time=1-00:00:00
#SBATCH --account=your_slurm_account
#SBATCH --partition=batch
#SBATCH --output=./log/%x_%A_%a.out    # %A = array job id, %a = array task id
#SBATCH --error=./log/%x_%A_%a.err

set -o pipefail

# -----------------------------------------------------------
# Environment
# -----------------------------------------------------------
source ../config.sh
module load adapterremoval/2.3.1
source ~/miniconda3/etc/profile.d/conda.sh
conda activate TubaSeq_Ultra
# Ensure bartender_single_com is on PATH (edit to your install location). A bare
# `sbatch` from a login shell usually inherits it via --export=ALL, but setting it
# here avoids silent failures when the job is submitted non-interactively.
export PATH="$HOME/bin:$PATH"

# -----------------------------------------------------------
# Paths
# -----------------------------------------------------------
working_dir="$PROJECT_DIR/01_data_collection"
input_data_info_address="$working_dir/data/NGS_address"
guide_ref="$working_dir/data/guide_reference-GS_single_guide.csv"
python_script_dir="$working_dir/python_scripts"

Step1_address="$working_dir/data/Merging"
Step2_address="$working_dir/data/Bartender"
Step5_address="$working_dir/data/Processed_data"
mkdir -p "$Step1_address" "$Step2_address" "$Step5_address"

# -----------------------------------------------------------
# Select THIS task's sample = line number SLURM_ARRAY_TASK_ID of NGS_address
# -----------------------------------------------------------
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
   echo "ERROR: SLURM_ARRAY_TASK_ID is not set. Submit via submit_02_parallel.sh or add --array." >&2
   exit 1
fi

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$input_data_info_address")
if [ -z "$line" ]; then
   echo "No sample on line ${SLURM_ARRAY_TASK_ID} of NGS_address; nothing to do."
   exit 0
fi

r1=$(echo "$line" | cut -d',' -f1)
r2=$(echo "$line" | cut -d',' -f2)
sampleID=$(echo "$line" | cut -d',' -f3)
echo "Array task ${SLURM_ARRAY_TASK_ID}: processing sample ${sampleID}"

# --- Step 1: Sequence trimming and merging (AdapterRemoval) ---
temp_folder1="$Step1_address/$sampleID"
mkdir -p "$temp_folder1"
AdapterRemoval --file1 "$r1" --file2 "$r2" \
   --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
   --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
   --basename "$temp_folder1/Merged" --collapse --gzip
echo "For sample ${sampleID}, sequence merging is finished."

# --- Step 2: Generate Bartender input (extract sgRNA + clonal barcode) ---
temp_folder2="$Step2_address/$sampleID"
mkdir -p "$temp_folder2"
python3 "$python_script_dir/single_guide_parsing.py" --a "$temp_folder1/Merged.collapsed.gz" --o "$temp_folder2"

# --- Step 3: Cluster sgRNA and map to reference ---
bartender_single_com -z -1 -d 2 -l 5 -f "$temp_folder2/gRNA.bartender" -o "$temp_folder2/gRNA"
temp_folder3="$Step2_address/$sampleID/Clonal_barcode"
mkdir -p "$temp_folder3"
python3 "$python_script_dir/single_guide_aggregate_sgRNA.py" --a1 "$temp_folder2/gRNA_barcode.csv" \
   --a2 "$temp_folder2/gRNA_cluster.csv" --a3 "$guide_ref" \
   --a5 "$temp_folder2/gRNA.bartender" --a6 "$temp_folder2/clonalbarcode.bartender" \
   --o "$temp_folder2/"

# --- Step 4: Cluster clonal barcodes for each sgRNA ---
while read -r line2; do
   new_name=${line2/.bartender/}
   bartender_single_com -z -1 -d 1 -l 5 -f "$line2" -o "$new_name"
done < "$temp_folder2/Bartender_input_address"

# --- Step 5: Combine data for this sample ---
temp_folder4="$Step5_address/$sampleID"
mkdir -p "$temp_folder4"
python3 "$python_script_dir/single_guide_aggregate_barcode.py" --a "$temp_folder2" --o "$temp_folder4/"

# Clean up per-sample intermediate to free disk space.
rm -r "$temp_folder3"
echo "Sample ${sampleID} finished."

# Per-task job statistics.
sacct --format=JobID,JobName,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j "${SLURM_JOB_ID}"
