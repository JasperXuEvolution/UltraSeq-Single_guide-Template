#!/bin/bash
# 02b-aggregate_samples.bash
#
# Final cross-sample merge for the PARALLEL workflow. Concatenates every
# per-sample Combined_deduplexed_df.csv (produced by the 02a array tasks) into
# data/Processed_data/gRNA_clonalbarcode_combined.csv.
#
# Submitted automatically by submit_02_parallel.sh with a dependency so it runs
# only after all 02a array tasks finish successfully.

#SBATCH --job-name=02b-aggregate_samples
#SBATCH --mail-user=your.email@example.com
#SBATCH --mail-type=END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=20g
#SBATCH --time=02:00:00
#SBATCH --account=your_slurm_account
#SBATCH --partition=batch
#SBATCH --output=./log/%x_%j.out
#SBATCH --error=./log/%x_%j.err

source ../config.sh
source ~/miniconda3/etc/profile.d/conda.sh
conda activate TubaSeq_Ultra
# Ensure bartender_single_com is on PATH (edit to your install location). A bare
# `sbatch` from a login shell usually inherits it via --export=ALL, but setting it
# here avoids silent failures when the job is submitted non-interactively.
export PATH="$HOME/bin:$PATH"

working_dir="$PROJECT_DIR/01_data_collection"
python_script_dir="$working_dir/python_scripts"
Step5_address="$working_dir/data/Processed_data"

python3 "$python_script_dir/single_guide_aggregate_sample.py" --o "$Step5_address/"

sacct --format=JobID,JobName,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j "${SLURM_JOB_ID}"
