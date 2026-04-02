#!/bin/bash
#
# Adaptive bootstrapping (non-plasmid). Set paths in ../../config.sh before running.
# SLURM: replace mail-user and account with your cluster values.
#

#SBATCH --job-name=01-BT_Adaptive
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.com
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=500g
#SBATCH --time=10-24:00:00
#SBATCH --account=your_slurm_account
#SBATCH --partition=batch

source ../../config.sh
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate TubaSeq_Ultra

working_dir="$PROJECT_DIR/03_bootstrapping"
BATCH_NAME="${BATCH_NAME:-UltraSeq_example}"
out_prefix="$working_dir/results/$BATCH_NAME"
mkdir -p "$working_dir/results"

input_data="$working_dir/data/GS_BTHC_KT_tumor.parquet"

python3 "$working_dir/python_scripts/TubaSeq_Ultra_Boostrapping.py" \
  --a0 "$input_data" \
  --a2 200 --a3 100 --a4 1000 --a5 BTHC --a6 KT --a7 100 \
  --o1 "$out_prefix" \
  --o2 "$out_prefix" \
  --l1 50 60 70 80 90 95 \
  --m 'N' --c 'Yes'

sacct --format=JobID,JobName,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j "$SLURM_JOBID"
