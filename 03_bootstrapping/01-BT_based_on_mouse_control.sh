#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=01-BT_based_on_mouse_control
#SBATCH --mail-user=xhq@stanford.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100g 
#SBATCH --time=10-24:00:00
#SBATCH --account=mwinslow
#SBATCH --partition=batch

# general input and output address
source ../config.sh
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate UltraSeq
working_dir="$PROJECT_DIR/04_bootstrapping"
input_data="/labs/mwinslow/Haiqing/Bootstrapping_analysis/GS_BT-All_Data-241230/Input_data/Batch1T4_data-Tumor.parquet"

# command
python3 "$working_dir/main_code/UltraSeq_Boostrapping_GW.py" \
--a0 "$input_data" \
--a2 200 --a3 100 --a4 1000 --a5 KTHC --a6 KT \
--o1 "$working_dir/data/GS_BT_AllData" \
--o2 "$working_dir/data/GS_BT_AllData" \
--l1 50 60 70 80 90 95 \
--m 'N' --c 'No' 

sacct --format=JobID,JobName,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID
