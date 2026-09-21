#!/bin/bash
# submit_02_parallel.sh
#
# Submit step 02 in PARALLEL: one SLURM array task per sample (02a), then a
# single aggregation job (02b) that runs only after every array task succeeds.
#
# Usage (run from 01_data_collection/):
#     bash submit_02_parallel.sh
#
# MAX_CONCURRENT caps how many samples run at once — raise/lower to fit your
# cluster limits (queue policy, memory, or Bartender license seats).

set -euo pipefail
cd "$(dirname "$0")"
mkdir -p ./log

source ../config.sh
input_data_info_address="$PROJECT_DIR/01_data_collection/data/NGS_address"

MAX_CONCURRENT=20

# Number of samples = number of non-empty lines in NGS_address.
N=$(grep -cve '^[[:space:]]*$' "$input_data_info_address")
if [ "$N" -lt 1 ]; then
   echo "No samples found in $input_data_info_address" >&2
   exit 1
fi
echo "Submitting $N samples (up to $MAX_CONCURRENT running at once)."

# 1) Per-sample array job.
array_jid=$(sbatch --parsable --array=1-"${N}"%"${MAX_CONCURRENT}" 02a-info_extraction_per_sample.bash)
echo "Array job:       $array_jid  (tasks 1-$N)"

# 2) Aggregation job — afterok waits for ALL array tasks to succeed.
#    Use afterany instead if you want to merge whatever completed even on failures.
agg_jid=$(sbatch --parsable --dependency=afterok:"${array_jid}" 02b-aggregate_samples.bash)
echo "Aggregation job: $agg_jid  (starts after $array_jid finishes OK)"
