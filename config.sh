#!/bin/bash
# config.sh
#
# Central path configuration for the pipeline. Edit these before running any
# SLURM job. All step scripts source this file via ../config.sh or ../../config.sh.

# Path to this repository (project root). Used to resolve step directories.
PROJECT_DIR="/path/to/UltraSeq-Single_guide"

# Directory where step 01 downloads raw NGS data (must NOT already exist for
# 01-data_download.bash; the script exits if it does).
NGS_DIR="/path/to/NGS_Raw_data/<run_name>"
