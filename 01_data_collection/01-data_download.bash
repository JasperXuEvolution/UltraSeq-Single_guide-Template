#!/bin/bash
# 01-data_download.bash
# This script downloads NGS data from a remote FTP server using lftp.

# SLURM directives for job submission
#SBATCH --job-name=01-data_download           # Job name
#SBATCH --mail-user=your.email@example.com    # Email for notifications
#SBATCH --cpus-per-task=2                     # Number of CPU cores per task
#SBATCH --nodes=1                             # Number of nodes
#SBATCH --ntasks-per-node=1                   # Number of tasks per node
#SBATCH --mem-per-cpu=2g                      # Memory per CPU core
#SBATCH --time=2-00:00:00                     # Maximum runtime (2 days)
#SBATCH --account=your_slurm_account          # Account name for billing
#SBATCH --output=./log/%x_%j.out              # Standard output log file
#SBATCH --error=./log/%x_%j.err               # Standard error log file
#SBATCH --partition=batch                     # Partition to submit the job

# Ensure the log directory exists
mkdir -p ./log

# Source the configuration file for environment variables
source ../config.sh

# Check if the NGS_DIR already exists to prevent overwriting
if [ -d "$NGS_DIR" ]; then
    echo "Error: Directory $NGS_DIR already exists. Please update config.sh to use a different directory."
    exit 1
fi

# Create the directory for NGS data
mkdir -p $NGS_DIR

# FTP connection parameters (provided by your sequencing vendor)
# NOTE: Do not commit real credentials. Fill these in locally only.
FTP_ACCOUNT="your_ftp_username"        # FTP account username
FTP_PASSWORD="your_ftp_password"       # FTP account password
FTP_HOST="ftp.your-vendor.com"         # FTP server hostname
FTP_PORT="3022"                        # FTP server port

# Download data using lftp
# - `mirror` recursively downloads files from the remote server to the local directory
# - `--verbose` enables detailed output
# - `--use-pget-n=8` enables parallel downloading with 8 connections
# - `-c` continues partially downloaded files
lftp -c "open sftp://$FTP_ACCOUNT:$FTP_PASSWORD@$FTP_HOST:$FTP_PORT; mirror --verbose --use-pget-n=8 -c . $NGS_DIR"

# Log job statistics using sacct
# - Provides details such as JobID, JobName, Submit time, Start time, End time, State, Partition, CPU usage, memory usage, and node list
sacct --format=JobID,JobName,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID
