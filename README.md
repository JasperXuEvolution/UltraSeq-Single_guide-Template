# UltraSeq Single Guide

Pipeline for single-guide UltraSeq: raw NGS processing (step 01), notebook-based cleaning and QC (step 02), then statistical bootstrapping on the processed tables (step 03). Configure paths in **`config.sh`** at the repository root before running SLURM jobs.

```
UltraSeq-Single_guide/
├── config.sh                 # PROJECT_DIR, NGS_DIR (edit for your system)
├── TubaSeq_Ultra.yml         # Conda environment spec (reference)
├── 01_data_collection/       # Step 01 — FASTQ download (optional) + extraction
│   ├── 01-data_download.bash
│   ├── 02-info_extraction_single_guide.bash   # serial: all samples in one job
│   ├── 02a-info_extraction_per_sample.bash    # parallel: SLURM array, one task/sample
│   ├── 02b-aggregate_samples.bash             # parallel: final cross-sample merge
│   ├── submit_02_parallel.sh                  # submits 02a array + 02b (afterok)
│   ├── data/                 # NGS_address, guide reference, intermediates
│   ├── python_scripts/       # Parsing and aggregation Python entrypoints
│   ├── auxiliary_code/       # Optional notebooks (address lists, QC)
│   └── README.md             # Full documentation for step 01
├── 02_data_cleaning_and_QC/  # Step 02 — QC notebooks → processed parquet for step 03
│   ├── *-QC_part1-*.ipynb    # Annotate + spike-in QC
│   ├── *-QC_part2-*.ipynb    # Sample QC + emit bootstrapping input
│   ├── UltraSeq-SampleSpecificAnalysisFunction-*.py  # Shared helpers
│   ├── data/                 # Inputs from step 01 + generated tables (git-ignored)
│   ├── figs/                 # QC figures (git-ignored)
│   └── README.md             # Full documentation for step 02
├── 03_bootstrapping/         # Step 03 — bootstrapping on parquet inputs
│   ├── bash/                 # SLURM wrappers (adaptive/normal × plasmid or not)
│   ├── data/                 # Input parquet files from step 02 (+ plasmid for some jobs)
│   ├── results/              # Output prefix root (created by scripts)
│   ├── python_scripts/       # TubaSeq_Ultra_Boostrapping.py
│   └── README.md             # Full documentation for step 03
└── README.md
```

## Requirements

- **SLURM** cluster with `sbatch` (steps 01 and 03 are SLURM jobs).
- **Conda env `TubaSeq_Ultra`** (see `TubaSeq_Ultra.yml`) with `pandas`, `regex`, `numpy`, etc. Adjust the `conda activate` line in the scripts if your env name differs.
- **AdapterRemoval** — loaded via `module load adapterremoval/2.3.1` (edit for your cluster).
- **Bartender** — `bartender_single_com` must be on `PATH` when the extraction job runs. The scripts add `export PATH="$HOME/bin:$PATH"` after `conda activate`; edit that to point at your Bartender install. (Without it, a job submitted non-interactively can silently skip clustering — the `.bartender` files are written but no `*_cluster.csv`, and aggregation then fails with "No objects to concatenate".)

## Step 01 — Data collection

- **Download:** `01_data_collection/01-data_download.bash` mirrors vendor SFTP data into `NGS_DIR` from `config.sh`.
- **Extraction (serial):** `02-info_extraction_single_guide.bash` runs AdapterRemoval → parsing → Bartender clustering → aggregation for every sample in `data/NGS_address`, in a single job.
- **Extraction (parallel):** `bash submit_02_parallel.sh` runs the same pipeline as a SLURM **job array** — one task per sample (`02a`), then a single merge job (`02b`) that starts only after all tasks succeed. Same output as the serial script, much faster for many samples.

Details: [01_data_collection/README.md](01_data_collection/README.md).

## Step 02 — Data cleaning and QC

- **Notebooks:** run `02_data_cleaning_and_QC/*-QC_part1-*.ipynb` then `*-QC_part2-*.ipynb` from that directory.
- **Input:** `gRNA_clonalbarcode_combined.csv` from step 01, plus guide/sample metadata (see step-02 README).
- **Output:** a processed `.parquet` table (e.g. `GS_BTHC_KT_tumor.parquet`) that step 03 consumes.

Details: [02_data_cleaning_and_QC/README.md](02_data_cleaning_and_QC/README.md).

## Step 03 — Bootstrapping

- **Scripts:** `03_bootstrapping/bash/` — `01`/`02` (non-plasmid) and `03`/`04` (plasmid) variants for adaptive vs normal bootstrapping.
- **Outputs:** under `03_bootstrapping/results/<BATCH_NAME>/` (default batch name is set in each shell script).

Details: [03_bootstrapping/README.md](03_bootstrapping/README.md).

## Configuration

Edit **`config.sh`**:

- **`PROJECT_DIR`** — path to this repository (or your project root used in paths).
- **`NGS_DIR`** — where step 01 download places raw NGS data.

Bootstrapping jobs source `config.sh` from `03_bootstrapping/bash/` via `../../config.sh`.
