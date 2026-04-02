# UltraSeq Single Guide

Pipeline for single-guide UltraSeq: raw NGS processing (step 01), then statistical bootstrapping on processed tables (step 03). Configure paths in **`config.sh`** at the repository root before running SLURM jobs.

```
UltraSeq-Single_guide/
├── config.sh                 # PROJECT_DIR, NGS_DIR (edit for your system)
├── TubaSeq_Ultra.yml         # Conda environment spec (reference)
├── 01_data_collection/       # Step 01 — FASTQ download (optional) + extraction
│   ├── 01-data_download.bash
│   ├── 02-info_extraction_single_guide.bash
│   ├── data/                 # NGS_address, guide reference, intermediates
│   ├── python_scripts/       # Parsing and aggregation Python entrypoints
│   ├── auxiliary_code/       # Optional notebooks (address lists, QC)
│   └── README.md             # Full documentation for step 01
├── 03_bootstrapping/         # Step 03 — bootstrapping on parquet inputs
│   ├── bash/                 # SLURM wrappers (adaptive/normal × plasmid or not)
│   ├── data/                 # Input parquet files (tumor; plasmid for some jobs)
│   ├── results/              # Output prefix root (created by scripts)
│   ├── python_scripts/       # TubaSeq_Ultra_Boostrapping.py
│   └── README.md             # Full documentation for step 03
└── README.md
```

## Step 01 — Data collection

- **Download:** `01_data_collection/01-data_download.bash` mirrors vendor SFTP data into `NGS_DIR` from `config.sh`.
- **Pipeline:** `01_data_collection/02-info_extraction_single_guide.bash` runs AdapterRemoval, parsing, Bartender clustering, and aggregation using sample lists in `data/NGS_address`.

Details: [01_data_collection/README.md](01_data_collection/README.md).

## Step 03 — Bootstrapping

- **Scripts:** `03_bootstrapping/bash/` — `01`/`02` (non-plasmid) and `03`/`04` (plasmid) variants for adaptive vs normal bootstrapping.
- **Outputs:** under `03_bootstrapping/results/<BATCH_NAME>/` (default batch name is set in each shell script).

Details: [03_bootstrapping/README.md](03_bootstrapping/README.md).

## Configuration

Edit **`config.sh`**:

- **`PROJECT_DIR`** — path to this repository (or your project root used in paths).
- **`NGS_DIR`** — where step 01 download places raw NGS data.

Bootstrapping jobs source `config.sh` from `03_bootstrapping/bash/` via `../../config.sh`.
