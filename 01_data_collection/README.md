# Data collection (step 01)

This step optionally downloads raw NGS data from a vendor SFTP mirror, then processes paired-end FASTQ files into merged reads, Bartender clustering, and per-sample plus combined tables for later UltraSeq steps.

## Layout

| Path | Role |
|------|------|
| `01-data_download.bash` | SLURM job: `lftp` mirror from the vendor into `NGS_DIR` (see repo `config.sh`). |
| `02-info_extraction_single_guide.bash` | SLURM job: AdapterRemoval → parsing → Bartender → aggregation for every line in `data/NGS_address`. |
| `data/NGS_address` | One sample per line: `R1_path,R2_path,sample_id` (comma-separated). |
| `data/NGS_address_test` | Small list for dry runs. |
| `data/guide_reference-GS_single_guide.csv` | Guide library with columns `Gene` and `gRNA`. |
| `data/reads_distribution_summary.csv` | Optional QC summary (paths must match your FASTQ layout). |
| `main_code/` | Python scripts invoked by `02-info_extraction_single_guide.bash`. |
| `auxiliary_code/` | Notebooks to build `NGS_address`-style files and read-level QC (optional). |

The pipeline creates working directories under `data/`, including `Merging/`, `Bartender/`, and `Processed_data/`.

## Configuration

Edit **`../config.sh`** (repository root):

- **`PROJECT_DIR`** — project root used to resolve `01_data_collection` paths.
- **`NGS_DIR`** — target directory for `01-data_download.bash` (must not already exist; the script exits if it does).

Adjust **`#SBATCH`** lines in each bash driver (mail, account, partition, memory, wall time) for your cluster.

## Download (`01-data_download.bash`)

Runs `lftp`/`mirror` into `NGS_DIR`. FTP host, port, and credentials are set inside the script; treat them as sensitive and avoid sharing a public copy with real passwords.

Submit from `01_data_collection/` (or fix paths to `./log`):

```bash
sbatch 01-data_download.bash
```

## Extraction pipeline (`02-info_extraction_single_guide.bash`)

**Environment:** `module load adapterremoval/2.3.1`, then Conda `~/miniconda3` with env **`TubaSeq_Ultra`** (change the `conda activate` line if your env name differs).

**Inputs:** `data/NGS_address` must list absolute paths to your `R1`/`R2` FASTQ files and a short `sample_id` per line.

**Outputs:** Per-sample folders under `data/Processed_data/<sample_id>/`, plus `gRNA_clonalbarcode_combined.csv` under `data/Processed_data/` after the final aggregation step.

```bash
sbatch 02-info_extraction_single_guide.bash
```

**Logs:** SLURM stdout/stderr go to `./log/<jobname>_<jobid>.out` / `.err`. The script also runs a `sed` filter assuming a file named `slurm-<jobid>.out` in the submission directory; if your cluster only writes `./log/...`, adjust that block or ignore the extra “clean” step.

## Python scripts (`main_code/`)

| Script | Purpose |
|--------|---------|
| `single_guide_parsing.py` | Extract gRNA and clonal barcode reads into Bartender-style inputs. |
| `single_guide_aggregate_sgRNA.py` | Map sgRNA clusters to the reference; emit per–sgRNA barcode Bartender inputs. |
| `single_guide_aggregate_barcode.py` | Merge barcode clustering with sgRNA assignments for one sample. |
| `single_guide_aggregate_sample.py` | Concatenate all samples into `gRNA_clonalbarcode_combined.csv`. |

## Auxiliary notebooks (`auxiliary_code/`)

- **`A1-NGS_data_address_generation.ipynb`** — help generate `NGS_address`-compatible path lists.
- **`A2-NGS_reads_distribution_after_parsing.ipynb`** — summarize read distributions (optional QC).
