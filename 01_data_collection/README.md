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
| `python_scripts/` | Python scripts invoked by `02-info_extraction_single_guide.bash`. |
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

There are two equivalent ways to run this step — pick one:

### Option A — serial (single job)

Processes every sample in one job, one after another. Simplest; slowest for many samples.

```bash
sbatch 02-info_extraction_single_guide.bash
```

### Option B — parallel (SLURM job array, recommended for many samples)

Runs one array task **per sample** (all samples at once, up to a concurrency cap), then a single aggregation job that starts only after every task succeeds. Same result as Option A, much faster wall-clock. Files:

| File | Role |
|------|------|
| `02a-info_extraction_per_sample.bash` | Job array; task *i* processes line *i* of `NGS_address` (steps 1–5 for one sample). |
| `02b-aggregate_samples.bash` | Final merge into `gRNA_clonalbarcode_combined.csv`. |
| `submit_02_parallel.sh` | Counts samples, sizes the array, and submits `02b` with `--dependency=afterok` on the array. |

```bash
bash submit_02_parallel.sh
```

Tune `MAX_CONCURRENT` in `submit_02_parallel.sh` and the per-task `#SBATCH` resources in `02a-...bash` for your cluster. Because each sample writes only to its own `Merging/<id>`, `Bartender/<id>`, and `Processed_data/<id>` folders, the tasks are independent; only `02b` reads across samples.

**Logs:** SLURM stdout/stderr go to `./log/`. The serial script names them `<jobname>_<jobid>`; the array names them `<jobname>_<arrayjobid>_<taskid>`. (The serial script also runs a legacy `sed` "clean" step assuming a `slurm-<jobid>.out` file in the submission directory — adjust or ignore it if your cluster only writes `./log/...`.)

## Python scripts (`python_scripts/`)

| Script | Purpose |
|--------|---------|
| `single_guide_parsing.py` | Extract gRNA and clonal barcode reads into Bartender-style inputs. |
| `single_guide_aggregate_sgRNA.py` | Map sgRNA clusters to the reference; emit per–sgRNA barcode Bartender inputs. |
| `single_guide_aggregate_barcode.py` | Merge barcode clustering with sgRNA assignments for one sample. |
| `single_guide_aggregate_sample.py` | Concatenate all samples into `gRNA_clonalbarcode_combined.csv`. |

## Auxiliary notebooks (`auxiliary_code/`)

- **`A1-NGS_data_address_generation.ipynb`** — help generate `NGS_address`-compatible path lists.
- **`A2-NGS_reads_distribution_after_parsing.ipynb`** — summarize read distributions (optional QC).
