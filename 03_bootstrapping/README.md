# Bootstrapping (step 03)

This step runs statistical bootstrapping on processed tumor (and optionally plasmid) data using `python_scripts/TubaSeq_Ultra_Boostrapping.py`. SLURM wrapper scripts live in `bash/`.

## Layout

| Path | Purpose |
|------|---------|
| `data/` | Input parquet files (e.g. tumor table; plasmid table for scripts 03–04). |
| `results/` | Output root; CSV prefixes are `results/<BATCH_NAME>/…` (see scripts). |
| `python_scripts/` | Main Python entrypoint. |
| `bash/` | SLURM job scripts. |

## Prerequisites

1. **Repository `config.sh`** (two levels up from `bash/`): set `PROJECT_DIR` so `working_dir="$PROJECT_DIR/03_bootstrapping"` is correct.

2. **Conda**: scripts `source ~/miniconda3/etc/profile.d/conda.sh` and run `conda activate TubaSeq_Ultra`. Adjust if your Conda install lives elsewhere.

3. **SLURM**: edit `#SBATCH` in each script (`mail-user`, `account`, `partition`, memory, wall time) for your cluster.

## Job scripts (`bash/`)

| Script | Description |
|--------|-------------|
| `01-BT_Adaptive.sh` | Non-plasmid, adaptive (`--c Yes`). |
| `02-BT_Normal.sh` | Non-plasmid, normal (`--c No`). |
| `03-BT_Plamid_Adaptive.sh` | Plasmid mode, adaptive; adds `--p` plasmid parquet. |
| `04-BT_Plamid_Normal.sh` | Plasmid mode, normal. |

Default inputs under `data/`:

- Tumor: `GS_BTHC_KT_tumor.parquet`
- Plasmid (03–04): `GS_MTAP-Plasmid.parquet`

Change paths inside the scripts if your filenames differ.

## Batch name and outputs

Each script sets:

```bash
BATCH_NAME="${BATCH_NAME:-UltraSeq_example}"
out_prefix="$working_dir/results/$BATCH_NAME"
mkdir -p "$working_dir/results"
```

If `BATCH_NAME` is unset, outputs use the prefix `results/UltraSeq_example`. To use another run label (when your site exports the submission environment):

```bash
export BATCH_NAME=MyRun
sbatch bash/01-BT_Adaptive.sh
```

Or edit the default in the script. See `TubaSeq_Ultra_Boostrapping.py` for exact CSV suffixes written to `--o1` / `--o2`.

## Submitting

From the repository root:

```bash
sbatch 03_bootstrapping/bash/01-BT_Adaptive.sh
```

The scripts run `sacct` at the end when `SLURM_JOBID` is set (normal `sbatch` runs).
