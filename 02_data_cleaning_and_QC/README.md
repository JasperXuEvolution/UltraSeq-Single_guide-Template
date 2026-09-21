# Data cleaning and QC (step 02)

Notebook-based quality control. This step takes the combined table from step 01
(`gRNA_clonalbarcode_combined.csv`), annotates it with guide and sample metadata,
converts read counts to cell numbers using spike-ins, runs sample-level QC, and
writes the processed **`.parquet`** table that step 03 (bootstrapping) consumes.

Run the notebooks in order from this directory (so relative `data/` paths resolve):

1. **`GS_BTHC-Batch1-QC_part1-251130.ipynb`** — raw data assembly + spike-in QC → `*_annotated_df.parquet`
2. **`GS_BTHC-Batch1-QC_part2-251206.ipynb`** — sample QC + contamination checks → final `.parquet` for bootstrapping

> The `GS_BTHC` naming is the worked example that flows through the whole
> repository (step 03 reads `GS_BTHC_KT_tumor.parquet`). Rename freely for your
> own project and update the `project_prefix` variable in each notebook.

## Layout

| Path | Role |
|------|------|
| `GS_BTHC-Batch1-QC_part1-251130.ipynb` | Load reads/guides/metadata, spike-in QC, compute cells-per-read, emit `annotated_df`. |
| `GS_BTHC-Batch1-QC_part2-251206.ipynb` | Top-N gene/tumor plots, per-sample QC, contamination analysis, emit final bootstrapping input. |
| `UltraSeq-SampleSpecificAnalysisFunction-241203.py` | Shared helper functions (tumor-size metrics, normalization, FDR). |
| `data/` | Inputs and generated tables (git-ignored except `.gitkeep`). |
| `figs/` | QC figures written by the notebooks (git-ignored). |

## Required inputs (place in `data/`)

| File | Description |
|------|-------------|
| `gRNA_clonalbarcode_combined.csv` | Output of step 01 (`01_data_collection/data/Processed_data/`). Columns: `gRNA`, `Clonal_barcode`, `Sample_ID`, `Frequency`. |
| `gRNA_information.csv` | Guide reference. Columns: `Gene`, `gRNA` — the **same schema** as step 01's `guide_reference-*.csv`, so the same file can be reused. (The notebook's `sgRNA→gRNA` rename is a defensive no-op when the column is already `gRNA`.) Spike-ins are detected by `"Spike"` appearing in the gene name. |
| `<project_prefix>_mice_info_standardized.csv` | Sample metadata. Columns: `Sample_ID`, `Mouse_Ear_Tag`, `Mouse_Genotype`, `Sex`, `Virus_Titer`, `Time_After_Tumor_Initiation(wks)`, `Total_Lung_Weight(g)`, `Pooling_library_name`, `Tissue_type`. (The notebook renames these to lower-case/underscore forms; the `'Sample ID'→Sample_ID` rename is a no-op when the column is already `Sample_ID`.) |
| `sgInert.csv` | List of inert (control) genes; column `Gene`. Used to label `Type = Inert`. |
| `KT_reference_tumor.parquet` | *(optional)* Reference cohort merged in at the end of part 2 to build `GS_BTHC_KT_tumor.parquet`. Drop this merge if you have no reference cohort. |

## Outputs (written to `data/`)

| File | Consumed by |
|------|-------------|
| `<project_prefix>_annotated_df.parquet` | part 2 of this step |
| `<project_prefix>_sample_summary_df.csv` | QC record |
| `<project_prefix>_final_df.parquet` | processed experimental table |
| `GS_BTHC_KT_tumor.parquet` | **step 03** (`03_bootstrapping/data/`) |

## Notes

- Edit the paths block near the top of each notebook (`parent_address`,
  `project_prefix`) before running.
- The example/illustration cells (e.g. pairwise sample comparisons) use
  placeholder `Sample_ID` values from the original dataset — replace them with
  IDs present in your own data.
- Spike-in names are hard-coded as `tuba-seq-v2_Spike-in-1/2/3`; adjust if your
  library differs.
