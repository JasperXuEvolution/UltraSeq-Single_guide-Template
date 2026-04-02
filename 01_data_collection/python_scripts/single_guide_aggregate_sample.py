#!/usr/bin/env python
# coding: utf-8
"""
README / Documentation
----------------------
This script combines deduplicated sgRNA and clonal barcode data from multiple mouse samples into a single DataFrame.
It assumes that the specified output directory (provided via the --o argument) contains subdirectories for each mouse sample,
and that each subdirectory contains a file named 'Combined_deduplexed_df.csv'. The script aggregates these files,
concatenates them, and writes out a final CSV file containing the combined results.

Usage:
------
Run the script from the command line with the following argument:
    --o : Path to the directory that contains subdirectories for each mouse sample.

Example:
    python combine_across_mice.py --o path/to/output_directory/

Outputs:
--------
- gRNA_clonalbarcode_combined.csv — A CSV file containing the combined deduplicated data from all mouse samples.
"""

import os
import argparse
import pandas as pd

def main():
    """
    Execute the workflow to combine deduplicated sgRNA and clonal barcode data from multiple mouse samples.

    The function performs the following steps:
      1. Parse the command-line argument to obtain the base output directory.
      2. Scan the base directory for subdirectories corresponding to individual mouse samples.
      3. For each sample subdirectory, read the 'Combined_deduplexed_df.csv' file.
      4. Concatenate all sample DataFrames into a single DataFrame.
      5. Save the final combined DataFrame as 'gRNA_clonalbarcode_combined.csv' in the base directory.
    """
    parser = argparse.ArgumentParser(
        description='Combine deduplicated sgRNA and clonal barcode data from multiple mouse samples.'
    )
    parser.add_argument("--o", required=True, help="Path to the directory containing sample subdirectories")
    args = parser.parse_args()
    
    # The base directory that contains subdirectories for each mouse sample.
    base_output_dir = args.o
    
    # Find all subdirectories in the base output directory (each represents a mouse sample).
    sample_folder_paths = []
    for entry_name in os.listdir(base_output_dir):
        entry_path = os.path.join(base_output_dir, entry_name)
        if os.path.isdir(entry_path):
            sample_folder_paths.append(entry_path)
    
    # List to store DataFrames from each mouse sample.
    combined_df_list = []
    for sample_folder in sample_folder_paths:
        # The sample ID is assumed to be the name of the subdirectory.
        sample_id = os.path.basename(sample_folder)
        # Construct the path to the deduplicated file for this sample.
        dedup_file_path = os.path.join(base_output_dir, sample_id, 'Combined_deduplexed_df.csv')
        # Read the CSV file into a DataFrame.
        sample_df = pd.read_csv(dedup_file_path)
        combined_df_list.append(sample_df)
    
    # Concatenate all sample DataFrames into a single DataFrame.
    final_df = pd.concat(combined_df_list, ignore_index=True)
    
    # Save the combined DataFrame to a CSV file in the base directory.
    output_file = os.path.join(base_output_dir, 'gRNA_clonalbarcode_combined.csv')
    final_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    main()
