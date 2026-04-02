#!/usr/bin/env python
# coding: utf-8
"""
README / Documentation
----------------------
This script combines clonal barcode data from multiple bartender output files for a single mouse sample.
It recursively searches for clonal barcode files within the specified sample folder, merges them with the 
associated sgRNA cluster data, and then deduplicates the merged results. The final outputs include a raw 
combined DataFrame with frequency counts and a deduplicated DataFrame containing only perfect sgRNA matches.

Usage:
------
Run the script from the command line with the following arguments:
    --a : Path to the input sample folder.
    --o : Output file prefix (directory or path where results will be saved).

Example:
    python combine_barcodes.py --a path/to/sample_folder --o path/to/output_prefix

Outputs:
--------
- Combined_ND_df.csv         — Raw combined data with frequency counts.
- Combined_deduplexed_df.csv — Deduplicated data with perfect sgRNA matches.
"""

import pandas as pd
import argparse
import glob
import os

def Combine_sgRNA_barcode_from_the_Same_mouse(input_folder_address):
    """
    Combine clonal barcode data from multiple bartender output files for a single mouse sample.
    
    This function recursively searches for clonal barcode files within the specified folder, merges 
    the data from corresponding barcode and cluster files using a helper function, and then deduplicates 
    the merged results by grouping on key columns. Finally, it merges the combined data with the 
    intermediate sgRNA reference DataFrame.
    
    Parameters:
        input_folder_address (str): Path to the folder containing the sample data.
        
    Returns:
        tuple: A tuple containing:
            - temp_final_raw (pd.DataFrame): Raw combined data with frequency counts.
            - temp_final_completely_deduplexed (pd.DataFrame): Deduplicated data containing only perfect sgRNA matches.
    """
    # Pattern to recursively search for files ending with '_cluster.csv' in Clonal_barcode directories.
    temp_pattern = '/**/Clonal_barcode/*_cluster.csv'
    barcode_address_list = glob.glob(input_folder_address + temp_pattern, recursive=True)
    
    # Path to the intermediate sgRNA reference file.
    temp_ref_address = input_folder_address + '/Intermediate_df.csv'
    temp_bc_df_list = []  # List to store DataFrames for each sgRNA file.
    
    for temp_a in barcode_address_list:
        # Generate file paths:
        # t1: Path for the barcode CSV file (by replacing 'cluster' with 'barcode').
        t1 = temp_a.replace('cluster', 'barcode')
        # t2: Path for the cluster CSV file (unchanged).
        t2 = temp_a
        # t3: Path for the bartender input file (replace '_cluster.csv' with '.bartender').
        t3 = temp_a.replace('_cluster.csv', '.bartender')
        
        temp_df = merge_barcode_and_sgRNA_output(t1, t2, t3)
        temp_bc_df_list.append(temp_df)
    
    # Concatenate all individual DataFrames into a single DataFrame.
    temp_total_bc_df = pd.concat(temp_bc_df_list).reset_index(drop=True)
    temp_sgRNA_ref_df = pd.read_csv(temp_ref_address)
    
    # Merge the total barcode DataFrame with the sgRNA reference DataFrame using 'Read_ID' and 'Clonal_barcode' as keys.
    temp_final = temp_total_bc_df.merge(temp_sgRNA_ref_df, on=['Read_ID', 'Clonal_barcode'])
    
    # Deduplicate the merged DataFrame by grouping on key columns.
    temp_name_list = ['gRNA_center', 'Clonal_barcode_center', 'gRNA', 'Clonal_barcode', 'Sample_ID']
    temp_final_raw = temp_final.groupby(temp_name_list, as_index=False)['Read_ID'].count()
    temp_final_raw.rename(columns={'Read_ID': 'Frequency'}, inplace=True)
    
    # Filter for rows where the sgRNA center perfectly matches the sgRNA sequence.
    temp_final_true_sgRNA = temp_final[temp_final.gRNA_center == temp_final.gRNA]
    
    # Further deduplicate based on perfect sgRNA matches.
    temp_final_completely_deduplexed = temp_final_true_sgRNA.groupby(
        ['gRNA_center', 'Clonal_barcode_center', 'Sample_ID'], as_index=False
    )['Read_ID'].count()
    temp_final_completely_deduplexed.rename(
        columns={'Read_ID': 'Frequency', 'gRNA_center': 'gRNA', 'Clonal_barcode_center': 'Clonal_barcode'},
        inplace=True
    )
    
    return (temp_final_raw, temp_final_completely_deduplexed)

def merge_barcode_and_sgRNA_output(input_barcode_address, input_cluster_address, bartender_input):
    """
    Merge clonal barcode and sgRNA cluster output files along with bartender input data.
    
    This function reads three input files:
      - The barcode file (dropping the 'Frequency' column).
      - The cluster file (dropping 'Cluster.Score' and 'time_point_1' columns).
      - The bartender input file (with specified column names).
    The function merges the barcode and cluster data on 'Cluster.ID', renames certain columns,
    and then merges the result with the bartender input DataFrame on 'Clonal_barcode'.
    
    Parameters:
        input_barcode_address (str): Path to the barcode CSV file.
        input_cluster_address (str): Path to the cluster CSV file.
        bartender_input (str): Path to the bartender input file.
        
    Returns:
        pd.DataFrame: Merged DataFrame containing combined clonal barcode data.
    """
    temp_df1 = pd.read_csv(input_barcode_address).drop(columns=['Frequency'])
    temp_df2 = pd.read_csv(input_cluster_address).drop(columns=['Cluster.Score', 'time_point_1'])
    temp_df3 = pd.read_csv(bartender_input, sep=',', names=['Clonal_barcode', 'Read_ID'], header=None)
    
    temp_merge = pd.merge(
        temp_df1, temp_df2, how='inner', on=['Cluster.ID'],
        left_index=False, right_index=False, sort=True, copy=True, indicator=False, validate=None
    )
    
    temp_merge.rename(
        columns={'Unique.reads': 'Clonal_barcode', 'Center': 'Clonal_barcode_center'},
        inplace=True
    )
    
    temp_merge = temp_merge.drop(columns=['Cluster.ID']).merge(temp_df3, on='Clonal_barcode', how='right')
    return temp_merge

def main():
    """
    Execute the workflow to combine sgRNA and clonal barcode information.
    
    The function parses command-line arguments to obtain the input sample folder and the output prefix,
    then processes the data using the Combine_sgRNA_barcode_from_the_Same_mouse function. The resulting 
    DataFrames are saved as CSV files:
      - Combined_ND_df.csv         — Raw combined data with frequency counts.
      - Combined_deduplexed_df.csv — Deduplicated data with perfect sgRNA matches.
    """
    parser = argparse.ArgumentParser(description='Combine sgRNA and clonal barcode information for a sample.')
    parser.add_argument("--a", required=True, help="Path to the input sample folder")
    parser.add_argument("--o", required=True, help="Output file prefix (directory/path for saving results)")
    args = parser.parse_args()
    
    temp_output_prefix = args.o
    temp_folder = args.a

    # Generate output file names.
    temp_o1 = os.path.join(temp_output_prefix, 'Combined_ND_df.csv')
    temp_o2 = os.path.join(temp_output_prefix, 'Combined_deduplexed_df.csv')
    
    # Combine the data for the sample.
    temp_df1, temp_df2 = Combine_sgRNA_barcode_from_the_Same_mouse(temp_folder)
    
    # Save the resulting DataFrames to CSV files.
    temp_df1.to_csv(temp_o1, index=False)
    temp_df2.to_csv(temp_o2, index=False)

if __name__ == "__main__":
    main()
