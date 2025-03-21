#!/usr/bin/env python
# coding: utf-8
"""
README / Documentation
----------------------
This script maps sgRNA clustering results to a reference and separates clonal barcode data based on sgRNA.
It merges outputs from bartender analyses (barcode and cluster CSV files), verifies the sgRNA reference library
using a minimal Hamming distance criterion, and generates filtered data frames along with bartender input files
for each sgRNA cluster and associated clonal_barcode input files.

The primary purpose is to identify the true sgRNA center for each sgRNA. Many sgRNAs within a short Hamming distance
from the center may result from PCR errors. Once the true center is determined, each sgRNA is assigned to its cluster,
and the corresponding clonal barcodes are grouped.

I only retained sgRNAs whose cluster center is in the sgRNA reference file.

Usage:
------
Run the script from the command line with the following arguments:
    --a1 : Input file (bartender output barcode.csv for sgRNA)
    --a2 : Input file (bartender output cluster.csv for sgRNA)
    --a3 : Input file for sgRNA reference
    --a4 : (Optional) Minimal acceptable Hamming distance between sgRNAs 
           (e.g., to tolerate up to 2 PCR errors, set minimal distance to 4)
    --a5 : Input file (output from single_guide_parsing.py for sgRNA)
    --a6 : Input file (output from single_guide_parsing.py for clonal_barcode)
    --o  : Output file prefix

Example:
    python script.py --a1 path/to/barcode.csv --a2 path/to/cluster.csv --a3 path/to/sgRNA_reference.csv \
                     --a4 4 --a5 path/to/sgRNA_bartender.csv --a6 path/to/barcode_bartender.csv --o output_prefix_

Outputs:
--------
- An intermediate CSV file containing the filtered sgRNA and clonal barcode data.
- Bartender input files for each sgRNA cluster, saved in a subdirectory with filenames based on the sgRNA center.
- A text file listing the paths to each generated bartender file.
"""

import pandas as pd
from itertools import combinations
import sys
import argparse

def merge_bartender_output(input_barcode_address, input_cluster_address):
    """
    Merge two CSV files (barcode and cluster outputs) on the 'Cluster.ID' column.

    Parameters:
        input_barcode_address (str): File path for the barcode CSV file.
        input_cluster_address (str): File path for the cluster CSV file.

    Returns:
        pd.DataFrame: Merged DataFrame resulting from an inner join on 'Cluster.ID'.
    """
    temp_df1 = pd.read_csv(input_barcode_address)
    temp_df2 = pd.read_csv(input_cluster_address)
    temp_merge = pd.merge(
        temp_df1, temp_df2, how='inner', on=['Cluster.ID'],
        left_index=False, right_index=False, sort=True, copy=True, indicator=False,
        validate=None
    )
    return temp_merge

def unique_read_to_cluster_dic(input_df):
    """
    Create a dictionary mapping each unique read to its corresponding sgRNA cluster center.

    Parameters:
        input_df (pd.DataFrame): Merged DataFrame from the bartender output files.

    Returns:
        dict: Dictionary with keys as 'Unique.reads' and values as 'Center'.
    """
    temp_UniqueRead_To_ClusterSeq_dic = {}
    for x, y in zip(input_df['Unique.reads'].to_list(), input_df['Center'].to_list()):
        temp_UniqueRead_To_ClusterSeq_dic[x] = y
    return temp_UniqueRead_To_ClusterSeq_dic

def hamming_distance(s1, s2):
    """
    Compute the Hamming distance between two strings of equal length.

    Parameters:
        s1 (str): First string.
        s2 (str): Second string.

    Returns:
        int: The Hamming distance (the number of positions at which the corresponding characters differ).
    """
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def All_Pairwise_Hamming_Distannce_from_df(input_df):
    """
    Calculate the Hamming distance for all pairwise combinations of sgRNA sequences.

    Parameters:
        input_df (pd.DataFrame): DataFrame containing a column 'gRNA' with sgRNA sequences.

    Returns:
        dict: Dictionary where keys are tuple pairs of sgRNA sequences and values are their Hamming distances.
    """
    # Generate all pairwise combinations of sgRNA sequences
    temp1 = list(combinations(input_df.gRNA, 2))
    temp_dict = {}
    for x in temp1:
        temp_dict[x] = hamming_distance(x[0], x[1])
    return temp_dict

def Check_sgRNA_Library_Distance(input_sgRNA_df, input_minimal_distance):
    """
    Verify that all sgRNAs in the reference library are sufficiently distinct based on a specified minimal Hamming distance.

    Parameters:
        input_sgRNA_df (pd.DataFrame): DataFrame containing sgRNA sequences and their lengths.
        input_minimal_distance (int): Minimal acceptable Hamming distance between any two sgRNAs.

    Returns:
        dict: Dictionary grouping sgRNA pairs (by their length) that have Hamming distances less than or equal to the minimal threshold.
              An empty sub-dictionary indicates that all pairs meet the minimal distance criterion.
    """
    # Group sgRNAs by their length
    df_groups = input_sgRNA_df.groupby("gRNA_length")
    top_dict = {}  # Keys: sgRNA lengths; Values: sub-dictionaries of sgRNA pairs and their distances
    for key, sub_df in df_groups:
        temp_dict = All_Pairwise_Hamming_Distannce_from_df(sub_df)
        temp_filtered_dict = {}
        for temp_key, value in temp_dict.items():
            if value <= input_minimal_distance:
                temp_filtered_dict[temp_key] = value
        top_dict[key] = temp_filtered_dict
    return top_dict

def Generate_Filtered_df(input_sgRNA_bartender_address1, input_barcode_bartender_address, input_sgRNA_dic):
    """
    Generate a filtered DataFrame mapping sgRNA sequences (that match the reference) to their corresponding clonal barcodes and read IDs.
    I only retained sgRNAs whose cluster center is in the sgRNA reference file 

    The function processes two bartender output files (one for sgRNA and one for clonal barcodes) line by line,
    matches the sgRNA sequence to the reference dictionary, and constructs a DataFrame with the relevant data.

    Parameters:
        input_sgRNA_bartender_address1 (str): File path for the sgRNA bartender output.
        input_barcode_bartender_address (str): File path for the clonal barcode bartender output.
        input_sgRNA_dic (dict): Dictionary mapping sgRNA sequences (Unique.reads) to their cluster center.

    Returns:
        pd.DataFrame: DataFrame containing the following columns:
                      'gRNA', 'gRNA_center', 'Clonal_barcode', 'Read_ID', 'Sample_ID'.
    """
    # Extract sample ID from the file path (assumes the sample ID is the name of the parent directory)
    temp_sampleID = input_sgRNA_bartender_address1.split('/')[-2]
    
    with open(input_sgRNA_bartender_address1, 'r') as handler1, open(input_barcode_bartender_address, 'r') as handler2:
        # Read header lines from both files
        temp1 = handler1.readline().strip()            
        temp2 = handler2.readline().strip()
        # Lists to store the matched data
        temp_sgRNA_list, temp_sgRNA_center_list, temp_barcode_list, temp_readsID_list = [], [], [], []
        # Process both files line by line simultaneously
        while bool(temp1) & bool(temp2):
            temp_gRNA = temp1.split(',')[0]  # Extract sgRNA sequence from the sgRNA file
            # Check if the sgRNA exists in the reference dictionary
            if temp_gRNA in input_sgRNA_dic.keys():
                temp_sgRNA_center = input_sgRNA_dic.get(temp_gRNA)
                temp_readsID = temp1.split(',')[1].split(' ')[0]  # Extract Read_ID from the sgRNA file
                temp_barcode = temp2.split(',')[0]  # Extract barcode from the barcode file
                temp_sgRNA_list.append(temp_gRNA)
                temp_sgRNA_center_list.append(temp_sgRNA_center)
                temp_barcode_list.append(temp_barcode)
                temp_readsID_list.append(temp_readsID)
            # Read next lines from both files
            temp1 = handler1.readline().strip()            
            temp2 = handler2.readline().strip()
    # Create and return the DataFrame from the collected data
    temp_output_df = pd.DataFrame({
        'gRNA': temp_sgRNA_list,
        'gRNA_center': temp_sgRNA_center_list,
        'Clonal_barcode': temp_barcode_list,
        'Read_ID': temp_readsID_list,
        'Sample_ID': temp_sampleID
    })
    return temp_output_df

def main():
    """
    Execute the workflow:
      - Parse command-line arguments.
      - Merge the sgRNA bartender outputs.
      - Read and verify the sgRNA reference (if a minimal distance is provided).
      - Filter the merged DataFrame to include only sgRNAs matching the reference.
      - Construct a mapping dictionary and generate a filtered DataFrame linking sgRNAs to clonal barcodes and read IDs.
      - Output the filtered DataFrame and create bartender input files grouped by sgRNA center.
    """
    parser = argparse.ArgumentParser(
        description='Map sgRNA clustering results to a reference and separate clonal barcode data based on sgRNA.'
    )
    parser.add_argument("--a1", required=True, help="Input file: bartender output barcode.csv for sgRNA")
    parser.add_argument("--a2", required=True, help="Input file: bartender output cluster.csv for sgRNA")
    parser.add_argument("--a3", required=True, help="Input file for sgRNA reference")
    parser.add_argument("--a4", required=False, help="Expected minimal Hamming distance between sgRNAs")
    parser.add_argument("--a5", required=True, help="Input file: bartender output for sgRNA from single_guide_parsing.py")
    parser.add_argument("--a6", required=True, help="Input file: bartender output for clonal_barcode from single_guide_parsing.py")
    parser.add_argument("--o", required=True, help="Output file prefix")
    args = parser.parse_args()
    
    temp_output_prefix = args.o
    
    # Merge the sgRNA bartender outputs into a single DataFrame
    df1 = merge_bartender_output(args.a1, args.a2)  # Clustered sgRNA DataFrame
    
    # Read the sgRNA reference and compute sequence lengths
    gRNA_df = pd.read_csv(args.a3)
    gRNA_df['gRNA_length'] = gRNA_df['gRNA'].apply(lambda x: len(x))
    gRNA_df = gRNA_df.drop_duplicates(keep='first')  # Remove duplicate sgRNAs
    
    # If a minimal distance is provided, verify the sgRNA reference quality
    if args.a4 is None:
        sample_to_exclude = []
        print("The sgRNA reference distance is not verified")
    else:
        # Check distances among reference sgRNAs
        sgRNA_distance_dict = Check_sgRNA_Library_Distance(gRNA_df, int(args.a4))
        for key, value in sgRNA_distance_dict.items():
            if value != {}:
                print('sgRNA distance error')
                sys.exit("Error message")
        print("sgRNA are good")
    
    # Update the 'Center' column in df1:
    # If 'Unique.reads' is in the reference, assign it directly; otherwise, retain the original 'Center' value.
    df1['Center'] = df1.apply(
        lambda x: x['Unique.reads'] if (x['Unique.reads'] in gRNA_df.gRNA.values) else x['Center'],
        axis=1
    )
    # Filter df1 to include only sgRNAs present in the reference
    df1_final = df1[df1.Center.isin(gRNA_df.gRNA)].copy()
    
    # Create a dictionary mapping unique reads to their corresponding cluster centers
    gRNA_ref_dict = unique_read_to_cluster_dic(df1_final)
    
    # Generate the filtered DataFrame with matched sgRNA and clonal barcode data
    Final_df = Generate_Filtered_df(args.a5, args.a6, gRNA_ref_dict)
    
    # Extract sample ID from the filtered DataFrame (assumes consistency across rows)
    temp_sample_ID = Final_df.Sample_ID.iloc[0]
    
    # Calculate total reads from the original and filtered DataFrames
    count1 = df1.Frequency.sum()
    count2 = df1_final.Frequency.sum()
    print(f"Sample {temp_sample_ID:s} has totally {count1:d} reads for sgRNA. "
          f"Among them {count2:d} match the reference sgRNA, which is about {count2/count1:.3f}")
    
    # Output the intermediate DataFrame to CSV
    Final_df.to_csv(temp_output_prefix + 'Intermediate_df.csv', index=False)
    
    # Group the filtered data by sgRNA center and prepare bartender input files for clonal barcodes
    sgRNA_groups = Final_df.groupby("gRNA_center")
    # Write the paths to the generated bartender files to an address file
    temp_address_file = temp_output_prefix + 'Bartender_input_address'
    with open(temp_address_file, 'w') as file_a:
        for groups in sgRNA_groups:
            g, value = groups
            temp_name = temp_output_prefix + 'Clonal_barcode/' + g + '.bartender'
            file_a.write(temp_name + '\n')
            # Write each group's clonal barcode and Read_ID data to its respective file
            value[['Clonal_barcode', 'Read_ID']].to_csv(temp_name, sep=',', header=False, index=False)
    
if __name__ == "__main__":
    main()
