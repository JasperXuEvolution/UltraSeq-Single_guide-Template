#!/usr/bin/env python
# coding: utf-8

"""
README / Documentation
----------------------
This script extracts single guide RNA (gRNA) and clonal barcode sequences from a merged FASTQ.gz file.
The extracted sequences are formatted for compatibility with the `bartender` software, enabling further processing.

Usage:
------
Execute the script from the command line with the following required arguments:
    --a : Path to the input FASTQ.gz file.
    --o : Path to the output directory where the extracted sequences will be saved.

Example:
    python single_guide_parsing.py --a path/to/input.fastq.gz --o path/to/output_directory

Outputs:
--------
- `gRNA.bartender` — Contains the extracted gRNA sequences.
- `clonalbarcode.bartender` — Contains the extracted clonal barcode sequences.

Each output file is formatted with two columns:

    <extracted_sequence>,<read_ID>

Summary Report:
---------------
Upon completion, the script prints a summary in the following format:

    Sample <sample_ID> has <total_reads> reads passing QC.
    Among them, <extracted_reads> contain both barcode and sgRNA, accounting for <fraction> of the total reads.

Regular Expression Pattern:
---------------------------
The script employs the following regular expression pattern to identify target sequences:

    (TAGTT){e<2}(.{16})TATGG(.{16,21})GTT(TAAGA){e<2}

Explanation:
- `(TAGTT){e<2}` — Matches the prefix with up to 1 errors.
- `(.{16})` — Extracts the 16-bp clonal barcode.
- `TATGG` — Fixed sequence element.
- `(.{16,21})` — Extracts the 16-21 bp gRNA sequence.
- `GTT(TAAGA){e<2}` — Matches the suffix with up to 1 errors.
"""

import gzip
import regex
import argparse

def main():
    """
    Parses a merged FASTQ gzipped file to extract gRNA and clonal barcode sequences.
    Outputs extracted sequences in bartender-compatible format.
    """
    parser = argparse.ArgumentParser(description='Extract gRNA and clonal barcode from a merged FASTQ gz file.')
    parser.add_argument("--a", required=True, help="Input FASTQ gz file path")
    parser.add_argument("--o", required=True, help="Output directory path")
    args = parser.parse_args()
    
    fastqgz_input_address = args.a
    output_dir = args.o
    
    # Define output file paths
    gRNA_output_address = f"{output_dir}/gRNA.bartender"
    clonal_barcode_output_address = f"{output_dir}/clonalbarcode.bartender"
    
    # Open output files for writing
    file_a = open(gRNA_output_address, 'wt')
    file_b = open(clonal_barcode_output_address, 'wt')
    
    temp_total_read, temp_extracted_read = 0, 0
    temp_sample_ID = output_dir.rstrip('/').split('/')[-1]  # Extract sample ID from output directory
    
    with gzip.open(fastqgz_input_address, 'rt') as handler:
        while True:
            temp_readID = handler.readline().rstrip()  # Read ID
            temp_sequence = handler.readline().rstrip()  # Sequence
            handler.readline()  # Skip plus line
            handler.readline()  # Skip quality score line
            
            if not temp_readID:
                break  # Exit loop if end of file is reached
            
            temp_total_read += 1
            
            # Regular expression pattern to extract barcode and gRNA
            temp_pattern = regex.compile(r'(TAGTT){e<2}(.{16})TATGG(.{16,21})GTT(TAAGA){e<2}')
            temp_search_result = temp_pattern.search(temp_sequence)
            
            if temp_search_result:
                temp_extracted_read += 1
                temp_gRNA = temp_search_result.group(3)
                temp_clonal_barcode = temp_search_result.group(2)
                
                # Write extracted sequences to output files
                file_a.write(f"{temp_gRNA},{temp_readID}\n")
                file_b.write(f"{temp_clonal_barcode},{temp_readID}\n")
    
    # Close output files
    file_a.close()
    file_b.close()
    
    # Print summary statistics
    print(f"Sample {temp_sample_ID} has {temp_total_read} reads passing QC. Among them, {temp_extracted_read} contain barcode and sgRNA, accounting for {temp_extracted_read/temp_total_read:.3f} of total reads.")
    
if __name__ == "__main__":
    main()
