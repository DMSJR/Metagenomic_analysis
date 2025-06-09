import pandas as pd
import os
import argparse

def process_diamond_results(normalized_cds_file, diamond_file, output_file):
    # Extract the gene name from the Diamond file name
    gene_name = os.path.basename(diamond_file).split('_')[1].split('.')[0]

    # Load the normalized CDS data
    cds_df = pd.read_csv(normalized_cds_file, sep='\t')

    # Check if the file has the required columns
    required_columns = {'CDS', 'Raw_Count', 'Normalized_Count'}
    if not required_columns.issubset(cds_df.columns):
        raise ValueError(f"The input file '{normalized_cds_file}' must contain the following columns: {required_columns}")

    # Load the Diamond results
    diamond_df = pd.read_csv(diamond_file, sep='\t', header=None,
                             names=['qseqid', 'qlen', 'sallseqid', 'slen', 'evalue', 'length', 'piden'])

    # Calculate coverage for each alignment
    diamond_df['coverage'] = diamond_df['length'] / diamond_df['slen']

    # Filter by identity and coverage thresholds
    filtered_diamond = diamond_df[(diamond_df['piden'] >= 40) & (diamond_df['coverage'] >= 0.4)]

    # Get the CDS IDs that meet the criteria
    valid_cds = set(filtered_diamond['qseqid'])

    # Filter the CDSs present in the normalized CDS file and that passed the filters
    cds_df = cds_df[cds_df['CDS'].isin(valid_cds)]

    # Add the gene name as a new column
    cds_df['Gene'] = gene_name

    # Group by gene and sum normalized values
    gene_reads = cds_df.groupby('Gene')['Normalized_Count'].sum().reset_index()

    # Rename columns for the final output
    gene_reads.columns = ['Gene', 'TotalNormalizedReads']

    # Save the results to a file
    gene_reads.to_csv(output_file, sep='\t', index=False)

    print(f"Table generated: {output_file}")
    print(f"Processed gene: {gene_name}")

if __name__ == "__main__":
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Process Diamond results and normalized read counts.")
    parser.add_argument("normalized_cds_file", help="File with normalized read counts per CDS.")
    parser.add_argument("diamond_file", help="Diamond result file.")
    parser.add_argument("output_file", help="Output file for results grouped by gene.")

    args = parser.parse_args()

    # Call the main function
    process_diamond_results(args.normalized_cds_file, args.diamond_file, args.output_file)
