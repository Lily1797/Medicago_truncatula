import pandas as pd

# Load the BLAST results
blast_results_file = 'output_blast.txt'
blast_results = pd.read_csv(blast_results_file, sep='\t', header=None)

# Load the position file
position_file = 'positions.txt'
positions = pd.read_csv(position_file, sep='\t', header=None)

# Rename columns for clarity
blast_results.columns = ['query_id', 'subject_id', 'identity', 'match', 'mismatch', 'gap', 'query_start', 'query_end', 'subject_start', 'subject_end', 'evalue', 'bitscore', 'query_length', 'subject_length', 'query_coverage', 'subject_coverage']

positions.columns = ['protein_id', 'start', 'end', 'chrom', 'strand']

# Merge query information
query_info = positions.rename(columns={'protein_id': 'query_id', 'start': 'query_start_chr', 'end': 'query_end_chr', 'chrom': 'query_chrom', 'strand': 'query_strand'})
add_query = blast_results.merge(query_info, on='query_id', how='inner')

# Merge subject information
subject_info = positions.rename(columns={'protein_id': 'subject_id', 'start': 'subject_start_chr', 'end': 'subject_end_chr', 'chrom': 'subject_chrom', 'strand': 'subject_strand'})
update_blast_results = add_query.merge(subject_info, on='subject_id', how='inner')

# Select relevant columns for the final output
final_columns = ['query_id', 'subject_id', 'identity', 'match', 'mismatch', 'gap', 'query_start', 'query_end', 'subject_start', 'subject_end', 'evalue', 'bitscore', 'query_length', 'subject_length', 'query_coverage', 'subject_coverage', 'query_start_chr', 'query_end_chr', 'query_chrom', 'query_strand', 'subject_start_chr', 'subject_end_chr', 'subject_chrom', 'subject_strand']
final_results = update_blast_results[final_columns]

# Save to a new file
final_results.to_csv('update_blast_results.txt', sep='\t', index=False, header=False)
