from itertools import combinations

def extract_dup(input_file, blast_file, output_file):
    """
    Extracts gene pairs from an MCL result file and filters BLAST results
    to retain only matches involving these pairs.

    Parameters:
        input_file (str): Path to the MCL result file where each line represents a gene family.
        blast_file (str): Path to the BLAST result file.
        output_file (str): Path to the output file for filtered BLAST results.
    """
    gene_pairs = set()

    # Step 1: Extract gene pairs
    with open(input_file, 'r') as infile:
        for line in infile:
            genes = line.strip().split()
            gene_pairs.update(tuple(sorted(pair)) for pair in combinations(genes, 2))

    # Step 2: Filter BLAST results
    with open(blast_file, 'r') as blast_in, open(output_file, 'w') as blast_out:
        for line in blast_in:
            fields = line.strip().split()
            gene1, gene2 = fields[20], fields[25]
            if tuple(sorted([gene1, gene2])) in gene_pairs:
                blast_out.write(line)

# Example usage
input_low = "MCL_low.tabular" 
input_high = "MCL_high.tabular"      
blast_low = "homolog_low_unique.txt" 
blast_high = "homolog_high_unique.txt"  
output_low = "dup_pairs_low.txt"  
output_high = "dup_pairs_high.txt"  

extract_dup(input_low, blast_low, output_low)
extract_dup(input_high, blast_high, output_high)