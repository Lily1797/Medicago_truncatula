import csv
from collections import defaultdict

def load_gene_families(gene_family_file):
    family_dict = {}
    with open(gene_family_file, 'r') as f:
        for line in f:
            genes = line.strip().split('\t')
            for gene in genes:
                family_dict[gene] = genes  
    return family_dict

# Step 2: Find tandem array genes allowing at most one spacer gene
def find_tandem_array_genes(blast_results, family_dict):
    tandem_arrays = []
    chrom_genes = defaultdict(list)
    added_genes = set()

    # Organize genes by chromosome and avoid duplicates
    for row in blast_results:
        query_id, subject_id = row[20], row[25]
        query_start, query_end = int(row[16]), int(row[17])
        subject_start, subject_end = int(row[21]), int(row[22])
        chrom = row[18]
        strand = row[19]

        # Add both query and subject if not already added
        if (query_id, chrom) not in added_genes:
            chrom_genes[chrom].append((query_start, query_end, query_id, strand))
            added_genes.add((query_id, chrom))
        if (subject_id, chrom) not in added_genes:
            chrom_genes[chrom].append((subject_start, subject_end, subject_id, strand))
            added_genes.add((subject_id, chrom))

    # Process each chromosome
    for chrom, genes in chrom_genes.items():
        genes.sort(key=lambda x: x[0])  # Sort by start position
        used_genes = set()
        
        for i, (start, end, gene_id, strand) in enumerate(genes):
            if gene_id in used_genes:
                continue

            # Initialize a new TAG with the current gene
            current_tag = [(gene_id, start, end, strand)]
            used_genes.add(gene_id)
            spacer_found = False

            # Check for adjacent homologous genes within the same family, allowing at most one spacer
            for j in range(i + 1, len(genes)):
                next_start, next_end, next_gene_id, next_strand = genes[j]
                if next_gene_id in used_genes:
                    continue

                # Check if next gene belongs to the same family and on the same strand
                if next_gene_id in family_dict.get(gene_id, []): # and next_strand == strand:
                    current_tag.append((next_gene_id, next_start, next_end, next_strand))
                    used_genes.add(next_gene_id)
                    spacer_found = False  # Reset spacer check as we found a homologous gene
                    #end = next_end  # Update TAG boundary
                else:
                    # Encounter a potential spacer
                    if spacer_found:
                        break  # Exit if a second spacer is encountered
                    spacer_found = True

            # Add the completed TAG if it has more than one gene
            if len(current_tag) > 1:
                tandem_arrays.append((chrom, current_tag))

    return tandem_arrays

# Read the BLAST results file without headers
def read_blast_results(file_path):
    blast_results = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            blast_results.append(row)
    return blast_results

# Write the TAGs to output file
def write_tags_to_file(tags, output_file_path):
    with open(output_file_path, 'w') as out_file:
        for idx, (chrom, tag) in enumerate(tags, 1):
            genes_info = ','.join([f"{gene_id}:{start}:{end}:{strand}" for gene_id, start, end, strand in tag])
            out_file.write(f"TAG{idx} Chromosome:{chrom} {genes_info}\n")
    print(f"Tandem array genes (TAGs) with positions have been saved to {output_file_path}")

# Main function
def main(file_path, gene_family_file, output_file_path):
    blast_results = read_blast_results(file_path)
    family_dict = load_gene_families(gene_family_file) 
    tags = find_tandem_array_genes(blast_results, family_dict)
    write_tags_to_file(tags, output_file_path)

# Usage
low_path = 'same_chromosomes_low.txt'  # Path to input file
gene_family_low = 'MCL_low.tabular'
output_low = 'TAGs_low.txt'  # Path to the output file
high_path = 'same_chromosomes_high.txt'  # Path to input file
gene_family_high = 'MCL_high.tabular'
output_high = 'TAGs_high.txt'  # Path to the output file

main(low_path, gene_family_low, output_low)
main(high_path, gene_family_high, output_high)