import matplotlib.pyplot as plt
from collections import Counter

# Load TAG data from file
def parse_tags_file(filename):
    tags = []

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            
            if line:  # Process non-empty lines
                try:
                    # Split the line to extract the TAG ID, Chromosome part, and Genes part
                    parts = line.split(" ")
                    tag_id = parts[0]  # Extract TAG ID
                    chromosome = parts[1].split(":")[1]  # Extract chromosome number after "Chromosome:"

                    # Parse genes and their positions
                    genes = []
                    gene_entries = parts[2].split(",")
                    for entry in gene_entries:
                        gene_info = entry.split(":")
                        gene_id = gene_info[0]
                        start = int(gene_info[1])
                        end = int(gene_info[2])
                        #strand = int(gene_info[3])
                        genes.append({
                            "gene_id": gene_id,
                            "start": start,
                            "end": end
                            #"strand": strand
                        })

                    # Store each TAG information in a dictionary
                    tags.append({
                        "tag_id": tag_id,
                        "chromosome": chromosome,
                        "genes": genes
                    })

                except (IndexError, ValueError) as e:
                    print(f"Error parsing line: {line}\nError: {e}")

    return tags

# Define the function to plot gene distributions in two TAG categories
def plot_all_gene_distribution(tags1, tags2):
    # Count the number of genes in each TAG
    tag_gene_counts1 = [len(tag["genes"]) for tag in tags1]
    tag_gene_counts2 = [len(tag["genes"]) for tag in tags2]

    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6))  # Adjusted figsize for better aspect ratio

    # Plot histogram for the first TAG category (low stringency)
    ax1.hist(tag_gene_counts1, bins=range(1, max(tag_gene_counts1) + 2), align='left', edgecolor='black', color='skyblue')
    ax1.set_xlabel("Number of Genes in TAG")
    ax1.set_ylabel("Frequency")
    ax1.set_title("Low Stringency")
    ax1.grid(axis='y', linestyle='--', alpha=0.7)
    ax1.set_xticks(range(1, max(tag_gene_counts1) + 1))
    ax1.tick_params(axis='x', rotation=45)  # Rotate x-tick labels
    ax1.set_ylim(0, 2800)

    # Plot histogram for the second TAG category (high stringency)
    ax2.hist(tag_gene_counts2, bins=range(1, max(tag_gene_counts2) + 2), align='left', edgecolor='black', color='salmon')
    ax2.set_xlabel("Number of Genes in TAG")
    ax2.set_title("High Stringency")
    ax2.grid(axis='y', linestyle='--', alpha=0.7)
    ax2.set_xticks(range(1, max(tag_gene_counts2) + 1))
    ax2.tick_params(axis='x', rotation=45)  # Rotate x-tick labels
    ax2.set_ylim(0, 2800)

    # Adjust layout and show plot
    plt.tight_layout()
    plt.savefig("gene_distribution.png", format='png', dpi=600)
    plt.show()


# Define the function to plot gene distributions in two TAG categories
def plot_tags_per_chromosome(tags1, tags2):
    # Count the number of TAGs for each chromosome
    chrom_counts1 = Counter(tag["chromosome"] for tag in tags1)
    chrom_counts2 = Counter(tag["chromosome"] for tag in tags2)
    
    # Sort chromosomes by chromosome number for logical ordering
    sorted_chroms1 = sorted(chrom_counts1.keys(), key=lambda x: int(x))
    tag_counts_per_chrom1 = [chrom_counts1[chrom] for chrom in sorted_chroms1]
    sorted_chroms2 = sorted(chrom_counts2.keys(), key=lambda x: int(x))
    tag_counts_per_chrom2 = [chrom_counts2[chrom] for chrom in sorted_chroms2]

    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6)) 
    
    # Plot bar chart for low stringency
    ax1.bar(sorted_chroms1, tag_counts_per_chrom1, color='skyblue', edgecolor='black')
    ax1.set_xlabel("Chromosome")
    ax1.set_ylabel("Number of TAGs")
    ax1.set_title("Low Stringency")
    ax1.tick_params(axis='x', rotation=45)  # Rotate x-tick labels
    ax1.grid(axis='y', linestyle='--', alpha=0.7)
    ax1.set_ylim(0, 650)
    
    # Plot bar chart for high stringency
    ax2.bar(sorted_chroms2, tag_counts_per_chrom2, color='salmon', edgecolor='black')
    ax2.set_xlabel("Chromosome")
    ax2.set_title("High Stringency")
    ax2.tick_params(axis='x', rotation=45)  # Rotate x-tick labels
    ax2.grid(axis='y', linestyle='--', alpha=0.7)   
    ax2.set_ylim(0, 650)
    
    # Adjust layout and show plot
    plt.tight_layout()
    plt.savefig("tag_per_chrom.png", format='png', dpi=600)
    plt.show()

#Ideogram
chromosome_lengths = {
    '1': 56706830,
    '2': 51972579,
    '3': 58931556,
    '4': 64763011,
    '5': 44819618,
    '6': 42866092,
    '7': 56236587,
    '8': 49719271,
}

# Duplicated Genes Inside and Outside TAGs in Low vs High Stringency Datasets
datasets = ['Low Stringency', 'High Stringency']
tag_genes = [13003, 12221]  # Duplicated genes inside TAGs
non_tag_genes = [43300 - 13003, 38853 - 12221]  # Duplicated genes outside TAGs

fig, ax = plt.subplots(figsize=(8, 6))

bars1 = ax.bar(datasets, non_tag_genes, label='Duplicated Genes Outside TAGs', color='gray')
bars2 = ax.bar(datasets, tag_genes, bottom=non_tag_genes, label='Duplicated Genes Inside TAGs', color='lightblue')

ax.set_xlabel('Dataset')
ax.set_ylabel('Gene Count')
ax.set_title('Duplicated Genes Inside and Outside TAGs in Low vs High Stringency Datasets')
ax.legend()

for i, (bar1, bar2) in enumerate(zip(bars1, bars2)):
    # Non-TAG gene count
    ax.text(bar1.get_x() + bar1.get_width() / 2, bar1.get_height() / 2,
            f'{non_tag_genes[i]}', ha='center', va='center', color='black', fontsize=10)
    # TAG gene count
    ax.text(bar2.get_x() + bar2.get_width() / 2, bar1.get_height() + bar2.get_height() / 2,
            f'{tag_genes[i]}', ha='center', va='center', color='black', fontsize=10)
    
plt.savefig('tag_nontag.png', format='png', dpi=600)
plt.show()


tags_low = 'TAGs_low.txt'
tags_high = 'TAGs_high.txt'
tags_low = parse_tags_file(tags_low)
tags_high = parse_tags_file(tags_high)
plot_all_gene_distribution(tags_low, tags_high)
plot_tags_per_chromosome(tags_low, tags_high)