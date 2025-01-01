import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#Distribution of Number of Genes across Families
def plot_gene_distribution(family_low, family_high):
    gene_low_counts = []
    gene_high_counts = []

    with open(family_low, 'r') as l:
        for line in l:
            genes_low = line.strip().split('\t')  
            gene_low_counts.append(len(genes_low))
    with open(family_high, 'r') as h:
        for line in h:
            genes_high = line.strip().split('\t')  
            gene_high_counts.append(len(genes_high))  

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6))

    ax1.hist(gene_low_counts, bins=range(1, max(gene_low_counts) + 2), color='skyblue', align='left')
    ax1.set_xlabel('Number of Genes per Family')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Low stringency')
    ax1.grid(axis='y', linestyle='--', alpha=0.7)

    ax2.hist(gene_high_counts, bins=range(1, max(gene_high_counts) + 2), color='salmon', align='left')
    ax2.set_xlabel('Number of Genes per Family')
    ax2.set_ylabel('Frequency')
    ax2.set_title('High stringency')
    ax2.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()
    plt.savefig("cluster_distribution.png", format='png', dpi=600)
    plt.show()

family_low = 'MCL_low.tabular'
family_high = 'MCL_high.tabular'
plot_gene_distribution(family_low, family_high)


# Duplicates and Singletons in Low vs High Stringency Datasets
datasets = ['Low Stringency', 'High Stringency']
duplicates = [43300, 38853]  
singletons = [50894 - 43300, 50894 - 38853] 
total_genes = [50894, 50894]  # Total genes for each dataset

fig, ax = plt.subplots(figsize=(8, 6))

bars1 = ax.bar(datasets, singletons, label='Singletons', color='gray')
bars2 = ax.bar(datasets, duplicates, bottom=singletons, label='Duplicates', color='lightblue')

ax.set_xlabel('Dataset')
ax.set_ylabel('Gene Count')
ax.set_title('Duplicated and Singletons in Low vs High Stringency Datasets')
ax.legend()

for i, (bar1, bar2) in enumerate(zip(bars1, bars2)):
    # singletons count and proportion
    single_proportion = singletons[i] / total_genes[i] * 100
    ax.text(bar1.get_x() + bar1.get_width() / 2, bar1.get_height() / 2,
            f'{singletons[i]} ({single_proportion:.1f}%)', ha='center', va='center', color='black', fontsize=10)
    
    # Duplicate gene count and proportion
    dup_proportion = duplicates[i] / total_genes[i] * 100
    ax.text(bar2.get_x() + bar2.get_width() / 2, bar1.get_height() + bar2.get_height() / 2,
            f'{duplicates[i]} ({dup_proportion:.1f}%)', ha='center', va='center', color='black', fontsize=10)

plt.savefig('dup_single.png', format='png', dpi=600)
plt.show()