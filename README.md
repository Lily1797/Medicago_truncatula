# Comparative genomics project for *Medicago truncatula*
 For detailed analysis workflow, visit the [Github repository](https://github.com/Lily1797/Medicago_truncatula.git)

## Data
In this project, we analyzed the dupplicate genes in the genomes of the *Medicago truncatula* plant. The *M. truncatula* genome contains eight chromosomes with the genome size of ~412 Mb.
Number of genes: 51519;
Number of isoforms: ??.

Count number of sequences located on chromosomes and scaffolds
```
awk '/^>/ { if (/chromosome:/) chr_count++; else if (/:scaffold/) scaffold_count++; } END { print "Chromosome sequences:", chr_count; print "Scaffold sequences:", scaffold_count; }' Medicago_truncatula.MedtrA17_4.0.pep.all.fa
```
Chromosome sequences: 54929;
Scaffold sequences: 2656.

## Workflow
### I.	Determine the number of duplicat genes in my genome
#### 1. Update Blast result
A Python script (blast_extend.py) enriched the BLAST results with four additional columns: query length, subject length, query coverage, and subject coverage. The enriched results were saved to output_blast.txt.
```
python scipts/blast_extend.py
```
As no mitochondrial or chloroplast sequences were identified in our genome, this step was omitted. To streamline subsequent genome map analysis, we filtered the BLAST results to include only sequences located on chromosomes. Additionally, ten columns were added to the filtered BLAST results, specifying the positions of the sequences on the chromosomes (query start, query end, query chrom, query strand, query geneID, subject start, subject end, subject chrom, subject strand and subject geneID). The updated results were saved to update_blast_results.txt.
```
# Extract Positions, chromosome, strand, and geneID from the FASTA File
awk '/^>/ {if (/chromosome:/) { split($0, a, " "); seq_id = a[1]; seq_id = substr(seq_id, 2); split(a[3], b, ":"); start = b[4]; end = b[5]; chrom = b[3]; strand = b[6]; split(a[4], c, ":"); geneid = c[2]; print seq_id "\t" start "\t" end "\t" chrom "\t" strand "\t" geneid;}}' Medicago_truncatula.MedtrA17_4.0.pep.all.fa > positions.txt

# Filter BLAST results and add information
python scipts/add_position.py
```
We then filter the result based on percentage of identity, query coverage and subject coverage to generate two distinct datasets: one with low stringency and another with high stringency.
```
awk '$3 >= 30 && $17 >= 30 && $18 >= 30' update_blast_results.txt > homolog_low.txt
awk '$3 >= 50 && $17 >= 40 && $18 >= 40' update_blast_results.txt > homolog_high.txt
```
Keep only unique records
```
awk '{ if ($1 < $2) {  key = $1 "\t" $2 } else { key = $2 "\t" $1  } if (!seen[key]++) {  print } }' homolog_low.txt > homolog_low_unique.txt
awk '{ if ($1 < $2) {  key = $1 "\t" $2 } else { key = $2 "\t" $1  } if (!seen[key]++) {  print } }' homolog_high.txt > homolog_high_unique.txt
```

#### Clustering
Extract specific columns from the homolog files and created input files for the MCL clustering method.
```
awk '{print $21, $26, $12}' homolog_low_unique.txt > cluster_input_low.txt
awk '{print $21, $26, $12}' homolog_high_unique.txt > cluster_input_high.txt
```
We ran the clustering with the MCL method for the two datasets. We identified 4938 (43300 genes) and 7006 (38853 genes) families for low and high stringency dataset, respectively, and then visualized the clustering result.
```
python scipts/plot_cluster.py
```
Distribution of number of genes across families
![Distribution_clusters](./figures/cluster_distribution.png)
Number of duplicated genes and singletons
![Duplicates_singletons](./figures/dup_single.png)

#### TAGs
Extract homologous pairs on the same chromosome.
```
awk -F'\t' '$19 == $24' homolog_low_unique.txt > same_chromosomes_low.txt
awk -F'\t' '$19 == $24' homolog_high_unique.txt > same_chromosomes_high.txt
```
Extract TAGs 
```
python scipts/TAGs_finder.py
```
We identified 4499 (13003 genes) and 4191 (12221 genes) TAGs for low and high stringency dataset, respectively.
What are the biggest TAGs?
```
# For low stringency set
$awk -F' ' '{print $1, $2, gsub(",", ",", $3)+1}' TAGs_low.txt | sort -k3,3nr | head -n 10
TAG2099 Chromosome:4 30
TAG2996 Chromosome:6 27
TAG118 Chromosome:1 26
TAG4177 Chromosome:8 23
TAG2768 Chromosome:5 22
TAG3010 Chromosome:6 22
TAG3464 Chromosome:7 20
TAG935 Chromosome:2 20
TAG2527 Chromosome:5 19
TAG3992 Chromosome:8 19

#For high stringency set
$awk -F' ' '{print $1, $2, gsub(",", ",", $3)+1}' TAGs_high.txt | sort -k3,3nr | head -n 10
TAG2825 Chromosome:6 31
TAG1971 Chromosome:4 30
TAG2819 Chromosome:6 28
TAG3091 Chromosome:6 26
TAG3941 Chromosome:8 24
TAG2599 Chromosome:5 22
TAG3259 Chromosome:7 22
TAG2831 Chromosome:6 21
TAG872 Chromosome:2 21
TAG1316 Chromosome:3 18
```
Visualize the results
```
python scipts/plot_TAGs.py
```
Distribution of number of genes in TAGs
![Distribution_TAGs](./figures/gene_distribution.png)
Barplot of number of TAGs per Chromosome
![TAGs_per_Chrom](./figures/tag_per_chrom.png)
Number of duplicated genes inside and outside TAGs
![TAG_ NonTAG](./figures/tag_nontag.png)
Inside the TAGs, are gene oriented significantly the same way?


#### Ks value
We calculated Ks values for the ?? stringency dataset as it is better.
Script to calculate all gene pairs inside families and store list of pairs
```
# There is a Blast command to retrieve the name of all pairs of genes (last year course) 
```
Retrieve the cds sequences on Ensembl plants, only the longset isoform

Script to calculate Ks values for all pairs with PAML
```{bash}
# for loop for all pairs of genes

```
Selection of Ks values < threshold  or with small SD

Analyze the results
Are the TAGs pairs different in age with the other duplicate pairs?
