# Medicago_truncatula
 Comparative genomics project

## Data
In this project, we analyzed the dupplicate genes in the genomes of the *Medicago truncatula* plant. The *M. truncatula* genome contains eight chromosomes with the genome size of ?? Mb.
The dataset includes three files...

Count number of sequences located on chromosomes and scaffolds
```
awk '/^>/ { if (/chromosome:/) chr_count++; else if (/:scaffold/) scaffold_count++; } END { print "Chromosome sequences:", chr_count; print "Scaffold sequences:", scaffold_count; }' Medicago_truncatula.MedtrA17_4.0.pep.all.fa
```

## Workflow
### I.	Determine the number of duplicate genes in my genome
#### 1. Update Blast result
A Python script (blast_extend.py) enriched the BLAST results with four additional columns: query length, subject length, query coverage, and subject coverage. The enriched results were saved to output_blast.txt.
```
python3 blast_extend.py
```
As no mitochondrial or chloroplast sequences were identified, this step was omitted. To streamline subsequent genome map analysis, we filtered the BLAST results to include only sequences located on chromosomes. Additionally, six columns were added to the filtered BLAST results, specifying the positions of the sequences on the chromosomes (query start, query end, query chrom, subject start, subject end, and subject chrom). The updated results were saved to update_blast_results.txt.
```
# Extract Positions from the FASTA File
awk '/^>/ {if (/chromosome:/) { split($0, a, " "); seq_id = a[1]; seq_id = substr(seq_id, 2); split(a[3], b, ":"); start = b[4]; end = b[5]; chrom = b[3]; print seq_id "\t" start "\t" end "\t" chrom;}}' Medicago_truncatula.MedtrA17_4.0.pep.all.fa > positions.txt
# Filter BLAST results and add position information
python3 add_position.py
```
Using a Python script (filter.py), we then identified homolog pairs within the file. The script employed various parameters (pid, que_cov, and sub_cov) to generate two distinct datasets: one with low stringency and another with high stringency.
```
python filter.py 30 30 30 --output_file homolog_low.txt
python filter.py 50 40 40 --output_file homolog_high.txt
```
#### Clustering
A Python script (write_text.py) extracted specific columns from the homolog files and created input files for the FTAG Finder v3 clustering method.
```
python write_text.py homolog_low.txt extracted_homolog_low.txt
python write_text.py homolog_high.txt extracted_homolog_high.txt
```