from Bio import SeqIO

# File paths
fasta_file = "Medicago_truncatula.MedtrA17_4.0.pep.all.fa"  # Replace with your FASTA file
blast_file = "Medicago_truncatula_Blastp_longIsoforme"  # Replace with your BLAST output file
output_file = "output_blast.txt"

# Parse the FASTA file to create a dictionary with sequence lengths
sequence_lengths = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    sequence_lengths[record.id] = len(record.seq)

# Function to calculate query and subject coverage
def calculate_coverage(alignment_length, seq_length):
    return (alignment_length / seq_length) * 100 if seq_length else 0

# Process the BLAST results and add the new columns
with open(blast_file, "r") as blast, open(output_file, "w") as out:
    for line in blast:
        columns = line.strip().split("\t")
        
        # Extract relevant columns
        query_id = columns[0]
        subject_id = columns[1]
        alignment_length = int(columns[3])

        # Get lengths from dictionary
        query_length = sequence_lengths.get(query_id, None)
        subject_length = sequence_lengths.get(subject_id, None)

        # Skip if either length is missing
        if query_length is None or subject_length is None:
            print(f"Missing sequence length for {query_id} or {subject_id}")
            continue

        # Calculate coverage
        query_coverage = calculate_coverage(alignment_length, query_length)
        subject_coverage = calculate_coverage(alignment_length, subject_length)

        # Write the original columns along with the new data
        out.write("\t".join(columns) + f"\t{query_length}\t{subject_length}\t{query_coverage:.2f}\t{subject_coverage:.2f}\n")

print("New BLAST results with added columns written to", output_file)
