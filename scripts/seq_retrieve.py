import sys
from Bio import SeqIO

def extract_sequences(blast_result_line, protein_fasta, cds_fasta):
    # Parse BLAST result line
    columns = blast_result_line.strip().split("\t")
    protein1, protein2 = columns[0], columns[1]

    # Load protein and CDS sequences
    protein_seqs = SeqIO.to_dict(SeqIO.parse(protein_fasta, "fasta"))
    cds_seqs = SeqIO.to_dict(SeqIO.parse(cds_fasta, "fasta"))

    # Extract sequences for the pair
    prot_seqs = []
    cds_seqs_list = []

    if protein1 in protein_seqs and protein2 in protein_seqs:
        prot_seqs.append(protein_seqs[protein1])
        prot_seqs.append(protein_seqs[protein2])
    else:
        print(f"Error: Protein sequences for {protein1} or {protein2} not found.")
        sys.exit(1)

    if protein1 in cds_seqs and protein2 in cds_seqs:
        cds_seqs_list.append(cds_seqs[protein1])
        cds_seqs_list.append(cds_seqs[protein2])
    else:
        print(f"Error: CDS sequences for {protein1} or {protein2} not found.")
        sys.exit(1)

    # Write sequences to output files
    SeqIO.write(prot_seqs, "prot.fst", "fasta")
    SeqIO.write(cds_seqs_list, "cds.fst", "fasta")

if __name__ == "__main__":
    # Expecting arguments from the command line
    if len(sys.argv) != 4:
        print("Usage: python3 seq_retrieve.py <blast_result_line> <protein_fasta> <cds_fasta>")
        sys.exit(1)

    blast_result_line = sys.argv[1]
    protein_fasta = sys.argv[2]
    cds_fasta = sys.argv[3]

    extract_sequences(blast_result_line, protein_fasta, cds_fasta)
