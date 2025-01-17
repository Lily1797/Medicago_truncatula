#!/bin/bash

# Paths to necessary files and directories
protein_fasta="Medicago_truncatula.MedtrA17_4.0.pep.all.fa"
cds_fasta="Medicago_truncatula.MedtrA17_4.0.cds.all.fa"
pal2nal_script="./pal2nal.pl"
control_file_template="yn00.ctl_master.txt"
control_file="yn00.ctl"
yn00_executable="/home/huongpm/software/paml/bin/yn00"
seq_retrieve_script="./seq_retrieve.py"
ks_extract_script="./extract_ks.py"
blast_file="testks.txt" 
output_file="Ks_calculated.txt"
yn_output="2YN.dS"

# Process each line in the blast file
while IFS= read -r pair; do
    echo "Processing pair: $pair"
    
    # Step 1: Retrieve the protein and CDS sequences
    python3 "$seq_retrieve_script" "$pair" "$protein_fasta" "$cds_fasta"
    
    # Step 2: Align the protein sequences
    clustalw2 -quiet -align -infile=prot.fst -outfile=prot.ali.aln
    
    # Step 3: Align CDS sequences based on protein alignment
    "$pal2nal_script" prot.ali.aln cds.fst -output paml > cds.ali.phy
    
    # Step 4: Modify the control file with the current alignment file
    awk -v file=cds.ali.phy '{gsub("XXXXX",file); print $0}' "$control_file_template" > "$control_file"
    
    # Step 5: Run yn00
    "$yn00_executable"
    
    # Step 6: Extract calculated values and append to the final output file
    python3 "$ks_extract_script" "$yn_output" "$output_file"

    
    echo "Finished processing pair: $pair"
done < "$blast_file"

echo "Pipeline completed. Results are in $output_file"
