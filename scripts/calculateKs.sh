#!/bin/bash

# Paths to necessary files and directories
protein_fasta="/home/huongpm//data/comparative_genomics/Medicago_truncatula/Medicago_truncatula.MedtrA17_4.0.pep.all.fa"
cds_fasta="/home/huongpm//data/comparative_genomics/Medicago_truncatula/Medicago_truncatula.MedtrA17_4.0.cds.all.fa"
pal2nal_script="/home/huongpm//data/comparative_genomics/Medicago_truncatula/pal2nal.pl"
control_file_template="/home/huongpm//data/comparative_genomics/Medicago_truncatula/yn00.ctl_master.txt"
yn00_executable="/home/huongpm/software/paml/bin/yn00"
seq_retrieve_script="/home/huongpm//data/comparative_genomics/Medicago_truncatula/seq_retrieve.py"
ks_extract_script="/home/huongpm//data/comparative_genomics/Medicago_truncatula/extract_ks.py"
blast_file="/home/huongpm//data/comparative_genomics/Medicago_truncatula/dup_pairs_high.txt"
output_file="/home/huongpm//data/comparative_genomics/Medicago_truncatula/Ks_calculated.txt"
yn_output="2YN.dS"

# Temporary directory for parallel processing
tmp_dir=$(mktemp -d)
export tmp_dir

# Function to process each pair
process_pair() {
    pair=$1
    echo "Processing pair: $pair"
    
    # Create a temporary working directory for this pair
    pair_dir="$tmp_dir/$pair"
    mkdir -p "$pair_dir"
    cd "$pair_dir" || exit 1
    
    # Step 1: Retrieve the protein and CDS sequences
    python3 "$seq_retrieve_script" "$pair" "$protein_fasta" "$cds_fasta"
    
    # Step 2: Align the protein sequences
    clustalw2 -quiet -align -infile=prot.fst -outfile=prot.ali.aln
    
    # Step 3: Align CDS sequences based on protein alignment
    "$pal2nal_script" prot.ali.aln cds.fst -output paml > cds.ali.phy
    
    # Step 4: Modify the control file with the current alignment file
    awk -v file=cds.ali.phy '{gsub("XXXXX",file); print $0}' "$control_file_template" > yn00.ctl
    
    # Step 5: Run yn00
    "$yn00_executable"
    
    # Step 6: Extract calculated values and append to the final output file
    python3 "$ks_extract_script" "$yn_output" "$tmp_dir/$pair.ks"

    echo "Finished processing pair: $pair"
}

export -f process_pair
export protein_fasta cds_fasta pal2nal_script control_file_template yn00_executable seq_retrieve_script ks_extract_script yn_output

# Use GNU parallel to process pairs in parallel
cat "$blast_file" | parallel --joblog "$tmp_dir/joblog" --results "$tmp_dir/results" process_pair

# Combine all results into the final output file
find "$tmp_dir" -name "*.ks" -exec cat {} + > "$output_file"

# Clean up
rm -rf "$tmp_dir"

echo "Pipeline completed. Results are in $output_file"
