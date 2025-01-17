import pandas as pd

# Input file paths
ks_file_path = "Ks_final.txt"
tag_file_path = "TAGs_high.txt"
output_file_path = "Ks_with_tags.txt"

# Load the Ks file, ensuring the first line is treated as the header
ks_df = pd.read_csv(ks_file_path, sep="\t", header=0)  # `header=0` means the first row is the header
ks_df.columns = ["gene1", "gene2", "Ks"]  # Ensure column names match

# Extract gene pairs from TAG file
tag_gene_pairs = set()
with open(tag_file_path, "r") as tag_file:
    for line in tag_file:
        line = line.strip()  # Remove leading/trailing whitespace
        if not line:  # Skip empty lines
            continue
        parts = line.split(" ")  # Split by space
        if len(parts) < 3:  # Ensure the line has at least three parts
            print(f"Skipping malformed line: {line}")
            continue
        tag_info = parts[2]  # The third part contains gene info
        genes = tag_info.split(",")  # Split by commas to get individual gene entries
        gene_names = [entry.split(":")[0] for entry in genes]  # Extract gene names (before the first colon)
        for i in range(len(gene_names)):
            for j in range(i + 1, len(gene_names)):
                # Add both gene1-gene2 and gene2-gene1 to account for order
                tag_gene_pairs.add((gene_names[i], gene_names[j]))
                tag_gene_pairs.add((gene_names[j], gene_names[i]))

# Check if each gene pair in the Ks file is in the TAG gene pairs
ks_df["Tag_Status"] = ks_df.apply(
    lambda row: "tag" if (row["gene1"], row["gene2"]) in tag_gene_pairs else "nontag",
    axis=1,
)

# Save the updated Ks file
ks_df.to_csv(output_file_path, sep="\t", index=False)
print(f"Updated Ks file saved to {output_file_path}")