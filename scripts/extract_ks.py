import sys

def extract_ks_values(yn_ds_file, output_file):
    try:
        with open(yn_ds_file, "r") as file:
            lines = file.readlines()  # Read all lines at once

        with open(output_file, "a") as outfile:
            # Extract values from the 2YN.dS file
            if len(lines) >= 3:  # Ensure there are enough lines
                seq1 = lines[1].strip()  # Second line: First sequence
                seq2 = lines[2].split()[0].strip()  # Third line: Second sequence
                dS_value = lines[2].split()[1].strip()  # Third line: dS value
                
                # Write the extracted values to the output file
                outfile.write(f"{seq1}\t{seq2}\t{dS_value}\n")
            else:
                print(f"Error: File '{yn_ds_file}' does not have the expected format.")
                sys.exit(1)

    except FileNotFoundError:
        print(f"Error: File '{yn_ds_file}' not found.")
        sys.exit(1)

if __name__ == "__main__":
    # Expecting arguments from the command line
    if len(sys.argv) != 3:
        print("Usage: python3 extract_ks.py <2YN.dS_file> <output_file>")
        sys.exit(1)

    yn_ds_file = sys.argv[1]
    output_file = sys.argv[2]

    extract_ks_values(yn_ds_file, output_file)
