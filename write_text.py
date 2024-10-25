import argparse

def extract_columns(input_file, output_file):
    with open(input_file, 'r') as infile:
        with open(output_file, 'w') as outfile:
            for line in infile:
                columns = line.split()
                if len(columns) > 11:
                    output_line = f"{columns[0]}\t{columns[1]}\t{columns[11]}\n"
                    outfile.write(output_line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract specific columns from a homolog file.')
    parser.add_argument('input_file', type=str, help='Input homolog file')
    parser.add_argument('output_file', type=str, help='Output file for extracted columns')

    args = parser.parse_args()

    extract_columns(args.input_file, args.output_file)