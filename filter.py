import argparse

def filter_homologs(input_file, output_file, pid=30, que_cov=0.7, sub_cov=0.7):
    count = 0

    with open(input_file, 'r') as infile:
        with open(output_file, 'w') as outfile:
            for line in infile:
                columns = line.split()
                if len(columns) > 17:
                    if float(columns[2]) > pid and float(columns[16]) > que_cov and float(columns[17]) > sub_cov:
                        outfile.write(line)
                        count += 1

    print(f'Number of homologs: {count}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter homolog pairs from BLAST results.')
    parser.add_argument('pid', type=float, help='Percentage identity threshold')
    parser.add_argument('que_cov', type=float, help='Query coverage threshold')
    parser.add_argument('sub_cov', type=float, help='Subject coverage threshold')
    parser.add_argument('--input_file', type=str, default='update_blast_results.txt', help='Input BLAST results file')
    parser.add_argument('--output_file', type=str, default='homologs.txt', help='Output file for filtered homologs')

    args = parser.parse_args()

    filter_homologs(args.input_file, args.output_file, args.pid, args.que_cov, args.sub_cov)