import argparse
import sys

# Function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Extracts all codons with Nm[+1] in 5UTR')
    parser.add_argument('-i', '--input',
                        help='path to input file',
                        required=True
                        )
    parser.add_argument('-o', '--output',
                        help='output file name',
                        required=True
                        )
    return parser.parse_args(args)

# Retrieve command line arguments and assign to variables
args = check_arg(sys.argv[1:])
infile = args.input
outfile = args.output

# Function to extract rows with near-cognate codons and write with matched codon
def extract_near_cognate_codons(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        sections = {"+1 Nm": (5, 8, [])}

        for line in infile:
            columns = line.strip().split('\t')
            if len(columns) < 5:
                continue  # Skip malformed lines
            seq = columns[4]

            for label, (start, end, collector) in sections.items():
                if end <= len(seq):
                    codon = seq[start:end]
                    new_line = line.strip() + '\t' + codon + '\n'
                    collector.append(new_line)

        for label, (_, _, lines) in sections.items():
            outfile.write(label + '\n')
            for item in lines:
                outfile.write(item)

# Run the function
extract_near_cognate_codons(infile, outfile)

print(f"Rows with near-cognate codons have been extracted to {outfile}")