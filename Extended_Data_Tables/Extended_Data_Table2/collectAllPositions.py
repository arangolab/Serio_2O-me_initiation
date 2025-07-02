import argparse
import sys
import os
import glob
import re
import pandas as pd

# Function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Extract Nm near-cognate codons and output cleaned Excel/bed files.')
    parser.add_argument('-i', '--input', help='path to input files (can use wildcards like *.txt)', required=True)
    parser.add_argument('-o', '--output', help='output file prefix (e.g., output => output.xlsx, output.bed)', required=True)
    return parser.parse_args(args)

# Get arguments
args = check_arg(sys.argv[1:])
input_pattern = args.input
output_prefix = args.output
excel_output = output_prefix + ".xlsx"
bed_output = output_prefix + ".bed"

# Codon positions (start, end) for each Nm window
nm_windows = {
    "-3 Nm": (8, 11),
    "-2 Nm": (7, 10),
    "-1 Nm": (6, 9),
    "+1 Nm": (5, 8),
    "+2 Nm": (4, 7),
    "+3 Nm": (3, 6),
    "+4 Nm": (2, 5)
}

# Target codons
near_cognate_codons = {"ATG", "GTG", "TTG", "CTG", "ATA", "ATC", "ATT", "AGG", "ACG", "AAG"}

# Parse filename for metadata
def parse_filename(filename):
    base_name = filename.replace('Seq.bed', '')
    match = re.match(r'^([a-z]+)(.+)', base_name)
    if match:
        return match.group(1), match.group(2)
    return '', base_name

# Main function
def extract_codons_all_positions(input_pattern, nm_windows, excel_output, bed_output):
    input_files = glob.glob(input_pattern)
    if not input_files:
        print(f"No files found matching pattern: {input_pattern}")
        return

    writer = pd.ExcelWriter(excel_output, engine='xlsxwriter')
    all_rows_combined = []

    for label, (start, end) in nm_windows.items():
        rows = []

        for file in input_files:
            filename = os.path.basename(file)
            reference, cell_line = parse_filename(filename)

            with open(file, 'r') as f:
                for line in f:
                    columns = line.strip().split('\t')
                    if len(columns) < 5:
                        continue

                    seq = columns[4]
                    if end <= len(seq):
                        codon = seq[start:end]
                        if codon in near_cognate_codons:
                            row = columns[:5] + [codon, reference, cell_line]
                            rows.append(row)

        if rows:
            df = pd.DataFrame(rows, columns=["Chr", "Nm_Position", "Gene", "Strand", "Sequence", "Codon", "Reference", "Cell_Line"])

            # Post-processing steps
            df = df.drop(columns=["Sequence"])
            df["Reference"] = df["Reference"].replace({
                "mlm": "Nanopore",
                "nju": "NJU-seq",
                "nm": "Nm-seq",
                "nmmt": "Nm-mut-seq"
            })
            df["Cell_Line"] = df["Cell_Line"].replace({
                "HEK": "HEK293T"
            })
            df["Codon"] = df["Codon"].str.replace("T", "U", regex=False)

            # Write sheet
            sheet_name = label.replace(" ", "")
            df.to_excel(writer, sheet_name=sheet_name[:31], index=False)  # Excel sheet names max out at 31 chars
            all_rows_combined.append(df)

    # Combine all and write to flat files
    if all_rows_combined:
        final_df = pd.concat(all_rows_combined, ignore_index=True)
        final_df.to_csv(bed_output, sep="\t", index=False)

    writer.close()
    print(f"Completed. Files written:\n- {excel_output} (sheets)\n- {bed_output} (.bed)")

# Run it
extract_codons_all_positions(input_pattern, nm_windows, excel_output, bed_output)