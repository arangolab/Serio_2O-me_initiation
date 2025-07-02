import pyBigWig
import pandas as pd
import numpy as np
import sys

# Get the input and output file paths from command-line arguments
bigwig_file = sys.argv[1]
bed_file = sys.argv[2]
output_file = sys.argv[3]

# Open the BigWig file
try:
    bw = pyBigWig.open(bigwig_file)
except Exception as e:
    print(f"Error opening BigWig file: {e}")
    sys.exit(1)

# Read the BED file
try:
    bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["chrom", "start", "end", "gene", "strand"], dtype={"chrom": str, "start": int, "end": int, "gene": str, "strand": str})
except Exception as e:
    print(f"Error reading BED file: {e}")
    sys.exit(1)

# Create a list to store the results
results = []

# Check chromosome names in BigWig file
chroms_in_bw = set(bw.chroms().keys())

# Iterate over each coordinate in the BED file
for index, row in bed_df.iterrows():
    chrom = row["chrom"]
    start = row["start"]
    end = row["end"]
    strand = row["strand"]
    
    if chrom not in chroms_in_bw:
        print(f"Chromosome {chrom} not found in BigWig file.")
        continue

    try:
        if strand == '+':
            # For positive strand, extract RPKM values from -13 to +16 relative to start
            rpkm_values = bw.values(chrom, max(start - 13, 0), start + 17, numpy=True)
            col_labels = list(range(-13, 17))
        elif strand == '-':
            # For negative strand, extract RPKM values from -16 to +13 relative to end and reverse the order
            rpkm_values = bw.values(chrom, max(end - 16, 0), end + 14, numpy=True)[::-1]
            col_labels = list(range(16, -14, -1))
        else:
            raise ValueError(f"Unexpected strand value: {strand}")
        
        # Create a result row
        result_row = [chrom, start, end, row["gene"], strand] + list(rpkm_values)
        
        # Append the result row to the results list
        results.append(result_row)
    except Exception as e:
        print(f"Error processing {chrom}:{start}-{end} on strand {strand} - {e}")

# Convert the results list to a DataFrame
results_df = pd.DataFrame(results, columns=["chrom", "start", "end", "gene", "strand"] + col_labels)

# Save the results to a file
results_df.to_csv(output_file, index=False)