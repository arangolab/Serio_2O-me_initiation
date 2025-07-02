import pandas as pd

# Load the BED file
input_file = "pos1Nmgenes.bed"
df = pd.read_csv(input_file, sep="\t", header=None)

# Assign column names for clarity
df.columns = ["chr", "start", "end", "gene", "strand", "codon"]

# Remove duplicates based on all columns except the codon
df_unique = df.drop_duplicates()

# Group by codon and write to separate files
for codon, group in df_unique.groupby("codon"):
    output_file = f"pos1Nm{codon}.bed"
    group.drop(columns=["codon"]).to_csv(output_file, sep="\t", header=False, index=False)

print("Files created successfully.")
