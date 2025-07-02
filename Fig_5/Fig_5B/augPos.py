import os
from collections import defaultdict

# Set the directory containing the BED files
bed_dir = "/projects/b1042/Arangolab/2ometh/5utrSequence/nmPositionsGenome/allCodons"
output_gene_file = "pos1Nmgenes.txt"
output_bed_file = "pos1Nmgenes.bed"

# Use a set to store unique gene names
unique_genes = set()

# Open the output .bed file for writing
with open(output_bed_file, 'w') as bed_out:
    # Loop through all .bed files in the directory
    for filename in os.listdir(bed_dir):
        if filename.endswith("Codons.bed"):
            filepath = os.path.join(bed_dir, filename)
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    # Skip section headers like "+1 Nm"
                    if line.startswith('+'):
                        continue
                    fields = line.split('\t')
                    if len(fields) == 6 and fields[5] in ("ATG", "CTG", "GTG", "ATC"):
                        chrom = fields[0]
                        pos = fields[1]
                        gene = fields[2]
                        strand = fields[3]
                        codon = fields[5]
                        unique_genes.add(gene)
                        bed_out.write(f"{chrom}\t{pos}\t{pos}\t{gene}\t{strand}\t{codon}\n")

# Write unique gene names
with open(output_gene_file, 'w') as out:
    for gene in sorted(unique_genes):
        out.write(gene + '\n')

print("Finished processing.")
print(f"Unique genes with Nm[+1] in ATG/CTG/GTG/ATC codon written to: {output_gene_file}")
print(f"Matching .bed lines written to: {output_bed_file}")