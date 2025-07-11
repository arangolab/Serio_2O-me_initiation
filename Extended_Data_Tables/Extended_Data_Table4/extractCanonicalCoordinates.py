def process_gtf_file(gtf_file, pos1_atg_genes_file, output_pos1_file, output_canonical_file):
    # Read the list of pos1ATG genes into a set for efficient lookup
    pos1_atg_genes = set()
    try:
        with open(pos1_atg_genes_file, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene:
                    pos1_atg_genes.add(gene)
    except FileNotFoundError:
        print(f"Warning: Gene list file '{pos1_atg_genes_file}' not found. Creating an empty set.")

    # Keep track of genes we've already written
    written_genes = set()

    with open(output_pos1_file, 'w') as pos1_out, open(output_canonical_file, 'w') as canonical_out:
        with open(gtf_file, 'r') as gtf:
            for line in gtf:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) >= 8 and fields[2] == "start_codon" and "tag \"Ensembl_canonical\"" in line:
                    chrom = fields[0]
                    strand = fields[6]
                    coord = fields[3] if strand == "+" else fields[4] if strand == "-" else None
                    if coord is None:
                        continue

                    attributes = fields[8]
                    gene_name = "NA"
                    gene_name_start = attributes.find("gene_name \"")
                    if gene_name_start != -1:
                        gene_name_start += len("gene_name \"")
                        gene_name_end = attributes.find("\"", gene_name_start)
                        if gene_name_end != -1:
                            gene_name = attributes[gene_name_start:gene_name_end]

                    # Skip if this gene has already been written
                    if gene_name in written_genes:
                        continue

                    output_line = f"{chrom}\t{coord}\t{coord}\t{gene_name}\t{strand}\n"

                    if gene_name in pos1_atg_genes:
                        pos1_out.write(output_line)
                    else:
                        canonical_out.write(output_line)

                    written_genes.add(gene_name)


gtf_file = "/projects/b1042/Arangolab/2ometh/canonicalStartSites/Homo_sapiens.GRCh38.113.chr.gtf"
pos1_atg_genes_file = "/projects/b1042/Arangolab/2ometh/5utrSequence/nmPositionsGenome/nearCognateNmBoxplot/pos1Nmgenes.txt"
output_pos1_file = "canonicalStartNMgenes.bed"
output_canonical_file = "canonicalStart.bed"

process_gtf_file(gtf_file, pos1_atg_genes_file, output_pos1_file, output_canonical_file)