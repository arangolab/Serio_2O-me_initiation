from Bio import SeqIO

gene_to_best_record = {}

for record in SeqIO.parse("results5UTR.fasta", "fasta"):
    if str(record.seq).strip().lower() == "sequenceunavailable":
        continue

    header_fields = record.description.split("|")
    if len(header_fields) < 6:
        continue  # skip malformed headers

    gene_name = header_fields[-1]
    try:
        cdna_length = int(header_fields[-2])
    except ValueError:
        continue

    existing_record = gene_to_best_record.get(gene_name)
    if not existing_record or cdna_length > int(existing_record.description.split("|")[-2]): # choose longest transcript available
        # Determine strand safely
        if len(header_fields) >= 7:
            strand = header_fields[-3]
        else:
            strand = "+"  # default if unknown

        if strand == "-1":
            record.seq = record.seq.reverse_complement()

        # Preserve entire original header exactly
        record.id = record.description
        record.name = record.description
        record.description = record.description

        gene_to_best_record[gene_name] = record

with open("filtered5utr.fasta", "w") as out_f:
    SeqIO.write(gene_to_best_record.values(), out_f, "fasta")
