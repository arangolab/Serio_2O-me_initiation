from pathlib import Path
import re

def fasta_to_bed(fasta_path, output_bed_path, codons=('ATG', 'CTG', 'GTG', 'ATC')):
    with open(fasta_path) as f:
        lines = f.readlines()

    bed_entries = []
    header = ''
    seq = ''

    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if header and seq:
                bed_entries += process_fasta_entry(header, seq, codons)
            header = line[1:]
            seq = ''
        else:
            seq += line.upper()

    # process the last entry
    if header and seq:
        bed_entries += process_fasta_entry(header, seq, codons)

    # Write to BED file
    with open(output_bed_path, 'w') as out:
        for entry in bed_entries:
            out.write('\t'.join(map(str, entry)) + '\n')

def process_fasta_entry(header, seq, codons):
    fields = header.split('|')
    if len(fields) < 8:
        return []

    chrom = fields[1]
    start_genomic = int(fields[2])
    gene_name = fields[-1]
    strand = fields[-3]  # need to check this is correct
    if strand == "1" or strand == "+1":
        strand = "+"
    elif strand == "-1":
        strand = "-"

    entries = []

    for codon in codons:
        for match in re.finditer(codon, seq):
            seq_pos = match.start()
            codon_start = start_genomic + seq_pos
            codon_end = codon_start
            entries.append([chrom, codon_start, codon_end, gene_name, strand, codon])

    return entries

# Example usage:
fasta_path = 'filtered5utr.fasta'
output_bed_path = '5utrNearCognates.bed'
fasta_to_bed(fasta_path, output_bed_path)
