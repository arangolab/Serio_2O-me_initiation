import os

def parse_fasta(fasta_file):
    """Parse FASTA file and return sequences with their headers."""
    sequences = []
    current_header = ""
    current_seq = ""
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header and current_seq:
                    sequences.append((current_header, current_seq))
                current_header = line[1:]
                current_seq = ""
            else:
                current_seq += line
        if current_header and current_seq:
            sequences.append((current_header, current_seq))
    
    return sequences

def parse_header(header):
    """Parse the header to extract chromosome, start, end, gene, strand, and transcript ID."""
    parts = header.split('|')
    if len(parts) >= 6:
        chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        gene = parts[3]
        strand = parts[4]
        transcript = parts[5]
        return chrom, start, end, gene, strand, transcript
    else:
        raise ValueError(f"Header format is incorrect: {header}")

def find_codon_positions(sequence, codon):
    """Find all occurrences of the specified codon in the sequence."""
    positions = []
    for i in range(len(sequence) - 2):
        if sequence[i:i+3].upper() == codon:
            positions.append(i)
    return positions

def convert_to_genomic_coords(positions, chrom, start, end, strand):
    """Convert sequence positions to genomic coordinates (for the 'A' base in the codon)."""
    coords = []
    for pos in positions:
        if strand == '+':
            genomic_start = start + pos
            genomic_end = genomic_start + 1
        else:
            genomic_end = end - pos
            genomic_start = genomic_end - 1
        coords.append((chrom, genomic_start, genomic_end))
    return coords

def main():
    input_file = "/projects/b1042/Arangolab/2ometh/5utrSequence/nmPositionsGenome/ATGcodon/fiveUTRsequence.fasta"
    codons = ["ATG", "GTG", "CTG", "ATC"]

    # Initialize writers for each codon
    output_files = {codon: open(f"fiveUTR_{codon}positions.bed", 'w') for codon in codons}

    sequences = parse_fasta(input_file)
    for header, seq in sequences:
        chrom, start, end, gene, strand, transcript = parse_header(header)
        for codon in codons:
            codon_positions = find_codon_positions(seq, codon)
            genomic_coords = convert_to_genomic_coords(codon_positions, chrom, start, end, strand)
            for gc_chrom, gc_start, gc_end in genomic_coords:
                output_files[codon].write(f"{gc_chrom}\t{gc_start}\t{gc_end - 1}\t{gene}\t{strand}\n")

    for f in output_files.values():
        f.close()

    print("Codon positions written to individual .bed files.")

if __name__ == "__main__":
    main()
