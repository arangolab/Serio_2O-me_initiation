# HR-Ribo-seq Processing

## Samples
The following datasets were downloaded from GEO (GSE162043, Arango et al, 2022 Mol Cell)

GSM4932133	HR_RPF_HeLa_WT_rep1 \
GSM4932134	HR_RPF_HeLa_WT_rep2

## Pre-processing

Adapters were trimmed using cutadapt/4.2 and the following code:

```
# Run cutadapt in single-end read mode (to see what each parameter means, look at the cutadapt version 4.2 manual)
cutadapt --match-read-wildcards -e 0.1 -O 1 --quality-cutoff 10 -m 22 -a CTGTAGGCACCATCAAT -o $OUTPUT_R1 $R1 > $METRICS_FILE
```

## Alignment

Reads were aligned to Human Genome (hg38) using hisat2/2.1.0 and the following code:

```
hisat2 -q -x "$reference_genome" -U "$R1" -S "$OUTPUT_DIR/${SAMPLE_NAME}.sam" --summary-file $OUTPUT_DIR/${SAMPLE_NAME}.txt
```

SAM files were then converted to BAM format, sorted, and indexed using samtools/1.6:

```
for FILE in "${FILES[@]}"; do
    # Get the sample name without the extension
    SAMPLE_NAME=$(basename "$FILE" .sam)
    # Convert SAM file to BAM file, sort it, and index the sorted BAM file
    samtools view -bS "$FILE" | samtools sort -o "$OUTPUT_DIR/${SAMPLE_NAME}.bam" -
    samtools index "$OUTPUT_DIR/${SAMPLE_NAME}.bam"
done
```
## Normalization

BAM files were converted to BIGWIG format using a bin of 1 nucleotide and the RPKM normalization using deeptools/3.5.1

```
bamCoverage -b "$FILE" -o "$OUTPUT_DIR/${SAMPLE_NAME}.bw" --outFileFormat bigwig --normalizeUsing RPKM --binSize 1
```

# Fig. 5A

# Extract the coordinates of Canonical AUG, upstream AUG, CUG, GUG, and AUC

1. xml wget command which outputs results5UTR.fasta

```
wget -O result.txt 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "5utr" />
		<Attribute name = "chromosome_name" />
		<Attribute name = "start_position" />
		<Attribute name = "end_position" />
		<Attribute name = "ensembl_transcript_id" />
		<Attribute name = "strand" />
		<Attribute name = "transcript_length" />
		<Attribute name = "external_gene_name" />
	</Dataset>
</Query>'
```
2. Select for the longest transcript per gene
```
python parseFasta.py
```
outputs: filtered5utr.fasta

3. Coordinates for upstream instances of AUG, CUG, GUG, AUC
```
python fastaToBed.py
```
Inputs: filtered5utr.fasta \
Outputs: 5utrNearCognates.bed

4. Separate entries by codon
```
python parseNearCognates.py
```
Inputs: 5utrNearCognates.bed
Outputs: 5utrNearCognatesATG.bed 5utrNearCognatesCTG.bed 5utrNearCognatesGTG.bed 

6. Extract the location of the first nucleotide of the canonical AUG start from GTF file 
```
python extractCanonicalCoordinates.py
```
Outputs: canonicalStart.bed (genomic coordinates of the canonical start codon)

# Extract the ribosome density around canonical and upstream codons 

