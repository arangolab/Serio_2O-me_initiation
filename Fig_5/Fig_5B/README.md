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

# Fig. 5B

1. Extract the coordinates and genes for Nm[+1] for all ATG, CTG, GTG, ATC codons, and output a list of any genes containing Nm.

```
python augPos.py
```

2. Separate the file by near-cognate codons
```
python parseNearCognates.py
```

3. Extract coordinates for 'A' in canonical 'ATG' start codon for genes without Nm present and coordinates for 'A' in canonical 'ATG' start codon for genes with Nm[+1] in ATG.
```
python extractCanonicalCoordinates.py
```

4. Extract the coordinates of any AUG/GUG/CUG/AUC in the 5'UTR of non-methylated transcripts

```
python extractATGcoordinates.py
```

5. Now we have prepped bed files for ribosome density extraction
*canonicalStart.bed # genomic coordinates for 'A' in canonical 'ATG' start codon for genes without Nm present
*canonicalStartNMgenes.bed # genomic coordinates for 'A' in canonical 'ATG' start codon for genes with Nm[+1] in ATG
*fiveUTR_<CODON>positions.bed # genomic coordinates for ATG/CTG/GTG/ATC within the 5UTR for genes without Nm present
*pos1Nm<CODON>.bed # genomic coordinates for Nm[+1] within the 5UTR

# Extract the ribosome density
6. Run the following Python script:

```
# genomic coordinates for 'A' in canonical 'ATG' start codon for genes without Nm present
BED_FILE="~/ribo.density.py "$BIGWIG_FILE" "$BED_FILE" "$OUTPUT_DIR/${SAMPLE_NAME}_canonical.csv"

# genomic coordinates for 'A' in canonical 'ATG' start codon for genes with Nm[+1] in ATG
BED_FILE="~/ribo.density.py "$BIGWIG_FILE" "$BED_FILE" "$OUTPUT_DIR/${SAMPLE_NAME}_canonicalNm.csv"

# genomic coordinates for ATG/CTG/GTG/ATC within the 5UTR for genes without Nm present
BED_FILE="~/ribo.density.py "$BIGWIG_FILE" "$BED_FILE" "$OUTPUT_DIR/${SAMPLE_NAME}_ATG.csv"

# genomic coordinates for Nm[+1] within the 5UTR
BED_FILE="~/ribo.density.py "$BIGWIG_FILE" "$BED_FILE" "$OUTPUT_DIR/${SAMPLE_NAME}_ATGNm.csv"

# genomic coordinates for ATG/CTG/GTG/ATC within the 5UTR for genes without Nm present
BED_FILE="~/ribo.density.py "$BIGWIG_FILE" "$BED_FILE" "$OUTPUT_DIR/${SAMPLE_NAME}_CTG.csv"

# genomic coordinates for Nm[+1] within the 5UTR
BED_FILE="~/ribo.density.py "$BIGWIG_FILE" "$BED_FILE" "$OUTPUT_DIR/${SAMPLE_NAME}_CTGNm.csv"

# genomic coordinates for ATG/CTG/GTG/ATC within the 5UTR for genes without Nm present
BED_FILE="~/ribo.density.py "$BIGWIG_FILE" "$BED_FILE" "$OUTPUT_DIR/${SAMPLE_NAME}_GTG.csv"

# genomic coordinates for Nm[+1] within the 5UTR
BED_FILE="~/ribo.density.py "$BIGWIG_FILE" "$BED_FILE" "$OUTPUT_DIR/${SAMPLE_NAME}_GTGNm.csv"

# genomic coordinates for ATG/CTG/GTG/ATC within the 5UTR for genes without Nm present
BED_FILE="~/ribo.density.py "$BIGWIG_FILE" "$BED_FILE" "$OUTPUT_DIR/${SAMPLE_NAME}_ATC.csv"

# genomic coordinates for Nm[+1] within the 5UTR
BED_FILE="~/ribo.density.py "$BIGWIG_FILE" "$BED_FILE" "$OUTPUT_DIR/${SAMPLE_NAME}_ATCNm.csv"
```

# Rstudio analysis and plots
7. Run the R code DensityBoxplotN.Rmd. The markdown is DensityBoxplotN.html. The data sets needed to reproduce the plots are:

* MH85_ATC.csv
* MH85_ATCNm.csv
* MH85_ATG.csv
* MH85_ATGNm.csv
* MH85_canonical.csv
* MH85_canonicalNm.csv
* MH85_CTG.csv
* MH85_CTGNm.csv
* MH85_GTG.csv
* MH85_GTGNm.csv
* MH87_ATC.csv
* MH87_ATCNm.csv
* MH87_ATG.csv
* MH87_ATGNm.csv
* MH87_canonical.csv
* MH87_canonicalNm.csv
* MH87_CTG.csv
* MH87_CTGNm.csv
* MH87_GTG.csv
* MH87_GTGNm.csv
