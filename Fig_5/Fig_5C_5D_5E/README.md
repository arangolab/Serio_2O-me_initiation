# Fig. 5C, 5D, 5E

## File preparation per study
Since we are using 2'O-methylation maps from different publications, file preparation differs slightly based on the original data provided.
### Nanopore-DRS (Li et. al. 2024)
Table S5 (excel file) contains a list of Nm sites detected in HEK293T, HeLa, and HepG2 cell lines. Each sheet contains the results for each cell line. The dataset was missing strand and mRNA consequence (5UTR, CDS, or 3UTR) information. 

We used this awk code to remove any Nm sites with a ratio below 10%. 
```
awk -F'\t' 'NR==1 || $3 >= 0.10 { print $1, $2 }' OFS='\t' HEKallData.txt > HEK_rawFiltered.txt
awk -F'\t' 'NR==1 || $3 >= 0.10 { print $1, $2 }' OFS='\t' HeLaallData.txt > HeLa_rawFiltered.txt
awk -F'\t' 'NR==1 || $3 >= 0.10 { print $1, $2 }' OFS='\t' HepG2allData.txt > HepG2_rawFiltered.txt
```
Awk was used to reformat and sort BED files to prepare for downstream analysis.
```
awk -F '\t' '{
  split($1, a, "_"); # original format is "chr1_12345"
  chrom = substr(a[1], 4);  # remove "chr"
  start = a[2];
  end = start + 1;
  gene = $2; # extract gene name
  print chrom "\t" start "\t" end "\t" gene;
}' HEK_rawFiltered.txt | sort -k1,1n -k2,2n > HEKconsequence10.bed

awk -F '\t' '{
  split($1, a, "_"); # original format is "chr1_12345"
  chrom = substr(a[1], 4);  # remove "chr"
  start = a[2];
  end = start + 1;
  gene = $2; # extract gene name
  print chrom "\t" start "\t" end "\t" gene;
}' HepG2_rawFiltered.txt | sort -k1,1n -k2,2n > HepG2consequence10.bed

awk -F '\t' '{
  split($1, a, "_"); # original format is "chr1_12345"
  chrom = substr(a[1], 4);  # remove "chr"
  start = a[2];
  end = start + 1;
  gene = $2; # extract gene name
  print chrom "\t" start "\t" end "\t" gene;
}' HeLa_rawFiltered.txt | sort -k1,1n -k2,2n > HeLaconsequence10.bed
```
Since the raw data did not include strand, we need to merge the data with a GTF file that contains the strand.

A python script was used to convert Ensembl's Homo_sapiens.GRCh38.113.chr.gtf file into a BED file with genomic coordinates per mRNA consequence for each gene, genomeConsequence0based.bed
```
python gtfToBed.py
```
Then, bedtools intersect was used to overlay Nm site BED files and genomeConsequence0based.bed.
```
bedtools intersect -a HEKconsequence10.bed -b genomeConsequence0based.bed -wa -wb > HEKconsequencetmp0based.bed
bedtools intersect -a HepG2consequence10.bed -b genomeConsequence0based.bed -wa -wb > HepG2consequencetmp0based.bed
bedtools intersect -a HeLaconsequence10.bed -b genomeConsequence0based.bed -wa -wb > HeLaconsequencetmp0based.bed
```
After bedtools, we added the strand to Nm site files by retaining Nm sites where the gene reported at the Nm site matched the GTF file gene. 
```
awk '$4 == $9' HEKconsequencetmp0based.bed > HEKlocation0based.bed
awk '$4 == $9' HepG2consequencetmp0based.bed > HepG2location0based.bed
awk '$4 == $9' HeLaconsequencetmp0based.bed > HeLalocation0based.bed
```
Since different transcripts from the same gene may have different mRNA consequences at an Nm site, this script only includes the most common consequence per Nm site. 
```
python addConsequence.py -i HEKlocation0based.bed -o HEKlocation0based_filtered.bed
python addConsequence.py -i HepG2location0based.bed -o HepG2location0based_filtered.bed
python addConsequence.py -i HeLalocation0based.bed -o HeLalocation0based_filtered.bed
```
We have added the strand and mRNA consequence information, but still need to reformat the BED file for downstream analysis.
```
awk 'NR > 1 { print $1, $2, $3 - 1, $4, $8, $10 }' OFS="\t" HEKlocation0based_filtered.bed > HEK0based_allSites.bed
awk 'NR > 1 { print $1, $2, $3 - 1, $4, $8, $10 }' OFS="\t" HeLalocation0based_filtered.bed > HeLa0based_allSites.bed
awk 'NR > 1 { print $1, $2, $3 - 1, $4, $8, $10 }' OFS="\t" HepG2location0based_filtered.bed > HepG20based_allSites.bed
```
These BED files can now be analyzed in parallel with the datasets from the other publications.

### NJU-seq (Tang et. al. 2024)
Table S2 (excel file) contains Nm sites detected in HEK293T, HeLa, and A549 cell lines. Each sheet contains the results for each cell line. Columns 'Chr', 'Position', 'Strand', 'Nm', 'ID', 'Name', and 'Distribution' were extracted from the original excel file into ```A549_Filtered.txt, HEK_Filtered.txt, and HeLa_Filtered.txt```. These files needed to be reformatted for downstream analysis.
```
awk 'BEGIN{OFS="\t"} {if ($6=="Forward") $6="+"; else if ($6=="Reverse") $6="-"; print}' A549_Filtered.txt | sort -k1,1 -k2,2n > A549_allSites.bed
awk 'BEGIN{OFS="\t"} {if ($6=="Forward") $6="+"; else if ($6=="Reverse") $6="-"; print}' HEK_Filtered.txt | sort -k1,1 -k2,2n > HEK_allSites.bed
awk 'BEGIN{OFS="\t"} {if ($6=="Forward") $6="+"; else if ($6=="Reverse") $6="-"; print}' HeLa_Filtered.txt | sort -k1,1 -k2,2n > HeLa_allSites.bed
```
The data can now be analyzed in parallel with the datasets from the other publications.

### Nm-mut-seq (Chen et. al. 2023)
All supplementary tables were combined into one excel file. Each sheet contains the results for each cell line. Sheets "S5_HeLa_mRNA_WT" and "S6_HepG2_mRNA_WT" contain all necessary genomic information. Columns 'chr', 'position', 'position2', 'mRNA segment', 'name', and 'strand' were extracted from the original excel file into ```HepG2_allSites.txt and HeLa_allSites.txt```. Files were then sorted for downstream analysis.
```
sort -k1,1 -k2,2n HepG2_allSites.txt > HepG2_allSites.bed
sort -k1,1 -k2,2n HeLa_allSites.txt > HeLa_allSites.bed
```
Now these files can be analyzed in parallel with the datasets from the other publications.

### Nm-seq (Dai et. al. 2018)
Nm site maps were downloaded from GSE90164. These BED files contained strand information, but not mRNA consequence. Each BED file contains the results for each cell line (HeLa and HEK293T). 

Bedtools intersect was used to add all mRNA consequences per Nm site. (see Nanopore-DRS section for the generation of genomeConsequence0based.bed)
```
bedtools intersect -a GSE90164_HeLamRNA.Nm.genome_sorted.bed -b genomeConsequence0based.bed -wa -wb -s > HeLaNm.hg38noLiftover_annotated.bed
bedtools intersect -a GSE90164_HEKmRNA.Nm.genome_sorted.bed -b genomeConsequence0based.bed -wa -wb -s > HEKNm.hg38noLiftover_annotated.bed
```
Since different transcripts from the same gene may have different mRNA consequences at an Nm site, this script only includes the most common consequence per Nm site. 
```
python extractConsequence.py -i HeLaNm.hg38noLiftover_annotated.bed -o HeLaNm.hg38noLiftover_filtered.bed
python extractConsequence.py -i HEKNm.hg38noLiftover_annotated.bed -o HEKNm.hg38noLiftover_filtered.bed
```
Now that we have added mRNA consequence, the data from this study can be analyzed in parallel with the datasets from the other publications.

## Extracting the sequence flanking Nm sites
Nm maps from all studies and cell lines were renamed to distinguish studies and cell lines as ```mlmHEKSites.bed, mlmHeLaSites.bed, mlmHepG2Sites.bed, njuA549Sites.bed, njuHeLaSites.bed, njuHEKSites.bed, nmHEKSites.bed, nmHeLaSites.bed, nmmtHeLaSites.bed, and nmmtHepG2Sites.bed```. Then, the sequence flanking Nm sites in the 5'UTR were extracted using a python script.
```
python extractSeq.py -i "$BED" -o "$BED_SEQ" --genome "Homo_sapiens.GRCh38.113"
```
Once the sequence was extracted, another script was used to identify Nm sites around near-cognate codons and categorize them by Nm position relative to the near-cognate codon.
```
python cognates.py -i "$BED_SEQ" -o "$BED_COGNATES"
```
## Prepping Nm files for downstream ribosome density extraction
To separate Nm sites around near-cognate codons per position relative to the start codon, a script was used to parse these files and output a BED file per Nm position.
```
python collectAllNmInfo.py
```
This results in ```all_1.bed, all_2.bed, all_3.bed, all_4.bed, all_m1.bed, all_m2.bed, and all_m3.bed``` which corresponds to Nm +1, +2, +3, +4, -1, -2, -3 in the near-cognate codon. These files can now be used for downstream processing with HR-Ribo-seq.
## HR-Ribo-seq Processing

1. Samples
The following datasets were downloaded from GEO (GSE162043, Arango et al, 2022 Mol Cell)

GSM4932133	HR_RPF_HeLa_WT_rep1 \
GSM4932134	HR_RPF_HeLa_WT_rep2

2. Pre-processing

Adapters were trimmed using cutadapt/4.2 and the following code:

```
# Run cutadapt in single-end read mode (to see what each parameter means, look at the cutadapt version 4.2 manual)
cutadapt --match-read-wildcards -e 0.1 -O 1 --quality-cutoff 10 -m 22 -a CTGTAGGCACCATCAAT -o $OUTPUT_R1 $R1 > $METRICS_FILE
```

3. Alignment

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
4. Normalization

BAM files were converted to BIGWIG format using a bin of 1 nucleotide and the RPKM normalization using deeptools/3.5.1

```
bamCoverage -b "$FILE" -o "$OUTPUT_DIR/${SAMPLE_NAME}.bw" --outFileFormat bigwig --normalizeUsing RPKM --binSize 1
```
## Extract the ribosome density
Run the following Python script to extract ribosome density around Nm positions.
```
BED_FILE="all_1.bed"
python ribo.density.py "$BIGWIG_FILE" "$BED_FILE" "$OUTPUT_DIR/${SAMPLE_NAME}_pos1.csv"

BED_FILE="all_2.bed"
python ribo.density.py "$BIGWIG_FILE" "$BED_FILE" "$OUTPUT_DIR/${SAMPLE_NAME}_pos2.csv"

BED_FILE="all_3.bed"
python ribo.density.py "$BIGWIG_FILE" "$BED_FILE" "$OUTPUT_DIR/${SAMPLE_NAME}_pos3.csv"

BED_FILE="all_4.bed"
python ribo.density.py "$BIGWIG_FILE" "$BED_FILE" "$OUTPUT_DIR/${SAMPLE_NAME}_pos4.csv"

BED_FILE="all_m1.bed"
python ribo.density.py "$BIGWIG_FILE" "$BED_FILE" "$OUTPUT_DIR/${SAMPLE_NAME}_posm1.csv"

BED_FILE="all_m2.bed"
python ribo.density.py "$BIGWIG_FILE" "$BED_FILE" "$OUTPUT_DIR/${SAMPLE_NAME}_posm2.csv"

BED_FILE="all_m3.bed"
python ribo.density.py "$BIGWIG_FILE" "$BED_FILE" "$OUTPUT_DIR/${SAMPLE_NAME}_posm3.csv"

```

## Rstudio analysis and plots
Run the R code DensityN.Rmd. The markdown is DensityN.html. The data sets needed to reproduce the plots are:

* MH85_pos1.csv
* MH85_pos2.csv
* MH85_posm1.csv
* MH85_posm2.csv
* MH85_posm3.csv
* MH87_pos1.csv
* MH87_pos2.csv
* MH87_posm1.csv
* MH87_posm2.csv
* MH87_posm3.csv
