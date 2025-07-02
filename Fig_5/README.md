# 2’O-methylation inhibits upstream translation initiation

## Samples
The following datasets were downloaded from GEO (GSE162043)

GSM4932133	HR_RPF_HeLa_WT_rep1
GSM4932134	HR_RPF_HeLa_WT_rep2

## Pre-processing

Adapters were trimmed using cutadapt/4.2 and the following code:

```
# Run cutadapt in single-end read mode (to see what each parameter means, look at the cutadapt version 4.2 manual)
cutadapt --match-read-wildcards -e 0.1 -O 1 --quality-cutoff 10 -m 22 -a CTGTAGGCACCATCAAT -o $OUTPUT_R1 $R1 > $METRICS_FILE
```

