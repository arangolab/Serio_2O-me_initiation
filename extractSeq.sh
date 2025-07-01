#!/bin/bash

#SBATCH --account=b1042 # genomics allocation
#SBATCH --partition=genomics # genomics partition that has the correct amount resources for the script
#SBATCH --array=0-9
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --job-name=extractFlankSeq
#SBATCH --output=extractFlankSeq_%A_%a.out # a log of the job processing will be output in the .out file
#SBATCH --error=extractFlankSeq_%A_%a.err # any errors in your code will be printed to the .err file
#SBATCH --mail-user=hannah.serio@northwestern.edu # if you want email updates on your job, include lines 14 and 15 with your email
#SBATCH --mail-type=ALL # get notified when job starts, finishes, if it fails, etc.

# Need to specify the location of your input, output, gtf, and transcriptome file
INPUT_DIR="/projects/b1042/Arangolab/2ometh/5utrSequence/allSites"
OUTPUT_DIR="/projects/b1042/Arangolab/2ometh/5utrSequence/5utrGenome"
GENOME="/projects/b1042/Arangolab/2ometh/canonicalStartSites/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Get a list of all files in your input directory with the same ending into a variable called FILES
FILES=("$INPUT_DIR"/*Sites.bed)
# Get a specific bam file from the array into a variable
BED=${FILES["$SLURM_ARRAY_TASK_ID"]}

# Remove the previous extension of file name and add new suffix for deduplicated files
REMOVE_PART="Sites.bed"
OUTPUT="${OUTPUT_DIR}/$(basename "$BED" $REMOVE_PART)Seq.bed"

module load python-miniconda3/4.10.3

eval "$(conda shell.bash hook)"
conda activate /home/hba9267/.conda/envs/gffutils

python /projects/b1042/Arangolab/2ometh/5utrSequence/5utrGenome/extractSeq.py -i "$BED" -o "$OUTPUT" --genome "$GENOME"

conda deactivate

module purge