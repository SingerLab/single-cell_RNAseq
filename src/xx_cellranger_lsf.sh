#!/bin/bash
#BSUB -n 1 -R "rusage[mem=4]" -W 24:00 -R "span[hosts=1]"

#<usage>
[[ $# -eq 5 ]] || {
echo "Usage:"
echo "This script expects the required arcuments to run cellranger"
echo " in order:  output directory, path/to/fastq/ \"sample_name\", \"project name\", and transcriptome directory"
echo ""
echo ""
echo 'bsub -n 8 -R "rusage[mem=7]" -W 48:00 -R "span[hosts=1]" "./src/01_cellranger_lsf.sh ML1105_01 seqdata/Sample_ML1105_IGO_09649_1/ ML1105_IGO_09649_1  "ML1105" "cellranger" > log/ML1105_01/cellranger.log 2>&1"'
echo ""
exit 1; }
#</usage>
set -e -x -o pipefail -u

OUTDIR=$1            ## ML1105_01
FASTQ_DIR=$2         ## seqdata/Sample_ML1105_IGO_09649_1/
SAMPLE=$3            ## ML1105_IGO_09649_1
PROJECT_NAME=$4      ## "ML1105"
TRANSCRIPTOME=$5     ## cellranger directory name

## cell ranger version
CELLRANGER=/work/singer/opt/src/cellranger-7.1.0/cellranger

## will identify if counts are more precise at a later point
$CELLRANGER count --id=${OUTDIR} \
 	    --fastqs=${FASTQ_DIR} \
	    --sample=${SAMPLE} \
	    --transcriptome=${TRANSCRIPTOME} \
	    --project=${PROJECT_NAME} \
	    --jobmode lsf \
	    --mempercore 24 \
	    --maxjobs 200 \
	    --nopreflight \
	    --nosecondary 

