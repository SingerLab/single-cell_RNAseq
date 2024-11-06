#!/bin/bash
#BSUB -n 12 -R "rusage[mem=7]" -W 64:00 -R "span[hosts=1]"
MAX_MEM=$(( 5*LSB_MAX_NUM_PROCESSORS - 1 ))
source /work/singer/opt/miniconda3/bin/activate single-cell-rnaseq

#<usage>
[[ $# -eq 5 ]] || {
echo "Usage:"
echo "This script expects the required arcuments to run cellranger"
echo " in order:  output directory, path/to/fastq/ \"sample_name\", \"project name\", and transcriptome directory"
echo ""
echo ""
echo 'bsub -n 8 -R "rusage[mem=7]" -W 48:00 -R "span[hosts=1]" "./src/01_cellranger.sh ML1105_01 seqdata/Sample_ML1105_IGO_09649_1/ ML1105_IGO_09649_1  "ML1105" "cellranger" > log/ML1105_01/cellranger.log 2>&1"'
echo ""
exit 1; }
#</usage>
set -e -x -o pipefail -u

OUTDIR=$1            ## ML1105_01
FASTQ_DIR=$2         ## seqdata/Sample_ML1105_IGO_09649_1/
SAMPLE=$3            ## ML1105_IGO_09649_1
PROJECT_NAME=$4      ## "ML1105"
REFERENCE=$5     ## cellranger directory name

## cell ranger version
CELLRANGER=/work/singer/opt/src/cellranger-7.1.0/cellranger

## reference analysis
$CELLRANGER vdj --id=${OUTDIR} \
	    --reference=${REFERENCE} \
 	    --fastqs=${FASTQ_DIR} \
	    --sample=${SAMPLE} \
	    --project=${PROJECT_NAME} \
	    --localcores=${LSB_MAX_NUM_PROCESSORS} \
	    --localmem=${MAX_MEM} 

## denovo analysis
$CELLRANGER vdj --id=${OUTDIR}_DENOVO \
	    --reference=${REFERENCE} \
 	    --fastqs=${FASTQ_DIR} \
	    --sample=${SAMPLE} \
	    --project=${PROJECT_NAME} \
	    --denovo \
	    --localcores=${LSB_MAX_NUM_PROCESSORS} \
	    --localmem=${MAX_MEM} 


## [[ ! -d ${OUTDIR} ]] || rm -rf ${OUTDIR}/SC_RNA_COUNTER_CS/
## [[ ! -d ${OUTDIR} ]] || rm -rf ${OUTDIR}/analysis/


