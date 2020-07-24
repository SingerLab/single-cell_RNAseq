#!/bin/bash
#BSUB -n 8 -R "rusage[mem=5]" -W 48:00 -R "span[hosts=1]"
MAX_MEM=38
source /home/gularter/opt/miniconda3/bin/activate single-cell-rnaseq

#<usage>
[[ $# -eq 4 ]] || {
echo "Usage:"
echo "This script expects the required arcuments to run cellranger"
echo " in order:  output director, path/to/fastq/ \"sample_name\" and project name"
echo ""
echo ""
echo 'bsub -n 8 -R "rusage[mem=5]" -W 48:00 -R "span[hosts=1]" "./src/01_cellranger.sh ML1105_01 rawdata/Sample_ML1105_IGO_09649_1/ ML1105_IGO_09649_1  "ML1105" > log/ML1105_01/cellranger.log 2>&1"'
echo ""
exit 1; }
#</usage>
set -e -x -o pipefail -u

OUTDIR=$1            ## ML1105_01
FASTQ_DIR=$2         ## rawdata/Sample_ML1105_IGO_09649_1/
SAMPLE=$3            ## ML1105_IGO_09649_1
PROJECT_NAME=$4      ## "ML1105"

~/src/cellranger-3.0.2/cellranger count --id=${OUTDIR} \
				  --fastqs=${FASTQ_DIR} \
				  --sample=${SAMPLE} \
				  --transcriptome=cellranger \
				  --project=${PROJECT_NAME} \
				  --localcores=${LSB_MAX_NUM_PROCESSORS} \
				  --localmem=${MAX_MEM}
