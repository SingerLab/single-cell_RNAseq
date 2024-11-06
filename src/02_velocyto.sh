#!/bin/bash
#BSUB -n 8 -R "rusage[mem=5]" -W 48:00 -R "span[hosts=1]"
MAX_MEM=$(( 5*LSB_MAX_NUM_PROCESSORS - 1 ))
## source /home/gularter/opt/miniconda3/bin/activate single-cell-rnaseq

#<usage>
[[ $# -eq 2 ]] || {
echo "Usage:"
echo "This script expects the path to the cellranger genome resource, and output directory."
echo " The cellranger genome resource to be linked in the on the same"
echo " directory from where it's launched, and must have <cellranger>/genes/genes.gtf."
echo ""
echo 'bsub -n 8 -R "rusage[mem=5]" -W 48:00 -R "span[hosts=1]" "./src/02_velocyto.sh cellranger/ ML1105_01/"'
echo ""
exit 1; }
#</usage>
set -e -x -o pipefail -u

TRANSCRIPTOME=$1
SAMPLE_DIR=$2

/work/singer/opt/miniconda3/bin/velocyto run10x \
	 -m ${TRANSCRIPTOME}/genes/hg19_rmsk.gtf \
	 -@ $LSB_MAX_NUM_PROCESSORS \
	 ${SAMPLE_DIR} \
	 ${TRANSCRIPTOME}/genes/genes.gtf 
