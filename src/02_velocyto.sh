#!/bin/bash
#BSUB -n 8 -R "rusage[mem=5]" -W 48:00 -R "span[hosts=1]"
MAX_MEM=38
## source /home/gularter/opt/miniconda3/bin/activate single-cell-rnaseq


#<usage>
[[ $# -eq 1 ]] || {
echo "Usage:"
echo "This script expects the path to the cell ranger output as the only argument."
echo " It does need the cellranger/ genome resource to be linked in the on the same"
echo " directory from where it's launched, and must have  cellranger/genes/genes.gtf."
echo ""
echo 'bsub -n 8 -R "rusage[mem=5]" -W 48:00 -R "span[hosts=1]" "./src/02_velocyto.sh ML1105_01/"'
echo ""
exit 1; }
#</usage>
set -e -x -o pipefail -u

SAMPLE_DIR=$1

/work/singer/opt/miniconda3/bin/velocyto run10x \
	 -m cellranger/genes/hg19_rmsk.gtf \
	 -@ $LSB_MAX_NUM_PROCESSORS \
	 ${SAMPLE_DIR} \
	 cellranger/genes/genes.gtf 
