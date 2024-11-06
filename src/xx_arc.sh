#!/bin/bash
#BSUB -n 12 -R "rusage[mem=7]" -W 64:00 -R "span[hosts=1]"
MEM_PER_JOB=$(echo $(printf "%d\n" ${LSB_CG_MEMLIMIT} )/1024^3 | bc ) ## $( echo $LSB_SUB_RES_REQ | sed -e 's/.*=//' -e 's/]//' )
MAX_MEM_GB=$(( MEM_PER_JOB*LSB_MAX_NUM_PROCESSORS - 1 ))
source /work/singer/opt/miniconda3/bin/activate single-cell-rnaseq

#<usage>
[[ $# -eq 3 ]] || {
echo "Usage:"
echo "This script expects the required arcuments to run cellranger_arc"
echo " in order:  output directory, path/to/fastq/ libraries.csv, and transcriptome directory"
echo ""
echo ""
echo 'bsub -n 16 -R "rusage[mem=7]" -W 48:00 -R "span[hosts=1]" ./src/01_arc.sh HCC3 hcc3_libraries.csv ML1105_IGO_09649_1 hg38_arc/'
echo ""
exit 1; }
#</usage>
set -e -x -o pipefail -u

CELLRANGER_ARC=/work/singer/opt/src/cellranger-arc-2.0.2/cellranger-arc

OUTDIR=$1            ## HCC3
LCSV=$2              ## hcc3_libraries.csv
TRANSCRIPTOME=$3     ## cellranger directory name

$CELLRANGER_ARC count --id=$OUTDIR \
		--reference=$TRANSCRIPTOME \
		--libraries=$LCSV \
		--localcores=${LSB_MAX_NUM_PROCESSORS} \
		--localmem=${MAX_MEM_GB}

