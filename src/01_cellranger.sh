#!/bin/bash
#BSUB -n 12 -R "rusage[mem=7]" -W 64:00 -R "span[hosts=1]"
MEM_PER_JOB=$(echo $(printf "%d\n" ${LSB_CG_MEMLIMIT} )/1024^3 | bc ) ## $( echo $LSB_SUB_RES_REQ | sed -e 's/.*=//' -e 's/]//' )
MAX_MEM=$(( ${MEM_PER_JOB} * $LSB_MAX_NUM_PROCESSORS - 1 ))

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
TRANSCRIPTOME=$5     ## cellranger directory name

## cell ranger version
CELLRANGER=/work/singer/opt/src/cellranger-7.1.0/cellranger
#% ~/src/cellranger-3.0.2/cellranger
#% /work/singer/opt/src/cellranger-6.0.1/cellranger
## differences between 3.0.2 and 6.0.1 are 10-30 less cells in 6.0.1; As of 23.Apr.2021 we keep 3.0.2
## will identify if counts are more precise at a later point
$CELLRANGER count --id=${OUTDIR} \
 	    --fastqs=${FASTQ_DIR} \
	    --sample=${SAMPLE} \
	    --transcriptome=${TRANSCRIPTOME} \
	    --project=${PROJECT_NAME} \
	    --localcores=${LSB_MAX_NUM_PROCESSORS} \
	    --localmem=${MAX_MEM} \
	    --nosecondary

[[ ! -d ${OUTDIR} ]] || rm -rf ${OUTDIR}/SC_RNA_COUNTER_CS/
[[ ! -d ${OUTDIR} ]] || rm -rf ${OUTDIR}/analysis/

## velocyto script
VELOCYTO_SCRIPT=./src/02_velocyto.sh

bsub -n 12 -R "rusage[mem=5]" -W 64:00 -R "span[hosts=1]" "$VELOCYTO_SCRIPT ${TRANSCRIPTOME} ${OUTDIR}"

