#!/bin/bash
#BSUB -n 12 -R "rusage[mem=7]" -W 64:00 -R "span[hosts=1]"
MEM_PER_JOB=$(echo $(printf "%d\n" ${LSB_CG_MEMLIMIT} )/1024^3 | bc )
MAX_MEM=$(( ${MEM_PER_JOB} * $LSB_MAX_NUM_PROCESSORS - 1 ))

source /work/singer/opt/miniconda3/bin/activate single-cell-rnaseq

set -euox pipefail

usage() {
    echo "Usage: $0 -g <transcriptome_ref> -p <project_name> -o <output_dir> -s <sample_name> -i <fastq_path>"
    exit 1
}

while getopts ":g:p:o:s:i:" opt; do
    case ${opt} in
	g ) TRANSCRIPTOME=$OPTARG ;;  ## hg38
	p ) PROJECT_NAME=$OPTARG ;;   ## ML1105
	o ) OUTDIR=$OPTARG ;;         ## ML1105_01 
	s ) SAMPLE=$OPTARG ;;         ## ML1105_IGO_09649_1
	i ) FASTQ_DIR=$OPTARG ;;      ## seqdata/Sample_ML1105_IGO_09649_1/
	\? ) echo "Invalid option: -$OPTARG" >&2; usage ;;
	: ) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

if [ -z "${TRANSCRIPTOME:-}" ] || [ -z "${PROJECT_NAME:-}" ] || [ -z "${OUTDIR:-}" ] || [ -z "${SAMPLE:-}" ] || [ -z "${FASTQ_DIR:-}" ]; then
    usage
fi

echo "[INFO] Running Cell Ranger count for gene expression..."


## cell ranger version
##use v9.0.1 for TIL & MEK
CELLRANGER=/work/singer/opt/src/cellranger-9.0.1/bin/cellranger 
## /work/singer/opt/src/cellranger-7.1.0/cellranger    
#% ~/src/cellranger-3.0.2/cellranger
#% /work/singer/opt/src/cellranger-6.0.1/cellranger
## differences between 3.0.2 and 6.0.1 are 10-30 less cells in 6.0.1; As of 23.Apr.2021 we keep 3.0.2
## will identify if counts are more precise at a later point
$CELLRANGER count --create-bam true \
	    --id=${OUTDIR} \
	    --fastqs=${FASTQ_DIR} \
	    --sample=${SAMPLE} \
	    --transcriptome=${TRANSCRIPTOME} \
	    --project=${PROJECT_NAME} \
	    --localcores=${LSB_MAX_NUM_PROCESSORS} \
	    --localmem=${MAX_MEM} \
	    --nosecondary

echo "[INFO] Cell Ranger count completed for ${SAMPLE}."

[[ ! -d ${OUTDIR} ]] || rm -rf ${OUTDIR}/SC_RNA_COUNTER_CS/
[[ ! -d ${OUTDIR} ]] || rm -rf ${OUTDIR}/analysis/

## velocyto script
VELOCYTO_SCRIPT=./src/02_velocyto.sh

bsub -n 12 -R "rusage[mem=5]" -W 64:00 -R "span[hosts=1]" "$VELOCYTO_SCRIPT ${TRANSCRIPTOME} ${OUTDIR}"


##__EOF__
