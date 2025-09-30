#!/bin/bash
#BSUB -n 12 -R "rusage[mem=7]" -W 64:00 -R "span[hosts=1]"
MEM_PER_JOB=$(echo $(printf "%d\n" ${LSB_CG_MEMLIMIT} )/1024^3 | bc )
MAX_MEM=$(( ${MEM_PER_JOB} * $LSB_MAX_NUM_PROCESSORS - 1 ))

## source /work/singer/opt/miniconda3/bin/activate single-cell-rnaseq

#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 -g <vdj_reference> -p <project_name> -o <output_dir> -s <sample_name> -i <fastq_path>"
    exit 1
}

while getopts ":g:p:o:s:i:" opt; do
    case ${opt} in
	g ) REFERENCE=$OPTARG ;;      ## cellranger VDJ directory name 
	p ) PROJECT_NAME=$OPTARG ;;   ## ML1105
	o ) OUTDIR=$OPTARG ;;         ## ML1105_01
	s ) SAMPLE=$OPTARG ;;         ## ML1105_IGO_09649_1
	i ) FASTQ_DIR=$OPTARG ;;      ## seqdata/Sample_ML1105_IGO_09649_1/
	\? ) echo "Invalid option: -$OPTARG" >&2; usage ;;
	: ) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

if [ -z "${REFERENCE:-}" ] || [ -z "${PROJECT_NAME:-}" ] || [ -z "${OUTDIR:-}" ] || [ -z "${SAMPLE:-}" ] || [ -z "${FASTQ_DIR:-}" ]; then
    usage
fi

echo "[INFO] Running Cell Ranger VDJ for immune profiling..."

## cell ranger version
CELLRANGER=/work/singer/opt/src/cellranger-7.1.0/cellranger

## check if VDJ reference or convert
if [[ "$REFERENCE" != *_vdj ]]; then
    REFERENCE="${REFERENCE}_vdj"
fi

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

echo "[INFO] Cell Ranger VDJ completed for ${SAMPLE}."

## __EOF__

