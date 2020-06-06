#!/bin/bash
#BSUB -n 8 -R "rusage[mem=6] -W 48:00 -R "span[hosts=1]"
MAX_MEM=48
source /home/gularter/opt/miniconda3/bin/activate single-cell-rnaseq

#<usage>
[[ $# -eq 7 ]] || {
echo "Usage:"
echo "This script expects the required arcuments to run cellranger"
echo " in order:  output director, path/to/fastq/ \"sample_name\", project name, and 'chemistry'"
echo ""
echo ""
echo 'bsub -n 8 -R "rusage[mem=6]" -W 48:00 -R "span[hosts=1]" "./src/01_gex_vdj.sh '
echo ""
exit 1; }
#</usage>
## setting pipeline parameters
set -e -x -o pipefail -u

OUTDIR=$1
VDJ_FASTQ_DIR=$2
GEX_FASTQ_DIR=$3
VDJ_SAMPLE=$4
GEX_SAMPLE=$5
PROJECT_NAME=$6
CHEM=$7

~/src/cellranger-3.0.2/cellranger count --id=GEX_VDJ_${OUTDIR}_${CHEM}_GEX \
				  --fastqs=${GEX_FASTQ_DIR} \
				  --sample=${GEX_SAMPLE} \
				  --transcriptome=cellranger \
				  --project=${PROJECT_NAME} \
				  --chemistry=${CHEM} \
				  --localcores=${LSB_MAX_NUM_PROCESSORS} \
				  --localmem=${MAX_MEM}

## parsing out the number of cells to force onto vdj
NCELLS=$( Rscript src/ncells.R --cellranger_metrics=GEX_VDJ_${OUTDIR}_${CHEM}_GEX/outs/metrics_summary.csv )

~/src/cellranger-3.0.2/cellranger vdj --id=GEX_VDJ_${OUTDIR}_${CHEM}_Im \
				  --fastqs=${VDJ_FASTQ_DIR} \
				  --sample=${VDJ_SAMPLE} \
				  --force-cells=${NCELLS} \
				  --reference=vdj --chain='all'\
				  --project=${PROJECT_NAME} \
				  --localcores=${LSB_MAX_NUM_PROCESSORS} \
				  --localmem=${MAX_MEM}


## B-cell specific
## ~/src/cellranger-3.0.2/cellranger vdj --id=GEX_VDJ_${OUTDIR}_B \
##				  --fastqs=${FASTQ_DIR} \
##				  --sample=${SAMPLE} \
##				  --reference=vdj/ \
##				  --localcores=${LSB_MAX_NUM_PROCESSORS} \
##				  --localmem=${MAX_MEM}

