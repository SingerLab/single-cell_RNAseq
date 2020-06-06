#!/bin/bash
#BSUB -n 8 -M 34 -W 48:00 -R "span[hosts=1]"

#<usage>
[[ $# -eq 2 ]] || {
    echo "Usage:"
    echo "This script runs kallisto bus.  Arguments are passed in the follwoing order:"
    echo " output director, path/to/fastq/ \"sample_name\" and project name"
    echo "Run as : "
    echo 'bsub -n 8 -M 34 -W 89 -R "span[hosts=1]" "./src/01a_kallisto.sh ML1105_02k rawdata/Sample_ML1105_IGO_09649_1/ > log/ML1105_01/kallisto.log 2>&1"'
    exit 1; }
#</usage>

source /home/gularter/opt/miniconda3/bin/activate single-cell-rnaseq
set -e -x -o pipefail -u

## KALLISTO_IDX=/ifs/depot/pi/resources/genomes/GRCh37/kallisto_index/gencode.v19
KALLISTO_IDX=/work/singer/gularter/genomes/homo_sapiens/Ensembl/GRCh37.p13/Sequence/kallisto_index/gencode.v30
OUTDIR=$1
FASTQ_FILES=( $( find $2 -name "*R*gz" | sort ) )

kallisto bus -i ${KALLISTO_IDX} -o ${OUTDIR} -x 10xv3 -t $LSB_MAX_NUM_PROCESSORS ${FASTQ_FILES[@]}

