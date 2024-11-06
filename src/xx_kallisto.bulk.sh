#!/bin/bash
#BSUB -n 3 -R "rusage[mem=4]" -W 89 -R "span[hosts=1]"

#<usage>
[[ $# -eq 2 ]] || {
    echo "Usage:"
    echo "This script runs kallisto quant.  Arguments are passed in the follwoing order:"
    echo " output director, path/to/fastq/ \"sample_name\" and project name"
    echo "Run as : "
    echo 'bsub -n 8 -M 34 -W 89 -R "span[hosts=1]" "./src/01a_kallisto.sh ML1105_02k rawdata/Sample_ML1105_IGO_09649_1/ > log/ML1105_01/kallisto.log 2>&1"'
    exit 1; }
#</usage>

source /home/gularter/opt/miniconda3/bin/activate rnaseq
set -e -x -o pipefail -u

## KALLISTO_IDX=/ifs/depot/pi/resources/genomes/GRCh37/kallisto_index/gencode.v19
KALLISTO_HOME=/work/singer/gularter/genomes/homo_sapiens/Ensembl/GRCh37.p13/Sequence/kallisto_index/
TX_FASTA=${KALLISTO_HOME}/gencode.v30lift37.transcripts.dd.fa 
KALLISTO_IDX=${KALLISTO_HOME}/gencode.v30
GENCODE_V30=${KALLISTO_HOME}/gencode.v30lift37.annotation.gtf

R1=( $( find -L $1 -name "*R1*gz" | sort ) )
R2=( $( echo ${R1[$LSB_JOBINDEX]} | sed 's/_R1_/_R2_/' ) )
	    
OUTDIR=$2/
subdir=$( basename ${R1[$LSB_JOBINDEX]} .fastq.gz | sed 's/_IGO.*//' )

[[ -d ${OUTDIR}/${subdir} ]] || mkdir -p ${OUTDIR}/${subdir}

kallisto quant --fusion --genomebam --gtf ${GENCODE_V30}  --bias -b 100 -l -s --seed=81 -t $LSB_MAX_NUM_PROCESSORS  -i ${KALLISTO_IDX}  -o ${OUTDIR}/${subdir}  ${R1[$LSB_JOBINDEX]} ${R2}

FRG_LENGTH=$( pizzly_get_fragment_length.py ${OUTDIR}/${subdir}/abundance.h5 )

## optional flatten_json.py 

pizzly -k 31 --align-score 2 --output ${OUTDIR}/${subdir}/${subdir} --gtf ${GENCODE_V30} --fasta ${TX_FASTA} ${OUTDIR}/${subdir}/fusion.txt



