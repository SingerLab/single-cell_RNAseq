#!/bin/bash
#BSUB -n 8 -R "rusage[mem=5]" -R "span[hosts=1]" -W 32:00 -J StarAln

#<usage>
[[ $# -gt 0 ]] || {
echo "Usage:"
echo "This script expects an input directory name as the first argument"
echo "and optionally an LSB_JOBINDEX as the second argument"
echo ""
echo "Examples to launch 100 jobs"
echo 'bsub -n 8 -R "rusage[mem=5]" -J "STAR[1-100]" -W 2:00 -R "span[hosts=1]" ./starFusion.sh ../path/to/fastq/'
echo 'for i in {0..100} ; do bsub -n 8 -R "rusage[mem=5]" -J "STAR[1-100]" -W 2:00 -R "span[hosts=1]" ./starFusion.sh ../path/to/fastq/ $i ; done'
echo 'for i in {0..100} ; do bsub -n 8 -R "rusage[mem=5]" -J "STAR[1-100]" -W 2:00 -R "span[hosts=1]" LSB_JOBINDEX=$i ./starFusion.sh ../path/to/fastq/ ; done'
echo ""
exit 1; }
#</usage>


## set to end bash if anything goes wrong
## set to print a copy of each command as script runss
set -e -x -o pipefail -u

## resources
CTAT19=
STARG=/ifs/depot/pi/resources/genomes/GRCh37/star_index/withoutJunctions/
GENCODE19=/ifs/depot/pi/resources/genomes/GRCh37/gencode/default/gencode.v19.annotation.gtf

## FASTQ directory
FASTQ_DIR=$1
EXTENSION=$2

## create fastq.gz arrays for R1 and R2 based on a path provided as a first argument
R1=(`ls ${FASTQ_DIR}/*R1_001.fastq.gz`)
R2=(`ls ${FASTQ_DIR}/*R3_001.fastq.gz`)

## get name of sample :  STRAIN = legacy from mouse strains
STRAIN=$( basename ${R1[$LSB_JOBINDEX]} ${EXTENSION} )

## output directory variable
OUTDIR=star/${STRAIN}

## print information 
echo "Output directory/ "  $OUTDIR
echo "Case Name         "  $STRAIN
echo "LSF Array Task  ID"  $LSB_JOBINDEX
echo "R1 fastq.gz       "  ${R1[$LSB_JOBINDEX]}
echo "R2 fastq.gz       "  ${R2[$LSB_JOBINDEX]}
echo "Shell Used        "  $SHELL

## making output directory if it doesn't exist
[[ -d $OUTDIR ]] ||  mkdir -p $OUTDIR

#STAR --genomeDir $STARG --readFilesCommand zcat --outReadsUnmapped Fastq \
#    --readFilesIn ${R1[$LSB_JOBINDEX]} ${R2[$LSB_JOBINDEX]} --runThreadN $NSLOTS \
#    --outSAMstrandField intronMotif \
#    --outSAMattrRGline ID:${NAME} SM:${STRAIN}_${NAME}_${POOL} LB:${LIBTYPE} CN:"MSKCC-IGO" DS:"SalivaryGlandACC" PL:Illumina PG:STAR \
#    --outFileNamePrefix $OUTDIR/ \
#    --outSAMtype BAM SortedByCoordinate --chimOutType WithinBAM


STAR --genomeDir ${STARG} --readFilesCommand zcat \
     --readFilesIn ${R1[$LSB_JOBINDEX]} ${R2[$LSB_JOBINDEX]} --runThreadN $LSB_MAX_NUM_PROCESSORS \
     --twopassMode Basic \
     --sjdbGTFfile $GENCODE19 \
     --outReadsUnmapped Fastq \
     --limitBAMsortRAM 31532137230 \
     --outFileNamePrefix $OUTDIR/ \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts

##     --outSAMattrRGline ID:${STRAIN} SM:${STRAIN}_01 LB:"RNAseq" CN:"MSKCC-IGO" DS:"SalivaryGlandACC" PL:Illumina PG:STAR \
