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
## CTAT19=
GTF_FILE=/work/singer/gularter/genomes/homo_sapiens/Ensembl/GRCh38.p5/Sequence/refdata-gex-GRCh38-and-mm10-2020-A/genes/genes.gtf
STARG=/work/singer/gularter/genomes/homo_sapiens/Ensembl/GRCh38.p5/Sequence/refdata-gex-GRCh38-and-mm10-2020-A/star/
STAR_TOOL=~/src/cellranger-3.0.2/STAR/5dda596/STAR 

## FASTQ directory
FASTQ_DIR=$1
EXTENSION=$2

## create fastq.gz arrays for R1 and R2 based on a path provided as a first argument
R1=($(ls ${FASTQ_DIR}/*_1.fastq.gz))
R2=($( echo ${R1[@]} | sed -e 's/_1.f/_2.f/' ))

## get name of sample :  STRAIN = legacy from mouse strains
STRAIN=$( basename ${R1[$LSB_JOBINDEX]} ${EXTENSION} )

## output directory variable
OUTDIR=star_out/${STRAIN}

## print information 
echo "Output directory/ "  $OUTDIR
echo "Case Name         "  $STRAIN
echo "LSF Array Task  ID"  $LSB_JOBINDEX
echo "R1 fastq.gz       "  ${R1[$LSB_JOBINDEX]}
echo "R2 fastq.gz       "  ${R2[$LSB_JOBINDEX]}
echo "Shell Used        "  $SHELL

## making output directory if it doesn't exist
[[ -d $OUTDIR ]] ||  mkdir -p $OUTDIR

$STAR_TOOL --genomeDir ${STARG} --readFilesCommand zcat \
     --readFilesIn ${R1[$LSB_JOBINDEX]} ${R2[$LSB_JOBINDEX]} \
     --runThreadN $LSB_MAX_NUM_PROCESSORS \
     --twopassMode Basic \
     --sjdbGTFfile $GTF_FILE \
     --outReadsUnmapped Fastq \
     --limitBAMsortRAM 31532137230 \
     --outFileNamePrefix $OUTDIR/ \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts

