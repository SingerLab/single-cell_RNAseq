#!/bin/bash
set -e -x -o pipefail -u

conda activate single-cell-rnaseq
BC_WHITELIST=$HOME/src/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/3M-february-2018.txt

~/src/STAR/bin/Linux_x86_64_static/STAR --genomeDir ~/cmo_genomes/GRCh37/star_index/ \
					--readFilesIn 
					--soloType Droplet --soloCBwhitelist $BC_WHITELIST \
					--runThreadN 8 --readFilesCommand zcat \
					--outSAMtype BAM SortedByCoordinate \
					--outFileNamePrefix star_solo/ML1105_

