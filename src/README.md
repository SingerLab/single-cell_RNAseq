# src/:

This directory contains the single-cell RNA-seq pipeline deployment scripts.  We implemented a common script to use across all 10X libraries which runs cellranger and velocyto.

## scripts
src/_init.sh : Generates a 00_cellranger_submit.sh to generate bsub submissions for the libraries in a project

src/01_cellranger.sh : runs cell ranger count

src/01_vdj.sh : runs cell ranger vdj

src/02_velocyto.sh : runs velocyto

src/03_cnv.pp.sh : deploys src/03_cnv.pp.R
src/03_cnv.pp.R : runs SCEVAN pipeline in a seurat object

src/ncells.R : pulls the number of cells from a metrics_summary.csv (removes commans and gets arund the \")

src/README.md : this file

## experimental scripts
src/xx_arc.sh
src/xx_cellranger_lsf.sh
src/xx_kallisto.bulk.sh
src/xx_kallisto.sc.sh
src/xx_star.map.sh
src/xx_STARsolo.sh
