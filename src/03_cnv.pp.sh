#!/bin/bash
#BSUB -n 8 -M 64 -W 24:00

source /work/singer/opt/miniconda3/bin/activate single-cell-rnaseq


#<usage>
[[ $# -eq 3 ]] || {
echo "Usage:"
echo "  Required arguments:"
echo " database : a Seurat Object `database` housing all the data, "
echo "            must contain column name bioID"
echo " bioID : sample ID to subset and predict CNA"
echo " downsample : logical TRUE/FALSE, weather to downsample to 1000 cells"
echo ""
echo "Example:"
echo "bsub -n 8 -M 64 -W 24:00 ./src/03_cnv.pp.R lpsS DD4388 FALSE"
exit 1; }
#</usage>

set -e -x -o pipefail -u

DATABASE=$1
BIOID=$2
DOWNSAMPLE=$3

RSCRIPT=/work/singer/opt/miniconda3/envs/single-cell-rnaseq/bin/Rscript

$RSCRIPT ./src/03_cnv.pp.R --database=$DATABASE --sample.name=$BIOID --downsample=$DOWNSAMPLE

## rsync -avz output/${BIOID}* tmp/03.cnv/${BIOID}/
## rsync -avz output/*${BIOID}* tmp/03.cnv/${BIOID}/
