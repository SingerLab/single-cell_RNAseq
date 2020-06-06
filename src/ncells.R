#!/home/gularter/opt/miniconda3/bin/Rscript

## parse the number of cells
args <- commandArgs(trailingOnly = TRUE)
## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      The R Script
 
      Arguments:
      --cellranger_metrics=path/to/outs/metrics_summary.csv       - character, path to varbin genome reference files
      --help                                                      - print this text
 
      Example:
      Rscript ncells.R --cellranger_metrics=path/to/outs/metrics_summary.csv
 ")
  q(save="no")
}

parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

if(! file.exists(argsL$cellranger_metrics)) {
    stop("Cell Ranger Metrics file ", argsL$cellranger_metrics, " not found")
    q(save = "no")
}


mtx <- read.csv(argsL$cellranger_metrics)

cat(as.numeric(gsub(",", "", as.character(mtx[,1]))))
