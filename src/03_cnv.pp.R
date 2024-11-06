#!/work/singer/opt/miniconda3/envs/single-cell-rnaseq/bin/Rscript

## Collect arguments
args <- commandArgs(trailingOnly = TRUE)
#' for debugging:
#' args <- c("--database=lpsS", "--sample.name=LPS7819", "--downsample=TRUE")

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      The R/Varbin Tables Script
 
      Arguments:
      --database=lpsS        - character, name of data set containg 'sample name'
      --sample.name=bioID    - character, sample name
      --downsample=FALSE     - logical, sample name")
  q(save="no")
}


## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

argsL


## retrieving data into matrices
## Run parameters
database <- argsL$database
sample.name <- argsL$sample.name


## ----setup--------------------------------------------------------------------
suppressPackageStartupMessages({
    library(Seurat)
    library(SCEVAN)
})


np <- as.integer(system("echo $LSB_MAX_NUM_PROCESSORS", intern = TRUE))
if(np == "" | is.na(np)) {
    np <- 1
}

infDir <- "tmp/03.cnv"
if(! dir.exists(infDir)) {
    dir.create(infDir, recursive = TRUE)
}

## output directories, etc
out.path <- file.path(infDir, sample.name)
if(! dir.exists(out.path)) {
    dir.create(out.path, recursive = TRUE)
}


## ----data---------------------------------------------------------------------
## load complete data
load(file.path("data", paste0(database, '.rda')))

## subset for sample of interest
sample.so <- subset(eval(parse(text = database)),
                    bioID == sample.name)

print(sample.so)

rm(list = database)
gc(reset = TRUE)

if(argsL$downsample) {
    if(ncol(sample.so) > 1000) {
        dso <- sample(1:ncol(sample.so), 1000)
        sample.so <- sample.so[, dso]
        sample.name <- paste0(sample.name, "_downsampled_to_1000_cells")
    } else {
        sample.so <- sample.so
    }
}

## extract read counts
sample.rna.counts = as.matrix(GetAssayData(object = sample.so, slot = "counts"))

## set seed
set.seed(2023)
## run pipeline
rx1 <- pipelineCNA(sample.rna.counts,
                   sample = sample.name, 
                   organism = "human",
                   par_cores = np,
                   SUBCLONES = TRUE,
                   plotTree = TRUE)

## save predictions in tmp/dir as rds
out.file.name <- paste0(sample.name, "_scevan_predictions.rds")
saveRDS(rx1, file = file.path(out.path, out.file.name))

## files in default output/ dir
out.files <- list.files("./output/", pattern = sample.name,
                        full.names = TRUE)

## copy output files to tmp/dir/
for(i in out.files) {
    print(i)
    system(paste0("cp ", '"', i, '" ', out.path))
}


## __EOF__
