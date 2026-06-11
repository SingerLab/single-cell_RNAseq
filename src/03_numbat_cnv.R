#!/usr/bin/env Rscript
#' run_numbat.R - Command line tool for Numbat CNV analysis
#'
#' Runs Numbat copy number variation analysis on a single 10X sample.
#' Supports custom reference expression profiles via RDS files.
#'
#' @usage
#' Rscript run_numbat.R --sample.dir <dir> [--ref <path>] [--out_dir <dir>] [--min_LLR <num>] [--ncores <int>]
#'
#' @param sample.dir  Path to sample directory containing outs/filtered_feature_bc_matrix.h5
#' @param ref         Path to reference expression RDS file. Default: numbat::ref_hca
#' @param out_dir     Output directory. Default: <sample.dir>/cnv_<ref_name>
#' @param min_LLR     Minimum log-likelihood ratio for CNV calling. Default: 3
#' @param ncores      Number of cores for parallel processing. Default: 4
#'
#' @examples
#' # Using default ref_hca
#' Rscript run_numbat.R --sample.dir P14_C0_D0
#'
#' # Using custom sarcoma reference
#' Rscript run_numbat.R --sample.dir P14_C0_D0 --ref refs/ref_sarc.rds
#'
#' # Custom parameters
#' Rscript run_numbat.R --sample.dir P14_C0_D0 --ref refs/ref_sarc.rds --min_LLR 5 --ncores 12

## ---- Parse command line arguments -------------------------------------------
library(optparse)

option_list <- list(
    make_option("--sample.dir", type = "character", default = NULL,
                help = "Sample directory with 10X output [required]"),
    make_option("--ref", type = "character", default = NULL,
                help = "Path to reference expression RDS file [default: numbat ref_hca]"),
    make_option("--out_dir", type = "character", default = NULL,
                help = "Output directory [default: <sample.dir>/cnv_<ref_name>]"),
    make_option("--min_LLR", type = "numeric", default = 3,
                help = "Minimum log-likelihood ratio [default: %default]"),
    make_option("--ncores", type = "integer", default = 4,
                help = "Number of cores [default: %default]")
)

opt <- parse_args(OptionParser(
    option_list = option_list,
    description = "Run Numbat CNV analysis on a 10X single-cell sample",
    epilogue = paste("Output: Saves Numbat result object as .rds in out_dir.",
                     "Default reference is numbat::ref_hca.",
                     "Custom references should be gene x cell type matrices saved as .rds")))

## ---- Validate required arguments --------------------------------------------
if (is.null(opt$sample.dir)) stop("--sample.dir is required. Use --help for usage.")
if (!dir.exists(opt$sample.dir)) stop("Sample directory not found: ", opt$sample.dir)

## ---- Load reference ---------------------------------------------------------
if (is.null(opt$ref)) {
    ## Default: use numbat's built-in ref_hca
    ref_name <- "hca"
    message("Using default reference: numbat::ref_hca")
    suppressPackageStartupMessages(library(numbat))
    ref_expression <- ref_hca
} else {
    ## Custom reference from RDS
    if (!file.exists(opt$ref)) stop("Reference file not found: ", opt$ref)
    ref_name <- gsub("^ref_|\\.rds$", "", basename(opt$ref), ignore.case = TRUE)
    message("Loading custom reference: ", opt$ref, " (name: ", ref_name, ")")
    suppressPackageStartupMessages(library(numbat))
    ref_expression <- readRDS(opt$ref)
}

## ---- Set output directory ---------------------------------------------------
if (is.null(opt$out_dir)) {
    opt$out_dir <- file.path(opt$sample.dir, paste0("cnv_", ref_name))
}
if (!dir.exists(opt$out_dir)) dir.create(opt$out_dir, recursive = TRUE)

## ---- Load remaining libraries -----------------------------------------------
suppressPackageStartupMessages(library(Seurat))

## ---- Construct input paths --------------------------------------------------
sample.name <- basename(opt$sample.dir)
h5.mat <- file.path(opt$sample.dir, "outs/filtered_feature_bc_matrix.h5")
allele.file <- file.path(opt$sample.dir, "numbat",
                         paste0(sample.name, "_allele_counts.tsv.gz"))

## ---- Validate inputs exist --------------------------------------------------
if (!file.exists(h5.mat)) stop("H5 matrix not found: ", h5.mat)
if (!file.exists(allele.file)) stop("Allele file not found: ", allele.file)

## ---- Run Numbat -------------------------------------------------------------
message(Sys.time(), " - Running Numbat for: ", sample.name)
message("  Reference: ", ref_name)
message("  Output: ", opt$out_dir)
message("  min_LLR: ", opt$min_LLR, " | ncores: ", opt$ncores)

count_mat <- Read10X_h5(h5.mat)
df_allele <- read.delim(gzfile(allele.file))

out <- run_numbat(
    count_mat,
    ref_expression,
    df_allele,
    genome = "hg38",
    min_LLR = opt$min_LLR,
    t = 1e-5,
    ncores = opt$ncores,
    plot = TRUE,
    out_dir = opt$out_dir
)

## ---- Save results -----------------------------------------------------------
rds.file <- file.path(opt$out_dir, paste0(sample.name, "_numbat.rds"))
saveRDS(out, rds.file)
message(Sys.time(), " - Results saved to: ", rds.file)


## ==== base analysis
##% ## libraries
##% library(numbat)
##% library(tidyverse)
##% library(Seurat)
##% 
##% 
##% sample.dirs <- list.files(pattern = "^P")
##% 
##% 
##% numbat_batch <- function(sample.dir) {
##% 
##%     h5.mat <- file.path(sample.dir, "outs/filtered_feature_bc_matrix.h5")
##%     allele.file = file.path(sample.dir, "numbat",
##%                             paste(sample.dir, "allele_counts.tsv.gz", sep ="_"))
##%     
##%     count_mat <- Read10X_h5(h5.mat)
##%     df_allele <- read.delim(gzfile(allele.file))
##%     
##%     ## run
##%     out = run_numbat(
##%         count_mat,   # gene x cell integer UMI count matrix 
##%         ref_hca,     # reference expression profile, a gene x cell type normalized expression level matrix
##%         df_allele,   # allele dataframe generated by pileup_and_phase script
##%         genome = "hg38",
##%         min_LLR = 3,
##%         t = 1e-5,
##%         ncores = 22,
##%         plot = TRUE,
##%         out_dir = file.path(sample.dir, "numbat", "cnv")
##%     )
##%     return(out)
##% }
##% 
##% numbat.cna <- lapply(sample.dirs, numbat_batch)
##% 
##% ##% count_mat_P43 <- Read10X_h5("P43_C0_D0_01/outs/filtered_feature_bc_matrix.h5")
##% ##% df_allele_P43 <- read.delim(gzfile("test_/P43_C0_D0_01_allele_counts.tsv.gz"))
##% ##%                          
##% ##% # run
##% ##% out = run_numbat(
##% ##%     count_mat_P43, # gene x cell integer UMI count matrix 
##% ##%     ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix
##% ##%     df_allele_P43, # allele dataframe generated by pileup_and_phase script
##% ##%     genome = "hg38",
##% ##%     min_LLR = 3,
##% ##%     t = 1e-5,
##% ##%     ncores = 10,
##% ##%     plot = TRUE,
##% ##%     out_dir = './test'
##% ##% )



