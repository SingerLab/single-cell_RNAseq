#' prepare cds from a list of files
#'
#' @export
prepare_cds <- function(file.list, sample.names) {

    cds.list <- lapply(1:length(file.list), function(i) {
        cdsx <- cds_with_sample_info(file.list[i], sample.names[i])
        return(cdsx)
    })
    names(cds.list) <- sample.names
    
    pre.cds <- combine_cds(cds.list)

    ## estimating library size
    pData(pre.cds)$Library_Size  <- colSums(pre.cds@assays@data$counts)
    selected.columns <- c("barcode", "Size_Factor",
                          "Library_Size", "sample", "patient.id",
                          "cycle", "day", "region", "replicate")
                                                           
    pData(pre.cds) <- pData(pre.cds)[, selected.columns]
    
    ## detect minimum expressing genes
    pre.cds <- detect_genes(pre.cds, min_expr = 0.1)
    expressed_genes <- row.names(subset(fData(pre.cds), num_cells_expressed >= 10))

    pre.cds <- pre.cds[expressed_genes, ]

    return(pre.cds)
} ## end prepare cds


#' construct cds and include sample information
#' 
cds_with_sample_info <- function(path, sample.name, ...) {

    cds <- load_cellranger_data(path, ...)
    
    pData(cds)$sample.id <- sample.name

    sAnnot <- data.frame(do.call(rbind, strsplit(pData(cds)$sample.id, "_")))
    
    if(ncol(sAnnot) == 3) {
        names(sAnnot) <- c("patient.id", "cycle", "day")
    } else {
        if(ncol(sAnnot) == 4) {
            names(sAnnot) <- c("patient.id", "cycle", "day", "replicate")
        } else {
            if(ncol(sAnnot) == 5) {
                names(sAnnot) <- c("patient.id", "cycle", "day",
                                   "region", "replicate")
            }
        }
    }
    
    pData(cds) <- cbind(pData(cds), sAnnot)
    
    return(cds)
    
} ## end cds_with_sample_info
