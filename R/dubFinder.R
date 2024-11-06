#' Doublet finding and subsequent removal
#' 
#' Removing duplicates with `DoubletFinder`.  Removal must be done on
#'   individual 10X libraries.  The `dubFinder()` goes through each
#'   10X library (mtx file) and runs both `Seurat`, and `DoubletFinder`
#'   complete pipelines.
#' 
#' Cells marked as _doublets_  can be picked up from the `dubFinder()`
#'   output named `s.data` on the copmlete using `pickup_dubs()`
#' 
#' @param a.data read count matrix from Read10X
#'
#' @param sample.name sample name
#'
#' @param pN.use  This defines the number of generated artificial doublets,
#'   expressed as a proportion of the merged real-artificial data. Default
#'   is set to 25%, based on observation that DoubletFinder performance is
#'   largely pN-invariant (McGinnis et. al, 2019, Cell Systems).
#'
#' @param PCs.use The number of statistically-significant principal components,
#'   specified as a range (e.g., PCs = 1:10)
#'
#' @param png.path file path to png figures
#'
#' @import ggplot2
#' @import dplyr
#' @import DoubletFinder
#' @import Seurat
#' 
#' @export
dubFinder <- function(a.data, sample.name,
                      pN.use = 0.24, PCs.use = 1:20,
                      min.genes = 200,
                      max.mt.pct = 20,
                      png.path = ".") {
    
    cat(paste(sample.name, "\n"))
    message(sample.name)
    
    ## creating seurat object
    seu_a <- CreateSeuratObject(counts = a.data,
                                project = gsub("_.*", "", sample.name))
    
    ## Variables to regress
    seu_a <- seu_a %>%
        PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt")

    ## Seurat pipeline
    if(is.null(max.mt.pct)) {
        max.mt.pct <- 100
    }
    
    seu_a <- seu_a %>%
        subset(subset = nFeature_RNA >= min.genes &
                   percent.mt <= max.mt.pct) %>%
        SCTransform(vars.to.regress = "percent.mt", verbose = TRUE) %>%
        FindVariableFeatures() %>%
        RunPCA() %>%
        RunUMAP(dims = 1:50) %>%
        FindNeighbors(dims = 1:50) %>%
        FindClusters()
    
    ## pK Identification (no ground-truth)
    sweep.res.list <- DoubletFinder::paramSweep_v3(seu_a, PCs = PCs.use, sct = TRUE)
    sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
    
    ## png of parameters
    pdf(file.path(png.path, paste0(sample.name, "_bcmvn_%03d.pdf")))
    bcmvn <- find.pK(sweep.stats)
    pK <- as.numeric(as.character(bcmvn$pK))
    plot(pK, bcmvn$BCmetric, pch = 16, type = "b",
         lty = 2, main = "The BCmvn distributions")
    pK.use <- pK[which.max(bcmvn$BCmetric)]
    abline(v = pK.use, lwd = 2, col = '#D40000', lty = 2)
    text(pK.use, max(bcmvn$BCmetric), as.character(pK.use),
         pos = 4, col = "#D40000")
    dev.off()
    
    ## Homotypic Doublet Proportion Estimate
    homotypic.prop <- DoubletFinder::modelHomotypic(seu_a@meta.data$seurat_clusters)
    nExp_poi <- round(0.08 * length(colnames(seu_a)))  ## Assuming 8% doublet
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    
    ## doublet finding initial run
    seu_a <- DoubletFinder::doubletFinder_v3(seu_a, PCs = PCs.use, pN = pN.use,
                              pK = pK.use, nExp = nExp_poi,
                              reuse.pANN = FALSE, sct = TRUE)

    ## rerun with adjusted poi
    ## run doublet finder
    pANN.a <- paste("pANN", pN.use, pK.use, nExp_poi, sep = "_")
    seu_a <- DoubletFinder::doubletFinder_v3(seu_a, PCs = PCs.use,
                                             pN = pN.use, pK = pK.use,
                              nExp = nExp_poi.adj,
                              reuse.pANN = pANN.a,
                              sct = TRUE)

    ## plot doublets and some features to asses doublet type
    DF.a <- paste("DF.classifications", pN.use, pK.use, nExp_poi, sep = "_")
    DF.b <- paste("DF.classifications", pN.use, pK.use, nExp_poi.adj, sep = "_")

    ## coalesce into single defined columns
    seu_a@meta.data$pANN.pN_pK_nExp <- paste(pN.use, pK.use, nExp_poi, sep = "_")
    seu_a@meta.data$pANN <- seu_a@meta.data[, pANN.a]
    seu_a@meta.data$df.params.high <- paste(pN.use, pK.use, nExp_poi, sep = "_")
    seu_a@meta.data$df.class.high <- seu_a@meta.data[, DF.a]
    seu_a@meta.data$df.params.low <- paste(pN.use, pK.use, nExp_poi.adj, sep = "_")
    seu_a@meta.data$df.class.low <- seu_a@meta.data[, DF.b]
    
    ## plot PNG of doublets
    pdf(file.path(png.path, paste0(sample.name, "_table_%03d.pdf")))
    gridExtra::grid.table(as.data.frame(table(seu_a@meta.data[, c(DF.a, DF.b)])))
    dev.off()
    
    pdf(file.path(png.path, paste0(sample.name, "_doublets_%03d.pdf")))
    Seurat::DimPlot(seu_a, group.by = DF.a) ##  + ggplot2::labs(title = paste(sample.name, DF.a))
    Seurat::DimPlot(seu_a, group.by = DF.b) ##  + ggplot2::labs(title = paste(sample.name, DF.b))
    dev.off()
    
    ## change orig.ident back to sm -- if not it exports only patient.id
    seu_a@meta.data$orig.ident <- sample.name
    
    return(seu_a)
} # end dubFinder

