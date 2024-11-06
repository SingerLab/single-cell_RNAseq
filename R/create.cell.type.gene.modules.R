#' create cell type gene modules object
#'
#' @export
create.cell.type.gene.modules <- function() {

    MacMono.markers <-     c("CD14", "CD68", "CD33", "TLR7",
                             "TLR9", "ITGAX", "FLT3", "CCL18")
    Monocyte.markers <-    c("CD14", "FCGR1A", "CD68", "S100A12")
    Macrophage <-          c("CD163")
    DC.markers <-          c("IL3RA", "CD1C", "BATF3", "THBD", "CD209")
    Bcell.markers <-       c("CD19", "MS4A1", "CD79A")
    tcell.markers <-       c("CD3D", "CD3E", "CD3G", "PTPRC")
    CD4.tcell.markers <-   c("CD4", "FOXP3", "IL2RA", "IL7R")
    CD8.tcell.markers <-   c("CD8A", "CD8B")
    NKcell.markers <-      c("NCAM1", "FCGR3A", "CD160", "ITGA2", "KLRB1")
    Endothelial.markers <- c("CDH5", "LYVE1", "TEK", "KDR", "RAMP2", "FLT1", "SELE")
    Fibroblast.markers <-  c("DES", "PDGFRA", "PDGFRB", "FAP")
    Granulocyte.markers <- c("CEACAM8", "FUT4", "SIGLEC8", "CCR3", "IL5RA")
    tumor.cell.markers <-  c("CDK4", "MDM2", "HMGA2", "JUN")
    
    cell.type.gene.modules <-
        list(MacMono.markers, Monocyte.markers, Macrophage,
             DC.markers, Bcell.markers,
             tcell.markers, CD4.tcell.markers,
             CD8.tcell.markers,
             NKcell.markers,
             Endothelial.markers,
             Fibroblast.markers,
             Granulocyte.markers)
    
    names(cell.type.gene.modules) <-
        paste0(c("MacMono", "Monocyte", "Macrophage",
                 "DC", "B.cell",
                 "T.cell", "CD4.tcell",
                 "CD8.tcell",
                 "NK.cell",
                 "Endothelial",
                 "Fibroblast",
                 "Granulocyte"), "_score")
    
    return(cell.type.gene.modules)
    
} # end create.cell.type.gene.modules

#' read hcdm markers and lps relevant genes
#'
#' @export
read_hcdm <- function() {
    hcdm.markers <- read.delim("../res/markers/hcdm.markers.txt", as.is = TRUE)
    cell.markers <- read.delim("../res/markers/cell.type.markers.txt", as.is = TRUE)
    immune.markers <- read.delim("../res/markers/immune.markers.txt", as.is = TRUE)
    ddls.genes <- read.delim("../res/markers/ddls.markers.txt", as.is = TRUE)
    fused.genes <- read.delim("../res/markers/fused.genes.txt", sep = "\t", as.is = TRUE)
    szabo.tcells <- as.list(read.delim("../res/sigs/szabo.tcell.txt",
                                       skip = 2, as.is = TRUE))
    
    azizi.mm <- sapply(c("../res/sigs/azizi.M1.txt", "../res/sigs/azizi.M2.txt"),
                       function(mm) readLines(mm)[-c(1:3)])
    names(azizi.mm) <- c("M1", "M2")                 
    
    azizi.tcells <- lapply(
        as.list(read.delim("../res/sigs/azizi.tcell.txt",
                           skip = 2, as.is = TRUE, na.strings = "")),
        function(ii) as.character(na.omit(ii)))

    azizi.cells <- c(azizi.mm, azizi.tcells)
                               
    usethis::use_data(hcdm.markers, cell.markers, immune.markers,
                      ddls.genes, fused.genes, szabo.tcells, azizi.cells,
                      overwrite = TRUE)

}
