#' split seurat object by cells of interst
#' and covert to a cds
split_cells_cds <- function(seurat_obj, keep.cells, dims = 1:20, min_expr = 0.1) {
    cds <- seurat_obj[, keep.cells] %>%
        DietSeurat() %>%
        SCTransform(vars.to.regress = c("percent.mt")) %>%
        FindVariableFeatures() %>% RunPCA() %>%
        FindNeighbors(dims = dims) %>%
        FindClusters() %>%
        RunUMAP(dims = dims) %>%
        as.cell_data_set() %>%
        cluster_cells() %>%
        learn_graph()
    cds <- detect_genes(cds, min_expr = min_expr)
    cds <- estimate_size_factors(cds, round_exprs = FALSE)
    rowData(cds)$gene_short_name <- rownames(rowData(cds))
    pData(cds)$cluster <- clusters(cds)
    pData(cds)$partition <- partitions(cds)
    return(cds)
}

#' split seurat object by cells of interst
#' and covert to a SCE
split_cells_sce <- function(seurat_obj, keep.cells, dims = 1:20, min_expr = 0.1) {
    cds <- seurat_obj[, keep.cells] %>%
        DietSeurat() %>%
        SCTransform(vars.to.regress = c("percent.mt")) %>%
        FindVariableFeatures() %>% RunPCA() %>%
        FindNeighbors(dims = dims) %>%
        FindClusters() %>%
        RunUMAP(dims = dims) %>%
        as.SingleCellExperiment()
    return(cds)
}
