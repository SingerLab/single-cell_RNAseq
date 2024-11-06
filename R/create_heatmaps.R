#' create pathway expression heatmap
#' @export
create_pathway_heatmaps <- function(cds, fit_deg, gsea, t2g, n.pathways = 10, ...) {
    
    hfx <- function(agg_mat, ...) {
        pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                           scale="column", clustering_method="ward.D2",
                           ...)
    }

    set.seed(2021)
    module_df <- t2g %>%
        filter(gs_name %in% (gsea %>% head(n = n.pathways) %>% rownames)) %>%
        mutate(id = gene_symbol,
               module = gs_name,
               supermodule = 1) %>%
        select(id, module, supermodule)
    
    cell_group_df <- tibble::tibble(cell = row.names(colData(cds)), 
                                    cell_group = cds$cluster)
    
    cell_part_df <- tibble::tibble(cell = row.names(colData(cds)),
                                   cell_group = cds$partition)
    
    cell_trt_df <- tibble::tibble(cell = row.names(colData(cds)),
                                  cell_group = cds$cycle)
    
    
    agg_mats <- lapply(list(cell_group_df, cell_part_df, cell_trt_df),
                       function(cg)
                           aggregate_gene_expression(cds, module_df, cg)
                       )
    xx <- lapply(agg_mats, hfx, fontsize = 10)
    
    return(xx)
    
}

#' create pathway expression heatmap
#' @export
create_gene_module_heatmaps <- function(cds, fit_deg, gsea,
                                        module_df = NULL, t2g, ...) {
    
    hfx <- function(agg_mat, ...) {
        pheatmap::pheatmap(agg_mat, cluster_rows = TRUE, cluster_cols = TRUE,
                           scale = "column", clustering_method="ward.D2",
                           ...)
    }

    set.seed(2021)
    if(is.null(module_df)) {
        module_df <- find_gene_modules(cds[fit_deg$gene_short_name,],
                                       resolution=1e-2)
        module_df$module <- paste("Module", module_df$module)
    } else {
        
        module_df <- module_df
    }
    
    cell_group_df <- tibble::tibble(cell = row.names(colData(cds)), 
                                    cell_group = cds$cluster)
    
    cell_part_df <- tibble::tibble(cell = row.names(colData(cds)),
                                   cell_group = cds$partition)
    
    cell_trt_df <- tibble::tibble(cell = row.names(colData(cds)),
                                  cell_group = cds$cycle)
    

    agg_mats <- lapply(list(cell_group_df, cell_part_df, cell_trt_df),
                       function(cg)
                           aggregate_gene_expression(cds, module_df, cg)
                       )
    xx <- lapply(agg_mats, hfx, fontsize = 10)

    return(xx)
        
}


#' create gene or pathway module matrices for biopsies
#'  the three matrices are built using three cell groupings
#'  
#' @export
create_gene_module_agg_matrices <- function(cds, fit_deg,
                                            module_df = NULL,  ...) {
    
    set.seed(2021)
    if(is.null(module_df)) {
        module_df <- find_gene_modules(cds[fit_deg$gene_short_name,],
                                       resolution=1e-2)
        module_df$module <- paste("Module", module_df$module)
        
    } else {

        module_df <- module_df
        
    }
    
    cell_group_df <- tibble::tibble(cell = row.names(colData(cds)), 
                                    cell_group = cds$cluster)
    
    cell_part_df <- tibble::tibble(cell = row.names(colData(cds)),
                                   cell_group = cds$partition)
    
    cell_trt_df <- tibble::tibble(cell = row.names(colData(cds)),
                                  cell_group = cds$cycle)
    

    agg_mats <- lapply(list(cell_group_df, cell_part_df, cell_trt_df),
                       function(cg)
                           aggregate_gene_expression(cds, module_df, cg)
                       )
    names(agg_mats) <- c("by.cluster", "by.partition", "by.trt")
    
    return(agg_mats)
        
}

module_matrix_from_pathway_genes <- function(gsea, t2g, n.pathways = 10, ...) {

    set.seed(2021)
    module_df <- t2g %>%
        filter(gs_name %in% (gsea %>% head(n = n.pathways) %>% rownames)) %>%
        mutate(id = gene_symbol,
               module = gs_name,
               supermodule = 1) %>%
        select(id, module, supermodule)

    return(module_df)
    
}
