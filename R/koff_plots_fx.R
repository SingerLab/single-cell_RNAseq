#' condition bubble plot
#' @param cds cds
#' @param gene_names gene of interst
#' @param group_cells_by grouping variable
#' @param normalize size factor normalization, default is TRUE
summarize_fq_fc <- function(cds, gene_names, group_cells_by, normalize = TRUE, ...) {

    if(!all(gene_names %in% rownames(cds))) {
        excl <- setdiff(gene_names, rownames(cds))
        gene_names <- intersect(gene_names, rownames(cds))
        warning(paste(excl, collapse = " "), " not available\n  using all other genes")
    } else {
        gene_names <- gene_names
    }
    
    cds_subset <- cds[gene_names, ]
    marker_exprs <- t(SingleCellExperiment::counts(cds_subset))
    if(normalize) {
        marker_exprs <- Matrix::t(Matrix::t(marker_exprs)/size_factors(cds_subset))
    } 
    grp <- pData(cds_subset)[, group_cells_by]
    
    pct.positive <- data.frame(grp, (as.matrix(marker_exprs) > 0)) %>%
        pivot_longer(cols = all_of(gene_names), names_to = "gene",
                     values_to = "count") %>%
        group_by(grp, gene) %>%
        summarise(mean = mean(count) * 100) %>%
        mutate(pct.cells.expressing = mean,
               grp = paste0(grp, ".fq")) %>%
        select(-mean) %>%
        pivot_wider(names_from = "grp", values_from = "pct.cells.expressing")

    gene.means <- data.frame(grp, as.matrix(marker_exprs)) %>%
        pivot_longer(cols = all_of(gene_names), names_to = "gene",
                     values_to = "count") %>%
        group_by(grp, gene) %>%
        mutate(grp = paste0(grp, ".mean")) %>%
        summarise(mean = mean(count)) %>%
        pivot_wider(names_from = "grp", values_from = "mean") 
    gene.means$fc <- (gene.means$Palbo.mean / gene.means$CTRL.mean)

    plot.dat <- left_join(pct.positive, gene.means)

    return(plot.dat)

}

koff_biomarker_plot <- function(cds, gene_names, group_cells_by, max.fc = 40, ...) {

    plot.dat <- summarize_fq_fc(cds = cds, gene_names = gene_names,
                                group_cells_by = group_cells_by, ...)

    plot.dat$fc[plot.dat$fc > max.fc] <- max.fc
    
    plot.dat %>%
        ggplot(aes(x = log10(CTRL.fq), y = log2(Palbo.fq),
                   size = fc, colour = fc,
                   label = gene)) +
        geom_point() +
        ggrepel::geom_text_repel() +
        viridis::scale_colour_viridis() +
        theme_light() +
        scale_x_continuous(name = bquote(~ log[10]~"(Percent Positive Untreated)"),
                           breaks = c(-1, 0, 1, 2),
                           labels = c(0.1, 1, 10, 100)) +
        scale_y_continuous(name = bquote(~ log[2]~"(Percent Positive Palbociclib)"),
                           breaks = c(1:6, log2(100)),
                           labels = c(2, 4, 8, 16, 32, 64, 100))
}


#' condition bubble plot
#' @param cds cds
#' @param gene_names gene of interst
#' @param group_cells_by grouping variable
#' @param normalize size factor normalization, default is TRUE
expression_matrix_from_cds <- function(cds, gene_names, normalize = TRUE, ...) {

    if(!all(gene_names %in% rownames(cds))) {
        excl <- setdiff(gene_names, rownames(cds))
        gene_names <- intersect(gene_names, rownames(cds))
        warning(paste(excl, collapse = " "), " not available\n  using all other genes")
    } else {
        gene_names <- gene_names
    }
    
    cds_subset <- cds[gene_names, ]
    marker_exprs <- t(SingleCellExperiment::counts(cds_subset))
    if(normalize) {
        marker_exprs <- Matrix::t(Matrix::t(marker_exprs)/size_factors(cds_subset))
    } 

    return(marker_exprs)
}


## ====seurat_koff_plot========================================================
#' condition bubble plot
#' @param seu seu
#' @param gene_names gene of interst
#' @param group_cells_by grouping variable
#' @param normalize size factor normalization, default is TRUE
summarize_fq_fc_seu <- function(seu, gene_names = c("TUBB", "HPRT1", "GAPDH"),
                                group_cells_by = "orig.ident",
                                normalize = TRUE, ...) {

    if(!all(gene_names %in% rownames(seu))) {
        excl <- setdiff(gene_names, rownames(seu))
        gene_names <- intersect(gene_names, rownames(seu))
        warning(paste(excl, collapse = " "), " not available\n  using all other genes")
    } else {
        gene_names <- gene_names
    }
    
    seu_subset <- seu[gene_names, ]
    marker_exprs <- t(SingleCellExperiment::counts(seu_subset))
    if(normalize) {
        marker_exprs <- Matrix::t(Matrix::t(marker_exprs)/size_factors(seu_subset))
    } 
    grp <- pData(seu_subset)[, group_cells_by]
    
    pct.positive <- data.frame(grp, (as.matrix(marker_exprs) > 0)) %>%
        pivot_longer(cols = all_of(gene_names), names_to = "gene",
                     values_to = "count") %>%
        group_by(grp, gene) %>%
        summarise(mean = mean(count) * 100) %>%
        mutate(pct.cells.expressing = mean,
               grp = paste0(grp, ".fq")) %>%
        select(-mean) %>%
        pivot_wider(names_from = "grp", values_from = "pct.cells.expressing")

    gene.means <- data.frame(grp, as.matrix(marker_exprs)) %>%
        pivot_longer(cols = all_of(gene_names), names_to = "gene",
                     values_to = "count") %>%
        group_by(grp, gene) %>%
        mutate(grp = paste0(grp, ".mean")) %>%
        summarise(mean = mean(count)) %>%
        pivot_wider(names_from = "grp", values_from = "mean") 
    ## gene.means$fc <- (gene.means$Palbo.mean / gene.means$CTRL.mean)

    plot.dat <- left_join(pct.positive, gene.means)

    return(plot.dat)

}

koff_biomarker_plot <- function(seu, gene_names, group_cells_by, max.fc = 40, ...) {

    plot.dat <- summarize_fq_fc(seu = seu, gene_names = gene_names,
                                group_cells_by = group_cells_by, ...)

    plot.dat$fc[plot.dat$fc > max.fc] <- max.fc
    
    plot.dat %>%
        ggplot(aes(x = log10(CTRL.fq), y = log2(Palbo.fq),
                   size = fc, colour = fc,
                   label = gene)) +
        geom_point() +
        ggrepel::geom_text_repel() +
        viridis::scale_colour_viridis() +
        theme_light() +
        scale_x_continuous(name = bquote(~ log[10]~"(Percent Positive Untreated)"),
                           breaks = c(-1, 0, 1, 2),
                           labels = c(0.1, 1, 10, 100)) +
        scale_y_continuous(name = bquote(~ log[2]~"(Percent Positive Palbociclib)"),
                           breaks = c(1:6, log2(100)),
                           labels = c(2, 4, 8, 16, 32, 64, 100))
}


#' condition bubble plot
#' @param seu seu
#' @param gene_names gene of interst
#' @param group_cells_by grouping variable
#' @param normalize size factor normalization, default is TRUE
expression_matrix_from_seu <- function(seu, gene_names, normalize = TRUE, ...) {

    if(!all(gene_names %in% rownames(seu))) {
        excl <- setdiff(gene_names, rownames(seu))
        gene_names <- intersect(gene_names, rownames(seu))
        warning(paste(excl, collapse = " "), " not available\n  using all other genes")
    } else {
        gene_names <- gene_names
    }
    
    seu_subset <- seu[gene_names, ]
    marker_exprs <- t(SingleCellExperiment::counts(seu_subset))
    if(normalize) {
        marker_exprs <- Matrix::t(Matrix::t(marker_exprs)/size_factors(seu_subset))
    } 

    return(marker_exprs)
}

