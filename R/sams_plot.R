#' ----plots_for_sam------------------------------------------------------------
sams_plot_s9 <- function(cds, seu, s1_trt = "trt", s1_colors = NULL,
                         s2_clust = "cluster", s2_colors = NULL,
                         s3_phase = "Phase", s3_colors = NULL,
                         s4_gene1 = "CDKN2A",
                         s5_gene2 = "ANGPTL4",
                         s6_score1 = "SASP_Score",
                         s7_pct.mt = "percent.mt") {

    if(is.null(s1_colors)) {
        s1_colors <- pals::alphabet()[1:length(unique(pData(cds)[,s1_trt]))]
        names(s1_colors) <- unique(pData(cds)[, s1_trt])
    }
    
    S1 <- simple_plot_cells(cds, color_cells_by = s1_trt,
                            show_trajectory_graph = FALSE) +
        scale_colour_manual(name = s1_trt, values = s1_colors) +
        theme_light() + labs(title = s1_trt) +
        theme(axis.title = element_blank())
    
    if(is.null(s2_colors)) {
        s2_colors <- pals::alphabet()[1:length(unique(pData(cds)[,s2_clust]))]
        names(s2_clust) <- unique(pData(cds)[, s2_clust])
    }

    S2 <- simple_plot_cells(cds, color_cells_by = s2_clust, aps = 0.9,
                            show_trajectory_graph = FALSE) +
        scale_colour_manual(name = s2_clust, values = s2_colors) +
        theme_light() + labs(title = s2_clust) +
        theme(axis.title = element_blank())

    if(is.null(s3_colors)) {
        s3_colors <- pals::alphabet()[1:length(unique(pData(cds)[,s3_phase]))]
        names(s3_colors) <- unique(pData(cds)[, s3_phase])
    }

    S3 <- simple_plot_cells(cds, color_cells_by = s3_phase, aps = 0.9,
                            show_trajectory_graph = FALSE) +
        scale_colour_manual(name = s3_phase, values = s3_colors) +
        theme_light() + labs(title = s3_phase) +
        theme(axis.title = element_blank())

    S4 <- simple_plot_cells(cds, genes = s4_gene1, aps = 0.9,
                            label_cell_groups  = FALSE,
                            show_trajectory_graph = FALSE) +
        labs(subtitle = s4_gene1) +
        theme_light() + theme(plot.title = element_blank(),
                              axis.title = element_blank()) + NoLegend() &
        scale_colour_gradientn(colours = RColorBrewer::brewer.pal(5, "YlGnBu"))
    
    S5 <- simple_plot_cells(cds, genes = s5_gene2, aps = 0.9,
                            label_cell_groups  = FALSE,
                            show_trajectory_graph = FALSE) +
        labs(subtitle = s5_gene2) +
        theme_light() + theme(plot.title = element_blank(),
                              axis.title = element_blank()) + NoLegend() &
        scale_colour_gradientn(colours = RColorBrewer::brewer.pal(5, "YlGnBu"))
    
    S6 <- simple_plot_cells(cds, color_cells_by = s6_score1, aps = 0.9,
                            label_cell_groups  = FALSE,
                            show_trajectory_graph = FALSE) +
        labs(subtitle = s6_score1) +
        theme_light() + theme(plot.title = element_blank(),
                              axis.title = element_blank()) &
        scale_colour_gradientn(colours = RColorBrewer::brewer.pal(5, "YlGnBu"))
    
    S7 <- simple_plot_cells(cds, color_cells_by = s7_pct.mt,
                            label_cell_groups  = FALSE,
                            show_trajectory_graph = FALSE) +
        theme_light() + labs(title = s7_pct.mt) +
        theme(axis.title = element_blank())
    
    phaseBarPlot <- pData(cds) %>% data.frame() %>%
        ggplot(aes(trt)) +
        geom_bar(aes(fill = Phase), position = "fill") +
        geom_text(size = 4, aes(label = ..count..), stat='count', angle = 0,
                  color = "black", position = position_fill(vjust = .05)) +
        scale_fill_manual(name = s3_phase, values = s3_colors) +
        ylim(0, 1)
    
    S8.1 <- phaseBarPlot + coord_flip() + theme_light() +
        theme(axis.title.y = element_blank()) + NoLegend() +
        labs(subtitle = "Phase")

    mitoColors <- c("<=24" = "#1E068F", ">24" = "#AC2793", ">60" = "#F79143")
    mitoBarPlot <- pData(cds) %>% data.frame() %>%
        ggplot(aes(trt)) +
        geom_bar(aes(fill = mt.cat), position = "fill") +
        geom_text(size = 4, aes(label = ..count..), stat='count', angle = 0,
                  color = "white", position = position_fill(vjust = .05)) +
        scale_fill_manual(name = "% MT", values = mitoColors) +
        labs(subtitle = "% MT") +
        ylim(0, 1)

    S8.2 <- mitoBarPlot + coord_flip() +
        theme_light() + theme(axis.title.y = element_blank()) +
        NoLegend() + labs(subtitle = "Mitochondria (%)")
    
    give.nm <- function(x) {
        return(c(y = max(x), label = length(x))) 
    }
    
    S9.1 <- plot_genes_violin(cds[s4_gene1,], group_cells_by = s2_clust,
                              ncol = 1) +
        stat_summary(fun.data = give.nm, geom = "text", col = "#D40000",
                     size = 1,
                     angle = 45, hjust = 1, vjust = 1) +
        scale_fill_manual(values = s2_colors) +
        labs(subtitle = s4_gene1) +
        theme_light() + NoLegend() + theme(axis.title.x = element_blank())
    
    S9.2 <- plot_genes_violin(cds[s5_gene2,], group_cells_by = s2_clust,
                              ncol = 1) +
        stat_summary(fun.data = give.nm, geom = "text", col = "#D40000",
                     size = 1,
                     angle = 45, hjust = 1, vjust = 1) +
        scale_fill_manual(values = s2_colors) +
        labs(subtitle = s5_gene2) +
        theme_light() + NoLegend() + 
        theme(axis.title.y = element_blank())
    
    ## plot it
    out <- ( S1 | S2 | S3 ) /
        ( S4 | S5 | S6 ) /
        ( S7 | (S8.1 / S8.2) | (S9.1 / S9.2) ) +
        plot_layout(guides = "collect")
    
    return(out)
    
}

#' ----plots_for_sam------------------------------------------------------------
sams_plot_s6 <- function(cds, seu, s1_trt = "trt", s1_colors = NULL,
                         s2_clust = "cluster", s2_colors = NULL,
                         s3_phase = "Phase", s3_colors = NULL,
                         s4_gene1 = "CDKN2A",
                         s5_gene2 = "ANGPTL4",
                         s6_score1 = "SASP_Score",
                         s7_pct.mt = "percent.mt") {
    
    if(is.null(s1_colors)) {
        s1_colors <- pals::alphabet()[1:length(unique(pData(cds)[,s1_trt]))]
        names(s1_colors) <- unique(pData(cds)[, s1_trt])
    }

    S1 <- simple_plot_cells(cds, color_cells_by = s1_trt,
                            show_trajectory_graph = FALSE) +
        scale_colour_manual(name = s1_trt, values = s1_colors) +
        theme_light() + labs(title = s1_trt) +
        theme(axis.title = element_blank())
    
    if(is.null(s2_colors)) {
        s2_colors <- pals::alphabet()[1:length(unique(pData(cds)[,s2_clust]))]
        names(s2_clust) <- unique(pData(cds)[, s2_clust])
    }

    S2 <- simple_plot_cells(cds, color_cells_by = s2_clust, aps = 0.9,
                            show_trajectory_graph = FALSE) +
        scale_colour_manual(name = s2_clust, values = s2_colors) +
        theme_light() + labs(title = s2_clust) +
        theme(axis.title = element_blank())

    if(is.null(s3_colors)) {
        s3_colors <- pals::alphabet()[1:length(unique(pData(cds)[,s3_phase]))]
        names(s3_colors) <- unique(pData(cds)[, s3_phase])
    }

    S3 <- simple_plot_cells(cds, color_cells_by = s3_phase, aps = 0.9,
                            show_trajectory_graph = FALSE) +
        scale_colour_manual(name = s3_phase, values = s3_colors) +
        theme_light() + labs(title = s3_phase) +
        theme(axis.title = element_blank())

    S4 <- simple_plot_cells(cds, genes = s4_gene1, aps = 0.9,
                            label_cell_groups  = FALSE,
                            show_trajectory_graph = FALSE) +
        labs(subtitle = s4_gene1) +
        theme_light() + theme(plot.title = element_blank(),
                              axis.title = element_blank()) + NoLegend() &
        scale_colour_gradientn(colours = RColorBrewer::brewer.pal(5, "YlGnBu"))
    
    S5 <- simple_plot_cells(cds, genes = s5_gene2, aps = 0.9,
                            label_cell_groups  = FALSE,
                            show_trajectory_graph = FALSE) +
        labs(subtitle = s5_gene2) +
        theme_light() + theme(plot.title = element_blank(),
                              axis.title = element_blank()) + NoLegend() &
        scale_colour_gradientn(colours = RColorBrewer::brewer.pal(5, "YlGnBu"))
    
    S6 <- simple_plot_cells(cds, color_cells_by = s6_score1, aps = 0.9,
                            label_cell_groups  = FALSE,
                            show_trajectory_graph = FALSE) +
        labs(subtitle = s6_score1) +
        theme_light() + theme(plot.title = element_blank(),
                              axis.title = element_blank()) &
        scale_colour_gradientn(colours = RColorBrewer::brewer.pal(5, "YlGnBu"))
    
    out <- (S1 | S2 | S3) /
        (S4 | S5 | S6) + plot_layout(guides = "collect")
    
    return(out)
}

#%     S4 <- plot_density(object = seu, features = s4_gene1) +
#%         labs(subtitle = s4_gene1) +
#%         theme_light() + theme(plot.title = element_blank(),
#%                               axis.title = element_blank()) + NoLegend() &
#%         scale_colour_gradientn(colours = RColorBrewer::brewer.pal(5, "YlGnBu"))
#% 
#%     S5 <- plot_density(object = seu, features = s5_gene2) +
#%         labs(subtitle = s5_gene2) +
#%         theme_light() + theme(plot.title = element_blank(),
#%                               axis.title = element_blank()) + NoLegend() &
#%         scale_colour_gradientn(colours = RColorBrewer::brewer.pal(5, "YlGnBu"))
#% 
#%     S6 <- plot_density(object = seu, features = s6_score1) +
#%         labs(subtitle = s6_score1) +
#%         theme_light() + theme(plot.title = element_blank(),
#%                               axis.title = element_blank()) &
#%         scale_colour_gradientn(colours = RColorBrewer::brewer.pal(5, "YlGnBu"))

