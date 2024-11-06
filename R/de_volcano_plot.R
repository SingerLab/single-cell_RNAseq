#' monocle3 differential expression volcano plot
#' 
#' @export
de_volcano_plot <- function(tbl, x.col = "estimate", y.col = "q_value",
                            effect.highlight.thrs = 2,
                            sig.value.highlight.thrs = 10^-8,
                            genes.of.interst = c("CDK4", "HMGA2", "MDM2")) {

    tbl.use <- tbl[, c("gene_short_name", x.col, y.col)]
    names(tbl.use) <- c("gene_short_name", "x", "y")
    
    volc <- tbl.use %>% ggplot(aes(x, -log10(y),
                         label = gene_short_name)) +
        geom_point(alpha = 0.4,
                   colour = ifelse(tbl.use$x >= effect.highlight.thrs &
                                   tbl.use$y <= sig.value.highlight.thrs,
                                   "#D40000",
                            ifelse(
                                tbl.use$x <= effect.highlight.thrs*-1 &
                                tbl.use$y <= sig.value.highlight.thrs,
                                "#0065BF",
                                "gray40"))) +
        ggrepel::geom_label_repel(
                     data = rbind(
                         tbl.use %>%
                         filter(gene_short_name %in% genes.of.interst),
                         tbl.use %>% top_n(-1, y),
                         tbl.use %>%
                         filter(x >= effect.highlight.thrs &
                                y <= sig.value.highlight.thrs |
                                x <= effect.highlight.thrs*-1 &
                                y <= sig.value.highlight.thrs)
                         ),
                     box.padding = .6, point.padding = .6, 
                 direction = "both", max.overlaps = 30) +
        theme_light()

    return(volc)

}
