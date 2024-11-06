#' plot_cells wrapper with predetermined label size, point size and alpha
#'
#' @param cds a cell_data_set object
#'
#' @param cps cell point size
#'
#' @param aps cell alpha value
#'
#' @param gls gene/cluster label size
#'
#' @export
simple_plot_cells <- function(cds, cps = 0.8, aps = 0.48, group_label_size = 4, ...) {
    plot_cells(cds,
               cell_size = cps,
               alpha = aps,
               ...)
}

