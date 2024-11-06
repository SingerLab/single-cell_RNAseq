#' function to rpocess all cell lines in the same way
#'
#' @param cds an unprocessed cds
#'
#' @param ndim number of pca dimensions
#'
#' @param niter number of iterations
#'
#' @import dplyr
#' 
#' @export
reprocess_cell_line <- function(pre.cds, ndim, niter, ...) {
    cds = pre.cds %>%
        preprocess_cds(num_dim = ndim, ...) %>%
        reduce_dimension(., ...) %>%
        cluster_cells(., num_iter = niter, ...) %>%
        learn_graph(., ...)
    cds
}
