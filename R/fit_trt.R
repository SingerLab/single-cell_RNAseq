#' fit model treatment vs control in individual cell lines
#'
#' @param cds a cds object
#' 
#' @export
fit_trt <- function(cds) {
    n.cores <- as.numeric(system("echo $LSB_MAX_NUM_PROCESSORS", intern = TRUE))
    
    trt_ctrl_fit <- cds %>%
        fit_models(model_formula_str = "~trt",
                   cores = n.cores, verbose = TRUE) %>%
        coefficient_table()
    ## remove model and model summary, they add a significant amount of weight
    ## to the table, if needed, re-create, it's faster than loading and saving
    trt_ctrl_fit  <-  trt_ctrl_fit %>%
        select(!(model:model_summary))
    
    trt_ctrl_fit
}
