# @url: https://stackoverflow.com/questions/15660829/how-to-add-a-number-of-observations-per-group-and-use-group-mean-in-ggplot2-boxp

#' function for number of observations
#' @examples
#' geom_boxplot(fill = "grey80", colour = "#3366FF") +
#'   stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
#'   stat_summary(fun.data = mean.n, geom = "text", fun.y = mean, colour = "red")
#' @export
give.n <- function(x) {
  return(c(y = -1, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

#' function for mean labels
#' @examples
#' ggplot(mtcars, aes(factor(cyl), mpg, label=rownames(mtcars))) +
#'   geom_boxplot(fill = "grey80", colour = "#3366FF") +
#'   stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
#'   stat_summary(fun.data = mean.n, geom = "text", fun.y = mean, colour = "red")
#' @export
give.mean <- function(x) {
  return(c(y = -1, label = round(mean(x),2))) 
  # experiment with the multiplier to find the perfect position
}

#' function for number of observations
#' @examples
#' geom_boxplot(fill = "grey80", colour = "#3366FF") +
#'   stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
#'   stat_summary(fun.data = mean.n, geom = "text", fun.y = mean, colour = "red")
#' @export
give.n0 <- function(x) {
  return(c(y = -0.02, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}
