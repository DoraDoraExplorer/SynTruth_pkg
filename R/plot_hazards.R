#' Plot hazard curves
#'
#' @param gep Synthetic data (first output of generate_data function)
#' @param shape Shape of Weibull distribution
#' @param mytitle A title (character)
#'
#' @return Hazard curve plot for 4 patient groups
#' @export
#' @import ggplot2
#'
#' @examples
#' plot_hazards(gep, shape = params$surv_d_params$shape, mytitle = 'My hazard curves')
#'
plot_hazards <- function(gep, shape, mytitle){
  surv <- gep$surv_time
  lambda <- gep$scale

  h_vec <- lambda * shape * surv^(shape-1)

  df <- data.frame(surv = surv, hazard = h_vec, group = gep$group)

  ggplot(df, aes(x=surv, y=hazard)) +
    geom_line(aes(colour = group)) +
    xlab('time') +
    ylab('hazard') +
    ggtitle(mytitle)
}
