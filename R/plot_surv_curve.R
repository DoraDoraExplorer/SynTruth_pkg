#' Plot survival (Kaplan-Meier) curve
#'
#' @param gep Synthetic data (First output of generate_data function)
#' @param output_hrs Estimated HRs (Second output of generate_data function)
#' @param mytitle A title
#'
#' @return A plot of survival curves for 4 patient groups
#' @export
#'
#' @import ggplot2
#' @importFrom survival survfit Surv
#' @examples
#' mytitle = "Best dataset"
#' plot_km_curve(gep, output_hrs, mytitle)
plot_surv_curve <- function(gep, output_hrs, mytitle){

  gep$group <- as.factor(gep$group)
  tx_fit <- survfit(Surv(surv_time, status) ~ group, data = gep)

  mytheme <- theme(plot.title = element_text(hjust = 0.5, size = 14),
                   legend.title = element_text(face="bold", size = 13),
                   axis.text.x = element_text(size = 13),
                   axis.text.y = element_text(size = 13),
                   legend.key.size = unit(1, 'cm'),
                   legend.text = element_text(size = 12))

  subtitle_hrs <- paste(c('HRs B0:', 'B1:', 'NB1:'), as.character(output_hrs))
  subtitle_hrs <- paste(subtitle_hrs, collapse = ', ')

  p1 <- autoplot(tx_fit) +
    labs(x = "Survival time (arbitrary measure)",
         y = "Cumulative survival (%)",
         title = mytitle,
         subtitle = subtitle_hrs) +
    mytheme


  p1
  return(p1)
}
