#' Plot proportional hazards
#'
#' @param gep Synthetic data (First output of generate_data function)
#'
#' @return a plot and table of Schoenfeld residuals
#' @export
#' @importFrom survival coxph Surv cox.zph
#' @examples
#' test_ph(gep)
test_ph <- function(gep){

  fit <- coxph(Surv(surv_time, status) ~ group,
               data = gep)

  # Plot Schoenfeld residuals

  sch_resid <- cox.zph(fit)
  print(sch_resid)

  plot(sch_resid,
       xlab = 'time',
       main = 'Schoenfeld residuals')

}
