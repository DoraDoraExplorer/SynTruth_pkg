#' Get best gep
#'
#' Runs generate data n_repeats times. Creates a plot of HRs. A dataset with HRs
#' under the given threshold is extracted as best_gep and its estimated HRs is the best_hrs.
#'
#' @param params About params: see generate_data function.
#' @param n_repeats Number of repeats
#' @param threshold The sum of the difference between the input HRs and the estimated HRs.
#' The first dataset with estimated HRs under this difference is extracted as best_gep.
#' If the threshold is 0.31, and the estimated HRs are  0.1, 0.1, 0.1, then their sum
#' is 0.3, which is under the threshold and this dataset is be extracted.
#'
#' @return list of:
#' 1) plot of estimated HRs vs input HRs.
#' 2) the best gep (synthetic data)
#' 3) the best HRs
#' @export
#' @import ggplot2
#'
#' @examples
#' get_best_gep(params, n_repeats = 100, threshold = 0.3)
#'
get_best_gep <- function(params, n_repeats, threshold){

  input_hrs <- data.frame(HR_B0_NB0 = params$HR_B0_NB0,
                          HR_B1_NB0 = params$HR_B1_NB0,
                          HR_NB1_NB0 = params$HR_NB1_NB0)

  pred_repeats_df <- data.frame(HR_B0_NB0 = rep(NA, n_repeats),
                                HR_B1_NB0 = rep(NA, n_repeats),
                                HR_NB1_NB0 = rep(NA, n_repeats))

  best_gep <- NULL
  best_output_hrs <- NULL
  for (i in 1:n_repeats){
    gep_hrs <- generate_data(params)
    output_hrs <- gep_hrs$hrs
    if (is.null(best_gep) & all(abs(input_hrs-output_hrs) < threshold)){
      best_gep <- gep_hrs$gep
      best_output_hrs <- output_hrs
    }

    pred_repeats_df$HR_B0_NB0[i] <- output_hrs['groupB0']
    pred_repeats_df$HR_B1_NB0[i] <- output_hrs['groupB1']
    pred_repeats_df$HR_NB1_NB0[i] <- output_hrs['groupNB1']
  }

  pred_repeats_df <- gather(pred_repeats_df, key = 'group')

  p <- ggplot(pred_repeats_df, aes(group, value)) +
    geom_boxplot(aes(colour = group)) +
    ylab("HR") +
    geom_point(gather(input_hrs, key = 'group'), mapping = aes(group, value), size = 3) +
    theme(text = element_text(size = 12),
          legend.text = element_text(size = 12),
          axis.text.y = element_text(size = 13),
          axis.title.x = element_blank()) +
    ggtitle('Input and estimated HRs from Cox-model fit on patient groups')

  return(list(p = p, best_gep = best_gep, best_output_hs = best_output_hrs))
}
