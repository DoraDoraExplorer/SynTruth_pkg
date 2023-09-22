#' Validate SynTruth in scenarios using the Cox-model
#'
#' For each scenario, fits the Cox-model on marker genes and treatment n_repeats times.
#' For nonmarker scenario, fits the Cox-model on 10 randomly chosen nonmarker genes n_repeats times.
#'
#'
#' @param scenario_names Names of scenarios (output of the make_scenarios function).
#' @param scenarios List where each element is a parameter list (Output of the make_scenarios function)
#' @param n_repeats Number of repeats
#'
#' @return A boxplot of the estimated HRs for each group and the input HRs.
#' @export
#'
#' @import ggplot2
#' @import survival
#' @importFrom tidyr gather
#'
#' @examples
#' cox_pred_hrs <- validate_cox_in_scenarios(scenario_names = scenario_names,
#' scenarios = scenario_params_list,
#' n_repeats = 100)
#'
#'
validate_cox_in_scenarios <- function(scenario_names, scenarios,
                                      n_repeats){

  scenario_df <- data.frame(key = c(), value = c(), scenario = c())

  for (i in 1:length(scenarios)){
    params <- scenarios[[i]]
    input_hrs <- data.frame(HR_B0_NB0 = params$HR_B0_NB0,
                            HR_B1_NB0 = params$HR_B1_NB0,
                            HR_NB1_NB0 = params$HR_NB1_NB0)

    pred_repeats_df <- data.frame(HR_B0_NB0 = rep(NA, n_repeats),
                                  HR_B1_NB0 = rep(NA, n_repeats),
                                  HR_NB1_NB0 = rep(NA, n_repeats))
    #browser()
    for (j in 1:n_repeats){
      gep <- generate_data(params)$gep
      pred_HRs <- fit_cox_on_markers_tx(gep = gep,
                                        params = params,
                                        output_format = 'HR')
      pred_repeats_df$HR_B0_NB0[j] <- pred_HRs['HR_B0_NB0']
      pred_repeats_df$HR_B1_NB0[j] <- pred_HRs['HR_B1_NB0']
      pred_repeats_df$HR_NB1_NB0[j] <- pred_HRs['HR_NB1_NB0']
    }

    pred_repeats_df <- gather(pred_repeats_df, key = 'group')
    pred_repeats_df$scenario <- scenario_names[i]

    scenario_df <- rbind(scenario_df, pred_repeats_df)

  }

  # fit on nonmarker genes.

  nm_pred_repeats_df <- data.frame(HR_B0_NB0 = rep(NA, n_repeats),
                                   HR_B1_NB0 = rep(NA, n_repeats),
                                   HR_NB1_NB0 = rep(NA, n_repeats))
  params <- scenarios[[1]]
  for (i in 1:n_repeats){
    gep <- generate_data(params)$gep
    pred_HRs <- fit_cox_on_nonmarkers_tx(gep = gep, params = params)

    nm_pred_repeats_df$HR_B0_NB0[i] <- pred_HRs['HR_B0_NB0']
    nm_pred_repeats_df$HR_B1_NB0[i] <- pred_HRs['HR_B1_NB0']
    nm_pred_repeats_df$HR_NB1_NB0[i] <- pred_HRs['HR_NB1_NB0']

  }

  nm_pred_repeats_df <- gather(nm_pred_repeats_df, key = 'group')
  nm_pred_repeats_df$scenario <- 'non-marker'
  scenario_df <- rbind(scenario_df, nm_pred_repeats_df)
  scenario_df <- scenario_df[order(scenario_df$scenario),]
  p <- ggplot(scenario_df, aes(group, value)) +
    geom_boxplot(aes(colour = group)) +
    facet_wrap(~ scenario) +
    ylab("HR predictions") +
    geom_point(gather(input_hrs, key = 'group'), mapping = aes(group, value), size = 3) +
    theme(legend.position = c(0.6, 0.02),
          legend.justification = c("right", "bottom"),
          strip.text = element_text(size = 12),
          legend.text = element_text(size = 12))
  #scale_color_manual(name = "", values = c("input HRs" = "black"))

  #print(paste('Number of repeats:', n_repeats))

  return(p)
}
