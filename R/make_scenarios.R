#' Make scenarios
#'
#' Makes a list where each element is a different param list.
#'
#' @param scenarios list of the parameters that are different in each list.
#' @param params see generate_data function.
#'
#' @return a list of:
#' 1) scenario list (list where each element is a different param list)
#' 2) a vector of names of the scenarios.
#' Each name consists of those parameters that are different per scenario.
#'
#' @export
#'
#' @examples
#' scenarios <- list(mu_diffs = c(1, 5),
#' gene_effects = c('AND', 'OR', "no_pattern"))
#'
#' params <- list(
#' marker_blocksizes = c(10),
#' mus_NB = c(0),
#' mu_diffs = c(2),
#' sds_B = c(0.4),
#' sds_NB = c(0.4),
#' rhos_B = c(0.8),
#' rhos_NB = c(0.8),
#' gene_effects = c("AND"),
#' mu_nonmarker = 0,
#' sds_nonmarker = c(0.4),
#' rhos_nonmarker = c(0.8),
#' n_pts = 1000,
#' n_genes = 100,
#' fraction_pts_benefit = 0.5,
#' fraction_tx_1 = 0.5,
#' fraction_censored = 0,
#' noise = list(type = "random", mean = 0, sd = 0),
#' surv_distribution = "Weibull",
#' surv_d_params = list(shape=1.5),
#' scale_NB0 = 10,
#' HR_B0_NB0 = 0.8,
#' HR_NB1_NB0 = 0.7,
#' HR_B1_NB0 = 0.5)
#'
#' scenario_params_names <- make_scenarios(scenarios, params)
#' scenario_params_list <- scenario_params_names$params_list
#' scenario_names <- scenario_params_names$sc_names
#' scenario_names
#'

make_scenarios <- function(scenarios, params){
  scs <- expand.grid(scenarios)
  n_scs <- dim(scs)[1]

  params_list <- as.list(rep(NA, n_scs))
  for (i in 1:n_scs){
    non_changing_params <- params[names(params) %in% names(scs) == FALSE]
    scenario <- append(non_changing_params, scs[i,])

    params_list[[i]] <- scenario
  }
  sc_names <- apply(scs, 1, function(x) paste(x, collapse = '-'))

  return(list(params_list = params_list, sc_names = sc_names))

}
