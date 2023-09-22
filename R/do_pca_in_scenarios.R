#' Do PCA in scenarios
#'
#' @param scenario_names Names of scenarios (Output of the make_scenarios function)
#' @param scenarios List where each element is a parameter list (Output of the make_scenarios function)
#' @param nrow Number of rows (for plotting)
#'
#' @return PCA scores plot with scenario subplots
#' @export
#'
#' @import plotly
#' @import ggfortify
#' @import factoExtra
#' @importFrom gridExtra arrangeGrob grid.arrange
#'
#' @examples
#' do_pca_in_scenarios(scenario_names,
#' scenarios = scenario_params_list,
#' nrow = 2)
#'
do_pca_in_scenarios <- function(scenario_names, scenarios, nrow){

  mytheme <- theme(plot.title = element_text(hjust = 0.5, size = 14),
                   legend.title = element_text(face="bold", size = 11),
                   axis.text.x = element_text(size = 9),
                   axis.text.y = element_text(size = 9),
                   legend.position = "none")

  mytheme_legend <- theme(plot.title = element_text(hjust = 0.5, size = 14),
                          legend.title = element_text(face="bold", size = 11),
                          axis.text.x = element_text(size = 9),
                          axis.text.y = element_text(size = 9),
                          legend.key.size = unit(0.8, 'cm'),
                          legend.text = element_text(size = 11),
                          legend.position = c(1.3, 1.1),
                          legend.justification = c("right", "top"))


  plot_list <- as.list(rep(NA, length(scenarios)))

  for (i in 1:length(scenarios)){
    params <- scenarios[[i]]
    gep <- generate_data(params)$gep
    gep_markers <- gep[,1:sum(params$marker_blocksizes)]
    nm_start_index <- sum(params$marker_blocksizes) + 1
    gep_nm <- gep[,nm_start_index:params$n_genes]

    pca_res <- prcomp(gep_markers, center = TRUE)

    gep$benefit <- as.factor(gep$benefit)

    p <- autoplot(pca_res, data = gep, colour = 'benefit') +
      ggtitle(scenario_names[[i]]) +
      mytheme
    plot_list[[i]] <- p
  }

  # Add nonmarker scenario

  gep_nm <- gep[,nm_start_index:params$n_genes]
  pca_nm_res <- prcomp(gep_nm, center = TRUE)
  nm_p <- autoplot(pca_nm_res, data = gep, colour = 'benefit') +
    ggtitle('Non-marker') +
    mytheme_legend

  plot_list <- append(plot_list, list(nm_p))

  # Make grid
  mygrid <- grid.arrange(grobs = plot_list, nrow = nrow)

  return(mygrid)
}
