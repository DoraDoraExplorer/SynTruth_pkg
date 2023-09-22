#' Do k-means clustering in scenarios
#'
#' @param scenario_names Names of scenarios (Output of the make_scenarios function)
#' @param scenarios List where each element is a parameter list (Output of the make_scenarios function)
#' @param n_repeats Number of repeats
#' @param nrows Number of rows (for plotting)
#'
#' @return A plot with subplots of the silhouette scores resulting from k-means clustering in each scenario.
#' @export
#'
#' @import factoextra
#' @importFrom cluster silhouette
#' @import ggplot2
#'
#' @examples
#' do_km_in_scenarios(scenario_names = scenario_names,
#' scenarios = scenario_params_list,
#' n_repeats = 100,
#' nrows = 3)
#'
do_km_in_scenarios <- function(scenario_names,
                               scenarios,
                               n_repeats,
                               nrows){

  scenario_df <- data.frame(stts = c(), scenario = c())

  for (i in 1:length(scenarios)){
    stts_vect <- c()
    params <- scenarios[[i]]
    for (j in 1:n_repeats){
      gep <- generate_data(params)$gep
      gep_genes <- gep[,1:params$n_genes]
      # Predictions
      pred <- kmeans(gep_genes, centers = 2, nstart = 25)$cluster
      stt <- silhouette(pred, dist(gep_genes))
      stts <- mean(stt[,3])
      stts_vect <- append(stts_vect, stts)
    }

    stts_df <- data.frame(stts = stts_vect, scenario = scenario_names[i])
    scenario_df <- rbind(scenario_df, stts_df)
  }

  # Add nonmarker scenario

  nm_stts_vect <- c()

  # nm_accuracies <- c()
  params <- scenarios[[1]]
  for (i in 1:n_repeats){
    gep <- generate_data(params)$gep
    start_nm_index <- sum(params$marker_blocksizes) + 1
    gep_nm_genes <- gep[,start_nm_index:params$n_genes]
    # Predictions
    pred <- kmeans(gep_nm_genes, centers = 2, nstart = 25)$cluster
    stt <- silhouette(pred, dist(gep_nm_genes))
    stts <- mean(stt[,3])
    nm_stts_vect <- append(nm_stts_vect, stts)
  }

  # append nm to scenarios
  nm_stts_df <- data.frame(stts = nm_stts_vect, scenario = 'non-marker')
  scenario_df <- rbind(scenario_df, nm_stts_df)
  scenario_df$stts <- round(scenario_df$stts, 3)

  scaleFUN <- function(x) sprintf("%.3f", x)

  p <- ggplot(scenario_df, aes(y = stts)) +
    geom_boxplot() +
    facet_wrap(~ scenario, nrow = nrows) + #scales = 'free'
    ylab("Silhouette score") +
    scale_y_continuous(labels = scaleFUN) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          #strip.text = element_text(size = 12),
          text=element_text(size=15))

  return(p)
}
