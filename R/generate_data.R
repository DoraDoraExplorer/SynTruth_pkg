#' Generate data
#' @description Generates a list of:
#'
#' 1) gep: synthetic data
#'
#' 2) hrs: HRs estimated from data. This is done by fitting a Cox-model on the group variable.
#'
#' @param params a list of all parameters to generate data.
#'
#' @return a list of gep and HRs
#' @export
#'
#' @examples
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
#' noisemean = 0,
#' noisesd = 1,
#' surv_distribution = "Weibull",
#' shape = 1.5,
#' scale_NB0 = 10,
#' HR_B0_NB0 = 0.8,
#' HR_NB1_NB0 = 0.7,
#' HR_B1_NB0 = 0.5)
#'
#' gep_hrs <- generate_data(params)
#'
#' gep <- gep_hrs$gep
#' hrs <- gep_hrs$hrs

generate_data <- function(params){


  marker_block_list <- generate_marker_blocks(
    marker_blocksizes = params$marker_blocksizes,
    mus_NB = params$mus_NB, mu_diffs = params$mu_diffs,
    sds_B = params$sds_B, sds_NB = params$sds_NB,
    rhos_B = params$rhos_B, rhos_NB =  params$rhos_NB,
    n_pts = params$n_pts,
    fraction_pts_benefit = params$fraction_pts_benefit,
    gene_effects = params$gene_effects)


  nonmarker_block_list <- generate_nonmarker_blocks(
    n_genes = params$n_genes,
    marker_blocksizes = params$marker_blocksizes,
    mu_nonmarker = params$mu_nonmarker,
    sds_nonmarker = params$sds_nonmarker,
    rhos_nonmarker = params$rhos_nonmarker,
    n_pts = params$n_pts,
    fraction_pts_benefit = params$fraction_pts_benefit)

  gepdata <- generate_gepdata(
    marker_block_list = marker_block_list,
    nonmarker_block_list = nonmarker_block_list,
    n_pts = params$n_pts,
    n_genes = params$n_genes,
    fraction_pts_benefit = params$fraction_pts_benefit,
    fraction_tx_1 = params$fraction_tx_1,
    fraction_censored = params$fraction_censored,
    noisemean = params$noisemean,
    noisesd = params$noisesd
    )


  gep_surv <- generate_gep_surv(gepdata = gepdata,
                                n_pts = params$n_pts,
                                surv_distribution = 'Weibull',
                                shape = params$shape,
                                scale_NB0 = params$scale_NB0,
                                HR_B0_NB0 = params$HR_B0_NB0,
                                HR_NB1_NB0 = params$HR_NB1_NB0,
                                HR_B1_NB0 = params$HR_B1_NB0)

  hrs <- get_hrs(gep_surv)

  return(list(gep = gep_surv, hrs = hrs))

}
