generate_gep_surv <- function(gepdata,
                              n_pts,
                              surv_distribution,
                              shape,
                              scale_NB0,
                              HR_B0_NB0,
                              HR_NB1_NB0,
                              HR_B1_NB0){


  # on.exit({
  #   if (!is.null(oldseed)) { .GlobalEnv$.Random.seed <- oldseed }
  #   else { rm(".Random.seed", envir = .GlobalEnv) }
  # }, add = TRUE)
  #
  # oldseed <- NULL
  # if (exists(".Random.seed", .GlobalEnv)) { oldseed <- .GlobalEnv$.Random.seed }
  #
  # set.seed(seed)


  #Determine Weibull scale parameters per group

  if (surv_distribution == "Weibull"){

    #shape <- surv_d_params$shape

    scale_B0 <-exp(-log(HR_B0_NB0 )/shape + log(scale_NB0))
    scale_NB1<-exp(-log(HR_NB1_NB0)/shape + log(scale_NB0))
    scale_B1<-exp(-log(HR_B1_NB0)/shape + log(scale_NB0))

    gepdata$scale <- ifelse(gepdata$group == 'B1', scale_B1,
                            ifelse(gepdata$group == 'B0', scale_B0,
                                   ifelse(gepdata$group == 'NB0', scale_NB0, scale_NB1)))

    gepdata$surv_time <- rweibull(n = nrow(gepdata),shape=shape, scale=as.numeric(gepdata$scale))
  }
  return(gepdata)
}


convert_cormat_to_covmat <- function(sd_vector, cormat){
  # https://stats.stackexchange.com/questions/62850/obtaining-covariance-matrix-from-correlation-matrix
  return(outer(sd_vector, sd_vector) * cormat)
}



convert_blocklist_to_bdiag <- function(block_list){
  #from individual blocks, bdiag generates a block diagonal matrix
  return(as.matrix(Matrix::bdiag(block_list)))
}



generate_marker_blocks <- function(marker_blocksizes,
                                   mus_NB,
                                   mu_diffs,
                                   sds_B, sds_NB,
                                   rhos_B, rhos_NB,
                                   n_pts,
                                   fraction_pts_benefit,
                                   gene_effects){

  marker_block_list = list()

  mu_B_vector = c()
  mu_NB_vector = c()
  sd_B_vector = c()
  sd_NB_vector = c()
  cormat_B_list <- as.list(rep(NA, length(marker_blocksizes)))
  cormat_NB_list <- as.list(rep(NA, length(marker_blocksizes)))
  gene_effect_layer_list = as.list(rep(NA, length(marker_blocksizes)))

  #browser()
  for (i in 1:length(marker_blocksizes)){
    marker_block <- block$new(blocktype = "marker",
                              size = marker_blocksizes[i],
                              mu_NB = mus_NB[i],
                              mu_diff = mu_diffs[i],
                              sd_B = sds_B[i],
                              sd_NB = sds_NB[i],
                              rho_B = rhos_B[i],
                              rho_NB = rhos_NB[i],
                              n_pts = n_pts,
                              fraction_pts_benefit = fraction_pts_benefit,
                              gene_effect = gene_effects[i],
                              mu_nonmarker = NULL,
                              sd_nonmarker = NULL,
                              rho_nonmarker = NULL
    )

    # 1. get blockvectors and block cormats
    mu_B_blockvector <- marker_block$mu_B_blockvector
    mu_NB_blockvector <- marker_block$mu_NB_blockvector
    sd_B_blockvector <- marker_block$sd_B_blockvector
    sd_NB_blockvector <- marker_block$sd_NB_blockvector
    block_cormat_B <- marker_block$block_cormat_B
    block_cormat_NB <- marker_block$block_cormat_NB

    # 2. from blockvectors, create vectors and from block cormats, create lists.
    mu_B_vector <- append(mu_B_vector, mu_B_blockvector)
    mu_NB_vector <- append(mu_NB_vector, mu_NB_blockvector)
    sd_B_vector <- append(sd_B_vector, sd_B_blockvector)
    sd_NB_vector <- append(sd_NB_vector, sd_NB_blockvector)
    cormat_B_list[[i]] <- block_cormat_B
    cormat_NB_list[[i]] <- block_cormat_NB

    # 3. get the gene effect (relationship) layers
    gene_effect_layer_list[[i]] <- marker_block$gene_effect_layer
  }

  # 4. make block diagonal cormat (B and NB)
  diag_cormat_B <- convert_blocklist_to_bdiag(cormat_B_list)
  diag_cormat_NB <- convert_blocklist_to_bdiag(cormat_NB_list)


  # 5. convert block diagonal cormat (B and NB) to covmat
  covmat_B <- convert_cormat_to_covmat(sd_B_vector, diag_cormat_B)
  covmat_NB <- convert_cormat_to_covmat(sd_NB_vector, diag_cormat_NB)

  # 6. create an overlay matrix from the layer blocks
  overlay_mat_marker <- do.call(cbind, gene_effect_layer_list)


  return(list(mu_B_vector, mu_NB_vector,
              covmat_B, covmat_NB,
              overlay_mat_marker))
}


set_nonmarker_blocksizes <- function(n_genes, marker_blocksizes, sds_nonmarker){
  n_nonmarker_genes <- n_genes-sum(marker_blocksizes)
  n_nonmarker_blocks <- length(sds_nonmarker)
  if (n_nonmarker_genes %% n_nonmarker_blocks == 0){
    return(rep(n_nonmarker_genes/n_nonmarker_blocks, n_nonmarker_blocks))
  }
  else {
    remainder <- n_nonmarker_genes %% n_nonmarker_blocks
    nonmarker_blocksizes <- c(rep(n_nonmarker_genes %/% n_nonmarker_blocks,
                                  n_nonmarker_blocks-1),
                              n_nonmarker_genes %/% n_nonmarker_blocks + remainder)
    return(nonmarker_blocksizes)
  }

}



generate_nonmarker_blocks <- function(n_genes,
                                      marker_blocksizes,
                                      mu_nonmarker,
                                      sds_nonmarker,
                                      rhos_nonmarker,
                                      n_pts,
                                      fraction_pts_benefit){

  # set nonmarker blocksizes
  nonmarker_blocksizes <- set_nonmarker_blocksizes(n_genes,
                                                   marker_blocksizes,
                                                   sds_nonmarker)

  mu_vector = c()
  sd_vector = c()
  cormat_list <- as.list(rep(NA, length(nonmarker_blocksizes)))
  gene_effect_layer_list = as.list(rep(NA, length(marker_blocksizes)))


  for (i in 1:length(nonmarker_blocksizes)){
    nonmarker_block = block$new(
      blocktype = "nonmarker",
      size = nonmarker_blocksizes[i],
      mu_nonmarker = mu_nonmarker,
      sd_nonmarker = sds_nonmarker[i],
      rho_nonmarker = rhos_nonmarker[i],
      n_pts = n_pts,
      fraction_pts_benefit = fraction_pts_benefit,
      gene_effect = NULL,
      mu_NB = NULL,
      mu_diff = NULL,
      sd_B = NULL,
      sd_NB = NULL,
      rho_B = NULL,
      rho_NB = NULL
    )


    mu_blockvector <- nonmarker_block$mu_nonmarker_blockvector
    sd_blockvector <- nonmarker_block$sd_nonmarker_blockvector
    block_cormat <- nonmarker_block$block_cormat_nonmarker

    # 2. from blockvectors, create vectors and from block cormats, create a list.
    mu_vector <- append(mu_vector, mu_blockvector)
    sd_vector <- append(sd_vector, sd_blockvector)
    cormat_list[[i]] <- block_cormat

    # 3. gene effect layer list
    gene_effect_layer_list[[i]] <- nonmarker_block$gene_effect_layer
  }

  # 4. make block diagonal cormat
  diag_cormat <- convert_blocklist_to_bdiag(cormat_list)

  # 5. convert block diagonal cormat to covmat
  covmat <- convert_cormat_to_covmat(sd_vector, diag_cormat)

  # 6. overlay matrix
  overlay_mat_nonmarker <- do.call(cbind, gene_effect_layer_list)

  return(list(mu_vector, covmat, overlay_mat_nonmarker))

}




generate_gepdata <- function(marker_block_list,
                             nonmarker_block_list,
                             n_pts,
                             n_genes,
                             fraction_pts_benefit,
                             fraction_tx_1,
                             fraction_censored,
                             noisemean,
                             noisesd){

  #browser()
  gep <- gepClass$new(marker_block_list = marker_block_list,
                      nonmarker_block_list = nonmarker_block_list,
                      n_genes = n_genes,
                      n_pts = n_pts,
                      fraction_pts_benefit = fraction_pts_benefit,
                      fraction_tx_1 = fraction_tx_1,
                      fraction_censored = fraction_censored,
                      noisemean = noisemean,
                      noisesd = noisesd)

  return(gep$gep)

}

get_hrs <- function(gep){
  gep$group <- as.factor(gep$group)
  gep$group <- relevel(gep$group, ref = "NB0")

  coxph <- coxph(Surv(surv_time, status) ~ group,
                 data = gep)

  hrs <- round(exp(coxph$coefficients), 2)

  return(hrs)

}
