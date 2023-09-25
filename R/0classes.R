#' @import R6
#' @importFrom mvtnorm rmvnorm
#' @importFrom rockchalk lazyCov
#'
#'
checker <- R6::R6Class('checker', lock_objects = FALSE)


# checker$set("public", "initialize", function(){
#   self$check_noise(x)
#   self$check_surv_effect(x)
#   self$check_surv_effect_sign(x)
#   self$check_surv_distribution(x)
#   self$check_numeric(x)
#   self$check_numeric_pos(x)
#   self$check_numeric_01(x)
#   self$is.wholenumber(x)
#   self$check_integer_pos(x)
#   self$check_logical(x)
#
#   return(self)
#
# })



checker$set("public", "check_noise", function(x,
                                              objnm = deparse(substitute(x))){
  params <- list("type", "mean", "sd")
  right_datatype <- c(x[[1]] == "random", is.numeric(x[[2]]), is.numeric(x[[3]]))
  if(!is.list(x)){
    stop(do.call(sprintf(fmt = "'%s' must be a list containing '%s'",
                         append(params, objnm, after = 0))))
  } else if(!all(names(x) %in% params)){
    stop(do.call(sprintf, c(fmt = "'%s' must contain all of the following parameters:
                            '%s'", append(params, objnm, after = 0))))
  } else if(!all(right_datatype)){
    stop(sprintf("In the '%s' parameter, the first element of the list should be
    'random' and the rest should be numeric" , objnm))
  } else{
    return(x)
  }
})


checker$set("public", "check_gene_effect", function(x,
                                                    objnm = deparse(substitute(x))){

  levels <- c("no_pattern", "AND", "OR")

  if(!x %in% levels){
    stop(sprintf("'%s' must be a either from the following:
    'no_pattern', 'AND', 'OR'", objnm))
  }
  else{
    return(x)
  }
})


checker$set("public", "check_surv_distribution", function(x,
                                                          objnm = deparse(substitute(x))){

  levels <- c("Weibull")

  if(!x %in% levels){
    stop(sprintf("'%s' must be a either from the following:
    'Weibull' ", objnm))
  }
  else{
    return(x)
  }
})


checker$set("public", "check_numeric", function(x,
                                                objnm = deparse(substitute(x))){
  if(!is.numeric(x)){
    stop(sprintf("'%s' must be numeric",objnm))
  }
  else{
    return(x)
  }
})


checker$set("public", "check_numeric_pos", function(x,
                                                    objnm = deparse(substitute(x))){
  if(!is.numeric(x)){
    stop(sprintf("'%s' must be numeric",objnm))
  } else if(x <= 0){
    stop(sprintf("'%s' must be > 0",objnm))
  } else{
    return(x)
  }
})


checker$set("public", "check_numeric_01", function(x,
                                                   objnm = deparse(substitute(x))){
  if(!is.numeric(x)){
    stop(sprintf("'%s' must be numeric and between 0 and 1",objnm))
  } else if(x < 0 | x > 1){
    stop(sprintf("'%s' must be >= 0 and <= 1",objnm))
  } else{
    return(x)
  }
})


checker$set("public", "is.wholenumber", function(x,
                                                 tol = .Machine$double.eps^0.5){
  return(abs(x - round(x)) < tol)
})


checker$set("public", "check_integer_pos", function(x,
                                                    objnm = deparse(substitute(x))){
  if(!is.numeric(x)){
    stop(sprintf("'%s' must be an integer > 0",objnm))
  } else if(!self$is.wholenumber(x)){
    stop(sprintf("'%s' must be an integer > 0",objnm))
  } else if(x <= 0){
    stop(sprintf("'%s' must be > 0",objnm))
  } else{
    return(x)
  }
})


checker$set("public", "check_logical", function(x,
                                                objnm = deparse(substitute(x))){
  if(!is.logical(x)){
    stop(sprintf("'%s' must be logical",objnm))
  }
  else{
    return(x)
  }
})


# checker$set("public", "check_surv_d_parameters", function(x,
#                                                           objnm = deparse(substitute(x))){
#
#   params <- list("shape")
#
#   if(!is.list(x)){
#     stop(sprintf("'%s' must be a list", objnm))
#   } else if(!names(x) %in% params){
#     stop(do.call(sprintf, c(fmt = "'%s' must contain a parameter from: '%s'",
#                             append(params, objnm, after = 0))))
#   } else if(!is.numeric(x[[1]])){
#     stop(sprintf("'%s' must be a list with numeric element > 0", objnm))
#   } else if (x[[1]] <=0) {
#     stop(sprintf("'%s' must be a list with numeric element > 0", objnm))
#   } else{
#     return(x)
#   }
# })


block <- R6::R6Class("block", list(blocktype = NULL,
                                   size = NULL,
                                   mu_NB = NULL,
                                   mu_diff = NULL,
                                   mu_nonmarker = NULL,
                                   sd_B = NULL,
                                   sd_NB = NULL,
                                   sd_nonmarker = NULL,
                                   rho_B = NULL,
                                   rho_NB = NULL,
                                   rho_nonmarker = NULL,
                                   n_pts = NULL,
                                   fraction_pts_benefit = NULL,
                                   gene_effect = NULL),
                     lock_objects = FALSE,
                     private = list(check = checker$new())
)



block$set("public", "initialize", function(blocktype = NULL,
                                           size = NULL,
                                           mu_NB = NULL,
                                           mu_diff = NULL,
                                           mu_nonmarker = NULL,
                                           sd_B = NULL,
                                           sd_NB = NULL,
                                           sd_nonmarker = NULL,
                                           rho_B = NULL,
                                           rho_NB = NULL,
                                           rho_nonmarker = NULL,
                                           n_pts = NULL,
                                           fraction_pts_benefit = NULL,
                                           gene_effect = NULL){

  # Common variables between marker and nonmarker blocks

  self$blocktype <- blocktype

  if (missing(size)) {self$size <- 5}
  else {self$size <- private$check$check_integer_pos(size)}

  if (missing(n_pts)) {self$n_pts <- 1000}
  else {self$n_pts <- private$check$check_integer_pos(n_pts)}

  if (missing(fraction_pts_benefit)) {self$fraction_pts_benefit <- 0.5}
  else {self$fraction_pts_benefit <- private$check$check_numeric_01(fraction_pts_benefit)}


  if (self$blocktype == "marker"){
    # marker-block specific variables are initialized:
    if (missing(mu_NB)) {self$mu_NB <- 0}
    else {self$mu_NB <- private$check$check_numeric(mu_NB)}

    if (missing(mu_diff)) {self$mu_diff <- 1}
    else {self$mu_diff <- private$check$check_numeric(mu_diff)}

    if (missing(sd_B)) {self$sd_B <- 0.4}
    else {self$sd_B <- private$check$check_numeric_pos(sd_B)}

    if (missing(sd_NB)) {self$sd_NB <- 0.4}
    else {self$sd_NB <- private$check$check_numeric_pos(sd_NB)}

    if (missing(rho_B)) {self$rho_B <- 0.8}
    else {self$rho_B <- private$check$check_numeric_01(rho_B)}

    if (missing(rho_NB)) {self$rho_NB <- 0.8}
    else {self$rho_NB <- private$check$check_numeric_01(rho_NB)}

    if (missing(gene_effect)) {self$gene_effect <- "ADD"}
    else {self$gene_effect <- private$check$check_gene_effect(gene_effect)}


    # nonmarker-block specific variables are initialized as NULL:
    if (missing(mu_nonmarker)) {self$mu_nonmarker <- NULL}
    else {self$mu_nonmarker <- mu_nonmarker}

    if (missing(sd_nonmarker)) {self$sd_nonmarker <- NULL}
    else {self$sd_nonmarker <- sd_nonmarker}

    if (missing(rho_nonmarker)) {self$rho_nonmarker <- NULL}
    else {self$rho_nonmarker <- rho_nonmarker}
  }

  else{
    # nonmarker-block specific variables are initialized:
    if (missing(mu_nonmarker)) {self$mu_nonmarker <- 0}
    else {self$mu_nonmarker <- private$check$check_numeric(mu_nonmarker)}

    if (missing(sd_nonmarker)) {self$sd_nonmarker <- 0.4}
    else {self$sd_nonmarker <- private$check$check_numeric_pos(sd_nonmarker)}

    if (missing(rho_nonmarker)) {self$rho_nonmarker <- 0.8}
    else {self$rho_nonmarker <- private$check$check_numeric_01(rho_nonmarker)}


    # marker-block specific are initialized as NULL:
    if (missing(mu_NB)) {self$mu_NB <- NULL}
    else {self$mu_NB <- mu_NB}

    if (missing(mu_diff)) {self$mu_diff <- NULL}
    else {self$mu_diff <- mu_diff}

    if (missing(sd_B)) {self$sd_B <- NULL}
    else {self$sd_B <- sd_B}

    if (missing(sd_NB)) {self$sd_NB <- NULL}
    else {self$sd_NB <- sd_NB}

    if (missing(rho_B)) {self$rho_B <- NULL}
    else {self$rho_B <- rho_B}

    if (missing(rho_NB)) {self$rho_NB <- NULL}
    else {self$rho_NB <- rho_NB}

    self$gene_effect <- NULL
  }

  return(self)
})


block$set("active", "mu_B", function(value){
  return(self$mu_NB + self$mu_diff)
})


block$set("public", "generate_blockvector", function(sd){
  return(rep(sd, self$size))
})


block$set("active", "mu_B_blockvector", function(value){
  return(self$generate_blockvector(self$mu_B))
})


block$set("active", "mu_NB_blockvector", function(value){
  return(self$generate_blockvector(self$mu_NB))
})


block$set("active", "mu_nonmarker_blockvector", function(value){
  return(self$generate_blockvector(self$mu_nonmarker))
})



block$set("active", "sd_B_blockvector", function(value){
  return(self$generate_blockvector(self$sd_B))
})


block$set("active", "sd_NB_blockvector", function(value){
  return(self$generate_blockvector(self$sd_NB))
})


block$set("active", "sd_nonmarker_blockvector", function(value){
  return(self$generate_blockvector(self$sd_nonmarker))
})



block$set("public", "generate_block_cormat", function(rho){
  # generates cormat for the block from the rho and size.
  return(rockchalk::lazyCov(rho, Sd = 1, d = self$size))
})


block$set("active", "block_cormat_B", function(value){
  return(self$generate_block_cormat(self$rho_B))
})


block$set("active", "block_cormat_NB", function(value){
  return(self$generate_block_cormat(self$rho_NB))
})


block$set("active", "block_cormat_nonmarker", function(value){
  return(self$generate_block_cormat(self$rho_nonmarker))
})


block$set("active", "n_pt_B", function(value){
  # sets size of Benefit patient group (patient group with Benefit variant of marker genes)
  return(round(self$n_pts * self$fraction_pts_benefit, 0))
})


block$set("active", "n_pt_NB", function(value){
  # sets size of NB patient group.
  return(self$n_pts - self$n_pt_B)
})


block$set("active", "gene_effect_layer", function(value){
  #' This function creates an overlay matrix that is 0 by default.
  #' Depending on gene effects, for each patient belonging to one of the B or NB group,
  #' min 0, max blocksize - 1 number of genes will be assigned the negative mu-diff value.
  #' The overlay matrix is added later to the gene expression. Thus the
  #' given genes will have a value as if they were sampled from the opposite distribution.

  # 1. Base overlay matrix
  mat <- matrix(0, nrow = self$n_pts, ncol = self$size)
  # If nonmarker block
  if (is.null(self$gene_effect)){
    return(mat)
  }
  else if (self$gene_effect == 'AND'){
    # for each NB patient, draw a random number that is the number of genes
    # to assign the mu_diff to.
    rowrand <- round(runif(self$n_pt_NB, 0, self$size-1))
    for (i in 1:self$n_pt_NB){
      # for each patient in the group, draw as many rand.int. as the rowrand,
      # to get the index of the gene to be affected.
      colrand <- round(runif(n = rowrand[i], min = 0, max = self$size))
      # assign negative mu_diff to the given genes.
      mat[self$n_pt_B+i, colrand] <- self$mu_diff
    }
    return(mat)
  }
  else if (self$gene_effect == 'OR'){
    # for each B patient, draw a random number that is the number of genes to
    # assign the negative mu_diff to.
    rowrand <- round(runif(self$n_pt_B, 0, self$size-1))
    for (i in 1:self$n_pt_B){
      colrand <- round(runif(rowrand[i], 0, self$size))
      mat[i, colrand] <- self$mu_diff * -1
    }
    return(mat)
  }
  else{
    # if no_pattern the overlay has only 0s
    return(mat)
  }
})


gepClass <- R6::R6Class("gep",
                        public = list(
                          marker_block_list = NULL,
                          nonmarker_block_list = NULL,
                          n_genes = NA,
                          n_pts = NA,
                          fraction_pts_benefit = NA,
                          fraction_tx_1 = NA,
                          fraction_censored = NA,
                          noise = NULL),
                        lock_objects = FALSE,
                        private = list(check = checker$new()))


gepClass$set("public", "initialize", function(marker_block_list = NULL,
                                              nonmarker_block_list = NULL,
                                              n_genes = NA,
                                              n_pts = NA,
                                              fraction_pts_benefit = NA,
                                              fraction_tx_1 = NA,
                                              fraction_censored = NA,
                                              noise = NULL){

  self$marker_block_list <- marker_block_list
  self$nonmarker_block_list <- nonmarker_block_list

  if (missing(n_genes)) {self$n_genes <- 100}
  else {self$n_genes <- private$check$check_integer_pos(n_genes)}

  if (missing(n_pts)) {self$n_pts <- 1000}
  else {self$n_pts <- private$check$check_integer_pos(n_pts)}

  if (missing(fraction_pts_benefit)) {self$fraction_pts_benefit <- 0.5}
  else {self$fraction_pts_benefit <- private$check$check_numeric_01(fraction_pts_benefit)}

  if (missing(fraction_tx_1)) {self$fraction_tx_1 <- 0.5}
  else {self$fraction_tx_1 <- private$check$check_numeric_01(fraction_tx_1)}

  if (missing(fraction_censored)) {self$fraction_censored <- 0.5}
  else {self$fraction_censored <- private$check$check_numeric_01(fraction_censored)}

  if (missing(noise)) {self$noise <- "no"}
  else {self$noise <- private$check$check_noise(noise)}

  return(self)

})


gepClass$set("active", "n_pt_B", function(value){
  # sets size of Benefit patient group (patient group with Benefit variant of marker genes)
  return(round(self$n_pts * self$fraction_pts_benefit, 0))
})


gepClass$set("active", "n_pt_NB", function(value){
  # sets size of NB patient group.
  return(self$n_pts - self$n_pt_B)
})


gepClass$set("active", "mu_B_vector", function(value){
  return(self$marker_block_list[[1]])
})


gepClass$set("active", "mu_NB_vector", function(value){
  return(self$marker_block_list[[2]])
})


gepClass$set("active", "covmat_B", function(value){
  return(self$marker_block_list[[3]])
})


gepClass$set("active", "covmat_NB", function(value){
  return(self$marker_block_list[[4]])
})


gepClass$set("active", "overlay_mat_marker", function(value){
  return(self$marker_block_list[[5]])
})


gepClass$set("active", "mu_nonmarker_vector", function(value){
  return(self$nonmarker_block_list[[1]])
})


gepClass$set("active", "covmat_nonmarker", function(value){
  return(self$nonmarker_block_list[[2]])
})


gepClass$set("active", "overlay_mat_nonmarker", function(value){
  return(self$nonmarker_block_list[[3]])
})


gepClass$set("active", "gep_B", function(value){
  return(mvtnorm::rmvnorm(self$n_pt_B, self$mu_B_vector, self$covmat_B))
})


gepClass$set("active", "gep_NB", function(value){
  return(mvtnorm::rmvnorm(self$n_pt_NB, self$mu_NB_vector, self$covmat_B))
})


gepClass$set("active", "gep_nonmarker", function(value){
  return(mvtnorm::rmvnorm(self$n_pts, self$mu_nonmarker_vector, self$covmat_nonmarker))
})


gepClass$set("active", "gep_marker", function(value){
  return(rbind(self$gep_B, self$gep_NB))
})


gepClass$set("active", "overlay_mat", function(value){
  return(cbind(self$overlay_mat_marker, self$overlay_mat_nonmarker))
})


gepClass$set("active", "gep_m_nm", function(value){
  # Assembles gene expression.

  # 1. binds marker (m), nonmarker (nm) columns together.
  gep_m_nm <- cbind(self$gep_marker, self$gep_nonmarker)
  # 2. Add overlay to gep
  gep_m_nm <- gep_m_nm + self$overlay_mat
  # 3. Add names to dimensions.
  dimnames(gep_m_nm) <- list(paste("pt",seq_len(nrow(gep_m_nm)),sep="_"), paste("gene",seq_len(ncol(gep_m_nm)),sep="_"))
  # 4. Add noise
  if (self$noise$type == "random"){
    return(gep_m_nm + matrix(rnorm(self$n_pts*self$n_genes, mean = self$noise$mean,
                                   sd = self$noise$sd),
                             nrow = self$n_pts))
  }
})



gepClass$set("active", "gep", function(value){
  # Assembles gep with benefit, tx, status (event)
  benefit <- c(rep(1, self$n_pt_B), rep(0, self$n_pt_NB))
  tx <- sample(c(1,0),prob=c(self$fraction_tx_1, 1-self$fraction_tx_1), self$n_pts, replace=TRUE)
  group <- ifelse(benefit=='1', paste("B", tx, sep = ""), paste("NB", tx, sep=''))
  status <- sample(c(0,1),prob=c(self$fraction_censored, 1-self$fraction_censored), self$n_pts, replace=TRUE)
  gep <- cbind(as.data.frame(self$gep_m_nm), benefit, tx, group, status)

  return(gep)
})




geneEffectBlockClass <- R6::R6Class("geneEffectBlockClass",
                                    list(gep = NULL,
                                         gene_indices = NULL,
                                         gene_effect = NULL,
                                         params = NULL),
                                    lock_objects = FALSE,
                                    private = list(check = checker$new())
)


geneEffectBlockClass$set("public", "initialize",
                         function(gep = NULL,
                                  gene_indices = NULL,
                                  gene_effect = NULL,
                                  params = NULL){


                           self$gep <- gep
                           self$gene_indices <- gene_indices
                           self$params <- params
                           self$gene_effect <- private$check$check_gene_effect(gene_effect)

                           return(self)
                         }
)


geneEffectBlockClass$set("active", "genes", function(value){
  return(abs(self$gep[, self$gene_indices[[1]]:self$gene_indices[[2]]]))
})


geneEffectBlockClass$set("active", "median_per_gene", function(value){
  return(apply(self$genes,2,median,na.rm=TRUE))
})


geneEffectBlockClass$set("active", "mu_NB", function(value){
  return(self$params$mus_NB)
})

geneEffectBlockClass$set("active", "mu_B", function(value){
  return(self$mu_NB + self$params$mu_diffs)
})


geneEffectBlockClass$set("public", "opt_func",function(par){
  ##Count number of genes are exceeding the threshold per patient
  y_hat<-colSums(t(self$genes) > par)

  if (self$gene_effect == "OR") {y_hat <- y_hat > 0}
  if (self$gene_effect == "AND") {y_hat <- y_hat >= dim(self$genes)[2]-1}
  ##Return the number of correct assignments (negative to minimize)
  return(-1*sum(self$gep$benefit == y_hat))
})




geneEffectBlockClass$set("active", "th", function(value){
  if (self$gene_effect == 'OR'){
    #par = mean(c(max(self$genes), min(self$genes))) #self$median_per_gene
    par = mean(c(self$mu_B, self$mu_NB))
    th_per_gene <- optim(par = par, fn = self$opt_func)$par
  } else if (self$gene_effect == 'AND'){
    #par = mean(c(max(self$genes), min(self$genes)))
    par = mean(c(self$mu_B, self$mu_NB))
    th_per_gene <- optim(par = par, fn = self$opt_func)$par
  }

  return(matrix(rep(th_per_gene, each = dim(self$genes)[1]), nrow = dim(self$genes)[1]))
})


geneEffectBlockClass$set("active", "is_higher_than_th", function(value){
  return(abs(self$genes) > self$th)
})


geneEffectBlockClass$set("public", "generate_and_block", function(value){
  # are all genes / patient higher than th? returns a vector with 1 bool value / patient
  return(as.integer(apply(self$is_higher_than_th, 1, all)))
})


geneEffectBlockClass$set("public", "generate_or_block", function(value){
  return(as.integer(apply(self$is_higher_than_th, 1, any)))
})


geneEffectBlockClass$set("public", "generate_add_block", function(value){
  # are all genes / patient higher than th? returns a vector with 1 bool value / patient
  return(as.integer(self$is_higher_than_th))
})



geneEffectBlockClass$set("active", "gene_effect_block", function(value){
  if (self$gene_effect == "AND"){
    return(self$generate_and_block())
  } else if (self$gene_effect == "OR"){
    return(self$generate_or_block())
  } else {
    return(self$generate_add_block())
  }
})



