convert_blocksizes_to_indices <- function(marker_blocksizes){
  block_indices <- as.list(rep(NA, length(marker_blocksizes)))
  block_indices[[1]] <- c(1, marker_blocksizes[1])

  if (length(marker_blocksizes) > 1){
    counter <- marker_blocksizes[1]
    for (i in 2:length(marker_blocksizes)){
      index_prev <- counter+1
      index <- index_prev + marker_blocksizes[i]-1
      block_indices[[i]] <- c(index_prev, index)
      counter <- block_indices[[i]][2]
    }
  }

  return(block_indices)
}




fit_cox_on_nonmarkers_tx <- function(gep, params){
  # take a random sample from the nonmarkers with size = max(marker_blocksizes)
  sample_size <- max(params$marker_blocksizes)
  first_nm <- sum(params$marker_blocksizes) + 1
  nonmarker_genes <- gep[,first_nm:params$n_genes]

  sample_index_1 <- ceiling(runif(1, min = first_nm, max = dim(nonmarker_genes)[2]))
  sample_index_2 <- sample_index_1 + sample_size - 1
  sample_indices <- c(sample_index_1, sample_index_2)

  # make formula
  colnames <- colnames(gep)[sample_indices[1]:sample_indices[2]]

  formula_gene <- paste(colnames, collapse = ' + ')
  formula_tx_int <- paste(paste(colnames, collapse = " * tx + "), '* tx')
  formula_all <- paste(formula_gene, formula_tx_int, 'tx', sep = ' + ')
  coxformula <- as.formula(paste('Surv(surv_time, status) ~', formula_all))

  rel_risk <- predict(coxph(coxformula, gep), type = 'risk')
  rel_risk_table <- data.frame(rel_risk = rel_risk,
                               group = gep$group)

  HRs <- get_hrs_from_cox_rel_risk(rel_risk_table)

  return(HRs)


}



get_hrs_from_cox_rel_risk <- function(rel_risk_table){
  risk_means <- aggregate(rel_risk_table$rel_risk,
                          list(rel_risk_table$group),
                          mean)
  colnames(risk_means) <- c('group', 'mean')


  mean_risk_B1 <- risk_means$mean[risk_means$group == 'B1']
  mean_risk_B0 <- risk_means$mean[risk_means$group == 'B0']
  mean_risk_NB1 <- risk_means$mean[risk_means$group == 'NB1']
  mean_risk_NB0 <- risk_means$mean[risk_means$group == 'NB0']
  HRs <- c('HR_B0_NB0' = exp(mean_risk_B0)/exp(mean_risk_NB0),
           'HR_B1_NB0' = exp(mean_risk_B1)/exp(mean_risk_NB0),
           'HR_NB1_NB0' = exp(mean_risk_NB1)/exp(mean_risk_NB0))
  HRs <- round(HRs, 2)
  return(HRs)
}



fit_cox_on_markers_tx <- function(gep, params, output_format = 'HR'){
  marker_blocksizes <- params$marker_blocksizes
  block_indices <- convert_blocksizes_to_indices(marker_blocksizes)
  gene_effects <- params$gene_effects

  # 1. Generate cox formula according to marker genes and gep.
  # 2. Generate a new variable in gep with binarized values of genes.

  block_formula_genes <- c()
  block_formula_tx_ints <- c()

  for (i in 1:length(marker_blocksizes)){
    if (gene_effects[i] == "no_pattern"){
      colnames <- colnames(gep)[block_indices[[i]][1]:block_indices[[i]][2]]

      block_formula_gene <- paste(colnames, collapse = ' + ')
      block_formula_tx_int <- paste(colnames, collapse = " * tx + ")


    } else {
      # 1 make formula
      block_formula_gene <- ifelse(gene_effects[i] == 'AND',
                                   paste('ANDblock_', as.character(i), sep = ''),
                                   paste('ORblock_', as.character(i), sep = ''))
      block_formula_tx_int <- block_formula_gene

      # 2. add new var to gep
      gene_effect_block_inst <- geneEffectBlockClass$new(gep = gep,
                                                         gene_indices = block_indices[[i]],
                                                         gene_effect = gene_effects[i],
                                                         params = params)
      #browser()
      gep$tempcolname <- gene_effect_block_inst$gene_effect_block
      names(gep)[names(gep) == "tempcolname"] <- block_formula_gene


    }
    block_formula_genes <- append(block_formula_genes, block_formula_gene)
    block_formula_tx_ints <- append(block_formula_tx_ints, block_formula_tx_int)

  }
  #browser()
  formula_genes <- paste(block_formula_genes, collapse = ' + ')
  formula_tx_int <- paste(paste(block_formula_tx_ints, collapse = ' + '), "* tx")

  formula_all <- paste(formula_genes, formula_tx_int, 'tx', sep = ' + ')

  coxformula <- as.formula(paste('Surv(surv_time, status) ~', formula_all))
  #print(paste('Coxformula:', paste('Surv(surv_time, status) ~', formula_all)))
  #browser()
  rel_risk <- predict(coxph(coxformula, gep), type = "risk")

  rel_risk_table <- data.frame(rel_risk = rel_risk,
                               group = gep$group)

  HRs <- get_hrs_from_cox_rel_risk(rel_risk_table)

  if (output_format == 'HR'){
    return(HRs)
  } else if (output_format == 'plot') {
    gep$group <- as.factor(gep$group)
    risk_plot <- plot(rel_risk ~ gep$group,
                      ylab = "Risk relative to sample mean",
                      xlab = "Group",
                      main = 'Risks relative to sample mean')

    return(risk_plot)
  }

}
