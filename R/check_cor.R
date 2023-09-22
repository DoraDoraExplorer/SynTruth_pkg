#' Check correlation
#'
#' @param gep Gene expression data
#' @param mytitle A title (character)
#'
#' @return Correlation plot
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#' @examples
#' marker_genes <- best_gep[,1:sum(params$marker_blocksizes)]
#' nonmarker_genes <- best_gep[,sum(params$marker_blocksizes):params$n_genes]
#' check_cor(nonmarker_genes, mytitle = "My correlation plot")

check_cor <- function(gep, mytitle){
  #browser()
  corr_mat <- round(cor(gep),2)
  melted_corr_mat <- reshape2::melt(corr_mat)
  ggplot(data = melted_corr_mat, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    ggtitle(mytitle) +
    ylab("Genes") +
    xlab("Genes") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                         limit = c(-1,1), space = "Lab",name="Pearson\nCorrelation")
}
