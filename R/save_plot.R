#' Save plot
#'
#' @param myfilename Filename (character)
#' @param myplot The plot object.
#' @param width Width of plot
#' @param height Height of plot.
#' @param dir directory  (character)
#'
#' @return saves a plot to
#' @export
#' @importFrom ggplot2 ggsave
#' @examples
#' myfilename <- 'My plot'
#' save_plot(myfilename, myplot = survplot, dir = './plots/', width = 12, height = 10)
save_plot <- function(myfilename, myplot, dir, width = 12, height = 8){
  ggsave(paste(dir, myfilename, ".png", sep = ''), myplot,
         width = width, height = height, units = "cm")
}
