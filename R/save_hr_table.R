#' Save HR table
#'
#' @param hr_table HRs (output of get_hr_table function)
#' @param myfilename A filename (character)
#' @param dir A directory for saving (character)
#'
#' @return Saves the HR table
#' @export
#' @importFrom gt gtsave
#'
#' @examples
#' save_hr_table(my_hr_table, dir = './plot/hrs/', myfilename = "my_filename")
save_hr_table <- function(hr_table, dir, myfilename){
  gt::gtsave(as_gt(hr_table), paste(dir, myfilename, '.png', sep=''))
}
