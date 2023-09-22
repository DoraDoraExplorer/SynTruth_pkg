#' Get difference between input and output HRs
#' Output HRs are those estimated from data by fitting a Cox-model using the group as covariate.
#'
#' @param input_hrs User-defined HRs.
#' @param output_hrs estimated HRs.
#' @param output_format ["table", "sum", "diff"]
#'
#' @return Difference between input and output HRs in 3 different formats: table, the sum of differences or a vector of differences.
#' @export
#'
#' @examples
#' diff_io_hrs <- get_diff_io_hrs(input_hrs = c(params$HR_B0_NB0, params$HR_B1_NB0,
#' params$HR_NB1_NB0),
#' output_hrs = best_output_hrs,
#' output_format = 'table')

get_diff_io_hrs <- function(input_hrs, output_hrs, output_format){

  diff_df <- data.frame(input_hrs, output_hrs)
  diff_df <- round(diff_df, 2)
  diff_df['diff'] <- round(diff_df$input_hrs-diff_df$output_hrs, 2)
  diff_df <- rbind(diff_df, data.frame(input_hrs = '-',
                                       output_hrs = '-',
                                       diff = sum(abs(diff_df$diff)),
                                       row.names = 'sum_abs'))

  if(output_format == 'table'){
    return(diff_df)
  } else if (output_format == 'sum'){
    return(diff_df$diff[rownames(diff_df) == 'sum_abs'])
  } else if (output_format == 'diff'){
    diffvect <- diff_df$diff[1:3]
    names(diffvect) <- rownames(diff_df)[1:3]
    return(diffvect)
  }
}
