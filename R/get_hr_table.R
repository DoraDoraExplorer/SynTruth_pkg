#' Get HR table
#'
#' @param gep Synthetic data (First output of generate_data_function)
#' @param output_format ["per_benefit", "per_group"]
#' @param mytitle A title (character)
#'
#' @return a HR table
#' @export
#'
#' @importFrom survival coxph Surv
#' @importFrom magrittr %>%
#' @importFrom gtsummary tbl_regression modify_caption tbl_stack
#'
#' @examples
#' get_hr_table(gep, output_format = 'per_benefit', mytitle = 'My title')
get_hr_table <- function(gep, output_format = "per_benefit", mytitle){
  gep$group <- as.factor(gep$group)

  if (output_format == 'per_benefit'){
    gep_B <- gep[gep$benefit == 1,]
    gep_NB <- gep[gep$benefit == 0,]

    gep_B$group <- droplevels(gep_B$group)
    gep_NB$group <- droplevels(gep_NB$group)
    gep_B$group <- relevel(gep_B$group, ref = "B0")
    gep_NB$group <- relevel(gep_NB$group, ref = "NB0")

    coxph_B <- coxph(Surv(surv_time, status) ~ group,
                     data = gep_B) %>%
      gtsummary::tbl_regression(exp = TRUE) %>% modify_caption(mytitle)

    coxph_NB <- coxph(Surv(surv_time, status) ~ group,
                      data = gep_NB) %>%
      gtsummary::tbl_regression(exp = TRUE) %>% modify_caption(mytitle)

    merged <- tbl_stack(list(coxph_B, coxph_NB),
                        group_header = list('Benefit group', 'Non-Benefit group'))

    return(merged)

  } else if (output_format == 'per_group') {
    gep$group <- relevel(gep$group, ref = "NB0")

    coxph <- coxph(Surv(surv_time, status) ~ group,
                   data = gep) %>%
      gtsummary::tbl_regression(exp = TRUE) %>%
      # remove_row_type('NB0', type = "header") %>%
      modify_caption(mytitle)

    return(coxph)

  }
}
