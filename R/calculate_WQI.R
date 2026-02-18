#' Calculate Water Quality Index (WQI)
#'
#' @param data Data frame containing physicochemical parameters
#' @param limits Named list of threshold functions
#' @return Data frame with WQI values
#' @export
calculate_WQI <- function(data, limits) {

  library(dplyr)

  # Aplicar límites
  data_flagged <- data %>%
    mutate(across(all_of(names(limits)),
                  ~ limits[[cur_column()]](.),
                  .names = "fail_{.col}"))

  # F1
  F1 <- data_flagged %>%
    summarise(
      failed_vars = sum(colSums(across(starts_with("fail_"))) > 0),
      total_vars = length(limits)
    ) %>%
    mutate(F1 = (failed_vars / total_vars) * 100)

  return(F1)
}
