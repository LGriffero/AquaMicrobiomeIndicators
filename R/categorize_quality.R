#' Categorize Water Quality Index
#'
#' Categorizes WQI or AQI values into qualitative classes.
#'
#' @param data Data frame containing index values.
#' @param index_column Name of column containing WQI/AQI values.
#' @param breaks Numeric vector of breakpoints (default = CCME scheme).
#' @param labels Character vector of category names.
#'
#' @return Data frame with added quality category column.
#' @export
categorize_quality <- function(data,
                               index_column = "AQI",
                               breaks = c(0, 45, 65, 80, 95, 100),
                               labels = c("Poor",
                                          "Marginal",
                                          "Fair",
                                          "Good",
                                          "Excellent")) {

  data$Quality_Category <- cut(
    data[[index_column]],
    breaks = breaks,
    labels = labels,
    include.lowest = TRUE,
    right = FALSE
  )

  return(data)
}
