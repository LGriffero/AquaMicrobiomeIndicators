#' Calculate Aquatic Quality Index (AQI)
#'
#' Calculates AQI following the CCME-WQI conceptual framework
#' using contaminant concentration thresholds.
#'
#' @param data Data frame containing grouping columns and contaminant concentrations.
#' @param thresholds Data frame with columns: contaminant and limit.
#' @param group_cols Character vector with grouping column names.
#'
#' @return Data frame with F1, F2, F3 and final AQI per group.
#' @export
calculate_AQI <- function(data, thresholds, group_cols) {

  # -------------------
  # Input validation
  # -------------------
  if (missing(group_cols)) {
    stop("You must provide grouping columns via 'group_cols'.")
  }

  if (!all(group_cols %in% names(data))) {
    stop("Some grouping columns are not present in the data.")
  }

  if (!all(c("contaminant", "limit") %in% names(thresholds))) {
    stop("Thresholds must contain columns: 'contaminant' and 'limit'.")
  }

  # -------------------
  # Pivot to long format
  # -------------------
  # Convert to long format (ONLY numeric contaminant columns)
  data_long <- data |>
    dplyr::select(
      dplyr::all_of(group_cols),
      dplyr::where(is.numeric)
    ) |>
    tidyr::pivot_longer(
      cols = -dplyr::all_of(group_cols),
      names_to = "contaminant",
      values_to = "concentration"
    )

  # -------------------
  # Join thresholds
  # -------------------
  data_eval <- dplyr::left_join(
    data_long,
    thresholds,
    by = "contaminant"
  )

  # -------------------
  # Determine exceedances
  # -------------------
  data_eval <- dplyr::mutate(
    data_eval,
    exceed = .data$concentration > .data$limit
  )

  # -------------------
  # F1 – Scope
  # -------------------
  fails_by_group <- data_eval |>
    dplyr::group_by(
      dplyr::across(dplyr::all_of(group_cols)),
      .data$contaminant
    ) |>
    dplyr::summarise(
      failed = any(.data$exceed, na.rm = TRUE),
      .groups = "drop"
    )

  F1_results <- fails_by_group |>
    dplyr::group_by(
      dplyr::across(dplyr::all_of(group_cols))
    ) |>
    dplyr::summarise(
      n_parameters = dplyr::n(),
      n_failed = sum(.data$failed),
      F1 = (.data$n_failed / .data$n_parameters) * 100,
      .groups = "drop"
    )

  # -------------------
  # F2 – Frequency
  # -------------------
  F2_results <- data_eval |>
    dplyr::group_by(
      dplyr::across(dplyr::all_of(group_cols))
    ) |>
    dplyr::summarise(
      n_tests = sum(!is.na(.data$concentration)),
      n_failed_tests = sum(.data$exceed, na.rm = TRUE),
      F2 = (.data$n_failed_tests / .data$n_tests) * 100,
      .groups = "drop"
    )

  # -------------------
  # F3 – Amplitude
  # -------------------
  data_excursion <- data_eval |>
    dplyr::mutate(
      excursion = dplyr::if_else(
        .data$exceed,
        (.data$concentration / .data$limit) - 1,
        0
      )
    )

  NSE_results <- data_excursion |>
    dplyr::group_by(
      dplyr::across(dplyr::all_of(group_cols))
    ) |>
    dplyr::summarise(
      NSE = sum(.data$excursion, na.rm = TRUE) /
        sum(!is.na(.data$concentration)),
      .groups = "drop"
    )

  F3_results <- NSE_results |>
    dplyr::mutate(
      F3 = .data$NSE / (0.01 * .data$NSE + 0.01)
    )

  # -------------------
  # Final AQI
  # -------------------
  AQI_components <- F1_results |>
    dplyr::left_join(F2_results, by = group_cols) |>
    dplyr::left_join(F3_results, by = group_cols) |>
    dplyr::mutate(
      AQI = 100 - (
        base::sqrt(.data$F1^2 +
                     .data$F2^2 +
                     .data$F3^2) / 1.732
      )
    )

  return(AQI_components)
}
