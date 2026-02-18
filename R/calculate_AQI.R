#' Calculate Aquatic Quality Index (AQI)
#'
#' Calculates AQI following the CCME-WQI conceptual framework
#' using contaminant concentration thresholds (e.g., PNEC values).
#'
#' @param data Data frame containing Site column and contaminant concentrations.
#' @param thresholds Data frame with columns: contaminant and limit.
#'
#' @return Data frame with F1, F2, F3 and final AQI per site.
#' @export
calculate_AQI <- function(data, thresholds) {

  # Convert to long format
  data_long <- tidyr::pivot_longer(
    data,
    cols = -Site,
    names_to = "contaminant",
    values_to = "concentration"
  )

  # Join thresholds
  data_eval <- dplyr::left_join(
    data_long,
    thresholds,
    by = "contaminant"
  )

  # Determine exceedances
  data_eval <- dplyr::mutate(
    data_eval,
    exceed = concentration > limit
  )

  # -------------------
  # F1 – Scope
  # -------------------
  fails_by_site <- data_eval %>%
    dplyr::group_by(Site, contaminant) %>%
    dplyr::summarise(
      failed = any(exceed, na.rm = TRUE),
      .groups = "drop"
    )

  F1_results <- fails_by_site %>%
    dplyr::group_by(Site) %>%
    dplyr::summarise(
      n_parameters = dplyr::n(),
      n_failed = sum(failed),
      F1 = (n_failed / n_parameters) * 100,
      .groups = "drop"
    )

  # -------------------
  # F2 – Frequency
  # -------------------
  F2_results <- data_eval %>%
    dplyr::group_by(Site) %>%
    dplyr::summarise(
      n_tests = sum(!is.na(concentration)),
      n_failed_tests = sum(exceed, na.rm = TRUE),
      F2 = (n_failed_tests / n_tests) * 100,
      .groups = "drop"
    )

  # -------------------
  # F3 – Amplitude
  # -------------------
  data_excursion <- dplyr::mutate(
    data_eval,
    excursion = dplyr::if_else(
      exceed,
      (concentration / limit) - 1,
      0
    )
  )

  NSE_results <- data_excursion %>%
    dplyr::group_by(Site) %>%
    dplyr::summarise(
      NSE = sum(excursion, na.rm = TRUE) / sum(!is.na(concentration)),
      .groups = "drop"
    )

  F3_results <- dplyr::mutate(
    NSE_results,
    F3 = NSE / (0.01 * NSE + 0.01)
  )

  # -------------------
  # Final AQI
  # -------------------
  AQI_components <- F1_results %>%
    dplyr::left_join(F2_results, by = "Site") %>%
    dplyr::left_join(F3_results, by = "Site") %>%
    dplyr::mutate(
      AQI = 100 - (sqrt(F1^2 + F2^2 + F3^2) / 1.732)
    )

  return(AQI_components)
}
