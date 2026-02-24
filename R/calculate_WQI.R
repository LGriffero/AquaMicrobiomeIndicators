#' Calculate Water Quality Index (CCME-WQI)
#'
#' Calculates the Water Quality Index (WQI) following the official
#' Canadian Council of Ministers of the Environment (CCME) framework.
#'
#' Limits must be provided as a named list with structure:
#'
#'   list(
#'     Variable = list(type = "max", limit = value),
#'     Variable = list(type = "min", limit = value),
#'     Variable = list(type = "range", min = value, max = value)
#'   )
#'
#' @param data Data frame containing water quality variables.
#' @param limits Named list defining threshold type and values.
#' @param group_cols Character vector of grouping columns.
#'
#' @return Tibble with F1, F2, F3 and WQI.
#' @export
calculate_WQI <- function(data,
                          limits,
                          group_cols) {

  # -------------------------
  # 0. Validation
  # -------------------------

  if (!all(group_cols %in% colnames(data))) {
    stop("Some grouping columns not found in data.")
  }

  if (!all(names(limits) %in% colnames(data))) {
    stop("Some variables in 'limits' not found in data.")
  }

  # -------------------------
  # 1. Internal CCME excursion function
  # -------------------------

  excursion_fun <- function(value, limit_info) {

    type <- limit_info$type

    if (type == "max") {

      limit <- limit_info$limit

      fail <- value > limit
      excursion <- ifelse(fail, (value / limit) - 1, 0)

    } else if (type == "min") {

      limit <- limit_info$limit

      fail <- value < limit
      excursion <- ifelse(fail, (limit / value) - 1, 0)

    } else if (type == "range") {

      min_limit <- limit_info$min
      max_limit <- limit_info$max

      fail <- value < min_limit | value > max_limit

      excursion <- ifelse(
        value < min_limit,
        (min_limit / value) - 1,
        ifelse(
          value > max_limit,
          (value / max_limit) - 1,
          0
        )
      )

    } else {

      stop("limit type must be 'max', 'min', or 'range'")
    }

    list(fail = fail, excursion = excursion)
  }

  # -------------------------
  # 2. Apply to each variable
  # -------------------------

  data_processed <- data

  fail_matrix <- list()
  excursion_matrix <- list()

  for (var in names(limits)) {

    result <- excursion_fun(data[[var]], limits[[var]])

    fail_matrix[[var]] <- result$fail
    excursion_matrix[[var]] <- result$excursion
  }

  fail_df <- as.data.frame(fail_matrix)
  colnames(fail_df) <- paste0("fail_", names(fail_df))

  excursion_df <- as.data.frame(excursion_matrix)
  colnames(excursion_df) <- paste0("exc_", names(excursion_df))

  # -------------------------
  # 3. Bind and group
  # -------------------------

  data_combined <- cbind(
    data[group_cols],
    fail_df,
    excursion_df
  )

  results <- dplyr::summarise(
    dplyr::group_by(
      data_combined,
      dplyr::across(dplyr::all_of(group_cols))
    ),

    total_vars = length(limits),

    failed_vars = sum(
      colSums(
        dplyr::across(dplyr::all_of(names(fail_df))),
        na.rm = TRUE
      ) > 0
    ),

    total_tests = dplyr::n() * length(limits),

    failed_tests = sum(
      as.matrix(
        dplyr::across(dplyr::all_of(names(fail_df)))
      ),
      na.rm = TRUE
    ),

    sum_excursions = sum(
      as.matrix(
        dplyr::across(dplyr::all_of(names(excursion_df)))
      ),
      na.rm = TRUE
    ),

    NSE = sum_excursions / total_tests,

    .groups = "drop"
  )

  # -------------------------
  # 4. Final CCME calculation
  # -------------------------

  results <- dplyr::mutate(
    results,

    F1 = (failed_vars / total_vars) * 100,
    F2 = (failed_tests / total_tests) * 100,
    F3 = NSE / (0.01 * NSE + 0.01),

    WQI = 100 - (
      base::sqrt(F1^2 + F2^2 + F3^2) / 1.732
    )
  )
  results <- dplyr::select(
    results,
    dplyr::all_of(group_cols),
    F1, F2, F3, WQI
  )

  return(results)
}
