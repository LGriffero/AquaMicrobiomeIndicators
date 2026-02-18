#' Run Indicator Species Analysis
#'
#' Performs indicator species analysis using IndVal, rphi or r.g
#' from the indicspecies package.
#'
#' @param asv_matrix Matrix or data frame (sites × taxa).
#' @param groups Vector of group assignments per site.
#' @param method Association method ("IndVal.g", "r", "r.g").
#' @param nperm Number of permutations.
#' @param presence_absence Logical; convert abundances to PA.
#'
#' @return multipatt object.
#' @export
run_indval_analysis <- function(asv_matrix,
                                groups,
                                method = "IndVal.g",
                                nperm = 999,
                                presence_absence = TRUE) {

  # Convert to matrix
  asv_matrix <- as.matrix(asv_matrix)

  # Convert to presence-absence if requested
  if (presence_absence) {
    asv_matrix <- ifelse(asv_matrix > 0, 1, 0)
  }

  # Run multipatt
  result <- indicspecies::multipatt(
    asv_matrix,
    cluster = groups,
    func = method,
    control = permute::how(nperm = nperm)
  )

  return(result)
}
