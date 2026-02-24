#' Run indicator species analysis with adaptive pruning
#'
#' Performs indicator species analysis using the IndVal framework
#' with adaptive pruning of specificity (At) and sqrtIV thresholds.
#'
#' @param asv_matrix Matrix or data.frame of taxa abundances
#'   (rows = sites, columns = taxa).
#' @param groups Vector of group assignments (same length as rows of asv_matrix).
#' @param B_threshold Numeric vector of fidelity (B) thresholds per group.
#'   If shorter than the number of groups, the first value is reused.
#' @param At_start Initial specificity (At) threshold.
#' @param At_min Minimum specificity threshold allowed during adaptive reduction.
#' @param At_step Step size for reducing At threshold.
#' @param sqrtIVt_start Initial sqrtIV threshold.
#' @param sqrtIVt_min Minimum sqrtIV threshold allowed during adaptive reduction.
#' @param sqrtIVt_step Step size for reducing sqrtIV threshold.
#' @param max.order Maximum order of species combinations.
#' @param presence_absence Logical; if TRUE, converts abundances to presence/absence.
#'
#' @return A list with:
#' \describe{
#'   \item{indicators}{Data frame of selected indicator taxa per group.}
#'   \item{coverage}{Coverage values per group.}
#' }
#'
#' @export
run_indval_analysis <- function(asv_matrix,
                                groups,
                                B_threshold = 0.5,
                                At_start = 1,
                                At_min = 0.6,
                                At_step = 0.05,
                                sqrtIVt_start = 0.8,
                                sqrtIVt_min = 0.4,
                                sqrtIVt_step = 0.05,
                                max.order = 3,
                                presence_absence = TRUE) {

  groups <- as.factor(groups)
  asv_matrix <- as.matrix(asv_matrix)

  if (presence_absence) {
    asv_matrix <- ifelse(asv_matrix > 0, 1, 0)
  }

  if (nrow(asv_matrix) != length(groups)) {
    stop("Mismatch between number of samples and group vector")
  }

  B <- indicspecies::strassoc(asv_matrix,
                              groups,
                              func = "B")

  groups_unique <- sort(unique(groups))

  indicator_list <- list()
  coverage_list <- list()

  for (g in groups_unique) {

    g_char <- as.character(g)

    if (!g_char %in% colnames(B)) next

    # 🔹 aplicar B_threshold específico del grupo
    g_index <- which(groups_unique == g)
    B_thr_g <- if (length(B_threshold) >= g_index) {
      B_threshold[g_index]
    } else {
      B_threshold[1]
    }

    sel_taxa <- which(B[, g_char] > B_thr_g)
    if (length(sel_taxa) == 0) next

    sc <- indicspecies::indicators(
      X = asv_matrix[, sel_taxa, drop = FALSE],
      cluster = groups,
      group = g,
      max.order = max.order
    )

    # --- Adaptive pruning doble ---
    At_current <- At_start
    sc_pruned <- NULL

    while (At_current >= At_min) {

      sqrtIVt_current <- sqrtIVt_start

      while (sqrtIVt_current >= sqrtIVt_min) {

        sc_try <- indicspecies::pruneindicators(
          sc,
          At = At_current,
          sqrtIVt = sqrtIVt_current
        )

        if (length(sc_try$sign) > 0) {
          sc_pruned <- sc_try
          break
        }

        sqrtIVt_current <- sqrtIVt_current - sqrtIVt_step
      }

      if (!is.null(sc_pruned)) break

      At_current <- At_current - At_step
    }

    if (is.null(sc_pruned)) next

    taxa_selected <- names(sc_pruned$sign)

    indicator_list[[g_char]] <- data.frame(
      taxon = taxa_selected,
      group = g,
      At_used = At_current,
      sqrtIVt_used = sqrtIVt_current,
      B_threshold = B_thr_g,
      stringsAsFactors = FALSE
    )

    coverage_list[[g_char]] <-
      indicspecies::coverage(sc_pruned)
  }

  if (length(indicator_list) == 0) {
    stop("No indicators detected even after adaptive pruning.")
  }

  indicator_table <- do.call(rbind, indicator_list)
  rownames(indicator_table) <- NULL

  return(list(
    indicators = indicator_table,
    coverage = coverage_list
  ))
}
