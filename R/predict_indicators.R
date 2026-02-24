
#' Predict ecological quality group using indicator taxa
#'
#' @param abundance_matrix Matrix/data.frame of taxa abundances (samples x taxa)
#' @param indicator_table Output from run_indval_analysis()
#' @param true_group Optional vector with true ecological group
#'
#' @return A list with predictions and performance metrics
#' @export
predict_indicators <- function(abundance_matrix,
                               indicator_table,
                               true_group = NULL) {

  abundance_matrix <- as.data.frame(abundance_matrix)

  # Taxa comunes
  common_taxa <- intersect(indicator_table$taxon,
                           colnames(abundance_matrix))

  if (length(common_taxa) == 0) {
    stop("No indicator taxa found in abundance_matrix")
  }

  indicator_table <- indicator_table[
    indicator_table$taxon %in% common_taxa, ]

  # Score por grupo
  scores <- sapply(unique(indicator_table$group), function(g) {

    taxa_g <- indicator_table$taxon[
      indicator_table$group == g]

    rowSums(abundance_matrix[, taxa_g,
                             drop = FALSE],
            na.rm = TRUE)
  })

  if (is.vector(scores)) {
    scores <- matrix(scores, ncol = 1)
    colnames(scores) <- unique(indicator_table$group)
  }

  predicted_group <- colnames(scores)[
    max.col(scores, ties.method = "first")]

  result <- list(
    predicted_group = predicted_group,
    scores = scores
  )

  # ----------------------------
  # Performance metrics
  # ----------------------------
  if (!is.null(true_group)) {

    if (length(true_group) != nrow(abundance_matrix)) {
      stop("Length of true_group must match number of samples")
    }

    confusion <- table(
      Real = true_group,
      Predicted = predicted_group
    )

    accuracy <- sum(diag(confusion)) / sum(confusion)

    sensitivity <- sapply(rownames(confusion), function(g) {
      TP <- confusion[g, g]
      FN <- sum(confusion[g, ]) - TP
      if ((TP + FN) == 0) return(NA)
      TP / (TP + FN)
    })

    total <- sum(confusion)
    row_marginals <- rowSums(confusion)
    col_marginals <- colSums(confusion)

    Pe <- sum(row_marginals * col_marginals) / (total^2)
    Po <- accuracy

    kappa <- (Po - Pe) / (1 - Pe)

    result$true_group <- true_group
    result$confusion_matrix <- confusion
    result$accuracy <- accuracy
    result$sensitivity <- sensitivity
    result$kappa <- kappa
  }

  return(result)
}
