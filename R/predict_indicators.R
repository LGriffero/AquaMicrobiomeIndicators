if (!is.null(true_group)) {

  if (length(true_group) != nrow(abundance_matrix)) {
    stop("Length of true_group must match number of samples")
  }

  confusion <- table(
    Real = true_group,
    Predicted = predicted_group
  )

  # Accuracy
  accuracy <- sum(diag(confusion)) / sum(confusion)

  # Sensitivity por grupo
  sensitivity <- sapply(rownames(confusion), function(g) {
    TP <- confusion[g, g]
    FN <- sum(confusion[g, ]) - TP
    if ((TP + FN) == 0) return(NA)
    TP / (TP + FN)
  })

  # Kappa
  total <- sum(confusion)
  Po <- accuracy

  row_marginals <- rowSums(confusion)
  col_marginals <- colSums(confusion)

  Pe <- sum(row_marginals * col_marginals) / (total^2)

  kappa <- (Po - Pe) / (1 - Pe)

  result$true_group <- true_group
  result$confusion_matrix <- confusion
  result$accuracy <- accuracy
  result$sensitivity <- sensitivity
  result$kappa <- kappa
}


