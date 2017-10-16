compiler::setCompilerOptions(optimize = 3)

BiasVariance <- compiler::cmpfun(function(predictions, posteriors) {
  n <- nrow(predictions)
  N <- ncol(predictions)
  labels <- sort(unique(c(predictions)))
  K <- length(labels)
  phats <- matrix(0, n, K)
  for (i in 1:K) {
    is.Y <- sapply(1:N, function(x) predictions[, x] == labels[i])
    phats[, i] <- apply(is.Y, 1L, mean)
  }
  prediction.idx <- max.col(phats)
  phat.max <- sapply(1:n, function(x) phats[x, prediction.idx[x]])
  Yhat <- labels[prediction.idx]
  bayes.idx <- max.col(posteriors)
  posterior.max <- sapply(1:n, function(x) posteriors[x, bayes.idx[x]])
  Ybayes <- labels[bayes.idx]
  B <- mean(Yhat != Ybayes)
  V <- mean(1 - phat.max)
  SE <- mean(posterior.max - sapply(1:n, function(x) posteriors[prediction.idx[x]]))
  VE <- mean(sapply(1:n, function(x) posteriors[prediction.idx[x]]) - apply(posteriors*phats, 1L, sum))
  BE <- mean(1 - posterior.max)
  
  return(bv.decomp <- list(B = B, V = V, SE = SE, VE = VE, BE = BE))
})