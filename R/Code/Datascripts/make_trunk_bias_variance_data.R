# make Trunk bias variance data
rm(list = ls())
library(MASS)
library(mvtnorm)
set.seed(123L)

p <- 10L
ns <- c(10L, 100L, 1000L, 10000L)
n.test <- 10000L
num.trials <- 100L

# define class-conditional distributions
mu <- cbind(-1/sqrt(1:p), 1/sqrt(1:p))
sig <- diag(p)
uY <- 0:1

# test set
Y <- sample(uY, n.test, replace = T)
X <- matrix(0, n.test, p)
for (k in uY) {
  isY <- Y == k
  X[isY, ] <- mvrnorm(sum(isY), mu = mu[, k + 1L], Sigma = sig)
}
posteriors <- matrix(0, n.test, 2L)
ccdens <- cbind(dmvnorm(X, mean = mu[, 1L], sigma = sig), dmvnorm(X, mean = mu[, 2L], sigma = sig))
posteriors[, 1L] <- ccdens[, 1L]/apply(ccdens, 1, sum)
posteriors[, 2L] <- 1 - posteriors[, 1L]

write.table(cbind(X, Y),
            file = paste0("~/R/Data/Trunk_bias_variance/Test/Trunk_bias_variance_test.csv"),
            sep = ",", row.names = F, col.names = F)

write.table(posteriors,
            file = paste0("~/R/Data/Trunk_bias_variance/Test/Trunk_bias_variance_test_posteriors.csv"),
            sep = ",", row.names = F, col.names = F)

for (i in 1:length(ns)) {
  n.train <- ns[i]
  print(paste0("n = ", n.train))
  for (trial in 1:num.trials) {
    # training set
    if (n.train == 10L) {
      Y <- c(rep(uY[1L], n.train/2), rep(uY[2L], n.train/2))
      X <- rbind(mvrnorm(n.train/2, mu = mu[, 1L], Sigma = sig), mvrnorm(n.train/2, mu = mu[, 2L], Sigma = sig))
    } else {
      Y <- sample(uY, n.train, replace = T)
      X <- matrix(0, n.train, p)
      for (k in uY) {
        isY <- Y == k
        X[isY, ] <- mvrnorm(sum(isY), mu = mu[, k + 1L], Sigma = sig)
      }
    }
    write.table(cbind(X, Y),
                file = paste0("~/R/Data/Trunk_bias_variance/Train/Trunk_bias_variance_train_n", n.train, "_trial", trial, ".csv"),
                sep = ",", row.names = F, col.names = F)
  }
}