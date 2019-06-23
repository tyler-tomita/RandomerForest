rm(list = ls())

ps <- c(2L, 8L)
ns <- c(4L, 10L, 100L, 1000L)
n.test <- 10000L
num.trials <- 100L
oversample.factor <- 1/c(pi/4, pi^4/24/(2^8))*4

for (j in 1:length(ps)) {
  p <- ps[j]
  print(paste0("p = ", p))
  # test set for orthogonal hyperspheres
  go <- T
  while (go) {
    X1 <- matrix(runif(round(n.test/2*oversample.factor[j])*p, -1/2, 1/2), round(n.test/2*oversample.factor[j]), p)
    X1.reject <- X1[apply(X1, 1, function(x) sqrt(sum(x^2)) <= 1/2), , drop = F]
    go <- nrow(X1.reject) < n.test/2
  }
  X1 <- t(t(X1.reject[1:(n.test/2), ]) - c(1/2, rep(0, p - 1L)))
  Y1 <- rep(0, n.test/2)
  go <- T
  while (go) {
    X2 <- matrix(runif(round(n.test/2*oversample.factor[j])*p, -1/2, 1/2), round(n.test/2*oversample.factor[j]), p)
    X2.reject <- X2[apply(X2, 1, function(x) sqrt(sum(x^2)) <= 1/2), , drop = F]
    go <- nrow(X2.reject) < n.test/2
  }
  X2 <- t(t(X2.reject[1:(n.test/2), ]) + c(1/2, rep(0, p - 1L)))
  Y2 <- rep(1, n.test/2)
  write.table(cbind(rbind(X1, X2), as.matrix(c(Y1, Y2))),
              file = paste0("~/R/Data/Hypersphere/Test/Hypersphere_orthogonal_test_p", p, ".csv"),
              sep = ",", row.names = F, col.names = F)
  
  # test set for oblique hyperspheres
  go <- T
  while (go) {
    X1 <- matrix(runif(round(n.test/2*oversample.factor[j])*p, -1/2, 1/2), round(n.test/2*oversample.factor[j]), p)
    X1.reject <- X1[apply(X1, 1, function(x) sqrt(sum(x^2)) <= 1/2), , drop = F]
    go <- nrow(X1.reject) < n.test/2
  }
  X1 <- t(t(X1.reject[1:(n.test/2), ]) - rep(1/sqrt(p)/2, p))
  Y1 <- rep(0, n.test/2)
  go <- T
  while (go) {
    X2 <- matrix(runif(round(n.test/2*oversample.factor[j])*p, -1/2, 1/2), round(n.test/2*oversample.factor[j]), p)
    X2.reject <- X2[apply(X2, 1, function(x) sqrt(sum(x^2)) <= 1/2), , drop = F]
    go <- nrow(X2.reject) < n.test/2
  }
  X2 <- t(t(X2.reject[1:(n.test/2), ]) + rep(1/sqrt(p)/2, p))
  Y2 <- rep(1, n.test/2)
  write.table(cbind(rbind(X1, X2), as.matrix(c(Y1, Y2))),
              file = paste0("~/R/Data/Hypersphere/Test/Hypersphere_oblique_test_p", p, ".csv"),
              sep = ",", row.names = F, col.names = F)
  
  for (i in 1:length(ns)) {
    n.train <- ns[i]
    print(paste0("ntrain = ", n.train))
    for (trial in 1:num.trials) {
      # training set for orthogonal hyperspheres
      while (go) {
        X1 <- matrix(runif(round(n.train/2*oversample.factor[j])*p, -1/2, 1/2), round(n.train/2*oversample.factor[j]), p)
        X1.reject <- X1[apply(X1, 1, function(x) sqrt(sum(x^2)) <= 1/2), , drop = F]
        go <- nrow(X1.reject) < n.train/2
      }
      X1 <- t(t(X1.reject[1:(n.train/2), ]) - c(1/2, rep(0, p - 1L)))
      Y1 <- rep(0, n.train/2)
      while (go) {
        X2 <- matrix(runif(round(n.train/2*oversample.factor[j])*p, -1/2, 1/2), round(n.train/2*oversample.factor[j]), p)
        X2.reject <- X2[apply(X2, 1, function(x) sqrt(sum(x^2)) <= 1/2), , drop = F]
        go <- nrow(X2.reject) < n.train/2
      }
      X2 <- t(t(X2.reject[1:(n.train/2), ]) + c(1/2, rep(0, p - 1L)))
      Y2 <- rep(1, n.train/2)
      write.table(cbind(rbind(X1, X2), as.matrix(c(Y1, Y2))),
                  file = paste0("~/R/Data/Hypersphere/Train/Hypersphere_orthogonal_train_p", p, "_n", n.train, "_trial", trial, ".csv"),
                  sep = ",", row.names = F, col.names = F)
      
      # training set for oblique hyperspheres
      while (go) {
        X1 <- matrix(runif(round(n.train/2*oversample.factor[j])*p, -1/2, 1/2), round(n.train/2*oversample.factor[j]), p)
        X1.reject <- X1[apply(X1, 1, function(x) sqrt(sum(x^2)) <= 1/2), , drop = F]
        go <- nrow(X1.reject) < n.train/2
      }
      X1 <- t(t(X1.reject[1:(n.train/2), ]) - rep(1/sqrt(p), p))
      Y1 <- rep(0, n.train/2)
      while (go) {
        X2 <- matrix(runif(round(n.train/2*oversample.factor[j])*p, -1/2, 1/2), round(n.train/2*oversample.factor[j]), p)
        X2.reject <- X2[apply(X2, 1, function(x) sqrt(sum(x^2)) <= 1/2), , drop = F]
        go <- nrow(X2.reject) < n.train/2
      }
      X2 <- t(t(X2.reject[1:(n.train/2), ]) + rep(1/sqrt(p), p))
      Y2 <- rep(1, n.train/2)
      write.table(cbind(rbind(X1, X2), as.matrix(c(Y1, Y2))),
                  file = paste0("~/R/Data/Hypersphere/Train/Hypersphere_oblique_train_p", p, "_n", n.train, "_trial", trial, ".csv"),
                  sep = ",", row.names = F, col.names = F)
    }
  }
}