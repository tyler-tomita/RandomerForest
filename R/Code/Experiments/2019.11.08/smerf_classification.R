rm(list = ls())
library(rerf)
library(caret)
library(MASS)

# function to compute mode
compute.mode = function(x){
  if (length(x) > 1L) {
    ta <- table(x)
    tam <- max(ta)
    if (is.numeric(x)) {
      mod <- as.numeric(names(ta)[ta == tam])
    } else {
      mod <- names(ta)[ta == tam]
    }
    if (length(mod) > 1L) {
      mod <- sample(mod, 1L)
    }
  } else {
    mod <- x
  }
  return(mod)
}

# gmm example
ntest <- 200L
ns <- c(20L, 80L, 200L)
p <- 20L
mean <- rbind(c(-1, rep(0, p-1L)), c(1, rep(0, p-1L)))
s <- diag(0.5, p, p)
num.trials <- 10L

Xtest <- rbind(mvrnorm(ntest/2, mean[1L, ], s), mvrnorm(ntest/2, mean[2L, ], s))
Ytest <- rep(c(1L, 2L), each=ntest/2)

acc.class <- matrix(0, nrow = num.trials, ncol = length(ns))
acc.sim <- matrix(0, nrow = num.trials, ncol = length(ns))

# RerF params
num.cores <- 1L
num.trees <- 500L
min.parent <- 2L
max.depth <- 0L
random.matrix <- RandMatBinary
d <- round(p^c(1/4, 1/2, 3/4, 1)) # p^d random features evaluated at each split node
sparsity <- 1/p
num.models <- length(d)

for (k in 1:length(ns)) {
  n <- ns[k]
  cat(paste0("n = ", n, "\n"))

  Y <- rep(c(1L, 2L), each=n/2)
  Ysim <- diag(1/2, n, n)
  for (j in 1:(n-1L)) {
    for (i in (j+1L):n) {
      Ysim[i, j] <- (Y[i] == Y[j])
    }
  }
  Ysim <- Ysim + t(Ysim)

  for (trial in 1:num.trials) {
    cat(paste0("trial ", trial, "\n"))

    X <- rbind(mvrnorm(n/2, mean[1L, ], s), mvrnorm(n/2, mean[2L, ], s))

    # train and tune smerf
    forests <- rep(list(NULL), num.models)
    oob.error <- rep(NA, num.models)
    for (i in seq.int(num.models)) {
      cat(paste0("d = ", d[i], "\n"))
      forests[[i]] <- RerF(X, Ysim, FUN = random.matrix, paramList = list(p = p, d = d[i], sparsity = sparsity, prob = 0.5), min.parent = min.parent, max.depth = max.depth, num.cores = num.cores, store.impurity = TRUE, store.oob = TRUE, task="similarity", eps=0)
      predictions <- Predict(X, forests[[i]], OOB = TRUE, num.cores = num.cores)
      oob.error[i] <- mean(abs(predictions[lower.tri(predictions)] - Ysim[lower.tri(Ysim)]))
    }
    # remove all classifiers except for the one with lowest oob error
    best.idx <- order(oob.error)[1L]
    forest <- forests[[best.idx]]
    rm("forests")
    gc()

    predictions <- Predict(Xtest, forest, OOB = FALSE, num.cores = num.cores, Xtrain = X)

    acc.sim[trial, k] <- mean(apply(predictions, 1, function(x) compute.mode(Y[x == max(x)])) == Ytest)

    # train and tune rerf
    forests <- rep(list(NULL), num.models)
    oob.error <- rep(NA, num.models)
    for (i in seq.int(num.models)) {
      cat(paste0("d = ", d[i], "\n"))
      forests[[i]] <- RerF(X, Y, FUN = random.matrix, paramList = list(p = p, d = d[i], sparsity = sparsity, prob = 0.5), min.parent = min.parent, max.depth = max.depth, num.cores = num.cores, store.impurity = TRUE, store.oob = TRUE, task="classification")
      predictions <- Predict(X, forests[[i]], OOB = TRUE, num.cores = num.cores)
      oob.error[i] <- mean(predictions != Y)
    }
    # remove all classifiers except for the one with lowest oob error
    best.idx <- order(oob.error)[1L]
    forest <- forests[[best.idx]]
    rm("forests")
    gc()

    predictions <- Predict(Xtest, forest, OOB = FALSE, num.cores = num.cores, Xtrain = X)

    acc.class[trial, k] <- mean(predictions == Ytest)
  }
}
