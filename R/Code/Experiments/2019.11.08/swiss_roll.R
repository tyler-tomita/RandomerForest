rm(list = ls())
gc()
library(MASS)
library(rerf)
library(ggplot2)
library(reshape2)

RFD.transform <- function(X, Q) {
  n <- nrow(X)
  all.pairs <- combn(n, 2)
  X.pairs <- t(apply(all.pairs, 2, function(x) cbind(X[x[2], ] - X[x[1], ], 1/2*(X[x[1], ] + X[x[2], ]))))
  Y.pairs <- Q[lower.tri(Q)]
  return(list(X = X.pairs, Y = Y.pairs))
}


# function to compute precision
compute.pr <- function(retrieved, relevant) {
  return(mean(retrieved %in% relevant))
}

# function to compute average precision (AP)
compute.ave.pr <- function(predicted.order, true.order) {
  n <- length(predicted.order)
  pr.sum <- 0
  for (k in 1:n) {
    if (predicted.order[k] %in% true.order[1:k]) {
      pr.sum <- pr.sum + compute.pr(predicted.order[1:k], true.order[1:k])
    }
  }
  return(pr.sum/n)
}

# function to compute mean average precision @k (MAP@k) from predicted and true pairwise similarities
compute.map <- function(Qhat, Q, k = nrow(Qhat)-1L) {
  n <- nrow(Q)
  # sort the indices from highest to lowest similarities
  Q.sort <- t(sapply(1:n, FUN = function(rw) order(Q[rw, -rw], decreasing = TRUE)))
  Qhat.sort <- t(sapply(1:n, FUN = function(rw) order(Qhat[rw, -rw], decreasing = TRUE)))
  mapr <- mean(sapply(1:n, FUN = function(x) compute.ave.pr(Qhat.sort[x, 1:k], Q.sort[x, 1:k])))
  return(mapr)
}

# function to compute mean spearman correlation from predicted and true pairwise similarities
compute.spearman <- function(Qhat, Q) {
  n <- nrow(Q)
  spearman <- mean(sapply(1:n, FUN = function(x) cor(Qhat[x, ], Q[x, ], method = "spearman")))
  return(spearman)
}

# function to rotate a point counterclockwise by angle phi
rotate2D <- function(x, phi) {
  cos.phi <- cos(phi)
  sin.phi <- sin(phi)
  R <- matrix(c(cos.phi, -sin.phi, sin.phi, cos.phi), nrow = 2L)
  return(x%*%R)
}

# Run Swiss Roll Simulation
seed <- 1234L
set.seed(seed)

ns <- c(25L, 50L, 100L)
ntest <- 250L
num.trials <- 5L
p <- 2L

phi_start <- 0
phi_end <- 8*pi
r_start <- 1
r_end <- 0.5

# guassian noise perpendicular to tangent
sgma <- 0.01

# SmerF params
num.cores <- 1L
num.trees <- 500L
min.parent <- 2L
max.depth <- 0L
random.matrix <- RandMatRF
d <- list("SmerF" = 1:p, "RFD" = 1:(2*p))
# d <- c(1L, 2L, 4L, 8L, 16L)
sparsity <- 1/p
# num.models <- length(d)
replacement <- TRUE
bagging <- 0.2
eps <- 0

mae <- list(SmerF=matrix(NA, num.trials, length(ns)), RFD=matrix(NA, num.trials, length(ns))) # Mean Absolute Error
map1 <- list(SmerF=matrix(NA, num.trials, length(ns)), RFD=matrix(NA, num.trials, length(ns))) # Mean Average Precision
map10 <- list(SmerF=matrix(NA, num.trials, length(ns)), RFD=matrix(NA, num.trials, length(ns))) # Mean Average Precision
map250 <- list(SmerF=matrix(NA, num.trials, length(ns)), RFD=matrix(NA, num.trials, length(ns))) # Mean Average Precision
spearman <- list(SmerF=matrix(NA, num.trials, length(ns)), RFD=matrix(NA, num.trials, length(ns))) # Spearman Rank Correlation
train.time <- list(SmerF=matrix(NA, num.trials, length(ns)), RFD=matrix(NA, num.trials, length(ns))) # Training Time
for (k in 1:length(ns)) {
    ntrain <- ns[k]
    cat(paste0("ntrain = ", ntrain, "\n"))
    n <- ntrain + ntest
    for (trial in 1:num.trials) {
      cat(paste0("trial ", trial, "\n"))

      train.idx <- sample.int(n = n, size = ntrain, replace = FALSE)

      x.manifold <- (runif(n))^(8/7) # position along swiss roll if unrolled into straight line
      r <- (1 - x.manifold)*(r_start - r_end) + r_end
      phi <- (1 - x.manifold)*(phi_start - phi_end) + phi_end

      # first generate gaussian noise points in 1D, which will then be rotated
      X <- cbind(r + rnorm(n, 0, sgma), rep(0, n))

      # now rotate each point by its corresponding angle
      X <- t(sapply(1:n, FUN = function(k) rotate2D(X[k, ], phi[k])))
      plot(X[, 1], X[, 2])

      # similarity matrix
      Q <- diag(1/2, n, n)
      for (j in 1:(n-1)) {
          for (i in (j+1):n) {
              Q[i, j] <- 1 - abs(x.manifold[i] - x.manifold[j])
          }
      }
      Q <- Q + t(Q)
      # Y <- r

      # Ysort <- sapply(1:ncol(Q[-train.idx, -train.idx]), FUN = function(x) {decOrder(Q[-train.idx, -train.idx][, x][-x])})

      RFD.train.data <- RFD.transform(X[train.idx, , drop = F], Q[train.idx, train.idx])
      RFD.test.data <- RFD.transform(X[-train.idx, , drop = F], Q[-train.idx, -train.idx])
      p.RFD <- p*2L

      # train and tune smerf
      num.models <- length(d$SmerF)
      forests <- rep(list(NULL), num.models)
      oob.mapr <- rep(NA, num.models)
      tr.time <- rep(NA, num.models)
      for (i in seq.int(num.models)) {
        cat(paste0("d = ", d$SmerF[i], "\n"))

        # train
        start.time <- proc.time()
        forests[[i]] <- RerF(X[train.idx, , drop = FALSE], Q[train.idx, train.idx], FUN = random.matrix,
                             paramList = list(p = p, d = d$SmerF[i], sparsity = sparsity, prob = 0.5),
                             min.parent = min.parent, max.depth = max.depth,
                             replacement = replacement, bagging = bagging,
                             num.cores = num.cores, store.impurity = TRUE,
                             store.oob = TRUE, task = "similarity", eps = eps,
                             honesty = FALSE)
        tr.time[i] <- (proc.time() - start.time)[[3L]]

        # evaluate on oob samples
        Qhat <- Predict(X[train.idx, , drop = FALSE], forests[[i]], OOB = TRUE, num.cores = num.cores)
        oob.mapr[i] <- compute.map(Qhat, Q[train.idx, train.idx])
      }

      # remove all forests except for the one with best oob performance
      best.idx <- order(oob.mapr, decreasing = TRUE)[1L]
      # forest <- forests[[best.idx]]
      rm("forests")
      gc()
      train.time$SmerF[trial, k] <- tr.time[best.idx]

      # Qhat <- Predict(X[-train.idx, , drop = FALSE], forest, OOB = FALSE, num.cores = num.cores)
      # mapr$SmerF[trial, k] <- compute.map(Qhat, Q[-train.idx, -train.idx])
      # mae$SmerF[trial, k] <- mean(abs(Qhat[lower.tri(Qhat)] - Q[-train.idx, -train.idx][lower.tri(Q[-train.idx, -train.idx])]))
      # spearman$SmerF[trial, k] <- compute.spearman(Qhat, Q[-train.idx, -train.idx])

      # retrain with leaves updated with OOB Samples
      forest <- RerF(X[train.idx, , drop = FALSE], Q[train.idx, train.idx], FUN = random.matrix,
                     paramList = list(p = p, d = d$SmerF[best.idx], sparsity = sparsity, prob = 0.5),
                     min.parent = min.parent, max.depth = max.depth,
                     replacement = replacement, bagging = bagging,
                     num.cores = num.cores, store.impurity = TRUE,
                     store.oob = TRUE, task = "similarity", eps = eps,
                     honesty = TRUE)

      Qhat <- Predict(X[-train.idx, , drop = FALSE], forest, OOB = FALSE, num.cores = num.cores)
      map1$SmerF[trial, k] <- compute.map(Qhat, Q[-train.idx, -train.idx], 1L)
      map10$SmerF[trial, k] <- compute.map(Qhat, Q[-train.idx, -train.idx], 10L)
      map250$SmerF[trial, k] <- compute.map(Qhat, Q[-train.idx, -train.idx])
      mae$SmerF[trial, k] <- mean(abs(Qhat[lower.tri(Qhat)] - Q[-train.idx, -train.idx][lower.tri(Q[-train.idx, -train.idx])]))
      spearman$SmerF[trial, k] <- compute.spearman(Qhat, Q[-train.idx, -train.idx])

      # train and tune RFD
      num.models <- length(d$RFD)
      forests <- rep(list(NULL), num.models)
      oob.mapr <- rep(NA, num.models)
      tr.time <- rep(NA, num.models)
      for (i in seq.int(num.models)) {
        cat(paste0("d = ", d$RFD[i], "\n"))

        # train
        start.time <- proc.time()
        forests[[i]] <- RerF(RFD.train.data$X, RFD.train.data$Y, FUN = random.matrix,
                             paramList = list(p = p.RFD, d = d$RFD[i], sparsity = sparsity, prob = 0.5),
                             min.parent = min.parent, max.depth = max.depth,
                             num.cores = num.cores, store.impurity = TRUE,
                             store.oob = TRUE, task="regression", eps=0)
        tr.time[i] <- (proc.time() - start.time)[[3L]]

        # evaluate on oob samples
        Yhat <- Predict(RFD.train.data$X, forests[[i]], OOB = TRUE, num.cores = num.cores)
        Qhat <- diag(1/2, ntrain, ntrain)
        Qhat[lower.tri(Qhat)] <- Yhat
        Qhat <- Qhat + t(Qhat)
        oob.mapr[i] <- compute.map(Qhat, Q[train.idx, train.idx])
      }

      # remove all forests except for the one with best oob performance
      best.idx <- order(oob.mapr, decreasing = TRUE)[1L]
      forest <- forests[[best.idx]]
      rm("forests")
      gc()
      train.time$RFD[trial, k] <- tr.time[best.idx]

      Yhat <- Predict(RFD.test.data$X, forest, OOB = FALSE, num.cores = num.cores)
      Qhat <- diag(1/2, ntest, ntest)
      Qhat[lower.tri(Qhat)] <- Yhat
      Qhat <- Qhat + t(Qhat)
      map1$RFD[trial, k] <- compute.map(Qhat, Q[-train.idx, -train.idx], 1L)
      map10$RFD[trial, k] <- compute.map(Qhat, Q[-train.idx, -train.idx], 10L)
      map250$RFD[trial, k] <- compute.map(Qhat, Q[-train.idx, -train.idx])
      mae$RFD[trial, k] <- mean(abs(Qhat[lower.tri(Qhat)] - Q[-train.idx, -train.idx][lower.tri(Q[-train.idx, -train.idx])]))
      spearman$RFD[trial, k] <- compute.spearman(Qhat, Q[-train.idx, -train.idx])


    }
    save(map1, map10, map250, mae, spearman, train.time, file = '/home/swiss_roll.RData')
}


# oobpred <- numeric(length(oobpairs))
# for (k in 1:nrow(oobpairs)) {
#   pr <- oobpairs[k, ]
#   oobpred[k] <- sum(Q[-pr, pr[1L]]*(Q[-pr, pr[2L]]/sum(Q[-pr, pr[2L]], na.rm=TRUE)))
# }

X <- as.matrix(c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5))
Q <- cbind(c(1, 0, 0, 0, 0, 0), c(0, 1, 0.5, 0, 0, 0), c(0, 0.5, 1, 0, 0, 0), c(0, 0, 0, 1, 0.5, 0), c(0, 0, 0, 0.5, 1, 0), c(0, 0, 0, 0, 0, 1))
forest <- RerF(X, Q, trees = 1L, FUN = RandMatRF,
               paramList = list(p = 1, d = 1, sparsity = 1, prob = 0.5),
               min.parent = 2L, max.depth = 0L,
               replacement = FALSE, bagging = 0,
               num.cores = 1L, store.impurity = TRUE,
               store.oob = TRUE, task = "similarity", eps = 0,
               honesty = FALSE)

Xtest <- as.matrix(c(1.25, 0.75, -0.75))
Qhat <- Predict(Xtest, forest, OOB = FALSE, num.cores = 1L)
Qhat
