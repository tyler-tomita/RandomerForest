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

# function to compute mean average precision (MAP) from predicted and true pairwise similarities
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

# Run Simulation
seed <- 1234L
set.seed(seed)

ns <- c(10, 20L, 40L, 80L, 120L)
ntest <- 160L
num.trials <- 10L
p <- 2L

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

mae <- list(SmerF=matrix(NA, num.trials, length(ns)),
            SmerF.h=matrix(NA, num.trials, length(ns)),
            RFD=matrix(NA, num.trials, length(ns)),
            RFD.h=matrix(NA, num.trials, length(ns))) # Mean Absolute Error
map1 <- list(SmerF=matrix(NA, num.trials, length(ns)),
             SmerF.h=matrix(NA, num.trials, length(ns)),
             RFD=matrix(NA, num.trials, length(ns)),
             RFD.h=matrix(NA, num.trials, length(ns))) # Mean Average Precision
map10 <- list(SmerF=matrix(NA, num.trials, length(ns)),
              SmerF.h=matrix(NA, num.trials, length(ns)),
              RFD=matrix(NA, num.trials, length(ns)),
              RFD.h=matrix(NA, num.trials, length(ns))) # Mean Average Precision
map200 <- list(SmerF=matrix(NA, num.trials, length(ns)),
               SmerF.h=matrix(NA, num.trials, length(ns)),
               RFD=matrix(NA, num.trials, length(ns)),
               RFD.h=matrix(NA, num.trials, length(ns))) # Mean Average Precision
spearman <- list(SmerF=matrix(NA, num.trials, length(ns)),
                 SmerF.h=matrix(NA, num.trials, length(ns)),
                 RFD=matrix(NA, num.trials, length(ns)),
                 RFD.h=matrix(NA, num.trials, length(ns))) # Spearman Rank Correlation
train.time <- list(SmerF=matrix(NA, num.trials, length(ns)),
                   RFD=matrix(NA, num.trials, length(ns))) # Training Time
for (k in 1:length(ns)) {
    ntrain <- ns[k]
    cat(paste0("ntrain = ", ntrain, "\n"))
    n <- ntrain + ntest
    for (trial in 1:num.trials) {
      cat(paste0("trial ", trial, "\n"))

      train.idx <- sample.int(n = n, size = ntrain, replace = FALSE)

      r <- sqrt(runif(n))
      r[] <- r[order(r)]
      theta <- runif(n, 0, 2*pi)

      X <- cbind(r*cos(theta), r*sin(theta))

      # similarity matrix
      Q <- diag(1/2, n, n)
      for (j in 1:(n-1)) {
          for (i in (j+1):n) {
              Q[i, j] <- 1 - abs(r[i] - r[j])
          }
      }
      Q <- Q + t(Q)

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
                             min.parent = min.parent, max.depth = max.depth, trees = num.trees,
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
      forest <- forests[[best.idx]]
      rm("forests")
      gc()
      train.time$SmerF[trial, k] <- tr.time[best.idx]

      Qhat <- Predict(X[-train.idx, , drop = FALSE], forest, OOB = FALSE, num.cores = num.cores)
      map1$SmerF[trial, k] <- compute.map(Qhat, Q[-train.idx, -train.idx], 1L)
      map10$SmerF[trial, k] <- compute.map(Qhat, Q[-train.idx, -train.idx], 10L)
      map200$SmerF[trial, k] <- compute.map(Qhat, Q[-train.idx, -train.idx])
      mae$SmerF[trial, k] <- mean(abs(Qhat[lower.tri(Qhat)] - Q[-train.idx, -train.idx][lower.tri(Q[-train.idx, -train.idx])]))
      spearman$SmerF[trial, k] <- compute.spearman(Qhat, Q[-train.idx, -train.idx])

      # retrain with leaves updated with OOB Samples
      forest <- RerF(X[train.idx, , drop = FALSE], Q[train.idx, train.idx], FUN = random.matrix,
                     paramList = list(p = p, d = d$SmerF[best.idx], sparsity = sparsity, prob = 0.5),
                     min.parent = min.parent, max.depth = max.depth, trees = num.trees,
                     replacement = replacement, bagging = bagging,
                     num.cores = num.cores, store.impurity = TRUE,
                     store.oob = TRUE, task = "similarity", eps = eps,
                     honesty = TRUE)

      Qhat <- Predict(X[-train.idx, , drop = FALSE], forest, OOB = FALSE, num.cores = num.cores)
      map1$SmerF.h[trial, k] <- compute.map(Qhat, Q[-train.idx, -train.idx], 1L)
      map10$SmerF.h[trial, k] <- compute.map(Qhat, Q[-train.idx, -train.idx], 10L)
      map200$SmerF.h[trial, k] <- compute.map(Qhat, Q[-train.idx, -train.idx])
      mae$SmerF.h[trial, k] <- mean(abs(Qhat[lower.tri(Qhat)] - Q[-train.idx, -train.idx][lower.tri(Q[-train.idx, -train.idx])]))
      spearman$SmerF.h[trial, k] <- compute.spearman(Qhat, Q[-train.idx, -train.idx])

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
                             min.parent = min.parent, max.depth = max.depth, trees = num.trees,
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
      map200$RFD[trial, k] <- compute.map(Qhat, Q[-train.idx, -train.idx])
      mae$RFD[trial, k] <- mean(abs(Qhat[lower.tri(Qhat)] - Q[-train.idx, -train.idx][lower.tri(Q[-train.idx, -train.idx])]))
      spearman$RFD[trial, k] <- compute.spearman(Qhat, Q[-train.idx, -train.idx])

      # retrain with leaves updated with OOB Samples
      forest <- RerF(RFD.train.data$X, RFD.train.data$Y, FUN = random.matrix,
                     paramList = list(p = p.RFD, d = d$RFD[best.idx], sparsity = sparsity, prob = 0.5),
                     min.parent = min.parent, max.depth = max.depth, trees = num.trees,
                     replacement = replacement, bagging = bagging,
                     num.cores = num.cores, store.impurity = TRUE,
                     store.oob = TRUE, task="regression", eps=0, honesty = TRUE)

      Yhat <- Predict(RFD.test.data$X, forest, OOB = FALSE, num.cores = num.cores)
      Qhat <- diag(1/2, ntest, ntest)
      Qhat[lower.tri(Qhat)] <- Yhat
      Qhat <- Qhat + t(Qhat)
      map1$RFD.h[trial, k] <- compute.map(Qhat, Q[-train.idx, -train.idx], 1L)
      map10$RFD.h[trial, k] <- compute.map(Qhat, Q[-train.idx, -train.idx], 10L)
      map200$RFD.h[trial, k] <- compute.map(Qhat, Q[-train.idx, -train.idx])
      mae$RFD.h[trial, k] <- mean(abs(Qhat[lower.tri(Qhat)] - Q[-train.idx, -train.idx][lower.tri(Q[-train.idx, -train.idx])]))
      spearman$RFD.h[trial, k] <- compute.spearman(Qhat, Q[-train.idx, -train.idx])
    }
    save(map1, map10, map200, mae, spearman, train.time, file = '/home/radial_similarity2.RData')
}

# forest <- RerF(RFD.train.data$X, RFD.train.data$Y, FUN = random.matrix,
#                paramList = list(p = p.RFD, d = d$RFD[i], sparsity = sparsity, prob = 0.5),
#                min.parent = min.parent, max.depth = max.depth, trees = num.trees,
#                replacement = replacement, bagging = bagging,
#                num.cores = num.cores, store.impurity = TRUE,
#                store.oob = TRUE, task="regression", eps=0, honesty = TRUE)
#
# any(sapply(1:num.trees, FUN = function(x) any(is.na(forest$trees[[x]]$leafPred))))

# sapply(map1, FUN = function(x) apply(x, 2, mean))
