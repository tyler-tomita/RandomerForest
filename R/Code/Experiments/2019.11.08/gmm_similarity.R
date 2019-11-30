library(MASS)
library(rerf)
library(ggplot2)
library(reshape2)
library(FastKNN)

# Define function for tuning number of sampled projections based on out-of-bag-error
tuneSimForest <- function(X, Ysim, projectionDist, numProjections) {
    oob.mae <- rep(NA, length(numProjections))
    n <- nrow(X)
    p <- ncol(X)

    for (i in 1:length(numProjections)) {
        forest <- RerFSim(X, Ysim,
                          bagging = 0.2, replacement = TRUE,
                          trees = 200L, min.parent = 2L, max.depth = 0,
                          mat.options = list(p = p, d = numProjections[i], random.matrix = projectionDist, rho = 1/p, 1/2),
                          store.oob = TRUE, seed = 123L, num.cores = 1L)

        Yoob <- OOBPredictSim(X, forest, num.cores = 1L)
        Yoob <- array(cbind(Yoob, t(Yoob)), dim = c(n, n, 2L))
        Yoob <- apply(Yoob, c(1, 2), function(a) sum(a, na.rm = TRUE)) - diag(n)
        oob.mae[i] <- mean(abs(Yoob[lower.tri(Yoob)] - Ysim[lower.tri(Yoob)]))
    }
    mn <- min(oob.mae)
    best.idx <- which(oob.mae == mn)
    if (length(best.idx) > 1L) {
        best.idx <- sample(best.idx, 1L)
    }
    return(list(oob.mae = oob.mae, best.idx = best.idx))
}

tuneRFD <- function(X, Y, Ysim, projectionDist, numProjections) {
    oob.mae <- rep(NA, length(numProjections))
    n <- nrow(X)
    p <- ncol(X)

    for (i in 1:length(numProjections)) {
        forest <- RFD(X, Y,
                         bagging = 0.2, replacement = TRUE,
                         trees = 200L, min.parent = 2L, max.depth = 0,
                         mat.options = list(p = ncol(X)*2, d = numProjections[i], random.matrix = projectionDist, rho = 1/ncol(X)/2, 1/2),
                         store.oob = TRUE, seed = 123L, num.cores = 2L)

        Yoob <- OOBPredictRFD(X, forest, num.cores = 1L)
        oob.mae[i] <- mean(abs(Yoob[lower.tri(Yoob)] - Ysim[lower.tri(Yoob)]))
    }
    mn <- min(oob.mae)
    best.idx <- which(oob.mae == mn)
    if (length(best.idx) > 1L) {
        best.idx <- sample(best.idx, 1L)
    }
    return(list(oob.mae = oob.mae, best.idx = best.idx))
}

# Define function for computing average precision at rank K
compute.avePk <- function(truth, estimate, k) {
    is.rel <- as.integer(estimate[1:k] %in% truth[1:k])
    num.rel <- sum(is.rel)
    Pk <- sapply(1:k, function(n) {mean(estimate[1:n] %in% truth[1:n])})
    avePk <- sum(Pk*is.rel)/num.rel
    return(avePk)
}

decOrder <- function(x) {order(x, decreasing = TRUE)}


# Run Simulation
seed <- 123L
set.seed(seed)

ns <- c(20L, 40L, 80L)
ntest <- 500L
num.trials <- 5L
p <- 2L

mae <- list(SmerF=matrix(NA, num.trials, length(ns)), RFD=matrix(NA, num.trials, length(ns))) # Mean Absolute Error
map1 <- list(SmerF=matrix(NA, num.trials, length(ns)), RFD=matrix(NA, num.trials, length(ns))) # Mean Average Precision
map10 <- list(SmerF=matrix(NA, num.trials, length(ns)), RFD=matrix(NA, num.trials, length(ns))) # Mean Average Precision
map100 <- list(SmerF=matrix(NA, num.trials, length(ns)), RFD=matrix(NA, num.trials, length(ns))) # Mean Average Precision
map500 <- list(SmerF=matrix(NA, num.trials, length(ns)), RFD=matrix(NA, num.trials, length(ns))) # Mean Average Precision
train.time <- list(SmerF=matrix(NA, num.trials, length(ns)), RFD=matrix(NA, num.trials, length(ns))) # Training Time
for (n.ind in 1:length(ns)) {
    ntrain <- ns[n.ind]
    n <- ntrain + ntest
    for (trial in 1:num.trials) {
        print(n.ind)
        print(trial)

        train.idx <- sample.int(n = n, size = ntrain, replace = FALSE)

        mu <- t(matrix(c(-2, -2, 2, 2), nrow = 2L))
        sg <- diag(c(1, 1))
        Y <- c(rep(0L, n/2), rep(1L, n/2))

        X <- matrix(nrow = n, ncol = p)
        Ysim <- diag(n)
        for (k in unique(Y)) {
            X[Y==k, ] <- mvrnorm(n = sum(Y==k), mu = mu[k+1L, ], Sigma = sg)
            Ysim[Y==k, Y==k] <- 1
        }

        Ysort <- sapply(1:ncol(Ysim[-train.idx, -train.idx]), FUN = function(x) {decOrder(Ysim[-train.idx, -train.idx][, x][-x])})

        # sim.pairs <- rbind(combn(1:(n/2), 2), combn((n/2+1):n, 2))
        # dissim.pairs <- expand.grid((n/2+1):n, 1:(n/2))[, c(2,1)]


        # tune SmerF
        ds <- c(1L, 2L, 4L, 8L)
        res <- tuneSimForest(X[train.idx, ], Ysim[train.idx, train.idx], 'binary', ds)

        # retrain SmerF with best hyperparameter
        start.time <- proc.time()

        forest <- RerFSim(X[train.idx, ], Ysim[train.idx, train.idx],
                          trees = 200L, min.parent = 2L, max.depth = 0,
                          mat.options = list(p = p, d = ds[res$best.idx], random.matrix = 'binary', rho = 1/p, 1/2),
                          store.oob = FALSE, seed = 123L, num.cores = 1L)

        train.time$SmerF[trial, n.ind] <- (proc.time() - start.time)[[3L]]

        # compute SmerF MAE on test set
        Yhats <- PredictSim(X = X[-train.idx, ], forest = forest, num.cores = 1L)
        Yhats <- Yhats + t(Yhats) - diag(ntest)
        mae$SmerF[trial, n.ind] <- mean(abs(Yhats[lower.tri(Yhats)] - Ysim[-train.idx, -train.idx][lower.tri(Yhats)]))
        YhatSort <- sapply(1:ncol(Yhats), FUN = function(x) {decOrder(Yhats[, x][-x])})
        map1$SmerF[trial, n.ind] <- mean(sapply(1:ntest, function(x) {compute.avePk(Ysort[, x], YhatSort[, x], 1L)}))
        map10$SmerF[trial, n.ind] <- mean(sapply(1:ntest, function(x) {compute.avePk(Ysort[, x], YhatSort[, x], 10L)}))
        map100$SmerF[trial, n.ind] <- mean(sapply(1:ntest, function(x) {compute.avePk(Ysort[, x], YhatSort[, x], 100L)}))
        map500$SmerF[trial, n.ind] <- mean(sapply(1:ntest, function(x) {compute.avePk(Ysort[, x], YhatSort[, x], 500L)}))

        # tune RFD
        ds <- c(1L, 2L, 3L, 4L)
        res <- tuneRFD(X[train.idx, ], Y[train.idx], Ysim[train.idx, train.idx], 'rf', ds)

        # retrain RFDReg with best hyperparameter
        start.time <- proc.time()

        forest <- RFD(X[train.idx, ], Y[train.idx],
                 bagging = 0.2, replacement = TRUE,
                 trees = 200L, min.parent = 2L, max.depth = 0,
                 mat.options = list(p = ncol(X)*2, d = ds[res$best.idx], random.matrix = 'rf', rho = 1/ncol(X)/2, 1/2),
                 store.oob = FALSE, seed = 123L, num.cores = 1L)

        train.time$RFD[trial, n.ind] <- (proc.time() - start.time)[[3L]]

        # compute RFDReg MAE on test set
        Yhats <- PredictRFD(X = X[-train.idx, ], forest = forest, num.cores = 1L)
        YhatSort <- sapply(1:ncol(Yhats), FUN = function(x) {decOrder(Yhats[, x][-x])})
        mae$RFD[trial, n.ind] <- mean(abs(Yhats[lower.tri(Yhats)] - Ysim[-train.idx, -train.idx][lower.tri(Yhats)]))
        map1$RFD[trial, n.ind] <- mean(sapply(1:ntest, function(x) {compute.avePk(Ysort[, x], YhatSort[, x], 1L)}))
        map10$RFD[trial, n.ind] <- mean(sapply(1:ntest, function(x) {compute.avePk(Ysort[, x], YhatSort[, x], 10L)}))
        map100$RFD[trial, n.ind] <- mean(sapply(1:ntest, function(x) {compute.avePk(Ysort[, x], YhatSort[, x], 100L)}))
        map500$RFD[trial, n.ind] <- mean(sapply(1:ntest, function(x) {compute.avePk(Ysort[, x], YhatSort[, x], 500L)}))
    }
    save(mae, map1, map10, map100, map500, train.time, file = '/Users/tyler/RandomerForest/R/Results/2018.12.03/gmm_similarity.RData')
}
