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


# Run Simulation
seed <- 1234L
set.seed(seed)

n <- 10L
p <- 2L

# SmerF params
num.cores <- 1L
num.trees <- 500L
min.parent <- 2L
max.depth <- 0L
random.matrix <- RandMatRF
d <- 1L
sparsity <- 1/p
# num.models <- length(d)
replacement <- TRUE
bagging <- 0.2
eps <- 0

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

RFD.data <- RFD.transform(X, Q)
p.RFD <- p*2L

forest <- RerF(RFD.data$X, RFD.data$Y, FUN = random.matrix,
               paramList = list(p = p.RFD, d = d, sparsity = sparsity, prob = 0.5),
               min.parent = min.parent, max.depth = max.depth, trees = num.trees,
               num.cores = num.cores, store.impurity = TRUE,
               store.oob = TRUE, task="regression", eps=0, honesty = TRUE)

# any(sapply(1:num.trees, FUN = function(x) any(is.na(forest$trees[[x]]$leafPred))))
