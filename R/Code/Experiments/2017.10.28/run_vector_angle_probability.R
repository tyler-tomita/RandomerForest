rm(list = ls())

source("~/RandomerForest/R/Code/Utils/RandMatDense.R")

VectorAngle <- function(u, v) {
  theta <- acos(u%*%v/sqrt(sum(u^2))/sqrt(sum(v^2)))*180/pi
}

IsClose <- function(u, v, target.theta) {
  theta <- VectorAngle(u, v)
  is.close <- (theta <= target.theta) || (theta >= (180 - target.theta))
}

ps <- 2^(1:6)
num.p <- length(ps)
deltas <- c(1, 10, 22.5, 45)
num.delta <- length(deltas)

num.trials <- 1e3L
random.matrix <- "continuous"

pclose1 <- list(sparse = matrix(0, length(deltas), length(ps)), dense = matrix(0, length(deltas), length(ps)))
sem1 <- list(sparse = matrix(0, length(deltas), length(ps)), dense = matrix(0, length(deltas), length(ps)))
for (j in 1:num.p) {
  p <- ps[j]
  d <- p
  target.vector <- c(1, rep(0, p - 1L))
  
  for (i in 1:num.delta) {
    delta <- deltas[i]
    is.close.sparse <- logical(num.trials)
    is.close.dense <- logical(num.trials)
    for (trial in 1:num.trials) {
      A.sparse <- RandMatDense(mat.options = list(p, d, random.matrix, 1/p))
      is.close.sparse[trial] <- any(apply(A.sparse, 2L, function(x) IsClose(target.vector, x, delta)), na.rm = T)
      A.dense <- RandMatDense(mat.options = list(p, d, random.matrix, 1))
      is.close.dense[trial] <- any(apply(A.dense, 2L, function(x) IsClose(target.vector, x, delta)), na.rm = T)
    }
    pclose1$sparse[i, j] <- mean(is.close.sparse)
    sem1$sparse[i, j] <- sd(is.close.sparse)/sqrt(num.trials)
    pclose1$dense[i, j] <- mean(is.close.dense)
    sem1$dense[i, j] <- sd(is.close.dense)/sqrt(num.trials)
  }
}

pclose2 <- list(sparse = matrix(0, length(deltas), length(ps)), dense = matrix(0, length(deltas), length(ps)))
sem2 <- list(sparse = matrix(0, length(deltas), length(ps)), dense = matrix(0, length(deltas), length(ps)))
for (j in 1:num.p) {
  p <- ps[j]
  target.vector <- rep(1, p)
  for (i in 1:num.delta) {
    delta <- deltas[i]
    is.close.sparse <- logical(num.trials)
    is.close.dense <- logical(num.trials)
    for (trial in 1:num.trials) {
      A.sparse <- RandMatDense(mat.options = list(p, d, random.matrix, 1/p))
      is.close.sparse[trial] <- any(apply(A.sparse, 2L, function(x) IsClose(target.vector, x, delta)), na.rm = T)
      A.dense <- RandMatDense(mat.options = list(p, d, random.matrix, 1))
      is.close.dense[trial] <- any(apply(A.dense, 2L, function(x) IsClose(target.vector, x, delta)), na.rm = T)
    }
    pclose2$sparse[i, j] <- mean(is.close.sparse)
    sem2$sparse[i, j] <- sd(is.close.sparse)/sqrt(num.trials)
    pclose2$dense[i, j] <- mean(is.close.dense)
    sem2$dense[i, j] <- sd(is.close.dense)/sqrt(num.trials)
  }
}

df <- data.frame(target.vector = c(rep("sparse", num.p*num.delta*2), rep("dense", num.p*num.delta*2)),
                 random.vector = rep(c(rep("sparse", num.p*num.delta), rep("dense", num.p*num.delta)), 2),
                 p = rep(c(sapply(ps, function(x) rep(x, num.delta))), 4),
                 theta = rep(deltas, num.p*4),
                 prob.success = c(pclose1$sparse, pclose1$dense, pclose2$sparse, pclose2$dense),
                 sem.success = c(sem1$sparse, sem1$dense, sem2$sparse, sem2$dense))


save(list = c("df", "d"), file = "~/RandomerForest/R/Results/2017.10.28/vector_angle_probability.RData")