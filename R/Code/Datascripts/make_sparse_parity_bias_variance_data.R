# make sparse parity bias variance data
rm(list = ls())

p <- 20L
ns <- c(500L, 1000L, 3000L, 5000L)
n.test <- 10000L
num.trials <- 100L

# test set
X <- matrix(runif(n.test*p, min = -1, max = 1), n.test, p)
Y <- apply(X[, 1:3], 1, function(x) sum(x > 0)%%2)

write.table(cbind(X, Y),
            file = paste0("~/R/Data/Sparse_parity_bias_variance/Test/Sparse_parity_bias_variance_test.csv"),
            sep = ",", row.names = F, col.names = F)

for (i in 1:length(ns)) {
  n.train <- ns[i]
  print(paste0("n = ", n.train))
  for (trial in 1:num.trials) {
    # training set
    X <- matrix(runif(n.train*p, min = -1, max = 1), n.train, p)
    Y <- apply(X[, 1:3], 1, function(x) sum(x > 0)%%2)
    
    write.table(cbind(X, Y),
                file = paste0("~/R/Data/Sparse_parity_bias_variance/Train/Sparse_parity_bias_variance_train_n", n.train, "_trial", trial, ".csv"),
                sep = ",", row.names = F, col.names = F)
  }
}