rm(list = ls())

n <- 1e5L
ps <- c(10L, 100L, 1000L)

for (p in ps) {
  X <- matrix(rnorm(n*p), nrow = n, ncol = p)

  write.table(X, file = paste0("~/work/tyler/Data/uci/gaussian_noise_p", p, ".csv"), sep = ",", row.names = F, col.names = F)
}
