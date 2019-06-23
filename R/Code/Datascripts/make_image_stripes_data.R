rm(list = ls())

ns <- c(10L, 20L, 50L)
imh <- 20L
imw <- 20L
p <- as.integer(imh*imw)
ntest <- 10000L
num.trials <- 10L

rmat0 <- function(...) {
  M <- matrix(0, imh, imw)
  M[sample.int(imh, 5L), ] <- 1
  return(M)
}
rmat1 <- function(...) {
  M <- matrix(0, imh, imw)
  M[, sample.int(imw, 5L)] <- 1
  return(M)
}

Xtest <- matrix(0, ntest, p)
# class 0 (horizontal stripes)
Xtest[1:(ntest/2), ] <- t(sapply(seq.int(ntest/2), FUN = rmat0))
# class 1 (vertical stripes)
Xtest[(ntest/2 + 1):ntest, ] <- t(sapply(seq.int(ntest/2), FUN = rmat1))
Ytest <- c(rep(0L, ntest/2), rep(1L, ntest/2))
write.table(cbind(Xtest, Ytest), file = "~/R/Data/Image_stripes/Test/Image_stripes_test_set.csv", quote = F, sep = ",", row.names = F, col.names = F)

for (ntrain in ns) {
  for (trial in 1:num.trials) {
    Xtrain <- matrix(0, ntrain, p)
    # class 0 (horizontal stripes)
    Xtrain[1:(ntrain/2), ] <- t(sapply(seq.int(ntrain/2), FUN = rmat0))
    # class 1 (vertical stripes)
    Xtrain[(ntrain/2 + 1):ntrain, ] <- t(sapply(seq.int(ntrain/2), FUN = rmat1))
    Ytrain <- c(rep(0L, ntrain/2), rep(1L, ntrain/2))
    write.table(cbind(Xtrain, Ytrain), file = paste0("~/R/Data/Image_stripes/Train/Image_stripes_train_set_n", ntrain, "_trial", trial, ".csv"), quote = F, sep = ",", row.names = F, col.names = F)
  }
}
