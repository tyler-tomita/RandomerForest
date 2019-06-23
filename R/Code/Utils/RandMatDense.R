# Create a Random Matrix 

library(RcppZiggurat)
library(compiler)

RandMatDense <- cmpfun(
  function(mat.options) {
    p <- mat.options[[1L]] # number of dimensions
    
    d <- mat.options[[2L]] # this determines the number of columns in the projection matrix.
    method <- mat.options[[3L]] # defines the distribution of the random projection matrix
    #Create the random matrix, a sparse matrix of 1's, -1's, and 0's.
    if (method == "binary") {
      rho <- mat.options[[4L]]
      nnzs <- round(p*d*rho)
      ind <- sort(sample.int((p*d), nnzs, replace = F))
      random.matrix <- matrix(0, p, d)
      random.matrix[ind] <- sample(c(1L, -1L), nnzs, replace = T)
    } else if (method == "continuous") {
      rho <- mat.options[[4L]]
      nnzs <- round(p*d*rho)
      ind <- sort(sample.int((p*d), nnzs, replace = F))
      random.matrix <- matrix(0, p, d)
      random.matrix[ind] <- zrnorm(nnzs)
    }
    return(random.matrix)
  }
)