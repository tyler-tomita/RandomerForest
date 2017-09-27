rm(list = ls())
source("~/work/tyler/RerF/Code/Classifiers/rfr_function.R")

nTrials <- 20L

filePath <- "~/work/tyler/Data/tests/"
outFile <- "~/work/tyler/RerF/Results/Tests/test_results_2017_08_22.RData"
# filePath <- "~/tests/"
# outFile <- "~/RerF/Results/Tests/test_results_2017_08_22.RData"

dataSets <- c("Sparse_parity", "Trunk", "Orthant", "mnist")

trainTime <- list()
OOBError <- list()
numNodes <- list()
testError <- list()
seed <- list()

for (ds in 1:length(dataSets)) {
# for (ds in 4L) {
  D <- dataSets[ds]
  print(paste("Dataset:", D))
  
  trainFile <- paste(filePath, D, "_train.csv", sep = "")
  testFile <- paste(filePath, D, "_test.csv", sep = "")
  
  Xtrain <- as.matrix(read.table(trainFile, header = F, sep = ",", quote = "", row.names = NULL))
  ntrain <- nrow(Xtrain)
  p <- ncol(Xtrain) - 1L
  Ytrain <- as.integer(Xtrain[, p + 1L]) + 1L
  Xtrain <- Xtrain[, -(p + 1L)]
  
  Xtest <- as.matrix(read.table(testFile, header = F, sep = ",", quote = "", row.names = NULL))
  ntest <- nrow(Xtest)
  Ytest <- as.integer(Xtest[, p + 1L] + 1L)
  Xtest <- Xtest[, -(p + 1L)] 
  
  labels <- unique(c(Ytrain,Ytest))
  nClasses <- length(labels)
  
  if (D == "mnist") {
    mtrys <- ceiling(c(p^(1/2), 10*p))
  } else {
    mtrys <- ceiling(p^c(1/2, 2))
  }
  # mtrys <-p^2
  # mtrys <- p
  # mtrys <- 10*p
  # mtrys <- ceiling(sqrt(p))
    
  nTrees <- 500L
  stratify <- T
  replacement <- T
  randomMatrix <- "binary"
  sparsity <- 1/p
  rotate <- F
  MinParent <- 2L
  MaxDepth <- "inf"
  COOB <- T
  CNS <- F
  NumCores <- 16L
  comp.mode <- "batch"
  rank.transform <- F
  
  train.time <- matrix(0, nTrials, length(mtrys))
  oob.error <- matrix(0, nTrials, length(mtrys))
  num.nodes <- matrix(0L, nTrials, length(mtrys))
  test.error <- matrix(0, nTrials, length(mtrys))
  seed[[D]] <- matrix(sample.int(1e4L, nTrials*length(mtrys)), nTrials, length(mtrys))
  
  for (trial in 1:nTrials) {
    print(paste("Trial", trial))
    for (m in 1:length(mtrys)) {
      # train the classifier
      cat("training\n")
      start.time <- proc.time()
      forest <- rerf(Xtrain, Ytrain, MinParent = MinParent, trees = nTrees, MaxDepth = MaxDepth, replacement = replacement,
                     stratify = stratify, FUN = randmat, options = list(p, mtrys[m], randomMatrix, sparsity), COOB = COOB,
                     CNS = CNS, NumCores = NumCores, seed = seed[[D]][trial, m], rotate = rotate, rank.transform = rank.transform)
      train.time[trial, m] <- (proc.time() - start.time)[[3]]
      cat("training complete\n")
      cat(paste0("elapsed time: ", train.time[trial, m], " sec\n"))
      
      # compute out-of-bag predictions
      cat("OOB\n")
      scores <- OOBpredict(Xtrain, forest, NumCores, comp.mode = comp.mode)
      oob.error[trial, m] <- mean(max.col(scores) != Ytrain)
      
      # compute average number of nodes per tree
      cat("node size\n")
      num.nodes[trial, m] <- mean(sapply(forest, FUN = function(tree) length(tree$treeMap)))
      
      # compute predictions on test set  
      cat("test predictions\n")
      scores <- predict(Xtest, forest, NumCores, comp.mode = comp.mode)
      test.error[trial, m] <- mean(max.col(scores) != Ytest)
      rm(forest)
      gc()
    }
  }
  trainTime[[D]] <- train.time
  OOBError[[D]] <- oob.error
  numNodes[[D]] <- num.nodes
  testError[[D]] <- test.error
}

save(file = outFile, list = c("trainTime", "OOBError", "numNodes", "testError", "seed"))
