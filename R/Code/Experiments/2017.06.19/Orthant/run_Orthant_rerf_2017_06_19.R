# evaluate RerF on Orthant

rm(list=ls())
options(scipen = 999)

rerfPath <- "~/"
dataPath <- "~/R/Data/Orthant/dat/Raw/"
source(paste(rerfPath, "RerF/Code/Classifiers/rfr_function.R", sep = ""))
source(paste(rerfPath, "RerF/Code/Classifiers/randmat.R", sep = ""))
source(paste(rerfPath, "RerF/Code/Classifiers/rerf_train.R", sep = ""))
source(paste(rerfPath, "RerF/Code/Classifiers/select_model.R", sep = ""))
source(paste(rerfPath, "RerF/Code/Classifiers/classprob.R", sep = ""))

# initialize arrays
ps <- c(2,4,6)
ns <- list(c(20,200,400), c(80,400,4000), c(400,2000,4000))
nTrials <- 10
testError <- array(vector(mode="double", length=length(ps)*length(ns[[1]])*nTrials), c(length(ns[[1]]),length(ps),nTrials))
trainTime <- array(vector(mode="double", length=length(ps)*length(ns[[1]])*nTrials), c(length(ns[[1]]),length(ps),nTrials))
best.idx <- array(vector(mode="double", length=length(ps)*length(ns[[1]])*nTrials), c(length(ns[[1]]),length(ps),nTrials))
labels <- c("neg","pos")

params <- list()
params$randomMatrix <- "binary"
params$MinParent <- 2
params$MaxDepth <- "inf"
params$bagging <- 0.2
params$COOB <- T
params$rescale <- "none"
params$replacement <- T
params$stratify <- T
params$NumCores <- 12

# loop over number of dimensions
for (j in 2) {
# for (j in 1:length(ps)) {
  p <- ps[j]
  print(paste("p = ", p, sep = ""))
  
  # read in test data for p dimensions
  D <- read.table(paste(dataPath, "Test/Orthant_test_p", p, ".dat", sep = ""), sep = ",")
  Xtest <- data.matrix(D[,1:ncol(D)-1])
  Ytest <- data.matrix(D[,ncol(D)])
  labels <- sort(unique(Ytest))
  ntest = nrow(Xtest)
  
  if (p <= 5) {
    params$d <- c(1:p, p^c(2, 3))
  } else {
    params$d <- ceiling(p^c(1/4, 1/2, 3/4, 1:3))
  }
  params$sparsity <- 1/p
  
  # loop over number of train observations
  for (i in 1:length(ns[[j]])) {
    ntrain <- ns[[j]][i]
    print(paste("n = ", ntrain, sep = ""))
    
    if (ntrain <= 1000) {
      params$trees <- 1000
    } else {
      params$trees <- 500
    }
    
    # loop over trials
    for (trial in 1:nTrials) {
      print(paste("trial = ",trial, sep = ""))
      
      # read in the nth trial of training data for p dimensions 
      D <- read.table(paste(dataPath, "Train/Orthant_train_p", ps[j], "_n", ntrain, "_trial", trial, ".dat", sep = ""), sep = ",")
      Xtrain <- data.matrix(D[,1:ncol(D)-1])
      Ytrain <- data.matrix(D[,ncol(D)])
      ntrain <- nrow(Xtrain)
      
      # train forests
      forests <- rerf_train(Xtrain, Ytrain, params)
      
      # select best forest
      forests.summary <- select_model(forests$forests, Xtrain, Ytrain)
      best.idx[i,j,trial] <- forests.summary$best.idx
      trainTime[i,j,trial] <- forests$trainTime[best.idx[i,j,trial]]
      
      # compute test error
      testError[i,j,trial] <- error_rate(Xtest, Ytest, forests$forests[[forests.summary$best.idx]])
      save(testError, trainTime, best.idx, params, file = paste(rerfPath, "RerF/Results/2017.06.19/Orthant/Orthant_rerf_2017_06_19.RData", sep = ""))
    }
  }
}
