# evaluate xgboost on benchmark datasets

rm(list=ls())
options(scipen = 999)

rerfPath <- "~/work/tyler/"
dataPath <- "~/work/tyler/Data/uci/processed/"
dataPath2 <- "~/work/tyler/Data/uci/gaussian_noise_p10.csv"
library(xgboost)
library(caret)
library(plyr)
source(paste0(rerfPath, "RandomerForest/R/Code/Utils/GetFolds.R"))

testError <- list()
colSample <- list()

dataSet <- "abalone"
fold <- GetFolds(paste0(dataPath, "cv_partitions/", dataSet, "_partitions.txt"))
nFolds <- length(fold)
X <- as.matrix(read.table(paste0(dataPath, "data/", dataSet, ".csv"), header = F, sep = ",", quote = "", row.names = NULL))
catMap <- NULL
p <- ncol(X) - 1L
Y <- as.integer(X[, p + 1L]) + 1L
Y <- paste0("Y",as.character(Y))
X <- X[, -(p + 1L)]
# remove columns with zero variance
X <- X[, apply(X, 2, function(x) any(as.logical(diff(x))))]
# mean-center and scale by sd
X <- scale(X)
n <- nrow(X)

# read in noisy features and append to X
Xnoise <- as.matrix(read.table(datapath2, header = F, sep = ",", quote = "", row.names = NULL, nrows = n))
X <- cbind(X, Xnoise)
p <- ncol(X)

# coerce to numeric so that xgboost will run properly
X[] <- apply(X, 2, as.numeric)

testError[[dataSet]] <- numeric(nFolds)
colSample[[dataSet]] <- integer(nFolds)

set.seed(02072018L)

# loop over folds
for (k in seq.int(nFolds)) {
  print(paste0("fold ", k))

  trainIdx <- unlist(fold[-k])
  testIdx <- fold[[k]]

  # evaluate models

  if (length(trainIdx) < 100) {
    subSample = 1
  } else {
    subSample = c(0.5,0.75,1)
  }
  xgb_grid <- expand.grid(
    nrounds = c(100, 1000),
    eta = c(0.001, 0.01),
    subsample = subSample,
    colsample_bytree = c(0.2, 0.4, 0.6, 0.8, 1),
    min_child_weight = 1,
    max_depth = c(4, 6, 8, 10, 100000),
    gamma = 0
  )

  nClasses <- length(unique(Y[trainIdx]))

  if (nClasses > 2) {
    # pack the training control parameters
    xgb_trcontrol <- trainControl(
      method = "cv",
      number = 5,
      verboseIter = FALSE,
      returnData = FALSE,
      returnResamp = "all",                                                        # save losses across all models
      classProbs = TRUE,                                                           # set to TRUE for AUC to be computed
      allowParallel = TRUE
    )
    obj <- "multi:softprob"
    Met <- "Accuracy"
    # train the model for each parameter combination in the grid, using CV to evaluate
    xgb_train <- train(
      x = X[trainIdx, ],
      y = as.factor(Y[trainIdx]),
      trControl = xgb_trcontrol,
      tuneGrid = xgb_grid,
      method = "xgbTree",
      objective = obj,
      num_class = nClasses,
      metric = Met,
      nthread = 24
    )
  } else {
    # pack the training control parameters
    xgb_trcontrol <- trainControl(
      method = "cv",
      number = 5,
      verboseIter = FALSE,
      returnData = FALSE,
      returnResamp = "all",                                                        # save losses across all models
      classProbs = TRUE,                                                           # set to TRUE for AUC to be computed
      summaryFunction = twoClassSummary,
      allowParallel = TRUE
    )
    obj <- "binary:logistic"
    Met <- "ROC"
    # train the model for each parameter combination in the grid, using CV to evaluate
    xgb_train <- train(
      x = X[trainIdx, ],
      y = as.factor(Y[trainIdx]),
      trControl = xgb_trcontrol,
      tuneGrid = xgb_grid,
      method = "xgbTree",
      objective = obj,
      metric = Met,
      nthread = 24
    )
  }


  scores <- predict(xgb_train$finalModel,X[testIdx, ])
  if(nClasses > 2) {
    scores <- matrix(scores, nrow = length(testIdx), ncol = nClasses, byrow = TRUE)
  } else {
    scores <- cbind(scores, 1 - scores)
  }
  predictions <- xgb_train$levels[max.col(scores)]
  testError[[dataSet]][k] <- sum(predictions != Y[testIdx])/length(testIdx)
  colSample[[dataSet]][k] <- xgb_train$finalModel$tuneValue$colsample_bytree
  save(testError, colSample, file = paste0(rerfPath, "RandomerForest/R/Results/2018.05.29/", dataSet, "_noise_D10_xgb_2019_05_29.RData"))
}
