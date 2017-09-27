select_model <- function(models, Xtrain, Ytrain, NumCores = 0L, seed = 1L) {
  library(dummies)
  library(AUC)
  set.seed(seed)
  
  forests <- models$forests
  
  ntrain <- length(Ytrain)
  nClasses <- length(unique(Ytrain))
  
  oobError <- vector(mode = "numeric", length = length(forests))
  oobAUC <- vector(mode = "numeric", length = length(forests))
  
  for (forest.idx in 1:length(forests)) {
    print(paste("evaluating forest ", forest.idx, sep = ""))
    nTrees <- length(forests[[forest.idx]])
    oobCounts <- vector(mode = "integer", length = ntrain)
    for (tree.idx in 1:nTrees) {
      # oob.idx <- forests[[forest.idx]]$OOBmat[[tree.idx]][,1]
      oobCounts[forests[[forest.idx]][[tree.idx]]$ind] <- oobCounts[forests[[forest.idx]][[tree.idx]]$ind] + 1
    }
    # oobPredictions <- max.col(oobScores)
    oobmat <- OOBpredict(Xtrain, Ytrain, forests[[forest.idx]], NumCores = NumCores, rank.transform = models$params$rank.transform)
    oobScores <- oobmat[[nTrees + 1]][, -c(1,2)]/oobCounts
    oobPredictions <- max.col(oobScores)
    oobError[forest.idx] <- sum(oobPredictions != Ytrain)/ntrain
    if (nClasses > 2) {
      Ybin <- as.factor(as.vector(dummy(Ytrain)))
      oobAUC[forest.idx] <- auc(roc(as.vector(oobScores), Ybin))
    } else {
      # Ytrain starts from 1, but here we need it to start from 0
      oobAUC[forest.idx] <- auc(roc(oobScores[,nClasses], as.factor(Ytrain - 1)))
    }
  }
  minError.idx <- which(oobError == min(oobError))
  if (length(minError.idx) > 1) {
    maxAUC.idx <- which(oobAUC[minError.idx] == max(oobAUC[minError.idx]))
    if (length(maxAUC.idx) > 1) {
      maxAUC.idx <- sample(maxAUC.idx, 1)
    }
    best.idx <- minError.idx[maxAUC.idx]  
  } else {
    best.idx <- minError.idx
  }
  return(list(oobError = oobError, oobAUC = oobAUC, best.idx = best.idx))
}