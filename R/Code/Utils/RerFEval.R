#' Rerf Evaluate 
#'
#' Evaluate Rerf models for several values of mat.options[[2]] and mat.options[[4]]
#'
#' @param Xtrain an n sample by d feature matrix (preferable) or data frame that will be used to train a forest.
#' @param Ytrain an n length vector of class labels
#' @param Xtest an n sample by d feature matrix (preferable) or data frame that will be used to test a forest.
#' @param Ytest an n length vector of class labels
#' @param params a list of parameters used in the various rerf functions
#' @param store.predictions ????? (store.predictions=FALSE)
#'
#' @return something ?????
#'
#' @author James and Tyler, jbrowne6@jhu.edu and
#' 
#' @importFrom AUC auc roc
#' @importFrom compiler setCompilerOptions cmpfun
#'

RerFEval <-
    function(Xtrain, Ytrain, Xtest, Ytest, params = list(trees = 500L, random.matrix = "binary", d = round(sqrt(ncol(Xtrain))), sparsity = 1/ncol(Xtrain), rotate = F, rank.transform = F, min.parent = 2L, max.depth = 0L, bagging = 1/exp(1), store.oob = T, store.ns = F, replacement = T, stratify = T, num.cores = 1L, seed = 1L, cat.map = NULL), store.predictions = F) {

        params.names <- names(params)
        
        nClasses <- length(unique(Ytrain))
        
        if (!("rotate" %in% params.names)) {
          params$rotate <- F
        }
        
        if (!("cat.map" %in% params.names)) {
          params$cat.map <- NULL
        }
      
        if (is.null(params$cat.map)) {
          p <- ncol(Xtrain)
        } else {
          if (!params$rotate) {
            p <- params$cat.map[[1L]][1L] - 1L + length(params$cat.map)
          } else{
            p <- ncol(Xtrain)
            params$cat.map <- NULL
          }
        }

        if (!("trees" %in% params.names)) {
            params$trees <- 500L
        }

        if (!("random.matrix" %in% params.names)) {
            params$random.matrix <- "binary"
        }

        if (!("d" %in% params.names)) {
            params$d <- round(sqrt(p))
        }

        if (!("sparsity" %in% params.names)) {
            if (params$random.matrix == "binary" || params$random.matrix == "continuous") {
                params$sparsity <- 1/p
            } else if (params$random.matrix == "frc" || params$random.matrix == "frcn") {
                params$sparsity <- 2
            } else if (params$random.matrix == "poisson") {
                params$sparsity <- 1
            } else if (params$random.matrix == "rf") {
                params$sparsity <- 1
            }
        }

        if (!("rank.transform" %in% params.names)) {
            params$rank.transform <- F
        }

        if (!("min.parent" %in% params.names)) {
            params$min.parent <- 2L
        }

        if (!("max.depth" %in% params.names)) {
            params$max.depth <- "inf"
        }

        if (!("bagging" %in% params.names)) {
            params$bagging <- 1/exp(1)
        }

        if (!("store.oob" %in% params.names)) {
            params$store.oob <- T
        }

        if (!("store.ns" %in% params.names)) {
            params$store.ns <- F
        }

        if (!("replacement" %in% params.names)) {
            params$replacement <- T
        }

        if (!("stratify" %in% params.names)) {
            params$stratify <- T
        }

        if (!("num.cores" %in% params.names)) {
            params$num.cores <- 1L
        }

        if (!("seed" %in% params.names)) {
            params$seed <- 1L
        }
        set.seed(params$seed)

        if (params$random.matrix == "binary" || params$random.matrix == "continuous" || params$random.matrix == "poisson" ||
            params$random.matrix == "frc" || params$random.matrix == "frcn") {
            nforest <- length(params$d)*length(params$sparsity)
            trainTime <- vector(mode = "numeric", length = nforest)
            oobTime <- vector(mode = "numeric", length = nforest)
            testTime <- vector(mode = "numeric", length = nforest)
            testError <- vector(mode = "numeric", length = nforest)
            oobError <- vector(mode = "numeric", length = nforest)
            oobAUC <- vector(mode = "numeric", length = nforest)
            treeStrength <- vector(mode = "numeric", length = nforest)
            treeCorrelation <- vector(mode = "numeric", length = nforest)
            numNodes <- vector(mode = "numeric", length = nforest)
            if (store.predictions) {
                Yhat <- matrix(0, nrow = Xtest, ncol = nforest)
            }
            for (i in 1:length(params$sparsity)) {
                for (j in 1:length(params$d)) {
                    mat.options <- list(p, params$d[j], params$random.matrix, params$sparsity[i], params$cat.map)
                    forest.idx <- (i - 1)*length(params$d) + j

                    print(paste("Evaluating forest ", as.character(forest.idx), " of ", as.character(nforest), sep = ""))

                    # train
                    print("training")
                    start.time <- proc.time()
                    forest <- RerF(Xtrain, Ytrain, trees = params$trees, mat.options = mat.options, rank.transform = params$rank.transform,
                                   min.parent = params$min.parent, max.depth = params$max.depth, bagging = params$bagging, store.oob = params$store.oob,
                                   store.ns = params$store.ns, replacement = params$replacement, stratify = params$stratify, num.cores = params$num.cores,
                                   seed = params$seed, cat.map = params$cat.map, rotate = params$rotate)
                    trainTime[forest.idx] <- (proc.time() - start.time)[[3L]]
                    print("training complete")
                    print(paste("elapsed time: ", trainTime[forest.idx], sep = ""))

                    # compute out-of-bag metrics
                    print("computing out-of-bag predictions")
                    start.time <- proc.time()
                    oobScores <<- OOBPredict(Xtrain, forest, num.cores = params$num.cores, output.scores = T)
                    oobTime[forest.idx] <- (proc.time() - start.time)[[3L]]
                    print("out-of-bag predictions complete")
                    print(paste("elapsed time: ", oobTime[forest.idx], sep = ""))
                    oobError[forest.idx] <- mean(forest$labels[max.col(oobScores)] != Ytrain)
                    if (nClasses > 2L) {
                        Ybin <- as.factor(as.vector(dummies::dummy(Ytrain)))
                        oobAUC[forest.idx] <- AUC::auc(AUC::roc(as.vector(oobScores), Ybin))
                    } else {
                        # Ytrain starts from 1, but here we need it to start from 0
                        oobAUC[forest.idx] <- AUC::auc(AUC::roc(oobScores[, 2L], as.factor(Ytrain - 1L)))
                    }

                    numNodes[forest.idx] <- mean(sapply(forest$trees, FUN = function(tree) length(tree$treeMap)))

                    # make predictions on test set
                    print("computing predictions on test set")
                    start.time <- proc.time()
                    if (store.predictions) {
                      Yhat[, forest.idx] <- Predict(Xtest, forest, num.cores = params$num.cores, Xtrain = Xtrain)
                      testTime[forest.idx] <- (proc.time() - start.time)[[3L]]
                      print("test set predictions complete")
                      print(paste("elapsed time: ", testTime[forest.idx], sep = ""))
                      testError[forest.idx] <- mean(Yhat[, forest.idx] != Ytest)
                    } else {
                      Yhat <- Predict(Xtest, forest, num.cores = params$num.cores, Xtrain = Xtrain)
                      testTime[forest.idx] <- (proc.time() - start.time)[[3L]]
                      print("test set predictions complete")
                      print(paste("elapsed time: ", testTime[forest.idx], sep = ""))
                      testError[forest.idx] <- mean(Yhat != Ytest)
                    }

                    # compute strength and correlation
                    print("computing tree strength and correlation")
                    preds <- Predict(Xtest, forest, num.cores = params$num.cores, Xtrain = Xtrain, aggregate.output = F)
                    sc <- StrCorr(preds, Ytest)
                    print("strength and correlation complete")
                    treeStrength[forest.idx] <- sc$s
                    treeCorrelation[forest.idx] <- sc$rho

                    # # save forest models
                    # save(forest, file = fileName)
                }
            }
        } else {
            params$d <- params$d[params$d <= p]
            nforest <- length(params$d)
            trainTime <- vector(mode = "numeric", length = nforest)
            oobTime <- vector(mode = "numeric", length = nforest)
            testTime <- vector(mode = "numeric", length = nforest)
            testError <- vector(mode = "numeric", length = nforest)
            oobError <- vector(mode = "numeric", length = nforest)
            oobAUC <- vector(mode = "numeric", length = nforest)
            treeStrength <- vector(mode = "numeric", length = nforest)
            treeCorrelation <- vector(mode = "numeric", length = nforest)
            numNodes <- vector(mode = "numeric", length = nforest)
            if (store.predictions) {
                Yhat <- matrix(0, nrow = Xtest, ncol = nforest)
            }
            for (forest.idx in 1:nforest) {
                mat.options <- list(p, params$d[forest.idx], params$random.matrix, NULL, params$cat.map)

                print(paste("Evaluating forest ", as.character(forest.idx), " of ", as.character(nforest), sep = ""))

                # train
                print("training")
                start.time <- proc.time()
                forest <- RerF(Xtrain, Ytrain, trees = params$trees, mat.options = mat.options, rank.transform = params$rank.transform,
                               min.parent = params$min.parent, max.depth = params$max.depth, bagging = params$bagging, store.oob = params$store.oob,
                               store.ns = params$store.ns, replacement = params$replacement, stratify = params$stratify, num.cores = params$num.cores,
                               seed = params$seed, cat.map = params$cat.map, rotate = params$rotate)
                trainTime[forest.idx] <- (proc.time() - start.time)[[3L]]
                print("training complete")
                print(paste("elapsed time: ", trainTime[forest.idx], sep = ""))

                # compute out-of-bag metrics
                print("computing out-of-bag predictions")
                start.time <- proc.time()
                oobScores <- OOBPredict(Xtrain, forest, num.cores = params$num.cores, output.scores = T)
                oobTime[forest.idx] <- (proc.time() - start.time)[[3L]]
                print("out-of-bag predictions complete")
                print(paste("elapsed time: ", oobTime[forest.idx], sep = ""))
                oobError[forest.idx] <- mean(forest$labels[max.col(oobScores)] != Ytrain)
                if (nClasses > 2L) {
                    Ybin <- as.factor(as.vector(dummies::dummy(Ytrain)))
                    oobAUC[forest.idx] <- AUC::auc(AUC::roc(as.vector(oobScores), Ybin))
                } else {
                    # Ytrain starts from 1, but here we need it to start from 0
                    oobAUC[forest.idx] <- AUC::auc(AUC::roc(oobScores[, 2L], as.factor(Ytrain - 1L)))
                }

                numNodes[forest.idx] <- mean(sapply(forest$trees, FUN = function(tree) length(tree$treeMap)))

                # make predictions on test set
                print("computing predictions on test set")
                start.time <- proc.time()
                if (store.predictions) {
                  Yhat[, forest.idx] <- Predict(Xtest, forest, num.cores = params$num.cores, Xtrain = Xtrain)
                  testTime[forest.idx] <- (proc.time() - start.time)[[3L]]
                  print("test set predictions complete")
                  print(paste("elapsed time: ", testTime[forest.idx], sep = ""))
                  testError[forest.idx] <- mean(Yhat[, forest.idx] != Ytest)
                } else {
                  Yhat <- Predict(Xtest, forest, num.cores = params$num.cores, Xtrain = Xtrain)
                  testTime[forest.idx] <- (proc.time() - start.time)[[3L]]
                  print("test set predictions complete")
                  print(paste("elapsed time: ", testTime[forest.idx], sep = ""))
                  testError[forest.idx] <- mean(Yhat != Ytest)
                }

                # compute strength and correlation
                print("computing tree strength and correlation")
                start.time <- proc.time()
                preds <- Predict(Xtest, forest, num.cores = params$num.cores, Xtrain = Xtrain, aggregate.output = F)
                sc <- StrCorr(preds, Ytest)
                print("strength and correlation complete")
                print(paste0("elapsed time: ", (proc.time() - start.time)[[3L]]))
                treeStrength[forest.idx] <- sc$s
                treeCorrelation[forest.idx] <- sc$rho
            }
        }

        # select best model
        minError.idx <- which(oobError == min(oobError))
        if (length(minError.idx) > 1L) {
            maxAUC.idx <- which(oobAUC[minError.idx] == max(oobAUC[minError.idx]))
            if (length(maxAUC.idx) > 1L) {
                maxAUC.idx <- sample(maxAUC.idx, 1L)
            }
            best.idx <- minError.idx[maxAUC.idx]  
        } else {
            best.idx <- minError.idx
        }

        if (store.predictions) {
            return(list(Yhat = Yhat[, best.idx], testError = testError, trainTime = trainTime, oobTime = oobTime, testTime = testTime, oobError = oobError, oobAUC = oobAUC, treeStrength = treeStrength, treeCorrelation = treeCorrelation, numNodes = numNodes, best.idx = best.idx, params = params))
        } else {
            return(list(testError = testError, trainTime = trainTime, oobTime = oobTime, testTime = testTime, oobError = oobError, oobAUC = oobAUC, treeStrength = treeStrength, treeCorrelation = treeCorrelation, numNodes = numNodes, best.idx = best.idx, params = params))
        }
    }
