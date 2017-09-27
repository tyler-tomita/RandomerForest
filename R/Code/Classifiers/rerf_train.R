rerf_train <- function(Xtrain, Ytrain, params = list(trees = 500L, randomMatrix = "binary", d = round(sqrt(ncol(Xtrain))), sparsity = 1/ncol(Xtrain), rotate = F, rank.transform = F, MinParent = 2, MaxDepth = "inf", bagging = 1/exp(1), COOB = T, rescale = "none", replacement = T, stratify = T, NumCores = 1, seed = 1L)) {
  
  p <- ncol(Xtrain)
  
  params.names <- names(params)
  
  if (!("trees" %in% params.names)) {
    params$trees <- 500L
  }
  
  if (!("randomMatrix" %in% params.names)) {
    params$randomMatrix <- "binary"
  }
  
  if (!("d" %in% params.names)) {
    params$d <- round(sqrt(p))
  }
  
  if (!("sparsity" %in% params.names)) {
    if (params$randomMatrix == "binary") {
      params$sparsity <- 1/p
    } else if (params$randomMatrix == "frc") {
      params$sparsity <- 2
    } else if (params$randomMatrix == "poisson") {
      params$sparsity <- 1
    } else if (params$randomMatrix == "rf") {
      params$sparsity <- 1
    }
  }
  
  if (!("rotate" %in% params.names)) {
    params$rotate <- F
  }
  
  if (!("rank.transform" %in% params.names)) {
    params$rank.transform <- F
  }
  
  if (!("MinParent" %in% params.names)) {
    params$MinParent <- 2
  }
  
  if (!("MaxDepth" %in% params.names)) {
    params$MaxDepth <- "inf"
  }
  
  if (!("bagging" %in% params.names)) {
    params$bagging <- 1/exp(1)
  }
  
  if (!("COOB" %in% params.names)) {
    params$COOB <- T
  }
  
  if (!("rescale" %in% params.names)) {
    params$rescale <- "none"
  }
  
  if (!("replacement" %in% params.names)) {
    params$replacement <- T
  }
  
  if (!("stratify" %in% params.names)) {
    params$stratify <- T
  }
  
  if (!("NumCores" %in% params.names)) {
    params$NumCores <- 1
  }
  
  if (!("seed" %in% params.names)) {
    params$seed <- 1L
  }
  
  # compile before the actual run so that training time of first iteration is consistent with the rest
  if (require(compiler)){
    if(!exists("comp_rfr")){
      setCompilerOptions("optimize"=3)
      comp_rfr <<- cmpfun(runrfr)
    }
  }
  
  if (params$randomMatrix == "binary" || params$randomMatrix == "continuous" || params$randomMatrix == "poisson" || params$randomMatrix == "frc") {
    nforest <- length(params$d)*length(params$sparsity)
    forests <- vector(mode = "list", nforest)
    trainTime <- vector(mode = "numeric", length = nforest)
    for (i in 1:length(params$sparsity)) {
      for (j in 1:length(params$d)) {
        options <- list(p, params$d[j], params$randomMatrix, params$sparsity[i])
        forest.idx <- (i - 1)*length(params$d) + j
        print(paste("Training forest ", as.character(forest.idx), " of ", as.character(nforest), sep = ""))
        start.time <- proc.time()
        forests[[forest.idx]] <- rfr(Xtrain, Ytrain, trees = params$trees, FUN = randmat, options = options, rotate = params$rotate, rank.transform = params$rank.transform, MinParent = params$MinParent, MaxDepth = params$MaxDepth, bagging = params$bagging, COOB = params$COOB, replacement = params$replacement, stratify = params$stratify, NumCores = params$NumCores, seed = seed)
        trainTime[forest.idx] <- (proc.time() - start.time)[[3]]
      }
    }
  } else {
    params$d <- params$d[params$d <= p]
    nforest <- length(params$d)
    forests <- vector(mode = "list", nforest)
    trainTime <- vector(mode = "numeric", length = nforest)
    for (j in 1:length(params$d)) {
      options <- list(p, params$d[j], params$randomMatrix, NULL)
      print(paste("Training forest ", as.character(j), " of ", as.character(nforest), sep = ""))
      start.time <- proc.time()
      forests[[j]] <- rfr(Xtrain, Ytrain, trees = params$trees, FUN = randmat, options = options, rotate = params$rotate, rank.transform = params$rank.transform, MinParent = params$MinParent, MaxDepth = params$MaxDepth, bagging = params$bagging, COOB = params$COOB, replacement = params$replacement, stratify = params$stratify, NumCores = params$NumCores, seed = seed)
      trainTime[j] <- (proc.time() - start.time)[[3]]
    }
  } 
  return(list(forests = forests, trainTime = trainTime, params = params))
}
