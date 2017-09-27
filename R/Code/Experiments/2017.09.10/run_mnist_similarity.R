rm(list = ls())
# source("~/work/tyler/RerF/Code/Classifiers/rfr_function.R")
source("~/R-RerF/rfr_function.R", chdir = T)
library(reshape2)
library(ggplot2)

# filePath <- "~/work/tyler/Data/tests/"
# outFile <- "~/work/tyler/RerF/Results/Tests/test_results_2017_08_22.RData"
filePath <- "~/tmp/big/processed/"
# outFile <- "~/RerF/Results/Tests/test_results_2017_08_22.RData"

# dataSets <- c("Sparse_parity", "Trunk", "Orthant", "mnist")
dataSets <- "mnist"

trainTime <- list()
OOBError <- list()
numNodes <- list()
testError <- list()
seed <- list()

D <- "mnist"

trainFile <- paste(filePath, D, ".train.csv", sep = "")
testFile <- paste(filePath, D, ".test.csv", sep = "")

Xtrain <- as.matrix(read.table(trainFile, header = F, sep = ",", quote = "", row.names = NULL))
ntrain <- nrow(Xtrain)
p <- ncol(Xtrain) - 1L
Ytrain <- as.integer(Xtrain[, p + 1L]) + 1L
Xtrain <- Xtrain[, -(p + 1L)]

Xtest <- as.matrix(read.table(testFile, header = F, sep = ",", quote = "", row.names = NULL))
ntest <- nrow(Xtest)
Ytest <- as.integer(Xtest[, p + 1L]) + 1L
Xtest <- Xtest[, -(p + 1L)]

nClasses <- length(unique(Ytrain))

mtry <- ceiling(sqrt(p))
nTrees <- 500L
stratify <- T
replacement <- T
randomMatrix <- "binary"
sparsity <- 1/p
rotate <- F
MinParent <- 20L
MaxDepth <- "inf"
COOB <- T
CNS <- F
NumCores <- 1L
comp.mode <- "batch"
rank.transform <- F

# train the classifier
cat("training\n")
start.time <- proc.time()
forest <- rerf(Xtrain, Ytrain, nClasses, MinParent = MinParent, trees = nTrees, MaxDepth = MaxDepth, replacement = replacement,
               stratify = stratify, FUN = randmat, options = list(p, mtry, randomMatrix, sparsity), COOB = COOB,
               CNS = CNS, NumCores = NumCores, seed = 123L, rotate = rotate, rank.transform = rank.transform)
cat("training complete\n")

# sortIdx <- order(Ytest)
sampleIdx <- integer(100L*10L)
for (i in seq.int(10L)) {
  sampleIdx[((i-1L)*100L + 1L):(i*100L)] <- sample(which(Ytest == i), 100L)
}

# compute similarity
similarity <- compute.similarity(Forest = forest, X = Xtest[sampleIdx, ], NumCores = 4L)

fmriu.plot.plot_graph <- function(mtx, title="",xlabel="ROI", ylabel="ROI", legend.name="metric", legend.show=TRUE, itype="sq",
                                 font.size=12, rem_diag=FALSE, include_diag=FALSE, limits=c(0, 1)) {
  if (itype == "ts") {
    mtx <- abs(cor(mtx))  # if a timeseries is passed in, correlate the features first
  }
  if (!include_diag) {
    diag(mtx) <- 0
  }
  dm <- melt(mtx)
  colnames(dm) <- c("x", "y", "value")
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  sqplot <- ggplot(dm, aes(x=x, y=y, fill=value)) +
    geom_tile() +
    scale_fill_gradientn(colours=jet.colors(7), name=legend.name, limits=limits) +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(title)
  if (legend.show) {
    sqplot <- sqplot +
      theme(text=element_text(size=font.size))
  } else {
    sqplot <- sqplot +
      theme(text=element_text(size=font.size, legend.position="none"))
  }
  return(sqplot)
}

fmriu.plot.plot_graph(mtx = similarity, xlabel = "", ylabel = "", title = "MNIST Test Set", include_diag = T, legend.name = "RerF\nsimilarity")
