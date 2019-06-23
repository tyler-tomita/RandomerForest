rm(list = ls())
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
source("~/RandomerForest/R/Code/Utils/GetCatMap.R")

filePath <- "~/RandomerForest/R/Results/2017.10.03/"

contents <- list.files(filePath)
contents <- contents[!grepl("frc", contents) & !grepl("rerfc", contents) & !grepl("rr-rf", contents)]
num.datasets <- length(contents)
catfiles <- list.files("~/tmp/uci/processed/categorical_map/")
load(paste0(filePath, contents[[1L]]))
fieldNames <- ls()
fieldNames <- fieldNames[(fieldNames != "filePath") & (fieldNames != "contents") & (fieldNames != "multiplot") &
                           (fieldNames != "catfiles") & (fieldNames != "GetCatMap")]
res <- vector("list", length(fieldNames))
names(res) <- fieldNames
classifiers <- vector("list", length(contents))

for (f in contents) {
  load(paste0(filePath, f))
  classifiers[[strsplit(f, "_2017_10_03.RData")[[1L]]]] <- names(get(fieldNames[[1L]])[[1L]])
  for (fn in fieldNames) {
    res[[fn]][[strsplit(f, "_2017_10_03.RData")[[1L]]]] <- get(fn)[[1L]]
  }
}

num.folds <- length(res$bestIdx[[1L]][[1L]])

# classifiers <- unique(unlist(classifiers))
classifiers <- c("rerfr", "rerfc", "rf", "frc")
num.classifiers <- length(classifiers)
class.map <- list(rerfr = "RerF(r)", rerfc = "RerF-C", rf = "RF", frc = "F-RC")
class.map2 <- class.map[names(class.map) != "rf"]

color.map <- c("#41ab5d", "#4292c6", "#f768a1")
line.width <- 1.5
marker.size <- 2

dataSets <- names(res$testError)
data.set <- vector("character", num.classifiers*num.datasets*num.folds)
classifier <- vector("character", num.classifiers*num.datasets*num.folds)
num.obs <- vector("integer", num.classifiers*num.datasets*num.folds)
num.dims <- vector("integer", num.classifiers*num.datasets*num.folds)
num.cat.dims <- vector("integer", num.classifiers*num.datasets*num.folds)
test.error <- vector("double", num.classifiers*num.datasets*num.folds)
chance.error <- vector("double", num.classifiers*num.datasets*num.folds)
train.time <- vector("double", num.classifiers*num.datasets*num.folds)
tree.strength <- vector("double", num.classifiers*num.datasets*num.folds)
tree.correlation <- vector("double", num.classifiers*num.datasets*num.folds)
sensitivity <- vector("double", num.classifiers*num.datasets*num.folds)
for (i in 1:num.classifiers) {
  cl <- classifiers[i]
  classifier[((i - 1L)*num.datasets*num.folds + 1L):(i*num.datasets*num.folds)] <- class.map[[cl]]
  for (j in 1:num.datasets) {
    ds <- dataSets[j]
    D <- as.matrix(read.table(paste0("~/tmp/uci/processed/data/", ds, ".csv"), header = F, sep = ",", quote = "", row.names = NULL))
    Y <- as.integer(D[, ncol(D)])
    data.set[((i - 1L)*num.datasets*num.folds + (j - 1L)*num.folds + 1L):(((i - 1L)*num.datasets + j)*num.folds)] <- ds
    num.obs[((i - 1L)*num.datasets*num.folds + (j - 1L)*num.folds + 1L):(((i - 1L)*num.datasets + j)*num.folds)] <- length(Y)
    if (paste0(ds, "_catmap.txt") %in% catfiles) {
      cat.map <- GetCatMap(paste0("~/tmp/uci/processed/categorical_map/", ds, "_catmap.txt"))
      num.cat.dims[((i - 1L)*num.datasets*num.folds + (j - 1L)*num.folds + 1L):(((i - 1L)*num.datasets + j)*num.folds)] <- length(cat.map)
      num.dims[((i - 1L)*num.datasets*num.folds + (j - 1L)*num.folds + 1L):(((i - 1L)*num.datasets + j)*num.folds)] <- cat.map[[1L]][1L] - 1L + length(cat.map)
    } else {
      num.dims[((i - 1L)*num.datasets*num.folds + (j - 1L)*num.folds + 1L):(((i - 1L)*num.datasets + j)*num.folds)] <- ncol(D) - 1L
    }
    chance.error[((i - 1L)*num.datasets*num.folds + (j - 1L)*num.folds + 1L):(((i - 1L)*num.datasets + j)*num.folds)] <- 1 - max(tabulate(Y))/length(Y)

    length.d <- length(res$params[[ds]][[cl]]$d)
    length.sparsity <- length(res$params[[ds]][[cl]]$sparsity)
    if ((cl == "rerfr") || (cl == "rerfc")) {
      if (length.sparsity > length(res$params[[ds]][["frc"]]$sparsity)) {
        ind <- 2:length.sparsity
        min.idx <- matrix(0, num.folds, length.sparsity - 1L)
      } else {
        ind <- 1:length.sparsity
        min.idx <- matrix(0, num.folds, length.sparsity)
      }
    } else {
      ind <- 1:length.sparsity
      min.idx <- matrix(0, num.folds, length.sparsity)
    }
    for (fold in 1:num.folds) {
      for (sp.idx in ind) {
        best.oob.error <- which(res$OOBError[[ds]][[cl]][fold, ((sp.idx - 1L)*length.d + 1L):(sp.idx*length.d)] == min(res$OOBError[[ds]][[cl]][fold, ((sp.idx - 1L)*length.d + 1L):(sp.idx*length.d)]))
        if (length(best.oob.error) > 1L) {
          best.oob.auc <- which(res$OOBAUC[[ds]][[cl]][fold, (sp.idx - 1L)*length.d + best.oob.error] == max(res$OOBAUC[[ds]][[cl]][fold, (sp.idx - 1L)*length.d + best.oob.error]))
          if (length(best.oob.auc) > 1L) {
            best.oob.auc <- sample(best.oob.auc, 1L)
          }
          best.oob.error <- best.oob.error[best.oob.auc]
        }
        if (cl == "frc") {
          min.idx[fold, sp.idx] <- (sp.idx - 1L)*length.d + best.oob.error
        } else {
          min.idx[fold, sp.idx - 1L] <- (sp.idx - 1L)*length.d + best.oob.error
        }
      }
    }
    test.error[((i - 1L)*num.datasets*num.folds + (j - 1L)*num.folds + 1L):(((i - 1L)*num.datasets + j)*num.folds)] <- sapply(1:num.folds, function(x) res$testError[[ds]][[cl]][x, res$bestIdx[[ds]][[cl]][x]])
    train.time[((i - 1L)*num.datasets*num.folds + (j - 1L)*num.folds + 1L):(((i - 1L)*num.datasets + j)*num.folds)] <- sapply(1:num.folds, function(x) res$trainTime[[ds]][[cl]][x, res$bestIdx[[ds]][[cl]][x]])
    tree.strength[((i - 1L)*num.datasets*num.folds + (j - 1L)*num.folds + 1L):(((i - 1L)*num.datasets + j)*num.folds)] <- sapply(1:num.folds, function(x) res$treeStrength[[ds]][[cl]][x, res$bestIdx[[ds]][[cl]][x]])
    tree.correlation[((i - 1L)*num.datasets*num.folds + (j - 1L)*num.folds + 1L):(((i - 1L)*num.datasets + j)*num.folds)] <- sapply(1:num.folds, function(x) res$treeCorr[[ds]][[cl]][x, res$bestIdx[[ds]][[cl]][x]])
    sensitivity[((i - 1L)*num.datasets*num.folds + (j - 1L)*num.folds + 1L):(((i - 1L)*num.datasets + j)*num.folds)] <- sapply(1:num.folds, function(x) sd(res$testError[[ds]][[cl]][x, min.idx[x, ]]))
  }
}

df <- data.frame(classifier = factor(classifier, levels = sapply(names(class.map), function(x) class.map[[x]])))
df$dataset <- data.set
df$p <- num.dims
df$n <- num.obs
df$error <- test.error
df$chance <- chance.error
df$train.time <- train.time
df$tree.strength <- tree.strength
df$tree.correlation <- tree.correlation
df$sensitivity <- sensitivity

# df2 <- data.frame(sensitivity.frc = df$sensitivity[classifier == "F-RC"], sensitivity.rerf = df$sensitivity[classifier == "RerF-C"])
# 
# p <- ggplot(df2, aes(x = sensitivity.rerf^(1/2), y = sensitivity.frc^(1/2))) +
#   geom_point(size = marker.size, color = "#4292c6") +
#   geom_abline(slope = 1, intercept = 0) +
#   coord_equal(xlim = c(0, 0.7), ylim = c(0, 0.7)) +
#   xlab("sparsity sensitivity (RerF)") +
#   ylab("sparsity sensitivity (F-RC)") +
#   ggtitle(paste0(num.datasets, " UCI Benchmark Datasets")) +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.title = element_text(face = "bold"))

df3 <- data.frame(sensitivity = (df$sensitivity[classifier == "F-RC"] - df$sensitivity[classifier == "RerF(r)"])/df$chance[classifier == "F-RC"])

p2 <- ggplot(df3, aes(x = sensitivity)) +
  geom_vline(xintercept = 0, size = 1, linetype = "dashed") +
  stat_density(adjust = 2, geom = "line", size = line.width, color = "#4292c6") +
  xlab("sensitivity to sparsity (F-RC - RerF(r))") +
  ylab("kernel density estimate") +
  ggtitle(paste0(num.datasets, " UCI Benchmark Datasets")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.title = element_text(face = "bold"))
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/sparsity_sensitivity_benchmarks.pdf", plot = p2, width = 11, height = 5, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/sparsity_sensitivity_benchmarks.png", plot = p2, width = 11, height = 5, units = "in")
