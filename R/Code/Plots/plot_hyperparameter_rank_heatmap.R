rm(list = ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(xtable)
library(gridExtra)
source("~/RandomerForest/R/Code/Utils/GetCatMap.R")
source("~/RandomerForest/R/Code/Utils/GetFolds.R")

# filePath <- "~/RandomerForest/R/Results/2017.11.28/"
# filePath2 <- "~/RandomerForest/R/Results/2018.01.31/"
filePath <- "~/RandomerForest/R/Results/2018.02.07/"
filePathCCF <- "~/RandomerForest/Results/2018.06.12/"
# filePath2 <- "~/RandomerForest/R/Results/2018.02.11/"
contents <- list.files(filePath)
contents <- contents[!grepl("frc", contents) & !grepl("rerfc", contents) & !grepl("rr-rf", contents) & !grepl("rr_rf", contents) & !grepl("xgb", contents)]
catfiles <- list.files("~/tmp/uci/processed/categorical_map/")
load(paste0(filePath, contents[[1L]]))
fieldNames <- ls()
fieldNames <- fieldNames[(fieldNames != "filePath") & (fieldNames != "contents") & (fieldNames != "multiplot") &
                           (fieldNames != "catfiles") & (fieldNames != "GetCatMap") & (fieldNames != "filePathCCF") &
                           (fieldNames != "GetFolds")]

res <- vector("list", length(fieldNames))
names(res) <- fieldNames

RankMatrix <-
  function(X, na.last = T, ties.method = "average") {
    if (is.matrix(X)) {
      X.rank <- apply(X, 2, FUN = function(x) rank(x, na.last = na.last, ties.method = ties.method))
    } else {
      X.rank <- rank(X, na.last = na.last, ties.method = ties.method)
    }
    return(X.rank)
  }

line.width <- 0.75
marker.size <- 1
color.map <- brewer.pal(5L, "Set2")

for (f in contents) {
  ds <- strsplit(f, "_2018_02_07.RData")[[1L]][1L]
  load(paste0(filePath, f))
  for (fn in fieldNames) {
    res[[fn]][[strsplit(f, "_2018_02_07.RData")[[1L]]]] <- get(fn)[[1L]]
  }
}

error.rates <- NULL
dataSets <- names(res$testError)
best <- list(d=NULL, sparsity=NULL)
rm.idx <- NULL
ct <- 0L
data.set <- NULL
k <- NULL
d <- NULL
sparsity <- NULL
rnk <- NULL
for (i in seq_along(dataSets)) {
  ct <- ct + 1L
  ds <- dataSets[i]
  if ((length(res$params[[ds]]$rerf$d) == 5L) && (length(res$params[[ds]]$rerf$sparsity) == 5L)) {
    p <- res$params[[ds]]$rerf$d[4L]
    D <- as.matrix(read.table(paste0("~/tmp/uci/processed/data/", ds, ".csv"), header = F, sep = ",", quote = "", row.names = NULL))
    Y <- as.integer(D[, ncol(D)])
    fold <- GetFolds(paste0("~/tmp/uci/processed/cv_partitions/", ds, "_partitions.txt"))
    chance.error <- sapply(1:5, function(k) 1 - max(tabulate(Y[fold[[k]]]))/length(Y[fold[[k]]]))
    # ex <- c(1/4, 1/2, 3/4, 1, 2)
    data.set <- c(data.set, rep(ds, 125))
    k <- c(k, c(apply(as.matrix(1:5), 1, function(x) rep(x,25))))
    d <- c(d, rep(1:5, 5))
    sparsity <- c(sparsity, c(apply(as.matrix(1:5), 1, function(x) rep(x,25))))
    rnk <- c(rnk, c(RankMatrix(t(res$testError[[ds]]$rerf))))
    # mn <- apply(res$testError[[ds]]$rerf, 1L, min)
    # mn.idx <- unlist(sapply(1:5, function(fold) sample(which(res$testError[[ds]]$rerf[fold, ] == mn[fold]), 1L, replace=FALSE)))
    # mn.idx <- unlist(sapply(1:5, function(fold) which(res$testError[[ds]]$rerf[fold, ] == mn[fold])))
    # best.params <- as.matrix(expand.grid(1:5, round(res$params[[ds]]$rerf$sparsity*p)))[mn.idx, , drop=FALSE]
    # best$d <- c(best$d, best.params[, 1L])
    # best$sparsity <- c(best$sparsity, best.params[, 2L])
    # error.rates <- cbind(error.rates, apply(apply(res$testError[[ds]]$rerf, 2, function(a) a/chance.error), 2, mean))
    # print(length(res$params[[ds]]$rerf$sparsity))
    # print(p)
  }
}

df <- data.frame(data.set=data.set, fold=k, d=d, sparsity=sparsity, rank=rnk)
df <- aggregate(rank ~ data.set+d+sparsity, df[-2], mean)
df <- aggregate(rank ~ d+sparsity, df[-1], median)

# best <- as.data.frame(best)

# p <- ggplot(best, aes(x=d, y=sparsity)) +
#   geom_count(aes(color=..n.., size=..n..)) +
#   # xlab("(# of projections, avg. nnzs per projection)") +
#   scale_x_continuous(breaks = 1:5, labels = c(bquote(p^{1/4}),bquote(p^{1/2}),bquote(p^{3/4}),bquote(p^1),bquote(p^2))) +
#   xlab("# of projections sampled") +
#   ylab("avg. # of nonzeros per projection") +
#   # ylab("frequency") +
#   ggtitle("Frequency of Best Hyperparameter Pairs") +
#   guides(color="legend") +
#   labs(color="frequency",size="frequency") +
#   theme(panel.grid.minor=element_blank())
# 
# ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/hyperparameter_frequency_scatter_plot.png", plot = p, width = 4, height = 4, units = "in")
# ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/hyperparameter_frequency_scatter_plot.pdf", plot = p, width = 4, height = 4, units = "in")

# df <- cbind(data.frame(err=apply(error.rates, 1L, mean), sem=apply(error.rates, 1L, function(a) sd(a)/sqrt(length(a)))), expand.grid(d=1:5, sparsity=1:5))
# df$label <- paste0(round(df$err,3), "\u00B1\n", round(df$sem,3))

p2 <- ggplot(df, aes(x=d, y=sparsity)) +
  geom_tile(aes(fill=rank)) +
  geom_text(aes(label=rank), color="white", size=3) +
  scale_x_continuous(breaks = 1:5, labels = c(bquote(p^{1/4}),bquote(p^{1/2}),bquote(p^{3/4}),bquote(p^1),bquote(p^2))) +
  xlab("No. of Projections Sampled") +
  ylab("Avg. No. of Nonzeros per Projection") +
  ggtitle("Median Rank of Hyperparameter Pairs") +
  labs(fill="Rank") +
  scale_fill_gradient(high="#132B43", low="#56B1F7") +
  guides(fill = guide_colorbar(reverse=TRUE))

ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/hyperparameter_rank_heatmap.png", plot = p2, width = 4, height = 4, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/hyperparameter_rank_heatmap.pdf", plot = p2, width = 4, height = 4, units = "in")
