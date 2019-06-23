rm(list = ls())
library(reshape2)
library(ggplot2)
library(ks)
library(meda)
source("~/RandomerForest/R/Code/Utils/multiplot.R")
source("~/RandomerForest/R/Code/Utils/GetCatMap.R")

filePath <- "~/RandomerForest/R/Results/2017.09.23/"

contents <- list.files(filePath)
contents <- contents[!grepl("frc", contents)]
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
  classifiers[[strsplit(f, "_2017_09_23.RData")[[1L]]]] <- names(get(fieldNames[[1L]])[[1L]])
  for (fn in fieldNames) {
    res[[fn]][[strsplit(f, "_2017_09_23.RData")[[1L]]]] <- get(fn)[[1L]]
  }
}

classifiers <- unique(unlist(classifiers))

dataSets <- names(res$testError)
mean.error <- matrix(NA, nrow = length(res$testError), ncol = length(classifiers))
rownames(mean.error) <- dataSets
colnames(mean.error) <- classifiers
train.time <- matrix(NA, nrow = length(res$testError), ncol = length(classifiers))
row.names(train.time) <- dataSets
colnames(train.time) <- classifiers
tree.strength <- matrix(NA, nrow = length(res$testError), ncol = length(classifiers))
row.names(tree.strength) <- dataSets
colnames(tree.strength) <- classifiers
tree.correlation <- matrix(NA, nrow = length(res$testError), ncol = length(classifiers))
row.names(tree.correlation) <- dataSets
colnames(tree.correlation) <- classifiers


chance.error <- double(length(dataSets))
n <- integer(length(dataSets))
p <- integer(length(dataSets))
pcat <- integer(length(dataSets))
for (i in seq.int(length(dataSets))) {
  ds <- dataSets[i]
  D <- as.matrix(read.table(paste0("~/tmp/uci/processed/data/", ds, ".csv"), header = F, sep = ",", quote = "", row.names = NULL))
  Y <- as.integer(D[, ncol(D)])
  n[i] <- length(Y)
  if (paste0(ds, "_catmap.txt") %in% catfiles) {
    cat.map <- GetCatMap(paste0("~/tmp/uci/processed/categorical_map/", ds, "_catmap.txt"))
    pcat[i] <- length(cat.map)
    p[i] <- cat.map[[1L]][1L] - 1L + pcat[i]
  } else {
    p[i] <- ncol(D) - 1L
  }
  chance.error[i] <- 1 - max(tabulate(Y))/length(Y)
  clNames <- names(res$testError[[ds]])
  for (j in seq.int(length(clNames))) {
    cl <- clNames[j]
    mean.error[i, j] <- mean(sapply(1:length(res$bestIdx[[ds]][[cl]]),
                                    function (x) res$testError[[ds]][[cl]][x, res$bestIdx[[ds]][[cl]][x]]))
    train.time[i, j] <- mean(sapply(1:length(res$bestIdx[[ds]][[cl]]),
                                    function (x) res$trainTime[[ds]][[cl]][x, res$bestIdx[[ds]][[cl]][x]]))
    tree.strength[i, j] <- mean(sapply(1:length(res$bestIdx[[ds]][[cl]]),
                                    function (x) res$treeStrength[[ds]][[cl]][x, res$bestIdx[[ds]][[cl]][x]]))
    tree.correlation[i, j] <- mean(sapply(1:length(res$bestIdx[[ds]][[cl]]),
                                    function (x) res$treeCorr[[ds]][[cl]][x, res$bestIdx[[ds]][[cl]][x]]))
  }
}

no.na <- apply(mean.error, 1, function(x) !any(is.na(x)))
chance.error <- chance.error[no.na]
mean.error <- mean.error[no.na, ]
train.time <- train.time[no.na, ]
tree.strength <- tree.strength[no.na, ]
tree.correlation <- tree.correlation[no.na, ]
n <- n[no.na]
p <- p[no.na]
pcat <- pcat[no.na]
dataSets <- dataSets[no.na]

results <- data.frame(data.set = melt(mean.error)[[1L]], classifier = melt(mean.error)[[2L]], n = rep(n, length(classifiers)),
                      p = rep(p, length(classifiers)), pcat = rep(pcat, length(classifiers)), mean.error = melt(mean.error)[[3L]],
                      train.time = melt(train.time)[[3L]], tree.strength = melt(tree.strength)[[3L]],
                      tree.correlation = melt(tree.correlation)[[3L]])

norm.error <- mean.error/chance.error

norm.rel.error <- (mean.error[, 2:ncol(mean.error)] - mean.error[, 1L])/chance.error
colnames(norm.rel.error) <- classifiers[2:length(classifiers)]
# norm.rel.error[norm.rel.error == 0] <- NA

results <- data.frame(data.set = melt(norm.error)[[1L]], classifier = melt(norm.error)[[2L]], num.points = rep(n, length(classifiers)),
                      num.dims = rep(p, length(classifiers)), num.cat.dims = rep(pcat, length(classifiers)), Lhat.normalized = melt(norm.error)[[3L]],
                      training.time = melt(train.time)[[3L]], strength = melt(tree.strength)[[3L]],
                      corr = melt(tree.correlation)[[3L]])

results.relative <- data.frame(data.set = melt(norm.rel.error)[[1L]], classifier = melt(norm.rel.error)[[2L]], num.points = rep(n, length(classifiers) - 1L),
                               num.dims = rep(p, length(classifiers) - 1L), num.cat.dims = rep(pcat, length(classifiers) - 1L), Lhat.relative = melt(norm.rel.error)[[3L]])

mn <- min(norm.rel.error, na.rm = T)
mx <- max(norm.rel.error, na.rm = T)
eval.points <- seq(mn, mx, (mx - mn)/100)
kscdf <- sapply(lapply(colnames(norm.rel.error), FUN = function(cname) kcde(x = norm.rel.error[!is.na(norm.rel.error[, cname]), cname], eval.points = eval.points)),
                FUN = function(x) x$estimate)
for (j in 1:ncol(kscdf)) {
  kscdf[kscdf[, j] == 0 & cumsum(kscdf[, j]) > 0, j] <- 1
}
colnames(kscdf) <- colnames(norm.rel.error)
# kscdf <- kcde(x = norm.rel.error[, 1], eval.points = eval.points)
ccol <- rainbow(length(classifiers) - 1L)

df <- melt(as.matrix(norm.rel.error))
names(df) <- c("ind","classifier","relative error")
p1 <- ggplot(df, aes(x=`relative error`)) +
  scale_color_manual(values=ccol) +
  geom_density(aes(group=classifier, colour=classifier)) +
  xlab("Normalized Relative Error") +
  ylab("Probability Density Estimate") +
  ggtitle(paste0(nrow(norm.rel.error), " UCI Benchmark Datasets"))
ggsave(filename = "~/RandomerForest/R/Figures/2017.09.23/Benchmark_kpde_2017_09_23.png", plot = p1)
ggsave(filename = "~/RandomerForest/R/Figures/2017.09.23/Benchmark_kpde_2017_09_23.pdf", plot = p1)

cdf.plot <- melt(kscdf)
names(cdf.plot) <- c("relative error", "classifier", "kcde")
cdf.plot$`relative error` <- eval.points[cdf.plot$`relative error`]
p2 <- ggplot(cdf.plot, aes(x = `relative error`, y = kcde)) +
  scale_color_manual(values = ccol) +
  geom_line(aes(group = classifier, colour = classifier)) +
  xlab("Normalized Relative Error") +
  ylab("Cumulative Distribution Estimate") +
  ggtitle(paste0(nrow(norm.rel.error), " UCI Benchmark Datasets"))
ggsave(filename = "~/RandomerForest/R/Figures/2017.09.23/Benchmark_kcde_2017_09_23.png", plot = p2)
ggsave(filename = "~/RandomerForest/R/Figures/2017.09.23/Benchmark_kcde_2017_09_23.pdf", plot = p2)

hmap <- d1heat(norm.rel.error, trunc = NA)
p3 <- plot(hmap, bincount = T) +
  xlab("Normalized Relative Error") +
  ylab("Classifier") +
  ggtitle(paste0(nrow(norm.rel.error), " UCI Benchmark Datasets"))
ggsave(filename = "~/RandomerForest/R/Figures/2017.09.23/Benchmark_heatmap_2017_09_23.png", plot = p3)
ggsave(filename = "~/RandomerForest/R/Figures/2017.09.23/Benchmark_heatmap_2017_09_23.pdf", plot = p3)

p4 <- ggplot(results.relative, aes(x = num.points, y = Lhat.relative)) +
  geom_point(aes(colour = classifier)) +
  scale_x_log10()
ggsave(filename = "~/RandomerForest/R/Figures/2017.09.23/L_vs_n_2017_09_23.png", plot = p4)
ggsave(filename = "~/RandomerForest/R/Figures/2017.09.23/L_vs_n_2017_09_23.pdf", plot = p4)
p5 <- ggplot(results.relative, aes(x = num.dims, y = Lhat.relative)) +
  geom_point(aes(colour = classifier)) +
  scale_x_log10()
ggsave(filename = "~/RandomerForest/R/Figures/2017.09.23/L_vs_p_2017_09_23.png", plot = p5)
ggsave(filename = "~/RandomerForest/R/Figures/2017.09.23/L_vs_p_2017_09_23.pdf", plot = p5)
p6 <- ggplot(results.relative, aes(x = num.cat.dims/num.dims, y = Lhat.relative)) +
  geom_point(aes(colour = classifier))
ggsave(filename = "~/RandomerForest/R/Figures/2017.09.23/L_vs_pcat_2017_09_23.png", plot = p6)
ggsave(filename = "~/RandomerForest/R/Figures/2017.09.23/L_vs_pcat_2017_09_23.pdf", plot = p6)
p7 <- ggplot(results.relative, aes(x = num.points/num.dims, y = Lhat.relative)) +
  geom_point(aes(colour = classifier)) +
  scale_x_log10()
ggsave(filename = "~/RandomerForest/R/Figures/2017.09.23/L_vs_n_over_p_2017_09_23.png", plot = p7)
ggsave(filename = "~/RandomerForest/R/Figures/2017.09.23/L_vs_n_over_p_2017_09_23.pdf", plot = p7)
p8 <- ggplot(results, aes(x = corr, y = strength)) +
  geom_point(aes(colour = classifier))
ggsave(filename = "~/RandomerForest/R/Figures/2017.09.23/strength_correlation_2017_09_23.png", plot = p8)
ggsave(filename = "~/RandomerForest/R/Figures/2017.09.23/strength_correlation_2017_09_23.pdf", plot = p8)

# multiplot(p1, p3, p3, cols = 1L)