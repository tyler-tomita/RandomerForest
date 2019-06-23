rm(list = ls())
library(reshape2)
library(ggplot2)
library(ks)
library(meda)
source("~/RandomerForest/R/Code/Utils/GetCatMap.R")

filePath <- "~/RandomerForest/R/Results/2017.10.03/"

contents <- list.files(filePath)
contents <- contents[!grepl("frc", contents) & !grepl("rerfc", contents) & !grepl("rr-rf", contents)]
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

# classifiers <- unique(unlist(classifiers))
classifiers <- c("rerfr", "rerfc", "rf", "frc")
class.map <- list(rerfr = "RerF(r)", rerfc = "RerF-C", rf = "RF", frc = "F-RC")
class.map2 <- class.map[names(class.map) != "rf"]

color.map <- c("#41ab5d", "#4292c6", "#f768a1")
line.width <- 1.5

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
  # clNames <- names(res$testError[[ds]])
  for (j in seq.int(length(classifiers))) {
    cl <- classifiers[j]
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
dataSets <- dataSets[no.na]

norm.error <- mean.error/chance.error

norm.rel.error <- (mean.error[, colnames(mean.error) != "rf"] - mean.error[, "rf"])/chance.error
colnames(norm.rel.error) <- classifiers[classifiers != "rf"]
# norm.rel.error[norm.rel.error == 0] <- NA

results <- data.frame(data.set = melt(norm.error)[[1L]], classifier = melt(norm.error)[[2L]], num.points = rep(n, length(classifiers)),
                      num.dims = rep(p, length(classifiers)), pcat = rep(pcat, length(classifiers)), Lhat.normalized = melt(norm.error)[[3L]],
                      training.time = melt(train.time)[[3L]], strength = melt(tree.strength)[[3L]],
                      corr = melt(tree.correlation)[[3L]])
results$classifier <- factor(unlist(class.map[results$classifier]), levels = sapply(names(class.map), function(x) class.map[[x]]))

results.relative <- data.frame(data.set = melt(norm.rel.error)[[1L]], classifier = melt(norm.rel.error)[[2L]], num.points = rep(n, length(classifiers) - 1L),
                               num.dims = rep(p, length(classifiers) - 1L), num.cat.dims = rep(pcat, length(classifiers) - 1L), Lhat.relative = melt(norm.rel.error)[[3L]])
results.relative$classifier <- factor(unlist(class.map2[results.relative$classifier]), levels = sapply(names(class.map2), function(x) class.map2[[x]]))

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

df <- melt(as.matrix(norm.rel.error))
names(df) <- c("ind","classifier","relative error")
df$classifier <- factor(unlist(class.map2[df$classifier]), levels = sapply(names(class.map2), function(x) class.map2[[x]]))
p1 <- ggplot(df, aes(x=`relative error`)) +
  geom_density(aes(group=classifier, colour=classifier), size = line.width) +
  xlab("Normalized Relative Error") +
  ylab("Probability Density Estimate") +
  ggtitle(paste0(nrow(norm.rel.error), " UCI Benchmark Datasets")) +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map)
ggsave(filename = "~/RandomerForest/R/Figures/2017.10.03/Benchmark_kpde_2017_10_03.png", plot = p1)
ggsave(filename = "~/RandomerForest/R/Figures/2017.10.03/Benchmark_kpde_2017_10_03.pdf", plot = p1)

cdf.plot <- melt(kscdf)
names(cdf.plot) <- c("relative error", "classifier", "kcde")
cdf.plot$classifier <- factor(unlist(class.map2[cdf.plot$classifier]), levels = sapply(names(class.map2), function(x) class.map2[[x]]))
cdf.plot$`relative error` <- eval.points[cdf.plot$`relative error`]
p2 <- ggplot(cdf.plot, aes(x = `relative error`, y = kcde)) +
  geom_line(aes(group = classifier, colour = classifier), size = line.width) +
  scale_x_continuous(limits = c(-0.25, 0.25)) +
  xlab("Normalized Relative Error") +
  ylab("Cumulative Distribution Estimate") +
  ggtitle(paste0(nrow(norm.rel.error), " UCI Benchmark Datasets")) +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map)
ggsave(filename = "~/RandomerForest/R/Figures/2017.10.03/Benchmark_kcde_2017_10_03.png", plot = p2)
ggsave(filename = "~/RandomerForest/R/Figures/2017.10.03/Benchmark_kcde_2017_10_03.pdf", plot = p2)

hmap <- d1heat(norm.rel.error, trunc = NA, breaks = seq(-1, 1, 0.25))
p3 <- plot(hmap, bincount = T) +
  xlab("Normalized Relative Error") +
  ylab("") +
  ggtitle(paste0(nrow(norm.rel.error), " UCI Benchmark Datasets")) +
  scale_y_discrete(labels = unlist(class.map))
ggsave(filename = "~/RandomerForest/R/Figures/2017.10.03/Benchmark_heatmap_2017_10_03.png", plot = p3)
ggsave(filename = "~/RandomerForest/R/Figures/2017.10.03/Benchmark_heatmap_2017_10_03.pdf", plot = p3)

p4 <- ggplot(results.relative, aes(x = num.points, y = Lhat.relative)) +
  geom_point(aes(colour = classifier), alpha = 0.3) +
  geom_smooth(se = F, aes(colour = classifier)) +
  scale_x_log10() +
  xlab("n") +
  ylab("normalized relative error") +
  ggtitle(paste0(nrow(norm.rel.error), " UCI Benchmark Datasets")) +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map)
ggsave(filename = "~/RandomerForest/R/Figures/2017.10.03/L_vs_n_2017_10_03.png", plot = p4)
ggsave(filename = "~/RandomerForest/R/Figures/2017.10.03/L_vs_n_2017_10_03.pdf", plot = p4)
p5 <- ggplot(results.relative, aes(x = num.dims, y = Lhat.relative)) +
  geom_point(aes(colour = classifier), alpha = 0.3) +
  geom_smooth(se = F, aes(colour = classifier)) +
  scale_x_log10() +
  xlab("p") +
  ylab("normalized relative error") +
  ggtitle(paste0(nrow(norm.rel.error), " UCI Benchmark Datasets")) +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map)
ggsave(filename = "~/RandomerForest/R/Figures/2017.10.03/L_vs_p_2017_10_03.png", plot = p5)
ggsave(filename = "~/RandomerForest/R/Figures/2017.10.03/L_vs_p_2017_10_03.pdf", plot = p5)
p6 <- ggplot(results.relative, aes(x = num.cat.dims/num.dims, y = Lhat.relative)) +
  geom_point(aes(colour = classifier), alpha = 0.3) +
  geom_smooth(se = F, aes(colour = classifier)) +
  xlab(expression("p"[cat]*"/p")) +
  ylab("normalized relative error") +
  ggtitle(paste0(nrow(norm.rel.error), " UCI Benchmark Datasets")) +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map)
ggsave(filename = "~/RandomerForest/R/Figures/2017.10.03/L_vs_pcat_2017_10_03.png", plot = p6)
ggsave(filename = "~/RandomerForest/R/Figures/2017.10.03/L_vs_pcat_2017_10_03.pdf", plot = p6)
p7 <- ggplot(results.relative, aes(x = num.points/num.dims, y = Lhat.relative)) +
  geom_point(aes(colour = classifier), alpha = 0.3) +
  geom_smooth(se = F, aes(colour = classifier)) +
  scale_x_log10() +
  xlab("n/p") +
  ylab("normalized relative error") +
  ggtitle(paste0(nrow(norm.rel.error), " UCI Benchmark Datasets")) +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map)
ggsave(filename = "~/RandomerForest/R/Figures/2017.10.03/L_vs_n_over_p_2017_10_03.png", plot = p7)
ggsave(filename = "~/RandomerForest/R/Figures/2017.10.03/L_vs_n_over_p_2017_10_03.pdf", plot = p7)
p8 <- ggplot(results[results$classifier != "RerF-C", ], aes(x = corr, y = strength)) +
  geom_point(aes(colour = classifier), alpha = 0.3) +
  geom_smooth(se = F, aes(colour = classifier)) +
  xlab("tree correlation") +
  ylab("tree strength") +
  ggtitle(paste0(nrow(norm.rel.error), " UCI Benchmark Datasets")) +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map)
ggsave(filename = "~/RandomerForest/R/Figures/2017.10.03/strength_correlation_2017_10_03.png", plot = p8)
ggsave(filename = "~/RandomerForest/R/Figures/2017.10.03/strength_correlation_2017_10_03.pdf", plot = p8)

# multiplot(p1, p3, p3, cols = 1L)

