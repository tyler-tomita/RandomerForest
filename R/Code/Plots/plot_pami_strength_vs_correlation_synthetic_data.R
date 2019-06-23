rm(list = ls())
library(ggplot2)
library(ggpubr)
library(scales)
library(RColorBrewer)

# Make plots for Sparse parity

load("~/RandomerForest/R/Results/2017.09.27/Sparse_parity_2017_09_27.RData")

length.n <- dim(testError[[1]])[1]
length.p <- dim(testError[[1]])[2]
num.trials <- dim(testError[[1]])[3]
classifiers <- c("rf", "rerf", "frc")
class.map <- list(rerf = "RerF", rf = "RF", frc = "F-RC")

ps <- c(3, 10, 20)
ns <- list(c(10, 100, 1000), c(100, 1000, 10000), c(1000, 5000, 10000))

# color.map <- c("#41ab5d", "#4292c6", "#f768a1")
# color.map <- brewer.pal(5L, "Set2")
color.map <- c("CCF" = "#6a51a3",
               "RF" = "#4a484c",
               "F-RC" = "#318dde",
               "RerF" = "#F41711",
               "XGBoost" = "#E78AC3")
line.width <- 1
marker.size <- 2

classifier <- vector("character", length(classifiers)*length.n*length.p*num.trials)
num.obs <- vector("integer", length(classifiers)*length.n*length.p*num.trials)
num.dims <- vector("integer", length(classifiers)*length.n*length.p*num.trials)
test.error <- vector("double", length(classifiers)*length.n*length.p*num.trials)
train.time <- vector("double", length(classifiers)*length.n*length.p*num.trials)
tree.strength <- vector("double", length(classifiers)*length.n*length.p*num.trials)
tree.correlation <- vector("double", length(classifiers)*length.n*length.p*num.trials)
for (i in 1:length(classifiers)) {
  cl <- classifiers[i]
  classifier[((i - 1L)*length.n*length.p*num.trials + 1L):(i*length.n*length.p*num.trials)] <- class.map[[cl]]
  for (j in 1:length.p) {
    num.dims[(((i - 1L)*length.p + (j - 1L))*length.n*num.trials + 1L):(((i - 1L)*length.p + j)*length.n*num.trials)] <- ps[j]
    for (k in 1:length.n) {
      num.obs[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- ns[[j]][k]
      test.error[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) testError[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
      train.time[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) trainTime[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
      tree.strength[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) treeStrength[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
      tree.correlation[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) treeCorr[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
    }
  }
}

df <- data.frame(classifier = factor(classifier, levels = sapply(names(class.map), function(x) class.map[[x]])))
df$p <- num.dims
df$n <- num.obs
df$error <- test.error
df$train.time <- train.time
df$tree.strength <- tree.strength
df$tree.correlation <- tree.correlation


panel.label <- textGrob("(A)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- min(df$tree.correlation[(df$p == 20L) & (df$n == 1000L)])
xmax <- max(df$tree.correlation[(df$p == 20L) & (df$n == 1000L)])
x.range <- xmax - xmin
xmin <- xmin - 0.1*x.range
xmax <- xmax + 0.1*x.range
x.range <- xmax - xmin
label.x <- xmin - 0.42*x.range
ymin <- min(df$tree.strength[(df$p == 20L) & (df$n == 1000L)])
ymax <- max(df$tree.strength[(df$p == 20L) & (df$n == 1000L)])
y.range <- ymax - ymin
ymin <- ymin - 0.1*y.range
ymax <- ymax + 0.1*y.range
label.y <- ymax + 0.48*y.range

p1 <- ggplot(df[(df$p == 20L) & (df$n == 1000L), ], aes(x = tree.correlation, y = tree.strength)) +
  geom_point(aes(colour = classifier), size = marker.size) +
  theme(axis.line = element_line(size = 1.5), text = element_text(size = 14)) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = seq(0.02, 0.04, 0.01)) +
  scale_y_continuous(limits = c(ymin, ymax)) +
  xlab("Tree Correlation") +
  ylab("Tree Strength") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  ggtitle("Sparse Parity", subtitle = "n = 1000") +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map[c(4, 2, 3)])

gt1 <- ggplot_gtable(ggplot_build(p1))
gt.legend <- gt1$grobs[[which(sapply(gt1$grobs, function(x) x$name == "guide-box"))]]
p1 <- p1 + guides(color = F)
gt1 <- ggplot_gtable(ggplot_build(p1))
gt1$layout$clip[] <- "off"

# Make plots for Orthant

load("~/RandomerForest/R/Results/2017.10.01/Orthant_2017_10_01.RData")

length.n <- dim(testError[[1]])[1]
length.p <- dim(testError[[1]])[2]
num.trials <- dim(testError[[1]])[3]
ps <- c(2, 4, 6)
ns <- list(c(20,200,400), c(80,400,4000), c(400,2000,4000))

classifier <- vector("character", length(classifiers)*length.n*length.p*num.trials)
num.obs <- vector("integer", length(classifiers)*length.n*length.p*num.trials)
num.dims <- vector("integer", length(classifiers)*length.n*length.p*num.trials)
test.error <- vector("double", length(classifiers)*length.n*length.p*num.trials)
train.time <- vector("double", length(classifiers)*length.n*length.p*num.trials)
tree.strength <- vector("double", length(classifiers)*length.n*length.p*num.trials)
tree.correlation <- vector("double", length(classifiers)*length.n*length.p*num.trials)
for (i in 1:length(classifiers)) {
  cl <- classifiers[i]
  classifier[((i - 1L)*length.n*length.p*num.trials + 1L):(i*length.n*length.p*num.trials)] <- class.map[[cl]]
  for (j in 1:length.p) {
    num.dims[(((i - 1L)*length.p + (j - 1L))*length.n*num.trials + 1L):(((i - 1L)*length.p + j)*length.n*num.trials)] <- ps[j]
    for (k in 1:length.n) {
      num.obs[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- ns[[j]][k]
      test.error[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) testError[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
      train.time[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) trainTime[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
      tree.strength[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) treeStrength[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
      tree.correlation[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) treeCorr[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
    }
  }
}

df <- data.frame(classifier = factor(classifier, levels = sapply(names(class.map), function(x) class.map[[x]])))
df$p <- num.dims
df$n <- num.obs
df$error <- test.error
df$train.time <- train.time
df$tree.strength <- tree.strength
df$tree.correlation <- tree.correlation

panel.label <- textGrob("(B)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- min(df$tree.correlation[(df$p == 6L) & (df$n == 400L)])
xmax <- max(df$tree.correlation[(df$p == 6L) & (df$n == 400L)])
x.range <- xmax - xmin
xmin <- xmin - 0.1*x.range
xmax <- xmax + 0.1*x.range
x.range <- xmax - xmin
label.x <- xmin - 0.42*x.range
ymin <- min(df$tree.strength[(df$p == 6L) & (df$n == 400L)])
ymax <- max(df$tree.strength[(df$p == 6L) & (df$n == 400L)])
y.range <- ymax - ymin
ymin <- ymin - 0.1*y.range
ymax <- ymax + 0.1*y.range
label.y <- ymax + 0.48*y.range

p2 <- ggplot(df[(df$p == 6L) & (df$n == 400L), ], aes(x = tree.correlation, y = tree.strength)) +
  geom_point(aes(colour = classifier), size = marker.size) +
  theme(axis.line = element_line(size = 1.5), text = element_text(size = 14)) +
  scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin, ymax)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Orthant", subtitle = "n = 400") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  guides(color=F) +
  theme_classic() +
  theme(plot.margin = unit(c(0.1, 0, 0.3, 0), "in"), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map[c(4, 2, 3)])

gt2 <- ggplot_gtable(ggplot_build(p2))
gt2$layout$clip[] <- "off"

# Make plots for Trunk

load("~/RandomerForest/R/Results/2017.10.07/Trunk_2017_10_07_aggregated2.RData")

length.n <- dim(testError[[1]])[1]
length.p <- dim(testError[[1]])[2]
num.trials <- dim(testError[[1]])[3]
ps <- c(10, 100, 1000)
ns <- list(c(10, 100, 1000, 10000), c(10, 100, 1000, 10000), c(10, 100, 1000, 10000))

classifier <- vector("character", length(classifiers)*length.n*length.p*num.trials)
num.obs <- vector("integer", length(classifiers)*length.n*length.p*num.trials)
num.dims <- vector("integer", length(classifiers)*length.n*length.p*num.trials)
test.error <- vector("double", length(classifiers)*length.n*length.p*num.trials)
train.time <- vector("double", length(classifiers)*length.n*length.p*num.trials)
tree.strength <- vector("double", length(classifiers)*length.n*length.p*num.trials)
tree.correlation <- vector("double", length(classifiers)*length.n*length.p*num.trials)
for (i in 1:length(classifiers)) {
  cl <- classifiers[i]
  classifier[((i - 1L)*length.n*length.p*num.trials + 1L):(i*length.n*length.p*num.trials)] <- class.map[[cl]]
  for (j in 1:length.p) {
    num.dims[(((i - 1L)*length.p + (j - 1L))*length.n*num.trials + 1L):(((i - 1L)*length.p + j)*length.n*num.trials)] <- ps[j]
    for (k in 1:length.n) {
      num.obs[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- ns[[j]][k]
      test.error[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, j, x] != 0, testError[[cl]][k, j, x, bestIdx[[cl]][k, j, x]], NA))
      train.time[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, j, x] != 0, trainTime[[cl]][k, j, x, bestIdx[[cl]][k, j, x]], NA))
      tree.strength[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, j, x] != 0, treeStrength[[cl]][k, j, x, bestIdx[[cl]][k, j, x]], NA))
      tree.correlation[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, j, x] != 0, treeCorr[[cl]][k, j, x, bestIdx[[cl]][k, j, x]], NA))
    }
  }
}

df <- data.frame(classifier = factor(classifier, levels = sapply(names(class.map), function(x) class.map[[x]])))
df$p <- num.dims
df$n <- num.obs
df$error <- test.error
df$train.time <- train.time
df$tree.strength <- tree.strength
df$tree.correlation <- tree.correlation

panel.label <- textGrob("(C)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- min(df$tree.correlation[(df$p == 10L) & (df$n == 10L)])
xmax <- max(df$tree.correlation[(df$p == 10L) & (df$n == 10L)])
x.range <- xmax - xmin
xmin <- xmin - 0.1*x.range
xmax <- xmax + 0.1*x.range
x.range <- xmax - xmin
label.x <- xmin - 0.42*x.range
ymin <- min(df$tree.strength[(df$p == 10L) & (df$n == 10L)])
ymax <- max(df$tree.strength[(df$p == 10L) & (df$n == 10L)])
y.range <- ymax - ymin
ymin <- ymin - 0.1*y.range
ymax <- ymax + 0.1*y.range
label.y <- ymax + 0.48*y.range

p3 <- ggplot(df[(df$p == 10L) & (df$n == 10L), ], aes(x = tree.correlation, y = tree.strength)) +
  geom_point(aes(colour = classifier), size = marker.size) +
  scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin, ymax)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Trunk", subtitle = "n = 10") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  guides(color=F) +
  theme_classic() +
  theme(plot.margin = unit(c(0.1, 0, 0.3, 0), "in"), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map[c(4, 2, 3)])

gt3 <- ggplot_gtable(ggplot_build(p3))
gt3$layout$clip[] <- "off"

panel.label <- textGrob("(D)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- min(df$tree.correlation[(df$p == 10L) & (df$n == 100L)])
xmax <- max(df$tree.correlation[(df$p == 10L) & (df$n == 10L)])
x.range <- xmax - xmin
xmin <- xmin - 0.1*x.range
xmax <- xmax + 0.1*x.range
x.range <- xmax - xmin
label.x <- xmin - 0.42*x.range
ymin <- min(df$tree.strength[(df$p == 10L) & (df$n == 100L)])
ymax <- max(df$tree.strength[(df$p == 10L) & (df$n == 100L)])
y.range <- ymax - ymin
ymin <- ymin - 0.1*y.range
ymax <- ymax + 0.1*y.range
label.y <- ymax + 0.48*y.range

p4 <- ggplot(df[(df$p == 10L) & (df$n == 100L), ], aes(x = tree.correlation, y = tree.strength)) +
  geom_point(aes(colour = classifier), size = marker.size) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = seq(0.15, 0.35, 0.1)) +
  scale_y_continuous(limits = c(ymin, ymax)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Trunk", subtitle = "n = 100") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  guides(color=F) +
  theme_classic() +
  theme(plot.margin = unit(c(0.1, 0, 0.3, 0), "in"), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map[c(4, 2, 3)])

gt4 <- ggplot_gtable(ggplot_build(p4))
gt4$layout$clip[] <- "off"

pall <- ggarrange(gt1, gt2, gt3, gt4, gt.legend, ncol = 5L, widths = c(1, 1, 1, 1, 0.5))
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/strength_vs_correlation_synthetic_data.pdf", plot = pall, width = 10.5, height = 2.5, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/strength_vs_correlation_synthetic_data.png", plot = pall, width = 10.5, height = 2.5, units = "in")
