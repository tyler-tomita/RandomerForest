rm(list = ls())
library(ggplot2)
library(ggpubr)
library(scales)
library(RColorBrewer)
library(grid)
library(gridExtra)

set.seed(123L)

# Make plots for Sparse parity

load("~/RandomerForest/R/Results/2017.09.27/Sparse_parity_2017_09_27.RData")

length.n <- dim(testError[[1]])[1]
length.p <- dim(testError[[1]])[2]
num.trials <- dim(testError[[1]])[3]
classifiers <- c("rf", "rerf", "frc", "ccf")
class.map <- list(rerf = "RerF", rf = "RF", frc = "F-RC", ccf = "CCF")

ps <- c(3, 10, 20)
ns <- list(c(10, 100, 1000), c(100, 1000, 10000), c(1000, 5000, 10000))

# color.map <- c("#41ab5d", "#4292c6", "#f768a1")
# color.map <- brewer.pal(5L, "Set2")
color.map2 <- brewer.pal(8L, "Dark2")
# color.map2 <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")[c(2L, 3L, 1L, 4L)]
color.map <- c("CCF" = "#6a51a3",
               "RF" = "#4a484c",
               "F-RC" = "#318dde",
               "RerF" = "#F41711",
               "XGBoost" = "#E78AC3")

line.width <- 0.75
marker.size <- 0.75

# generate Sparse parity data

nplot <- 400L
pplot <- 3L
X <- matrix(runif(nplot*pplot, min = -1, max = 1), nplot, pplot)
Y <- as.factor(apply(X, 1, function(x) sum(x > 0)%%2))
df.plot <- data.frame(x1 = X[, 1], x2 = X[, 2], x3 = X[, 3], class = Y)

panel.label <- textGrob("(A)", gp=gpar(fontsize=14, fontface="bold"))
label.x <- -1 - 0.5*2
label.y <- 0.31*2 + 1
panel.label2 <- textGrob("Sparse Parity", rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
label.x2 <- -1 - 0.75*2
label.y2 <- 0

p1 <- ggplot(df.plot[df.plot$x3 < 0, ], aes(x = x1, y = x2, color = class)) +
  geom_point(size = marker.size) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  scale_y_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  xlab(expression(X[1])) +
  # ylab(expression(atop(bold("Sparse Parity"), X[2]))) +
  ylab(expression(X[2])) +
  # ggtitle(expression(X[3] %in% (list(-1, 0))*","~X[4]*","~ldots*","~X[20] %in% (list(-1, 1)))) +
  ggtitle(expression(X[3] %in% (list(-1, 0)))) +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  annotation_custom(panel.label2, xmin = label.x2, xmax = label.x2, ymin = label.y2, ymax = label.y2) +
  guides(color = F) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.title = element_text(size = 14)) +
  scale_color_manual(values = color.map2[1:2])

gt1 <- ggplot_gtable(ggplot_build(p1))
gt1$layout$clip[] <- "off"

panel.label <- textGrob("(B)", gp=gpar(fontsize=14, fontface="bold"))
label.x <- -1 - 0.5*2
label.y <- 0.31*2 + 1

p2 <- ggplot(df.plot[df.plot$x3 > 0, ], aes(x = x1, y = x2, color = class)) +
  geom_point(size = marker.size) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  scale_y_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  xlab(expression(X[1])) +
  ylab(expression(X[2])) +
  # ggtitle(expression(X[3] %in% (list(0, 1))*","~X[4]*","~ldots*","~X[20] %in% (list(-1, 1)))) +
  ggtitle(expression(X[3] %in% (list(0, 1)))) +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  guides(color = F) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.title = element_text(size = 14)) +
  scale_color_manual(values = color.map2[1:2])

gt2 <- ggplot_gtable(ggplot_build(p2))
gt2$layout$clip[] <- "off"

# generate Orthant data

X <- matrix(runif(nplot*pplot, min = -1, max = 1), nplot, pplot)
Xbin <- X > 0
Y <- as.factor(apply(sapply(1:3, function(x) Xbin[, x]*2^(3 - x)), 1, sum))
df.plot <- data.frame(x1 = X[, 1], x2 = X[, 2], x3 = X[, 3], class = Y)

panel.label <- textGrob("(D)", gp=gpar(fontsize=14, fontface="bold"))
label.x <- -1 - 0.5*2
label.y <- 0.31*2 + 1
panel.label2 <- textGrob("Orthant", rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
label.x2 <- -1 - 0.75*2
label.y2 <- 0

p3 <- ggplot(df.plot[df.plot$x3 < 0, ], aes(x = x1, y = x2, color = class)) +
  geom_point(size = marker.size) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  scale_y_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  xlab(expression(X[1])) +
  # ylab(expression(atop(bold("Orthant"), X[2]))) +
  ylab(expression(X[2])) +
  # ggtitle(expression(X[3] %in% (list(-1, 0))*","~X[4]*","~ldots*","~X[6] %in% (list(-1, 0)))) +
  ggtitle(expression(X[3] %in% (list(-1, 0)))) +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  annotation_custom(panel.label2, xmin = label.x2, xmax = label.x2, ymin = label.y2, ymax = label.y2) +
  guides(color = F) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.title = element_text(size = 14)) +
  scale_color_manual(values = color.map2[1:4])

gt3 <- ggplot_gtable(ggplot_build(p3))
gt3$layout$clip[] <- "off"

panel.label <- textGrob("(E)", gp=gpar(fontsize=14, fontface="bold"))
label.x <- -1 - 0.5*2
label.y <- 0.31*2 + 1

p4 <- ggplot(df.plot[df.plot$x3 > 0, ], aes(x = x1, y = x2, color = class)) +
  geom_point(size = marker.size) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  scale_y_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  xlab(expression(X[1])) +
  ylab(expression(X[2])) +
  # ggtitle(expression(X[3] %in% (list(0, 1))*","~X[4]*","~ldots*","~X[6] %in% (list(-1, 0)))) +
  ggtitle(expression(X[3] %in% (list(0, 1)))) +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  guides(color = F) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.title = element_text(size = 14)) +
  scale_color_manual(values = color.map2[5:8])

gt4 <- ggplot_gtable(ggplot_build(p4))
gt4$layout$clip[] <- "off"

classifier <- vector("character", length(classifiers)*length.n*length.p*num.trials)
num.obs <- vector("integer", length(classifiers)*length.n*length.p*num.trials)
num.dims <- vector("integer", length(classifiers)*length.n*length.p*num.trials)
test.error <- vector("double", length(classifiers)*length.n*length.p*num.trials)
train.time <- vector("double", length(classifiers)*length.n*length.p*num.trials)
tree.strength <- vector("double", length(classifiers)*length.n*length.p*num.trials)
tree.correlation <- vector("double", length(classifiers)*length.n*length.p*num.trials)
error.ccf <- array(read.table("~/RandomerForest/Results/2018.07.02/Sparse_parity_ccf_2018_07_02.csv", header = FALSE, sep = ",", row.names = NULL)$V1, c(3,3,10))
for (i in 1:length(classifiers)) {
  cl <- classifiers[i]
  classifier[((i - 1L)*length.n*length.p*num.trials + 1L):(i*length.n*length.p*num.trials)] <- class.map[[cl]]
  for (j in 1:length.p) {
    num.dims[(((i - 1L)*length.p + (j - 1L))*length.n*num.trials + 1L):(((i - 1L)*length.p + j)*length.n*num.trials)] <- ps[j]
    for (k in 1:length.n) {
      if (cl != "ccf") {
        num.obs[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- ns[[j]][k]
        test.error[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) testError[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
        train.time[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) trainTime[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
        tree.strength[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) treeStrength[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
        tree.correlation[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) treeCorr[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
      } else {
        num.obs[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- ns[[j]][k]
        test.error[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- error.ccf[k, j, 1:num.trials]
        train.time[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- 0
        tree.strength[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- 0
        tree.correlation[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- 0
      }
    }
  }
}

df <- data.frame(classifier = factor(classifier, levels = sapply(names(class.map), function(x) class.map[[x]])))
df$p <- num.dims
df$n <- log10(num.obs)
df$error <- test.error
df$train.time <- train.time
df$tree.strength <- tree.strength
df$tree.correlation <- tree.correlation

df.mean <- aggregate(cbind(error,train.time,tree.strength,tree.correlation)~classifier+p+n, df, mean)
names(df.mean)[4:ncol(df.mean)] <- paste0(names(df)[4:ncol(df)], ".mean")
df.sem <- aggregate(cbind(error,train.time,tree.strength,tree.correlation)~classifier+p+n, df, function(x) sd(x)/sqrt(length(x)))
names(df.sem)[4:ncol(df.sem)] <- paste0(names(df)[4:ncol(df)], ".sem")
df.mean <- cbind(df.mean, df.sem[4:ncol(df.sem)])
df.mean$classifier <- factor(df.mean$classifier, levels = c("CCF", "RF", "F-RC", "RerF"))

# p1 <- ggplot(df.mean[df.mean$p == 3L, ], aes(x = n, y = error.mean)) +
#   geom_line(aes(colour = classifier), size = line.width) +
#   geom_errorbar(aes(ymin=error.mean-error.sem, ymax=error.mean+error.sem, colour = classifier), size = line.width, width = 0.1) +
#   scale_x_log10(breaks = ns[[1L]], labels = function(x) format(x, scientific = T)) +
#   xlab(expression("n"[train])) +
#   # annotation_custom(grob = textGrob("error rate", rot = 90), xmin = -0.1, xmax = -0.1, ymin = )
#   ylab(expression(atop(bold("Sparse Parity"),"error rate"))) +
#   ggtitle("p = 3") +
#   guides(color=guide_legend(title="Algorithm")) +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
#   scale_color_manual(values = color.map)
# 
# p2 <- ggplot(df.mean[df.mean$p == 10L, ], aes(x = n, y = error.mean)) +
#   geom_line(aes(colour = classifier), size = line.width) +
#   geom_errorbar(aes(ymin=error.mean-error.sem, ymax=error.mean+error.sem, colour = classifier), size = line.width, width = 0.1) +
#   scale_x_log10(breaks = ns[[2L]], labels = function(x) format(x, scientific = T)) +
#   xlab(expression("n"[train])) +
#   ylab(expression(atop(bold(" "),"error rate"))) +
#   ggtitle("p = 10") +
#   guides(color=guide_legend(title="Algorithm")) +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
#   scale_color_manual(values = color.map)

panel.label <- textGrob("(C)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- log10(ns[[3L]][1L])
xmax <- log10(ns[[3L]][3L])
x.range <- xmax - xmin
xmin <- xmin - 0.1*x.range
xmax <- xmax + 0.1*x.range
x.range <- xmax - xmin
label.x <- xmin - 0.5*x.range
label.y <- 1.31*0.45

p5 <- ggplot(df.mean[df.mean$p == 20L, ], aes(x = n, y = error.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  geom_errorbar(aes(ymin=error.mean-error.sem, ymax=error.mean+error.sem, colour = classifier), size = line.width, width = 0.1) +
  # scale_x_log10(breaks = ns[[3L]], labels = ns[[3L]]/1000L) +
  # scale_x_log10(breaks = ns[[3L]], labels = function(x) format(x, scientific = T)) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = log10(ns[[3L]]), labels = ns[[3L]]/1000L) +
  scale_y_continuous(limits = c(0, 0.45), breaks = c(0, 0.2, 0.4), labels = c("0.00", "0.20", "0.40")) +
  xlab("n (in thousands)") +
  ylab("Error Rate") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.margin = unit(c(0.43, 0, 0.06, 0), "in")) +
  scale_color_manual(values = color.map[c(4, 2, 3, 1)])

gt5 <- ggplot_gtable(ggplot_build(p5))
gt.legend <- gt5$grobs[[which(sapply(gt5$grobs, function(x) x$name == "guide-box"))]]
p5 <- p5 + guides(color = F)
gt5 <- ggplot_gtable(ggplot_build(p5))
gt5$layout$clip[] <- "off"

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
error.ccf <- array(read.table("~/RandomerForest/Results/2018.07.02/Orthant_ccf_2018_07_02.csv", header = FALSE, sep = ",", row.names = NULL)$V1, c(3,3,10))
for (i in 1:length(classifiers)) {
  cl <- classifiers[i]
  classifier[((i - 1L)*length.n*length.p*num.trials + 1L):(i*length.n*length.p*num.trials)] <- class.map[[cl]]
  for (j in 1:length.p) {
    num.dims[(((i - 1L)*length.p + (j - 1L))*length.n*num.trials + 1L):(((i - 1L)*length.p + j)*length.n*num.trials)] <- ps[j]
    for (k in 1:length.n) {
      if (cl != "ccf") {
        num.obs[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- ns[[j]][k]
        test.error[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) testError[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
        train.time[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) trainTime[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
        tree.strength[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) treeStrength[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
        tree.correlation[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) treeCorr[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
      } else {
        num.obs[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- ns[[j]][k]
        test.error[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- error.ccf[k, j, 1:num.trials]
        train.time[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- 0
        tree.strength[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- 0
        tree.correlation[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- 0
      }
    }
  }
}

df <- data.frame(classifier = factor(classifier, levels = sapply(names(class.map), function(x) class.map[[x]])))
df$p <- num.dims
df$n <- log10(num.obs)
df$error <- test.error
df$train.time <- train.time
df$tree.strength <- tree.strength
df$tree.correlation <- tree.correlation

df.mean <- aggregate(cbind(error,train.time,tree.strength,tree.correlation)~classifier+p+n, df, mean)
names(df.mean)[4:ncol(df.mean)] <- paste0(names(df)[4:ncol(df)], ".mean")
df.sem <- aggregate(cbind(error,train.time,tree.strength,tree.correlation)~classifier+p+n, df, function(x) sd(x)/sqrt(length(x)))
names(df.sem)[4:ncol(df.sem)] <- paste0(names(df)[4:ncol(df)], ".sem")
df.mean <- cbind(df.mean, df.sem[4:ncol(df.sem)])
df.mean$classifier <- factor(df.mean$classifier, levels = c("CCF", "RF", "F-RC", "RerF"))

# p4 <- ggplot(df.mean[df.mean$p == 2L, ], aes(x = n, y = error.mean)) +
#   geom_line(aes(colour = classifier), size = line.width) +
#   geom_errorbar(aes(ymin=error.mean-error.sem, ymax=error.mean+error.sem, colour = classifier), size = line.width, width = 0.1) +
#   scale_x_log10(breaks = ns[[1L]], labels = function(x) format(x, scientific = T)) +
#   xlab(expression("n"[train])) +
#   ylab(expression(atop(bold("Orthant"),"error rate"))) +
#   ggtitle("p = 2") +
#   guides(color=guide_legend(title="Algorithm")) +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
#   scale_color_manual(values = color.map)
# 
# p5 <- ggplot(df.mean[df.mean$p == 4L, ], aes(x = n, y = error.mean)) +
#   geom_line(aes(colour = classifier), size = line.width) +
#   geom_errorbar(aes(ymin=error.mean-error.sem, ymax=error.mean+error.sem, colour = classifier), size = line.width, width = 0.1) +
#   scale_x_log10(breaks = ns[[2L]], labels = function(x) format(x, scientific = T)) +
#   xlab(expression("n"[train])) +
#   ylab(expression(atop(bold(" "),"error rate"))) +
#   ggtitle("p = 4") +
#   guides(color=guide_legend(title="Algorithm")) +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
#   scale_color_manual(values = color.map)

panel.label <- textGrob("(F)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- log10(ns[[3L]][1L])
xmax <- log10(ns[[3L]][3L])
x.range <- xmax - xmin
xmin <- xmin - 0.1*x.range
xmax <- xmax + 0.1*x.range
x.range <- xmax - xmin
label.x <- xmin - 0.5*x.range
label.y <- 1.31*0.19

p6 <- ggplot(df.mean[df.mean$p == 6L, ], aes(x = n, y = error.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  geom_errorbar(aes(ymin=error.mean-error.sem, ymax=error.mean+error.sem, colour = classifier), size = line.width, width = 0.1) +
  # scale_x_log10(breaks = ns[[3L]], labels = ns[[3L]]/1000L) +
  # scale_x_log10(breaks = ns[[3L]], labels = function(x) format(x, scientific = T)) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = log10(ns[[3L]]), labels = ns[[3L]]/1000L) +
  scale_y_continuous(limits = c(0, 0.19), breaks = c(0, 0.05, 0.1, 0.15), labels = c("0.00", "0.05", "0.10", "0.15")) +
  xlab("n (in thousands)") +
  ylab("Error Rate") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  guides(color = F) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.margin = unit(c(0.43, 0, 0.06, 0), "in")) +
  scale_color_manual(values = color.map[c(4, 2, 3, 1)])

gt6 <- ggplot_gtable(ggplot_build(p6))
gt6$layout$clip[] <- "off"

# Make plots for Trunk

# load("~/RandomerForest/R/Results/2017.10.07/Trunk_2017_10_07_aggregated2.RData")
# 
# length.n <- dim(testError[[1]])[1]
# length.p <- dim(testError[[1]])[2]
# num.trials <- dim(testError[[1]])[3]
# 
# ps <- c(10, 100, 1000)
# ns <- list(c(10, 100, 1000, 10000), c(10, 100, 1000, 10000), c(10, 100, 1000, 10000))
# 
# color.map <- c("#41ab5d", "#4292c6", "#f768a1")
# line.width <- 1.5
# 
# classifier <- vector("character", length(classifiers)*length.n*length.p*num.trials)
# num.obs <- vector("integer", length(classifiers)*length.n*length.p*num.trials)
# num.dims <- vector("integer", length(classifiers)*length.n*length.p*num.trials)
# test.error <- vector("double", length(classifiers)*length.n*length.p*num.trials)
# train.time <- vector("double", length(classifiers)*length.n*length.p*num.trials)
# tree.strength <- vector("double", length(classifiers)*length.n*length.p*num.trials)
# tree.correlation <- vector("double", length(classifiers)*length.n*length.p*num.trials)
# for (i in 1:length(classifiers)) {
#   cl <- classifiers[i]
#   classifier[((i - 1L)*length.n*length.p*num.trials + 1L):(i*length.n*length.p*num.trials)] <- class.map[[cl]]
#   for (j in 1:length.p) {
#     num.dims[(((i - 1L)*length.p + (j - 1L))*length.n*num.trials + 1L):(((i - 1L)*length.p + j)*length.n*num.trials)] <- ps[j]
#     for (k in 1:length.n) {
#       num.obs[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- ns[[j]][k]
#       test.error[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, j, x] != 0, testError[[cl]][k, j, x, bestIdx[[cl]][k, j, x]], NA))
#       train.time[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, j, x] != 0, trainTime[[cl]][k, j, x, bestIdx[[cl]][k, j, x]], NA))
#       tree.strength[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, j, x] != 0, treeStrength[[cl]][k, j, x, bestIdx[[cl]][k, j, x]], NA))
#       tree.correlation[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, j, x] != 0, treeCorr[[cl]][k, j, x, bestIdx[[cl]][k, j, x]], NA))
#     }
#   }
# }
# 
# df <- data.frame(classifier = factor(classifier, levels = sapply(names(class.map), function(x) class.map[[x]])))
# df$p <- num.dims
# df$n <- log10(num.obs)
# df$error <- test.error
# df$train.time <- train.time
# df$tree.strength <- tree.strength
# df$tree.correlation <- tree.correlation
# 
# df.mean <- aggregate(cbind(error,train.time,tree.strength,tree.correlation)~classifier+p+n, df, function(x) mean(x, na.rm = T))
# names(df.mean)[4:ncol(df.mean)] <- paste0(names(df)[4:ncol(df)], ".mean")
# df.sem <- aggregate(cbind(error,train.time,tree.strength,tree.correlation)~classifier+p+n, df, function(x) sd(x, na.rm = T)/sqrt(length(x[!is.na(x)])))
# names(df.sem)[4:ncol(df.sem)] <- paste0(names(df)[4:ncol(df)], ".sem")
# df.mean <- cbind(df.mean, df.sem[4:ncol(df.sem)])
# # levels(df.mean$classifier) <- classifiers
# 
# p7 <- ggplot(df.mean[df.mean$p == 10L, ], aes(x = n, y = error.mean)) +
#   geom_line(aes(colour = classifier), size = line.width) +
#   geom_errorbar(aes(ymin=error.mean-error.sem, ymax=error.mean+error.sem, colour = classifier), size = line.width, width = 0.1) +
#   scale_x_log10(breaks = ns[[1L]], labels = function(x) format(x, scientific = T)) +
#   xlab(expression("n"[train])) +
#   # annotation_custom(grob = textGrob("error rate", rot = 90), xmin = -0.1, xmax = -0.1, ymin = )
#   ylab(expression(atop(bold("Trunk"),"error rate"))) +
#   ggtitle("p = 10") +
#   guides(color=guide_legend(title="Algorithm")) +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
#   scale_color_manual(values = color.map)
# 
# p8 <- ggplot(df.mean[df.mean$p == 100L, ], aes(x = n, y = error.mean)) +
#   geom_line(aes(colour = classifier), size = line.width) +
#   geom_errorbar(aes(ymin=error.mean-error.sem, ymax=error.mean+error.sem, colour = classifier), size = line.width, width = 0.1) +
#   scale_x_log10(breaks = ns[[2L]], labels = function(x) format(x, scientific = T)) +
#   xlab(expression("n"[train])) +
#   ylab(expression(atop(bold(" "),"error rate"))) +
#   ggtitle("p = 100") +
#   guides(color=guide_legend(title="Algorithm")) +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
#   scale_color_manual(values = color.map)
# 
# p9 <- ggplot(df.mean[df.mean$p == 1000L, ], aes(x = n, y = error.mean)) +
#   geom_line(aes(colour = classifier), size = line.width) +
#   geom_errorbar(aes(ymin=error.mean-error.sem, ymax=error.mean+error.sem, colour = classifier), size = line.width, width = 0.1) +
#   scale_x_log10(breaks = ns[[3L]], labels = function(x) format(x, scientific = T)) +
#   xlab(expression("n"[train])) +
#   ylab(expression(atop(bold(" "),"error rate"))) +
#   ggtitle("p = 1000") +
#   guides(color=guide_legend(title="Algorithm")) +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
#   scale_color_manual(values = color.map)

# pall1 <- ggarrange(p1, p2, ncol = 2L, labels = c("A", "B")) + theme(plot.margin = unit(c(0, 0.87,   0, 0.07), "in"))
# pall2 <- ggarrange(p3, p4, ncol = 2L, labels = c("C", "D"), common.legend = T, legend = "right")
# 
# pall <- ggarrange(pall1, pall2, ncol = 1L, nrow = 2L, heights = c(1, 0.9), common.legend = F)

# pall <- ggarrange(p1, p2, p5, p3, p4, p6, nrow = 2L, ncol = 3L, widths = c(0.72, 0.6, 0.95), labels = "AUTO")
# pall1 <- ggarrange(p1, p2, p3, p4, nrow = 2L, ncol = 2L, widths = c(1, 1), labels = c("A", "B", "D", "E"))
# pall2 <- ggarrange(p5, p6, nrow = 2L, legend = "right", common.legend = T, labels = c("C", "F"))
psub1 <- ggarrange(gt1, gt2, gt3, gt4, nrow = 2L, ncol = 2L, widths = c(1, 1))
psub2 <- ggarrange(gt5, gt6, nrow = 2L)
psub3 <- ggarrange(gt.legend)
pall <- ggarrange(psub1, psub2, psub3, ncol = 3L, widths = c(1, 0.5, 0.25))
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/error_rate_synthetic_data.pdf", plot = pall, width = 7.5, height = 4, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/error_rate_synthetic_data.png", plot = pall, width = 7.5, height = 4, units = "in")
