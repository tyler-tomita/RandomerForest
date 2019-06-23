rm(list = ls())
library(ggplot2)
library(ggpubr)
library(scales)
library(RColorBrewer)
library(grid)

# Make plots for Sparse parity

load("~/RandomerForest/R/Results/2017.09.27/Sparse_parity_2017_09_27.RData")

p <- 20L
n <- 5000L
num.trials <- dim(testError[[1]])[3]
# classifiers <- names(testError)
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
line.width <- 0.75
marker.size <- 0.75

cl <- "rf"
ind <- 1:length(params[[cl]][[8L]]$d)
length.d <- length(ind)
error.rate <- matrix(0, num.trials, length(ind))
for (trial in 1:num.trials) {
  for (d.idx in ind) {
    error.rate[trial, d.idx] <- testError[[cl]][2L, 3L, trial, d.idx]
  }
}
df <- data.frame(classifier = rep("RF", length(ind)), error.mean = apply(error.rate, 2, mean), error.sem = apply(error.rate, 2, function(x) sd(x)/sqrt(length(x))), d = params[[cl]][[8L]]$d)

cl <- "rerf"
ind <- 1:length(params[[cl]][[8L]]$d)
length.d <- length(ind)
length.sparsity <- length(params[[cl]][[8L]]$sparsity)
error.rate <- matrix(0, num.trials, length(ind))
for (trial in 1:num.trials) {
  for (d.idx in ind) {
    best.oob.error <- which(OOBError[[cl]][2L, 3L, trial, ((1:length.sparsity) - 1L)*length.d + d.idx] == min(OOBError[[cl]][2L, 3L, trial, ((1:length.sparsity) - 1L)*length.d + d.idx]))
    if (length(best.oob.error) > 1L) {
      best.oob.auc <- which(OOBAUC[[cl]][2L, 3L, trial, (best.oob.error - 1L)*length.d + d.idx] == max(OOBAUC[[cl]][2L, 3L, trial, (best.oob.error - 1L)*length.d + d.idx]))
      if (length(best.oob.auc) > 1L) {
        best.oob.auc <- sample(best.oob.auc, 1L)
      }
      best.oob.error <- best.oob.error[best.oob.auc]
    }
    error.rate[trial, d.idx] <- testError[[cl]][2L, 3L, trial, (best.oob.error - 1L)*length.d + d.idx]
  }
}

df2 <- data.frame(classifier = rep("RerF", length(ind)), error.mean = apply(error.rate, 2, mean), error.sem = apply(error.rate, 2, function(x) sd(x)/sqrt(length(x))), d = params[[cl]][[8L]]$d)

cl <- "frc"
ind <- 1:length(params[[cl]][[8L]]$d)
length.d <- length(ind)
length.sparsity <- length(params[[cl]][[8L]]$sparsity)
error.rate <- matrix(0, num.trials, length(ind))
for (trial in 1:num.trials) {
  for (d.idx in ind) {
    best.oob.error <- which(OOBError[[cl]][2L, 3L, trial, ((1:length.sparsity) - 1L)*length.d + d.idx] == min(OOBError[[cl]][2L, 3L, trial, ((1:length.sparsity) - 1L)*length.d + d.idx]))
    if (length(best.oob.error) > 1L) {
      best.oob.auc <- which(OOBAUC[[cl]][2L, 3L, trial, (best.oob.error - 1L)*length.d + d.idx] == max(OOBAUC[[cl]][2L, 3L, trial, (best.oob.error - 1L)*length.d + d.idx]))
      if (length(best.oob.auc) > 1L) {
        best.oob.auc <- sample(best.oob.auc, 1L)
      }
      best.oob.error <- best.oob.error[best.oob.auc]
    }
    error.rate[trial, d.idx] <- testError[[cl]][2L, 3L, trial, (best.oob.error - 1L)*length.d + d.idx]
  }
}

df3 <- data.frame(classifier = rep("F-RC", length(ind)), error.mean = apply(error.rate, 2, mean), error.sem = apply(error.rate, 2, function(x) sd(x)/sqrt(length(x))), d = params[[cl]][[8L]]$d)

df <- rbind(df, df2, df3)
df$classifier <- factor(df$classifier, levels = c("RF", "F-RC", "RerF"))
df$d <- log10(df$d)

panel.label <- textGrob("(C)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- min(df$d)
xmax <- max(df$d)
x.range <- xmax - xmin
xmin <- xmin - 0.1*x.range
xmax <- xmax + 0.1*x.range
x.range <- xmax - xmin
label.x <- xmin - 0.2*x.range
ymin <- min(df$error.mean)
ymax <- max(df$error.mean)
y.range <- ymax - ymin
ymin <- ymin - 0.1*y.range
ymax <- ymax + 0.1*y.range
# label.y <- ymax + 0.16*y.range
label.y <- ymax

p1 <- ggplot(df, aes(x = d, y = error.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  geom_errorbar(aes(ymin=error.mean-error.sem, ymax=error.mean+error.sem, colour = classifier), size = line.width, width = 0.1) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = log10(c(5, 20, 400)), labels = c("5", "20", "400")) +
  scale_y_continuous(limits = c(ymin, ymax)) +
  # scale_x_log10(breaks = c(5, 20, 400)) +
  xlab("d") +
  ylab("Error Rate") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  # ggtitle("Sparse Parity", subtitle = bquote(n[train] == .(n))) +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0.1, 0, 0.15), "in"),
        axis.line = element_line(size = line.width),
        text = element_text(size = 14),
        aspect.ratio = 1.5,
        title = element_text(size = 14),
        legend.position = "bottom") +
  scale_color_manual(values = color.map[2:4])

gt1 <- ggplot_gtable(ggplot_build(p1))
gt.legend <- gt1$grobs[[which(sapply(gt1$grobs, function(x) x$name == "guide-box"))]]
p1 <- p1 + guides(color = F)
gt1 <- ggplot_gtable(ggplot_build(p1))
gt1$layout$clip[] <- "off"

cl <- "rf"
ind <- 1:length(params[[cl]][[8L]]$d)
length.d <- length(ind)
train.time <- matrix(0, num.trials, length(ind))
for (trial in 1:num.trials) {
  for (d.idx in ind) {
    train.time[trial, d.idx] <- trainTime[[cl]][2L, 3L, trial, d.idx]
  }
}
df <- data.frame(classifier = rep("RF", length(ind)), time.mean = apply(log10(train.time), 2, mean), time.sem = apply(log10(train.time), 2, function(x) sd(x)/sqrt(length(x))), d = params[[cl]][[8L]]$d)

cl <- "rerf"
ind <- 1:length(params[[cl]][[8L]]$d)
length.d <- length(ind)
length.sparsity <- length(params[[cl]][[8L]]$sparsity)
train.time <- matrix(0, num.trials, length(ind))
for (trial in 1:num.trials) {
  for (d.idx in ind) {
    best.oob.error <- which(OOBError[[cl]][2L, 3L, trial, ((1:length.sparsity) - 1L)*length.d + d.idx] == min(OOBError[[cl]][2L, 3L, trial, ((1:length.sparsity) - 1L)*length.d + d.idx]))
    if (length(best.oob.error) > 1L) {
      best.oob.auc <- which(OOBAUC[[cl]][2L, 3L, trial, (best.oob.error - 1L)*length.d + d.idx] == max(OOBAUC[[cl]][2L, 3L, trial, (best.oob.error - 1L)*length.d + d.idx]))
      if (length(best.oob.auc) > 1L) {
        best.oob.auc <- sample(best.oob.auc, 1L)
      }
      best.oob.error <- best.oob.error[best.oob.auc]
    }
    train.time[trial, d.idx] <- trainTime[[cl]][2L, 3L, trial, (best.oob.error - 1L)*length.d + d.idx]
  }
}

df2 <- data.frame(classifier = rep("RerF", length(ind)), time.mean = apply(log10(train.time), 2, mean), time.sem = apply(log10(train.time), 2, function(x) sd(x)/sqrt(length(x))), d = params[[cl]][[8L]]$d)

cl <- "frc"
ind <- 1:length(params[[cl]][[8L]]$d)
length.d <- length(ind)
length.sparsity <- length(params[[cl]][[8L]]$sparsity)
train.time <- matrix(0, num.trials, length(ind))
for (trial in 1:num.trials) {
  for (d.idx in ind) {
    best.oob.error <- which(OOBError[[cl]][2L, 3L, trial, ((1:length.sparsity) - 1L)*length.d + d.idx] == min(OOBError[[cl]][2L, 3L, trial, ((1:length.sparsity) - 1L)*length.d + d.idx]))
    if (length(best.oob.error) > 1L) {
      best.oob.auc <- which(OOBAUC[[cl]][2L, 3L, trial, (best.oob.error - 1L)*length.d + d.idx] == max(OOBAUC[[cl]][2L, 3L, trial, (best.oob.error - 1L)*length.d + d.idx]))
      if (length(best.oob.auc) > 1L) {
        best.oob.auc <- sample(best.oob.auc, 1L)
      }
      best.oob.error <- best.oob.error[best.oob.auc]
    }
    train.time[trial, d.idx] <- trainTime[[cl]][2L, 3L, trial, (best.oob.error - 1L)*length.d + d.idx]
  }
}

df3 <- data.frame(classifier = rep("F-RC", length(ind)), time.mean = apply(log10(train.time), 2, mean), time.sem = apply(log10(train.time), 2, function(x) sd(x)/sqrt(length(x))), d = params[[cl]][[8L]]$d)

df <- rbind(df, df2, df3)
df$classifier <- factor(df$classifier, levels = c("RF", "F-RC", "RerF"))
df$d <- log10(df$d)

panel.label <- textGrob("(B)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- min(df$d)
xmax <- max(df$d)
x.range <- xmax - xmin
xmin <- xmin - 0.1*x.range
xmax <- xmax + 0.1*x.range
x.range <- xmax - xmin
label.x <- xmin - 0.2*x.range
ymin <- min(df$time.mean)
ymax <- max(df$time.mean)
y.range <- ymax - ymin
ymin <- ymin - 0.1*y.range
ymax <- ymax + 0.1*y.range
# label.y <- ymax + 0.16*y.range
label.y <- ymax

p2 <- ggplot(df, aes(x = d, y = time.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  geom_errorbar(aes(ymin=time.mean-time.sem, ymax=time.mean+time.sem, colour = classifier), size = line.width, width = 0.1) +
  # scale_x_log10(breaks = c(5, 20, 400)) +
  # scale_y_log10(breaks = c(5, 10^(1:2)), limits = c(5, 10^2.5)) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = log10(c(5, 20, 400)), labels = c("5", "20", "400")) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = log10(c(5, 10, 100)), labels = c("5", "10", "100")) +
  xlab("d") +
  ylab("Training Time (s)") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  guides(color=F) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1.5, title = element_text(size = 14)) +
  scale_color_manual(values = color.map[2:4])

gt2 <- ggplot_gtable(ggplot_build(p2))
gt2$layout$clip[] <- "off"

load("~/RandomerForest/R/Results/2017.09.27/Sparse_parity_2017_09_27.RData")

length.n <- dim(testError[[1]])[1]
length.p <- dim(testError[[1]])[2]
num.trials <- dim(testError[[1]])[3]
classifiers <- c("rf", "rerf", "frc")
class.map <- list(frc = "F-RC", rerf = "RerF", rf = "RF")

ps <- c(3, 10, 20)
ns <- list(c(10, 100, 1000), c(100, 1000, 10000), c(1000, 5000, 10000))

# color.map <- c("#41ab5d", "#4292c6", "#f768a1")
# color.map <- brewer.pal(5L, "Set2")
line.width <- 0.75

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

df.mean <- aggregate(cbind(error,train.time,tree.strength,tree.correlation)~classifier+p+n, df, mean)
names(df.mean)[4:ncol(df.mean)] <- paste0(names(df)[4:ncol(df)], ".mean")
df.sem <- aggregate(cbind(error,train.time,tree.strength,tree.correlation)~classifier+p+n, df, function(x) sd(x)/sqrt(length(x)))
names(df.sem)[4:ncol(df.sem)] <- paste0(names(df)[4:ncol(df)], ".sem")
df.mean <- cbind(df.mean, df.sem[4:ncol(df.sem)])
df.mean$n <- log10(df.mean$n)

panel.label <- textGrob("(A)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- min(df.mean$n[df.mean$p == 20L])
xmax <- max(df.mean$n[df.mean$p == 20L])
x.range <- xmax - xmin
xmin <- xmin - 0.1*x.range
xmax <- xmax + 0.1*x.range
x.range <- xmax - xmin
label.x <- xmin - 0.2*x.range
ymin <- min(df.mean$train.time.mean[df.mean$p == 20L])
ymax <- max(df.mean$train.time.mean[df.mean$p == 20L])
y.range <- ymax - ymin
ymin <- ymin - 0.1*y.range
ymax <- ymax + 0.1*y.range
# label.y <- ymax + 0.16*y.range
label.y <- ymax

p3 <- ggplot(df.mean[df.mean$p == 20L, ], aes(x = n, y = train.time.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  geom_errorbar(aes(ymin=train.time.mean-train.time.sem, ymax=train.time.mean+train.time.sem, colour = classifier), size = line.width, width = 0.1) +
  # scale_x_log10(breaks = ns[[3L]], labels = ns[[3L]]/1000L) +
  # scale_x_log10(breaks = ns[[3L]], labels = function(x) format(x, scientific = T)) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = log10(c(1000, 5000, 10000)), labels = c("1", "5", "10")) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = c(0, 200, 400), labels = c("0", "200", "400")) +
  xlab("n (in thousands)") +
  ylab("Training Time (s)") +
  # ggtitle("Sparse Parity") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  guides(color=F) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1.5) +
  scale_color_manual(values = color.map[2:4])

gt3 <- ggplot_gtable(ggplot_build(p3))
gt3$layout$clip[] <- "off"

# Make plots for Orthant

# load("~/RandomerForest/R/Results/2017.10.01/Orthant_2017_10_01.RData")
# 
# length.n <- dim(testError[[1]])[1]
# length.p <- dim(testError[[1]])[2]
# num.trials <- dim(testError[[1]])[3]
# ps <- c(2, 4, 6)
# ns <- list(c(20,200,400), c(80,400,4000), c(400,2000,4000))
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
#       test.error[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) testError[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
#       train.time[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) trainTime[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
#       tree.strength[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) treeStrength[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
#       tree.correlation[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) treeCorr[[cl]][k, j, x, bestIdx[[cl]][k, j, x]])
#     }
#   }
# }
# 
# df <- data.frame(classifier = factor(classifier, levels = sapply(names(class.map), function(x) class.map[[x]])))
# df$p <- num.dims
# df$n <- num.obs
# df$error <- test.error
# df$train.time <- train.time
# df$tree.strength <- tree.strength
# df$tree.correlation <- tree.correlation
# 
# df.mean <- aggregate(cbind(error,train.time,tree.strength,tree.correlation)~classifier+p+n, df, mean)
# names(df.mean)[4:ncol(df.mean)] <- paste0(names(df)[4:ncol(df)], ".mean")
# df.sem <- aggregate(cbind(error,train.time,tree.strength,tree.correlation)~classifier+p+n, df, function(x) sd(x)/sqrt(length(x)))
# names(df.sem)[4:ncol(df.sem)] <- paste0(names(df)[4:ncol(df)], ".sem")
# df.mean <- cbind(df.mean, df.sem[4:ncol(df.sem)])
# df.mean$n <- log10(df.mean$n)

# panel.label <- textGrob("(B)", gp=gpar(fontsize=14, fontface="bold"))
# xmin <- min(df.mean$n[df.mean$p == 6L])
# xmax <- max(df.mean$n[df.mean$p == 6L])
# x.range <- xmax - xmin
# xmin <- xmin - 0.1*x.range
# xmax <- xmax + 0.1*x.range
# x.range <- xmax - xmin
# label.x <- xmin - 0.15*x.range
# ymin <- min(df.mean$train.time.mean[df.mean$p == 6L])
# ymax <- max(df.mean$train.time.mean[df.mean$p == 6L])
# y.range <- ymax - ymin
# ymin <- ymin - 0.1*y.range
# ymax <- ymax + 0.1*y.range
# label.y <- ymax + 0.09*y.range
# 
# p4 <- ggplot(df.mean[df.mean$p == 6L, ], aes(x = n, y = train.time.mean)) +
#   geom_line(aes(colour = classifier), size = line.width) +
#   geom_errorbar(aes(ymin=train.time.mean-train.time.sem, ymax=train.time.mean+train.time.sem, colour = classifier), size = line.width, width = 0.1) +
#   # scale_x_log10(breaks = ns[[3L]], labels = ns[[3L]]/1000L) +
#   # scale_x_log10(breaks = ns[[3L]], labels = function(x) format(x, scientific = T)) +
#   scale_x_continuous(limits = c(xmin, xmax), breaks = log10(c(400, 2000, 4000)), labels = c("0.4", "2", "4")) +
#   xlab(expression(n[train]*"(in thousands)")) +
#   ylab("Training Time (s)") +
#   ggtitle("Orthant") +
#   annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
#   guides(color=F) +
#   theme_classic() +
#   theme(plot.margin = unit(c(0, 0.1, 0, 0.15), "in"), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1.5) +
#   scale_color_manual(values = color.map[c(3, 1, 2)])
# 
# gt4 <- ggplot_gtable(ggplot_build(p4))
# gt4$layout$clip[] <- "off"
# 
# ptext <- ggplot() + 
#   scale_x_continuous(limits = c(0, 1)) +
#   scale_y_continuous(limits = c(0, 1)) +
#   annotate("text", x = 0.5, y = 0.65, size=6, label = "Sparse Parity") + 
#   annotate("text", x = 0.5, y = 0.2, size=6, label = "n[train] == 5000", parse = T) + 
#   theme_classic() +
#   theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())

# pall <- ggarrange(p3, p4, ncol = 2L, common.legend = T, legend = "right", labels = "AUTO")

# psub1 <- ggarrange(gt3, gt4, ncol = 2L)
# psub2 <- ggarrange(gt2, gt1, ncol = 2L)
# psub3 <- ggarrange(psub1, ptext, psub2, nrow = 3L, heights = c(1, 0.35, 1))
# pall <- ggarrange(psub3, gt.legend, ncol = 2L, widths = c(1, 0.16))
psub <- ggarrange(gt3, gt2, gt1, ncol = 3L, widths = c(1, 1, 1))
pall <- ggarrange(psub, gt.legend, nrow = 2L, heights = c(1, 1/9))
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/d_gt_p.pdf", plot = pall, width = 5.75, height = 2.5, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/d_gt_p.png", plot = pall, width = 5.75, height = 2.5, units = "in")
