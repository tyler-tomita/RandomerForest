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
class.map <- list(frc = "F-RC", rerf = "RerF", rf = "RF")

ps <- c(3, 10, 20)
ns <- list(c(10, 100, 1000), c(100, 1000, 10000), c(1000, 5000, 10000))

# color.map <- c("#41ab5d", "#4292c6", "#f768a1")
color.map <- brewer.pal(5L, "Set2")
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
# levels(df.mean$classifier) <- classifiers

# p1 <- ggplot(df.mean[df.mean$p == 3L, ], aes(x = n, y = train.time.mean)) +
#   geom_line(aes(colour = classifier), size = line.width) +
#   geom_errorbar(aes(ymin=train.time.mean-train.time.sem, ymax=train.time.mean+train.time.sem, colour = classifier), size = line.width, width = 0.1) +
#   scale_x_log10(breaks = ns[[1L]], labels = function(x) format(x, scientific = T)) +
#   xlab(expression("n"[train])) +
#   # annotation_custom(grob = textGrob("Training Time (sec)", rot = 90), xmin = -0.1, xmax = -0.1, ymin = )
#   ylab(expression(atop(bold("Sparse Parity"),"Training Time (sec)"))) +
#   ggtitle("p = 3") +
#   guides(color=guide_legend(title="Algorithm")) +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
#   scale_color_manual(values = color.map)
# 
# p2 <- ggplot(df.mean[df.mean$p == 10L, ], aes(x = n, y = train.time.mean)) +
#   geom_line(aes(colour = classifier), size = line.width) +
#   geom_errorbar(aes(ymin=train.time.mean-train.time.sem, ymax=train.time.mean+train.time.sem, colour = classifier), size = line.width, width = 0.1) +
#   scale_x_log10(breaks = ns[[2L]], labels = function(x) format(x, scientific = T)) +
#   xlab(expression("n"[train])) +
#   ylab(expression(atop(bold(" "),"Training Time (sec)"))) +
#   ggtitle("p = 10") +
#   guides(color=guide_legend(title="Algorithm")) +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
#   scale_color_manual(values = color.map)

p3 <- ggplot(df.mean[df.mean$p == 20L, ], aes(x = n, y = train.time.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  geom_errorbar(aes(ymin=train.time.mean-train.time.sem, ymax=train.time.mean+train.time.sem, colour = classifier), size = line.width, width = 0.1) +
  scale_x_log10(breaks = ns[[3L]], labels = ns[[3L]]/1000L) +
  # scale_x_log10(breaks = ns[[3L]], labels = function(x) format(x, scientific = T)) +
  xlab(expression("n"[train]%*%10^-3)) +
  ylab("Training Time (sec)") +
  ggtitle("Sparse Parity") +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map[c(3, 1, 2)])

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

df.mean <- aggregate(cbind(error,train.time,tree.strength,tree.correlation)~classifier+p+n, df, mean)
names(df.mean)[4:ncol(df.mean)] <- paste0(names(df)[4:ncol(df)], ".mean")
df.sem <- aggregate(cbind(error,train.time,tree.strength,tree.correlation)~classifier+p+n, df, function(x) sd(x)/sqrt(length(x)))
names(df.sem)[4:ncol(df.sem)] <- paste0(names(df)[4:ncol(df)], ".sem")
df.mean <- cbind(df.mean, df.sem[4:ncol(df.sem)])
# levels(df.mean$classifier) <- classifiers

# p4 <- ggplot(df.mean[df.mean$p == 2L, ], aes(x = n, y = train.time.mean)) +
#   geom_line(aes(colour = classifier), size = line.width) +
#   geom_errorbar(aes(ymin=train.time.mean-train.time.sem, ymax=train.time.mean+train.time.sem, colour = classifier), size = line.width, width = 0.1) +
#   scale_x_log10(breaks = ns[[1L]], labels = function(x) format(x, scientific = T)) +
#   xlab(expression("n"[train])) +
#   ylab(expression(atop(bold("Orthant"),"Training Time (sec)"))) +
#   ggtitle("p = 2") +
#   guides(color=guide_legend(title="Algorithm")) +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
#   scale_color_manual(values = color.map)
# 
# p5 <- ggplot(df.mean[df.mean$p == 4L, ], aes(x = n, y = train.time.mean)) +
#   geom_line(aes(colour = classifier), size = line.width) +
#   geom_errorbar(aes(ymin=train.time.mean-train.time.sem, ymax=train.time.mean+train.time.sem, colour = classifier), size = line.width, width = 0.1) +
#   scale_x_log10(breaks = ns[[2L]], labels = function(x) format(x, scientific = T)) +
#   xlab(expression("n"[train])) +
#   ylab(expression(atop(bold(" "),"Training Time (sec)"))) +
#   ggtitle("p = 4") +
#   guides(color=guide_legend(title="Algorithm")) +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
#   scale_color_manual(values = color.map)

p6 <- ggplot(df.mean[df.mean$p == 6L, ], aes(x = n, y = train.time.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  geom_errorbar(aes(ymin=train.time.mean-train.time.sem, ymax=train.time.mean+train.time.sem, colour = classifier), size = line.width, width = 0.1) +
  scale_x_log10(breaks = ns[[3L]], labels = ns[[3L]]/1000L) +
  # scale_x_log10(breaks = ns[[3L]], labels = function(x) format(x, scientific = T)) +
  xlab(expression("n"[train]%*%10^-3)) +
  ylab("Training Time (sec)") +
  ggtitle("Orthant") +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map[c(3, 1, 2)])

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
# df$n <- num.obs
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
# p7 <- ggplot(df.mean[df.mean$p == 10L, ], aes(x = n, y = train.time.mean)) +
#   geom_line(aes(colour = classifier), size = line.width) +
#   geom_errorbar(aes(ymin=train.time.mean-train.time.sem, ymax=train.time.mean+train.time.sem, colour = classifier), size = line.width, width = 0.1) +
#   scale_x_log10(breaks = ns[[1L]], labels = function(x) format(x, scientific = T)) +
#   xlab(expression("n"[train])) +
#   # annotation_custom(grob = textGrob("Training Time (sec)", rot = 90), xmin = -0.1, xmax = -0.1, ymin = )
#   ylab(expression(atop(bold("Trunk"),"Training Time (sec)"))) +
#   ggtitle("p = 10") +
#   guides(color=guide_legend(title="Algorithm")) +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
#   scale_color_manual(values = color.map)
# 
# p8 <- ggplot(df.mean[df.mean$p == 100L, ], aes(x = n, y = train.time.mean)) +
#   geom_line(aes(colour = classifier), size = line.width) +
#   geom_errorbar(aes(ymin=train.time.mean-train.time.sem, ymax=train.time.mean+train.time.sem, colour = classifier), size = line.width, width = 0.1) +
#   scale_x_log10(breaks = ns[[2L]], labels = function(x) format(x, scientific = T)) +
#   xlab(expression("n"[train])) +
#   ylab(expression(atop(bold(" "),"Training Time (sec)"))) +
#   ggtitle("p = 100") +
#   guides(color=guide_legend(title="Algorithm")) +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
#   scale_color_manual(values = color.map)
# 
# p9 <- ggplot(df.mean[df.mean$p == 1000L, ], aes(x = n, y = train.time.mean)) +
#   geom_line(aes(colour = classifier), size = line.width) +
#   geom_errorbar(aes(ymin=train.time.mean-train.time.sem, ymax=train.time.mean+train.time.sem, colour = classifier), size = line.width, width = 0.1) +
#   scale_x_log10(breaks = ns[[3L]], labels = function(x) format(x, scientific = T)) +
#   xlab(expression("n"[train])) +
#   ylab(expression(atop(bold(" "),"Training Time (sec)"))) +
#   ggtitle("p = 1000") +
#   guides(color=guide_legend(title="Algorithm")) +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
#   scale_color_manual(values = color.map)

# pall <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3L, nrow = 3L, common.legend = T, widths = c(1, 0.95, 0.95), legend = "right", labels = "AUTO")
pall <- ggarrange(p3, p6, ncol = 2L, common.legend = T, legend = "right", labels = "AUTO")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/train_time_synthetic_data.pdf", plot = pall, width = 5.5, height = 2.5, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/train_time_synthetic_data.png", plot = pall, width = 5.5, height = 2.5, units = "in")
