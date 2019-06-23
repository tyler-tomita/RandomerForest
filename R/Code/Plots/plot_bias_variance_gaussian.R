rm(list = ls())
library(ggplot2)
library(ggpubr)
library(scales)
library(RColorBrewer)

# Make plots for Guassian Orthogonal

load("~/RandomerForest/R/Results/2017.10.15/Gaussian_orthogonal_2017_10_15_aggregated.RData")

length.n <- dim(testError[[1]])[1]
length.p <- dim(testError[[1]])[2]
num.trials <- dim(testError[[1]])[3]
# classifiers <- names(testError)
classifiers <- c("rf", "rerf", "frc")
class.map <- list(rerf = "RerF", rf = "RF", frc = "F-RC")

ps <- c(2L, 4L)
ns <- c(5L, 10L, 100L, 1000L)

color.map <- c("#41ab5d", "#4292c6", "#f768a1")
line.width <- 1.5

classifier <- vector("character", length(classifiers)*length.n*length.p*num.trials)
num.obs <- vector("integer", length(classifiers)*length.n*length.p*num.trials)
num.dims <- vector("integer", length(classifiers)*length.n*length.p*num.trials)
test.error <- vector("double", length(classifiers)*length.n*length.p*num.trials)
train.time <- vector("double", length(classifiers)*length.n*length.p*num.trials)
tree.strength <- vector("double", length(classifiers)*length.n*length.p*num.trials)
tree.correlation <- vector("double", length(classifiers)*length.n*length.p*num.trials)
bias.vec <- vector("double", length(classifiers)*length.n*length.p*num.trials)
var.vec <- vector("double", length(classifiers)*length.n*length.p*num.trials)
se.vec <- vector("double", length(classifiers)*length.n*length.p*num.trials)
ve.vec <- vector("double", length(classifiers)*length.n*length.p*num.trials)
be.vec <- vector("double", length(classifiers)*length.n*length.p*num.trials)
for (i in 1:length(classifiers)) {
  cl <- classifiers[i]
  classifier[((i - 1L)*length.n*length.p*num.trials + 1L):(i*length.n*length.p*num.trials)] <- class.map[[cl]]
  for (j in 1:length.p) {
    num.dims[(((i - 1L)*length.p + (j - 1L))*length.n*num.trials + 1L):(((i - 1L)*length.p + j)*length.n*num.trials)] <- ps[j]
    for (k in 1:length.n) {
      num.obs[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- ns[k]
      test.error[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, j, x] != 0, testError[[cl]][k, j, x, bestIdx[[cl]][k, j, x]], NA))
      train.time[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, j, x] != 0, trainTime[[cl]][k, j, x, bestIdx[[cl]][k, j, x]], NA))
      tree.strength[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, j, x] != 0, treeStrength[[cl]][k, j, x, bestIdx[[cl]][k, j, x]], NA))
      tree.correlation[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, j, x] != 0, treeCorr[[cl]][k, j, x, bestIdx[[cl]][k, j, x]], NA))
      bias.vec[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- B[[cl]][k, j]
      var.vec[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- V[[cl]][k, j]
      se.vec[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- SE[[cl]][k, j]
      ve.vec[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- VE[[cl]][k, j]
      be.vec[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- BE[[cl]][k, j]
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
df$bias <- bias.vec
df$var <- var.vec
df$se <- se.vec
df$ve <- ve.vec
df$be <- be.vec

df.mean <- aggregate(cbind(error,train.time,tree.strength,tree.correlation,bias,var,se,ve,be)~n+p+classifier, df, mean)
names(df.mean)[4:ncol(df.mean)] <- paste0(names(df)[4:ncol(df)], ".mean")
df.sem <- aggregate(cbind(error,train.time,tree.strength,tree.correlation,bias,var,se,ve,be)~n+p+classifier, df, function(x) sd(x)/sqrt(length(x)))
names(df.sem)[4:ncol(df.sem)] <- paste0(names(df)[4:ncol(df)], ".sem")
df.mean <- cbind(df.mean, df.sem[4:ncol(df.sem)])
# levels(df.mean$classifier) <- classifiers

p1 <- ggplot(df.mean[df.mean$p == 2L, ], aes(x = n, y = bias.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  scale_x_log10(breaks = ns[[1L]], labels = function(x) format(x, scientific = T)) +
  xlab(expression("n"[train])) +
  ylab("bias") +
  ggtitle("Orthogonal Gaussian") +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14)) +
  scale_color_manual(values = color.map)

p2 <- ggplot(df.mean[df.mean$p == 2L, ], aes(x = n, y = var.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  scale_x_log10(breaks = ns[[2L]], labels = function(x) format(x, scientific = T)) +
  xlab(expression("n"[train])) +
  ylab("variance") +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14)) +
  scale_color_manual(values = color.map)

p3 <- ggplot(df.mean[df.mean$p == 2L, ], aes(x = n, y = se.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  scale_x_log10(breaks = ns[[2L]], labels = function(x) format(x, scientific = T)) +
  xlab(expression("n"[train])) +
  ylab("systematic effect") +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14)) +
  scale_color_manual(values = color.map)

p4 <- ggplot(df.mean[df.mean$p == 2L, ], aes(x = n, y = ve.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  scale_x_log10(breaks = ns[[2L]], labels = function(x) format(x, scientific = T)) +
  xlab(expression("n"[train])) +
  ylab("variance effect") +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14)) +
  scale_color_manual(values = color.map)

df2 <- df.mean[df.mean$classifier == "RF", ]
levels(df2$classifier) <- c(levels(df2$classifier), "Bayes")
df2$classifier[] <- "Bayes"
df2$error.mean <- df2$be.mean
df2$error.sem <- df2$be.sem
df2 <- rbind(df.mean, df2)

p5 <- ggplot(df.mean[df.mean$p == 2L, ], aes(x = n, y = error.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  geom_errorbar(aes(ymin=error.mean-error.sem, ymax=error.mean+error.sem, colour = classifier), size = line.width, width = 0.1) +
  scale_x_log10(breaks = ns[[1L]], labels = function(x) format(x, scientific = T)) +
  xlab(expression("n"[train])) +
  ylab("error rate") +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14)) +
  scale_color_manual(values = c(color.map, "black"))

load("~/RandomerForest/R/Results/2017.10.15/Gaussian_oblique_2017_10_15_aggregated.RData")

classifier <- vector("character", length(classifiers)*length.n*length.p*num.trials)
num.obs <- vector("integer", length(classifiers)*length.n*length.p*num.trials)
num.dims <- vector("integer", length(classifiers)*length.n*length.p*num.trials)
test.error <- vector("double", length(classifiers)*length.n*length.p*num.trials)
train.time <- vector("double", length(classifiers)*length.n*length.p*num.trials)
tree.strength <- vector("double", length(classifiers)*length.n*length.p*num.trials)
tree.correlation <- vector("double", length(classifiers)*length.n*length.p*num.trials)
bias.vec <- vector("double", length(classifiers)*length.n*length.p*num.trials)
var.vec <- vector("double", length(classifiers)*length.n*length.p*num.trials)
se.vec <- vector("double", length(classifiers)*length.n*length.p*num.trials)
ve.vec <- vector("double", length(classifiers)*length.n*length.p*num.trials)
be.vec <- vector("double", length(classifiers)*length.n*length.p*num.trials)

for (i in 1:length(classifiers)) {
  cl <- classifiers[i]
  classifier[((i - 1L)*length.n*length.p*num.trials + 1L):(i*length.n*length.p*num.trials)] <- class.map[[cl]]
  for (j in 1:length.p) {
    num.dims[(((i - 1L)*length.p + (j - 1L))*length.n*num.trials + 1L):(((i - 1L)*length.p + j)*length.n*num.trials)] <- ps[j]
    for (k in 1:length.n) {
      num.obs[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- ns[k]
      test.error[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, j, x] != 0, testError[[cl]][k, j, x, bestIdx[[cl]][k, j, x]], NA))
      train.time[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, j, x] != 0, trainTime[[cl]][k, j, x, bestIdx[[cl]][k, j, x]], NA))
      tree.strength[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, j, x] != 0, treeStrength[[cl]][k, j, x, bestIdx[[cl]][k, j, x]], NA))
      tree.correlation[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, j, x] != 0, treeCorr[[cl]][k, j, x, bestIdx[[cl]][k, j, x]], NA))
      bias.vec[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- B[[cl]][k, j]
      var.vec[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- V[[cl]][k, j]
      se.vec[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- SE[[cl]][k, j]
      ve.vec[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- VE[[cl]][k, j]
      be.vec[((((i - 1L)*length.p + (j - 1L))*length.n + (k - 1L))*num.trials + 1L):((((i - 1L)*length.p + (j - 1L))*length.n + k)*num.trials)] <- BE[[cl]][k, j]
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
df$bias <- bias.vec
df$var <- var.vec
df$se <- se.vec
df$ve <- ve.vec
df$be <- be.vec

df.mean <- aggregate(cbind(error,train.time,tree.strength,tree.correlation,bias,var,se,ve,be)~n+p+classifier, df, mean)
names(df.mean)[4:ncol(df.mean)] <- paste0(names(df)[4:ncol(df)], ".mean")
df.sem <- aggregate(cbind(error,train.time,tree.strength,tree.correlation,bias,var,se,ve,be)~n+p+classifier, df, function(x) sd(x)/sqrt(length(x)))
names(df.sem)[4:ncol(df.sem)] <- paste0(names(df)[4:ncol(df)], ".sem")
df.mean <- cbind(df.mean, df.sem[4:ncol(df.sem)])
# levels(df.mean$classifier) <- classifiers

p6 <- ggplot(df.mean[df.mean$p == 2L, ], aes(x = n, y = bias.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  scale_x_log10(breaks = ns[[1L]], labels = function(x) format(x, scientific = T)) +
  xlab(expression("n"[train])) +
  ylab("bias") +
  ggtitle("Oblique Gaussian") +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14)) +
  scale_color_manual(values = color.map)

p7 <- ggplot(df.mean[df.mean$p == 2L, ], aes(x = n, y = var.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  scale_x_log10(breaks = ns[[2L]], labels = function(x) format(x, scientific = T)) +
  xlab(expression("n"[train])) +
  ylab("variance") +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14)) +
  scale_color_manual(values = color.map)

p8 <- ggplot(df.mean[df.mean$p == 2L, ], aes(x = n, y = se.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  scale_x_log10(breaks = ns[[2L]], labels = function(x) format(x, scientific = T)) +
  xlab(expression("n"[train])) +
  ylab("systematic effect") +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14)) +
  scale_color_manual(values = color.map)

p9 <- ggplot(df.mean[df.mean$p == 2L, ], aes(x = n, y = ve.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  scale_x_log10(breaks = ns[[2L]], labels = function(x) format(x, scientific = T)) +
  xlab(expression("n"[train])) +
  ylab("variance effect") +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14)) +
  scale_color_manual(values = color.map)

df2 <- df.mean[df.mean$classifier == "RF", ]
levels(df2$classifier) <- c(levels(df2$classifier), "Bayes")
df2$classifier[] <- "Bayes"
df2$error.mean <- df2$be.mean
df2$error.sem <- df2$be.sem
df2 <- rbind(df.mean, df2)

p10 <- ggplot(df.mean[df.mean$p == 2L, ], aes(x = n, y = error.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  geom_errorbar(aes(ymin=error.mean-error.sem, ymax=error.mean+error.sem, colour = classifier), size = line.width, width = 0.1) +
  scale_x_log10(breaks = ns[[1L]], labels = function(x) format(x, scientific = T)) +
  xlab(expression("n"[train])) +
  ylab("error rate") +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14)) +
  scale_color_manual(values = c(color.map, "black"))

pall <- ggarrange(p1, p6, p2, p7, p3, p8, p4, p9, p5, p10, ncol = 2L, nrow = 5L, common.legend = T, widths = c(1, 1), heights = c(1, 0.85, 0.85, 0.85, 0.85), legend = "right", labels = "AUTO")

ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/bias_variance_synthetic_data.pdf", plot = pall, width = 11, height = 25, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/bias_variance_synthetic_data.png", plot = pall, width = 11, height = 25, units = "in")
