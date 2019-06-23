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
marker.size <- 2

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

p1 <- ggplot(df[(df$p ==4L) & (df$n == 5L), ], aes(x = tree.correlation, y = tree.strength)) +
  geom_point(aes(colour = classifier), size = marker.size) +
  theme(axis.line = element_line(size = 1.5), text = element_text(size = 14)) +
  xlab("tree correlation") +
  ylab("tree strength") +
  ggtitle("Orthogonal Gaussian", subtitle = expression('p = 4, n'[train]*' = 5')) +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.title = element_text(face = "bold")) +
  scale_color_manual(values = color.map)

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

p2 <- ggplot(df[(df$p == 4L) & (df$n == 10L), ], aes(x = tree.correlation, y = tree.strength)) +
  geom_point(aes(colour = classifier), size = marker.size) +
  theme(axis.line = element_line(size = 1.5), text = element_text(size = 14)) +
  xlab("tree correlation") +
  ylab("tree strength") +
  ggtitle("Oblique Gaussian", subtitle = expression('p = 4, n'[train]*' = 10')) +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.title = element_text(face = "bold")) +
  scale_color_manual(values = color.map)

pall <- ggarrange(p1, p2, ncol = 2L, nrow = 1L, common.legend = T, widths = c(1, 1), legend = "right", labels = "AUTO")

ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/strength_vs_correlation_gaussian.pdf", plot = pall, width = 11, height = 5, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/strength_vs_correlation_gaussian.png", plot = pall, width = 11, height = 5, units = "in")
