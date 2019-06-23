rm(list = ls())
library(rerf)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)

VectorAngle <- function(u, v) {
  return(theta <- acos(u%*%v/sqrt(sum(u^2))/sqrt(sum(v^2)))*180/pi)
}

line.width <- 1
marker.size <- 0.75
# color.map <- c("#41ab5d", "#4292c6", "#f768a1")
color.map <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")

# sparse oblique hyperplanes

n <- 2000L
p <- 20L
d <- ceiling(p^(3))
num.trials <- 50L
sep <- 0.05
noise <- 0.01

# generate points for 2d plot
p.plot <- 2L
n.plot <- 100L

X1 <- matrix(0, n.plot/2, p.plot)
X1[, -p.plot] <- runif(n = n.plot/2*(p.plot - 1L), min = -1/2, max = 1/2)
X1[, p.plot] <- -apply(X1[, -p.plot, drop = F], 1, sum)
X1 <- X1 - t(sapply(rnorm(n.plot/2, mean = sep, sd = noise), function(x) rep(x/sqrt(p.plot), p.plot)))
X2 <- matrix(0, n.plot/2, p.plot)
X2[, -p.plot] <- runif(n = n.plot/2*(p.plot - 1L), min = -1/2, max = 1/2)
X2[, p.plot] <- -apply(X2[, -p.plot, drop = F], 1, sum)
X2 <- X2 + t(sapply(rnorm(n.plot/2, mean = sep, sd = noise), function(x) rep(x/sqrt(p.plot), p.plot)))
X <- rbind(X1, X2)
Y <- as.factor(c(rep(0, n.plot/2), rep(1, n.plot/2)))

df.plot <- data.frame(x1 = X[, 1L], x2 = X[, 2L], class = Y)


panel.label <- textGrob("(A)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- min(df.plot$x1)
xmax <- max(df.plot$x1)
x.range <- xmax - xmin
# xmin <- xmin - 0.1*x.range
# xmax <- xmax + 0.1*x.range
x.range <- xmax - xmin
label.x <- xmin - 0.3*x.range
ymin <- min(df.plot$x2)
ymax <- max(df.plot$x2)
y.range <- ymax - ymin
# ymin <- ymin - 0.1*y.range
# ymax <- ymax + 0.1*y.range
label.y <- ymax + 0.8*y.range

p1 <- ggplot(df.plot, aes(x = x1, y = x2, color = class)) +
  geom_point(size = marker.size) +
  scale_x_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, 0, 0.5)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, 0, 0.5)) +
  xlab(expression(X[1])) +
  ylab(expression(X[2])) +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  ggtitle("Parallel\nHyperplanes") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, title = element_text(size = 14)) +
  scale_color_manual(values = color.map[c(3, 4)]) +
  guides(color = F)

gt1 <- ggplot_gtable(ggplot_build(p1))
p1 <- p1 + guides(color = F)
gt1 <- ggplot_gtable(ggplot_build(p1))
gt1$layout$clip[] <- "off"

# generate points for experiment

X1 <- matrix(0, n/2, p)
X1[, -2L] <- runif(n = n/2*(p - 1L), min = -1/2, max = 1/2)
X1[, 2L] <- -X1[, 1L]
X1[, 1:2] <- X1[, 1:2] - t(sapply(rnorm(n/2, mean = sep, sd = noise), function(x) rep(x/sqrt(2), 2L)))
X2 <- matrix(0, n/2, p)
X2[, -2L] <- runif(n = n/2*(p - 1L), min = -1/2, max = 1/2)
X2[, 2L] <- -X2[, 1L]
X2[, 1:2] <- X2[, 1:2] + t(sapply(rnorm(n/2, mean = sep, sd = noise), function(x) rep(x/sqrt(2L), 2L)))
X <- rbind(X1, X2)
Y <- c(rep(0, n/2), rep(1, n/2))


# RerF
Xp <- matrix(0, n, num.trials)
for (trial in 1:num.trials) {
  forest <- RerF(X, Y, max.depth = 2L, trees = 1L, bagging = 0, mat.options = list(p=p, d=d, random.matrix="binary", 1/p, 1/2))
  indexHigh <- forest$trees[[1L]]$matAindex[2L]
  indexLow <- forest$trees[[1L]]$matAindex[1L] + 1L
  s <- (indexHigh - indexLow + 1L)/2L
  projection <- double(p)
  projection[forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L-1L]] <- forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L]
  if (VectorAngle(rep(1, p), projection) <= 90) {
    Xp[, trial] <- X[, forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L-1L], drop = F]%*%forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L]
  } else {
    Xp[, trial] <- X[, forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L-1L], drop = F]%*%-forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L]
  }
}
Xp <- apply(Xp, 2, function(x) (x - min(x))/(max(x) - min(x)))
df <- data.frame(projection = c(Xp), class = as.factor(rep(Y, num.trials)), random.matrix = rep("RerF", n*num.trials))

# RR-RF
Xp <- matrix(0, n, num.trials)
for (trial in 1:num.trials) {
  forest <- RerF(X, Y, max.depth = 2L, trees = 1L, bagging = 0, mat.options = list(p=p, d=d, random.matrix="continuous", 1))
  indexHigh <- forest$trees[[1L]]$matAindex[2L]
  indexLow <- forest$trees[[1L]]$matAindex[1L] + 1L
  s <- (indexHigh - indexLow + 1L)/2L
  # projection <- forest$trees[[1L]]$rotmat[, forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L-1L]]
  projection <- double(p)
  projection[forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L-1L]] <- forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L]
  if (VectorAngle(rep(1, p), projection) <= 90) {
    Xp[, trial] <- X%*%projection
  } else {
    Xp[, trial] <- X%*%-projection
  }
}
Xp <- apply(Xp, 2, function(x) (x - min(x))/(max(x) - min(x)))
df <- rbind(df, data.frame(projection = c(Xp), class = as.factor(rep(Y, num.trials)), random.matrix = rep("RR-RF", n*num.trials)))

# RF
Xp <- matrix(0, n, num.trials)
for (trial in 1:num.trials) {
  forest <- RerF(X, Y, max.depth = 2L, trees = 1L, bagging = 0, mat.options = list(p=p, d=p, random.matrix="rf"), rfPack = FALSE)
  indexHigh <- forest$trees[[1L]]$matAindex[2L]
  indexLow <- forest$trees[[1L]]$matAindex[1L] + 1L
  s <- (indexHigh - indexLow + 1L)/2L
  projection <- double(p)
  projection[forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L-1L]] <- forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L]
  if (VectorAngle(rep(1, p), projection) <= 90) {
    Xp[, trial] <- X[, forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L-1L], drop = F]%*%forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L]
  } else {
    Xp[, trial] <- X[, forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L-1L], drop = F]%*%-forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L]
  }
}
Xp <- apply(Xp, 2, function(x) (x - min(x))/(max(x) - min(x)))
df <- rbind(df, data.frame(projection = c(Xp), class = as.factor(rep(Y, num.trials)), random.matrix = rep("RF", n*num.trials)))

p2 <- ggplot(df[df$random.matrix == "RerF", ], aes(x = projection, color = class)) +
  # geom_density(adjust = 2, size = line.width) +
  stat_density(aes(fill = class), geom = "area", size = line.width, adjust = 3, position = "identity", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("RerF") +
  theme_classic() +
  theme(plot.margin = unit(c(0.36, 0, 0.35, 0), "in"), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, title = element_text(size = 14)) +
  scale_fill_manual(values = color.map[c(3, 4)]) +
  scale_color_manual(values = color.map[c(3, 4)]) +
  guides(color = F)

p3 <- ggplot(df[df$random.matrix == "RR-RF", ], aes(x = projection, color = class)) +
  # geom_density(adjust = 2, size = line.width) +
  stat_density(aes(fill = class), geom = "area", size = line.width, adjust = 3, position = "identity", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("RR-RF") +
  theme_classic() +
  theme(plot.margin = unit(c(0.36, 0, 0.35, 0), "in"), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, title = element_text(size = 14)) +
  scale_fill_manual(values = color.map[c(3, 4)]) +
  scale_color_manual(values = color.map[c(3, 4)]) +
  guides(color = F)

p4 <- ggplot(df[df$random.matrix == "RF", ], aes(x = projection, color = class)) +
  # geom_density(adjust = 2, size = line.width) +
  stat_density(aes(fill = class), geom = "area", size = line.width, adjust = 3, position = "identity", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  xlab("Best Projection") +
  ylab(expression(atop(bold("Sparse"), "KDE"))) +
  ggtitle("RF") +
  theme_classic() +
  theme(plot.margin = unit(c(0.36, 0, 0.1, 0), "in"), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, title = element_text(size = 14)) +
  scale_fill_manual(values = color.map[c(3, 4)]) +
  scale_color_manual(values = color.map[c(3, 4)]) +
  guides(color = F)

# dense oblique hyperplanes

n <- 2000L
p <- 20L
d <- ceiling(p^(3))
num.trials <- 50L
sep <- 0.05
noise <- 0.01

# generate points for 2d plot
# p.plot <- 3L
# n.plot <- 50L
# 
# X1 <- matrix(0, n.plot/2, p.plot)
# X1[, -p.plot] <- runif(n = n.plot/2*(p.plot - 1L), min = -1/2, max = 1/2)
# X1[, p.plot] <- -apply(X1[, -p.plot, drop = F], 1, sum)
# X1 <- X1 - t(sapply(rnorm(n.plot/2, mean = sep, sd = noise), function(x) rep(x/sqrt(p.plot), p.plot)))
# X2 <- matrix(0, n.plot/2, p.plot)
# X2[, -p.plot] <- runif(n = n.plot/2*(p.plot - 1L), min = -1/2, max = 1/2)
# X2[, p.plot] <- -apply(X2[, -p.plot, drop = F], 1, sum)
# X2 <- X2 + t(sapply(rnorm(n.plot/2, mean = sep, sd = noise), function(x) rep(x/sqrt(p.plot), p.plot)))
# X <- rbind(X1, X2)
# Y <- as.factor(c(rep(0, n.plot/2), rep(1, n.plot/2)))
# 
# df.plot <- data.frame(x1 = X[, 1L], x2 = X[, 2L], class = Y)
# 
# p1 <- ggplot(df.plot, aes(x = x1, y = x2, color = class)) +
#   geom_point(size = marker.size) +
#   scale_x_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, 0, 0.5)) +
#   scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, 0, 0.5)) +
#   xlab(expression(X[1])) +
#   ylab(expression(X[2])) +
#   ggtitle("Parallel\nHyperplanes") +
#   theme_classic() +
#   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, title = element_text(size = 14)) +
#   scale_color_manual(values = color.map[c(3, 4)]) 

# generate points for experiment
X1 <- matrix(0, n/2, p)
X1[, -p] <- runif(n = n/2*(p - 1L), min = -1/2, max = 1/2)
X1[, p] <- -apply(X1[, -p, drop = F], 1, sum)
X1 <- X1 - t(sapply(rnorm(n/2, mean = sep, sd = noise), function(x) rep(x/sqrt(p), p)))
X2 <- matrix(0, n/2, p)
X2[, -p] <- runif(n = n/2*(p - 1L), min = -1/2, max = 1/2)
X2[, p] <- -apply(X2[, -p, drop = F], 1, sum)
X2 <- X2 + t(sapply(rnorm(n/2, mean = sep, sd = noise), function(x) rep(x/sqrt(p), p)))
X <- rbind(X1, X2)
Y <- c(rep(0, n/2), rep(1, n/2))

# RerF
Xp <- matrix(0, n, num.trials)
for (trial in 1:num.trials) {
  forest <- RerF(X, Y, max.depth = 2L, trees = 1L, bagging = 0, mat.options = list(p=p, d=d, random.matrix="binary", 1/p, 1/2))
  indexHigh <- forest$trees[[1L]]$matAindex[2L]
  indexLow <- forest$trees[[1L]]$matAindex[1L] + 1L
  s <- (indexHigh - indexLow + 1L)/2L
  projection <- double(p)
  projection[forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L-1L]] <- forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L]
  if (VectorAngle(rep(1, p), projection) <= 90) {
    Xp[, trial] <- X[, forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L-1L], drop = F]%*%forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L]
  } else {
    Xp[, trial] <- X[, forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L-1L], drop = F]%*%-forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L]
  }
}
Xp <- apply(Xp, 2, function(x) (x - min(x))/(max(x) - min(x)))
df <- data.frame(projection = c(Xp), class = as.factor(rep(Y, num.trials)), random.matrix = rep("RerF", n*num.trials))

# RR-RF
Xp <- matrix(0, n, num.trials)
for (trial in 1:num.trials) {
  forest <- RerF(X, Y, max.depth = 2L, trees = 1L, bagging = 0, mat.options = list(p=p, d=d, random.matrix="continuous", 1))
  indexHigh <- forest$trees[[1L]]$matAindex[2L]
  indexLow <- forest$trees[[1L]]$matAindex[1L] + 1L
  s <- (indexHigh - indexLow + 1L)/2L
  # projection <- forest$trees[[1L]]$rotmat[, forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L-1L]]
  projection <- double(p)
  projection[forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L-1L]] <- forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L]
  if (VectorAngle(rep(1, p), projection) <= 90) {
    Xp[, trial] <- X%*%projection
  } else {
    Xp[, trial] <- X%*%-projection
  }
}
Xp <- apply(Xp, 2, function(x) (x - min(x))/(max(x) - min(x)))
df <- rbind(df, data.frame(projection = c(Xp), class = as.factor(rep(Y, num.trials)), random.matrix = rep("RR-RF", n*num.trials)))

Xp <- matrix(0, n, num.trials)
for (trial in 1:num.trials) {
  forest <- RerF(X, Y, max.depth = 2L, trees = 1L, bagging = 0, mat.options = list(p=p, d=p, random.matrix="rf"), rfPack = FALSE)
  indexHigh <- forest$trees[[1L]]$matAindex[2L]
  indexLow <- forest$trees[[1L]]$matAindex[1L] + 1L
  s <- (indexHigh - indexLow + 1L)/2L
  projection <- double(p)
  projection[forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L-1L]] <- forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L]
  if (VectorAngle(rep(1, p), projection) <= 90) {
    Xp[, trial] <- X[, forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L-1L], drop = F]%*%forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L]
  } else {
    Xp[, trial] <- X[, forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L-1L], drop = F]%*%-forest$trees[[1L]]$matAstore[indexLow:indexHigh][(1:s)*2L]
  }
}
Xp <- apply(Xp, 2, function(x) (x - min(x))/(max(x) - min(x)))
df <- rbind(df, data.frame(projection = c(Xp), class = as.factor(rep(Y, num.trials)), random.matrix = rep("RF", n*num.trials)))

p5 <- ggplot(df[df$random.matrix == "RerF", ], aes(x = projection, color = class)) +
  # geom_density(adjust = 2, size = line.width) +
  stat_density(aes(fill = class), geom = "area", size = line.width, adjust = 3, position = "identity", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("RerF") +
  theme_classic() +
  theme(plot.margin = unit(c(0.36, 0, 0.35, 0), "in"), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, title = element_text(size = 14)) +
  scale_fill_manual(values = color.map[c(3, 4)]) +
  scale_color_manual(values = color.map[c(3, 4)]) +
  guides(color = F)

p6 <- ggplot(df[df$random.matrix == "RR-RF", ], aes(x = projection, color = class)) +
  # geom_density(adjust = 2, size = line.width) +
  stat_density(aes(fill = class), geom = "area", size = line.width, adjust = 3, position = "identity", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("RR-RF") +
  theme_classic() +
  theme(plot.margin = unit(c(0.36, 0, 0.35, 0), "in"), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, title = element_text(size = 14)) +
  scale_fill_manual(values = color.map[c(3, 4)]) +
  scale_color_manual(values = color.map[c(3, 4)]) +
  guides(color = F)

p7 <- ggplot(df[df$random.matrix == "RF", ], aes(x = projection, color = class)) +
  # geom_density(adjust = 2, size = line.width) +
  stat_density(aes(fill = class), geom = "area", size = line.width, adjust = 3, position = "identity", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  xlab("Best Projection") +
  ylab(expression(atop(bold("Dense"), "KDE"))) +
  ggtitle("RF") +
  theme_classic() +
  theme(plot.margin = unit(c(0.36, 0, 0.1, 0), "in"), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, title = element_text(size = 14)) +
  scale_fill_manual(values = color.map[c(3, 4)]) +
  scale_color_manual(values = color.map[c(3, 4)]) +
  guides(color = F)

# p.scatter <- ggarrange(p1, labels = "(A)")
p.hist <- ggarrange(p4, p3, p2, p7, p6, p5, nrow = 2L, ncol = 3L, common.legend = T, widths = c(1, 0.9, 0.9), legend = "right", labels = paste0("(", LETTERS[2:7], ")"))

pall <- ggarrange(gt1, p.hist, ncol = 2L, common.legend = F, widths = c(0.25, 1))
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/pami_fig1_intuition2.pdf", plot = pall, width = 8, height = 4, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/pami_fig1_intuition2.png", plot = pall, width = 8, height = 4, units = "in")
