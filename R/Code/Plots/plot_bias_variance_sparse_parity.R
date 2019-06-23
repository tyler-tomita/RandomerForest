rm(list = ls())
library(ggplot2)
library(ggpubr)
library(scales)
library(RColorBrewer)
library(grid)

# Make plots for sparse parity

load("~/RandomerForest/R/Results/2017.10.31/Sparse_parity_bias_variance_2017_10_31.RData")

ns <- c(500L, 1000L, 3000L, 5000L)
length.n <- length(ns)
num.trials <- dim(testError[[1]])[2]
# classifiers <- names(testError)
classifiers <- c("rf", "rerf", "frc")
class.map <- list(rf = "RF", frc = "F-RC", rerf = "RerF")

# color.map <- c("#41ab5d", "#4292c6", "#f768a1")
# color.map <- brewer.pal(5L, "Set2")
color.map <- c("CCF" = "#6a51a3",
               "RF" = "#4a484c",
               "F-RC" = "#318dde",
               "RerF" = "#F41711",
               "XGBoost" = "#E78AC3")
line.width <- 1

classifier <- vector("character", length(classifiers)*length.n*num.trials)
num.obs <- vector("integer", length(classifiers)*length.n*num.trials)
num.dims <- vector("integer", length(classifiers)*length.n*num.trials)
test.error <- vector("double", length(classifiers)*length.n*num.trials)
train.time <- vector("double", length(classifiers)*length.n*num.trials)
tree.strength <- vector("double", length(classifiers)*length.n*num.trials)
tree.correlation <- vector("double", length(classifiers)*length.n*num.trials)
bias.vec <- vector("double", length(classifiers)*length.n*num.trials)
var.vec <- vector("double", length(classifiers)*length.n*num.trials)
se.vec <- vector("double", length(classifiers)*length.n*num.trials)
ve.vec <- vector("double", length(classifiers)*length.n*num.trials)
be.vec <- vector("double", length(classifiers)*length.n*num.trials)
for (i in 1:length(classifiers)) {
  cl <- classifiers[i]
  classifier[((i - 1L)*length.n*num.trials + 1L):(i*length.n*num.trials)] <- class.map[[cl]]
  for (k in 1:length.n) {
    num.obs[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- ns[k]
    test.error[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, x] != 0, testError[[cl]][k, x, bestIdx[[cl]][k, x]], NA))
    train.time[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, x] != 0, trainTime[[cl]][k, x, bestIdx[[cl]][k, x]], NA))
    tree.strength[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, x] != 0, treeStrength[[cl]][k, x, bestIdx[[cl]][k, x]], NA))
    tree.correlation[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, x] != 0, treeCorr[[cl]][k, x, bestIdx[[cl]][k, x]], NA))
    bias.vec[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- B[[cl]][k]
    var.vec[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- V[[cl]][k]
    se.vec[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- SE[[cl]][k]
    ve.vec[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- VE[[cl]][k]
    be.vec[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- BE[[cl]][k]
  }
}

df <- data.frame(classifier = factor(classifier, levels = sapply(names(class.map), function(x) class.map[[x]])))
df$n <- log10(num.obs)
df$error <- test.error
df$train.time <- train.time
df$tree.strength <- tree.strength
df$tree.correlation <- tree.correlation
df$bias <- bias.vec
df$var <- var.vec
df$se <- se.vec
df$ve <- ve.vec
df$be <- be.vec

df.mean <- aggregate(cbind(error,train.time,tree.strength,tree.correlation,bias,var,se,ve,be)~n+classifier, df, mean)
names(df.mean)[3:ncol(df.mean)] <- paste0(names(df)[3:ncol(df)], ".mean")
df.sem <- aggregate(cbind(error,train.time,tree.strength,tree.correlation,bias,var,se,ve,be)~n+classifier, df, function(x) sd(x)/sqrt(length(x)))
names(df.sem)[3:ncol(df.sem)] <- paste0(names(df)[3:ncol(df)], ".sem")
df.mean <- cbind(df.mean, df.sem[3:ncol(df.sem)])
# levels(df.mean$classifier) <- classifiers

panel.label <- textGrob("(A)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- log10(ns[1L])
xmax <- log10(ns[3L])
x.range <- xmax - xmin
xmin <- xmin - 0.1*x.range
xmax <- xmax + 0.1*x.range
x.range <- xmax - xmin
label.x <- xmin - 0.42*x.range
ymin <- min(df.mean$bias.mean)
ymax <- max(df.mean$bias.mean)
y.range <- ymax - ymin
# ymin <- ymin - 0.1*y.range
ymin <- 0
ymax <- ymax + 0.1*y.range
label.y <- ymax + 0.17*y.range

p1 <- ggplot(df.mean, aes(x = n, y = bias.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  # scale_x_log10(breaks = ns, labels = ns/1000L) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = log10(ns), labels = ns/1000L) +
  scale_y_continuous(limits = c(ymin, ymax)) +
  xlab("n (in thousands)") +
  ylab("Bias") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(plot.margin = unit(c(0.085 +0.1, 0, 0.065, 0), "in"), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map[c(2, 3, 4)])

gt1 <- ggplot_gtable(ggplot_build(p1))
gt.legend <- gt1$grobs[[which(sapply(gt1$grobs, function(x) x$name == "guide-box"))]]
p1 <- p1 + guides(color = F)
gt1 <- ggplot_gtable(ggplot_build(p1))
gt1$layout$clip[] <- "off"

panel.label <- textGrob("(B)", gp=gpar(fontsize=14, fontface="bold"))
label.x <- xmin - 0.42*x.range
ymin <- min(df.mean$var.mean)
ymax <- max(df.mean$var.mean)
y.range <- ymax - ymin
# ymin <- ymin - 0.1*y.range
ymin <- 0
ymax <- ymax + 0.1*y.range
label.y <- ymax + 0.17*y.range

p2 <- ggplot(df.mean, aes(x = n, y = var.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  # scale_x_log10(breaks = ns, labels = ns/1000L) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = log10(ns), labels = ns/1000L) +
  scale_y_continuous(limits = c(ymin, ymax)) +
  xlab(NULL) +
  ylab("Variance") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  guides(color=F) +
  theme_classic() +
  theme(plot.margin = unit(c(0.085 + 0.1, 0, 0.37, 0), "in"), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map[c(2, 3, 4)])

gt2 <- ggplot_gtable(ggplot_build(p2))
gt2$layout$clip[] <- "off"

p3 <- ggplot(df.mean, aes(x = n, y = se.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  scale_x_log10(breaks = ns, labels = ns/1000L) +
  xlab(NULL) +
  ylab("Sys. Effect") +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(plot.margin = unit(c(0.085, 0, 0.37, 0), "in"), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map[c(2, 3, 4)])

panel.label <- textGrob("(C)", gp=gpar(fontsize=14, fontface="bold"))
label.x <- xmin - 0.42*x.range
ymin <- min(df.mean$ve.mean)
ymax <- max(df.mean$ve.mean)
y.range <- ymax - ymin
# ymin <- ymin - 0.1*y.range
ymin <- 0
ymax <- ymax + 0.1*y.range
label.y <- ymax + 0.17*y.range

p4 <- ggplot(df.mean, aes(x = n, y = ve.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  # scale_x_log10(breaks = ns, labels = ns/1000L) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = log10(ns), labels = ns/1000L) +
  scale_y_continuous(limits = c(ymin, ymax)) +
  xlab(NULL) +
  ylab("Var. Effect") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  guides(color=F) +
  theme_classic() +
  theme(plot.margin = unit(c(0.085 + 0.1, 0, 0.37, 0), "in"), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map[c(2, 3, 4)])

gt4 <- ggplot_gtable(ggplot_build(p4))
gt4$layout$clip[] <- "off"

panel.label <- textGrob("(D)", gp=gpar(fontsize=14, fontface="bold"))
label.x <- xmin - 0.42*x.range
ymin <- min(df.mean$error.mean)
ymax <- max(df.mean$error.mean)
y.range <- ymax - ymin
# ymin <- ymin - 0.1*y.range
ymin <- 0
ymax <- ymax + 0.1*y.range
label.y <- ymax + 0.17*y.range

p5 <- ggplot(df.mean, aes(x = n, y = error.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  geom_errorbar(aes(ymin=error.mean-error.sem, ymax=error.mean+error.sem, colour = classifier), size = line.width, width = 0.1) +
  # scale_x_log10(breaks = ns, labels = ns/1000L) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = log10(ns), labels = ns/1000L) +
  scale_y_continuous(limits = c(ymin, ymax)) +
  xlab(NULL) +
  ylab("Error Rate") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  guides(color=F) +
  theme_classic() +
  theme(plot.margin = unit(c(0.085 + 0.1, 0, 0.37, 0), "in"), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map[c(2, 3, 4)])

gt5 <- ggplot_gtable(ggplot_build(p5))
gt5$layout$clip[] <- "off"


pall <- ggarrange(gt1, gt2, gt4, gt5, gt.legend, ncol = 5L, widths = c(1, 1, 1, 1, 0.5))
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/bias_variance_synthetic_data.pdf", plot = pall, width = 10.5, height = 2, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/bias_variance_synthetic_data.png", plot = pall, width = 10.5, height = 2, units = "in")
