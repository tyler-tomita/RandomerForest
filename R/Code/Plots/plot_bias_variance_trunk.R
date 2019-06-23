rm(list = ls())
library(ggplot2)
library(ggpubr)
library(scales)
library(RColorBrewer)
library(grid)

# Make plots for Trunk

load("~/RandomerForest/R/Results/2017.12.10/Trunk_bias_variance_2017_12_10.RData")

ns <- c(10L, 100L, 1000L, 10000L)
length.n <- length(ns)
num.trials <- dim(testError[[1]])[2]
# classifiers <- names(testError)
classifiers <- c("rf", "rerf", "frc", "lda", "md", "bayes")
class.map <- list(rf = "RF", frc = "F-RC", rerf = "RerF", lda = "LDA", md = "Mean-Diff", bayes = "Bayes")

# color.map <- c("#41ab5d", "#4292c6", "#f768a1")
color.map <- brewer.pal(5L, "Set2")
line.width <- 1

classifier <- vector("character", length(classifiers)*length.n*num.trials)
num.obs <- vector("integer", length(classifiers)*length.n*num.trials)
num.dims <- vector("integer", length(classifiers)*length.n*num.trials)
test.error <- vector("double", length(classifiers)*length.n*num.trials)
bias.vec <- vector("double", length(classifiers)*length.n*num.trials)
var.vec <- vector("double", length(classifiers)*length.n*num.trials)
se.vec <- vector("double", length(classifiers)*length.n*num.trials)
ve.vec <- vector("double", length(classifiers)*length.n*num.trials)
for (i in 1:length(classifiers)) {
  cl <- classifiers[i]
  classifier[((i - 1L)*length.n*num.trials + 1L):(i*length.n*num.trials)] <- class.map[[cl]]
  if (cl %in% c("rf", "rerf", "frc")) {
    for (k in 1:length.n) {
      num.obs[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- ns[k]
      test.error[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- sapply(1:num.trials, function(x) ifelse(bestIdx[[cl]][k, x] != 0, testError[[cl]][k, x, bestIdx[[cl]][k, x]], NA))
      bias.vec[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- B[[cl]][k]
      var.vec[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- V[[cl]][k]
      se.vec[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- SE[[cl]][k]
      ve.vec[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- VE[[cl]][k]
    }
  } else if (cl %in% c("lda", "md")) {
    for (k in 1:length.n) {
      num.obs[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- ns[k]
      test.error[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- testError[[cl]][k, ]
      bias.vec[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- B[[cl]][k]
      var.vec[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- V[[cl]][k]
      se.vec[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- SE[[cl]][k]
      ve.vec[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- VE[[cl]][k]
    }
  } else {
    for (k in 1:length.n) {
      num.obs[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- ns[k]
      test.error[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- BE$rf[1L]
      bias.vec[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- 0
      var.vec[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- 0
      se.vec[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- 0
      ve.vec[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- 0 
    }
  }
}

df <- data.frame(classifier = factor(classifier, levels = sapply(names(class.map), function(x) class.map[[x]])))
df$n <- log10(num.obs)
df$error <- test.error
df$bias <- bias.vec
df$var <- var.vec
df$se <- se.vec
df$ve <- ve.vec

df.mean <- aggregate(cbind(error,bias,var,se,ve)~n+classifier, df, mean)
names(df.mean)[3:ncol(df.mean)] <- paste0(names(df)[3:ncol(df)], ".mean")
df.sem <- aggregate(cbind(error,bias,var,se,ve)~n+classifier, df, function(x) sd(x)/sqrt(length(x)))
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
ymin <- ymin - 0.1*y.range
ymax <- ymax + 0.1*y.range
label.y <- ymax + 0.17*y.range

p1 <- ggplot(df.mean, aes(x = n, y = bias.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  # scale_x_log10(breaks = ns, labels = ns/1000L) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = log10(ns), labels = ns/1000L) +
  scale_y_continuous(limits = c(ymin, ymax)) +
  xlab(expression("n"[train]*" (in thousands)")) +
  ylab("Bias") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(plot.margin = unit(c(0.085 +0.1, 0, 0.065, 0), "in"), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = c(color.map[c(2, 3, 1, 4, 5)], "black"))

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
ymin <- ymin - 0.1*y.range
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
  scale_color_manual(values = c(color.map[c(2, 3, 1, 4, 5)], "black"))

gt2 <- ggplot_gtable(ggplot_build(p2))
gt2$layout$clip[] <- "off"

panel.label <- textGrob("(C)", gp=gpar(fontsize=14, fontface="bold"))
label.x <- xmin - 0.42*x.range
ymin <- min(df.mean$se.mean)
ymax <- max(df.mean$se.mean)
y.range <- ymax - ymin
ymin <- ymin - 0.1*y.range
ymax <- ymax + 0.1*y.range
label.y <- ymax + 0.17*y.range

p3 <- ggplot(df.mean, aes(x = n, y = se.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = log10(ns), labels = ns/1000L) +
  scale_y_continuous(limits = c(ymin, ymax)) +
  xlab(NULL) +
  ylab("Sys. Effect") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  guides(color=F) +
  theme_classic() +
  theme(plot.margin = unit(c(0.085 + 0.1, 0, 0.37, 0), "in"), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = c(color.map[c(2, 3, 1, 4, 5)], "black"))

gt3 <- ggplot_gtable(ggplot_build(p3))
gt3$layout$clip[] <- "off"

panel.label <- textGrob("(D)", gp=gpar(fontsize=14, fontface="bold"))
label.x <- xmin - 0.42*x.range
ymin <- min(df.mean$ve.mean)
ymax <- max(df.mean$ve.mean)
y.range <- ymax - ymin
ymin <- ymin - 0.1*y.range
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
  scale_color_manual(values = c(color.map[c(2, 3, 1, 4, 5)], "black"))

gt4 <- ggplot_gtable(ggplot_build(p4))
gt4$layout$clip[] <- "off"

panel.label <- textGrob("(E)", gp=gpar(fontsize=14, fontface="bold"))
label.x <- xmin - 0.42*x.range
ymin <- min(df.mean$error.mean)
ymax <- max(df.mean$error.mean)
y.range <- ymax - ymin
ymin <- ymin - 0.1*y.range
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
  scale_color_manual(values = c(color.map[c(2, 3, 1, 4, 5)], "black"))

gt5 <- ggplot_gtable(ggplot_build(p5))
gt5$layout$clip[] <- "off"


pall <- ggarrange(gt1, gt2, gt3, gt4, gt5, gt.legend, ncol = 6L, widths = c(1, 1, 1, 1, 1, 0.5))
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/bias_variance_trunk.pdf", plot = pall, width = 13, height = 2, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/bias_variance_trunk.png", plot = pall, width = 13, height = 2, units = "in")
