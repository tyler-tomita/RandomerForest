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

cl <- "rerf"
ind <- 2:length(params[[cl]][[8L]]$sparsity)
length.d <- length(params[[cl]][[8L]]$d)
error.rate <- matrix(0, num.trials, length(ind))
for (trial in 1:num.trials) {
  for (sp.idx in ind) {
    best.oob.error <- which(OOBError[[cl]][2L, 3L, trial, ((sp.idx - 1L)*length.d + 1L):(sp.idx*length.d)] == min(OOBError[[cl]][2L, 3L, trial, ((sp.idx - 1L)*length.d + 1L):(sp.idx*length.d)]))
    if (length(best.oob.error) > 1L) {
      best.oob.auc <- which(OOBAUC[[cl]][2L, 3L, trial, (sp.idx - 1L)*length.d + best.oob.error] == max(OOBAUC[[cl]][2L, 3L, trial, (sp.idx - 1L)*length.d + best.oob.error]))
      if (length(best.oob.auc) > 1L) {
        best.oob.auc <- sample(best.oob.auc, 1L)
      }
      best.oob.error <- best.oob.error[best.oob.auc]
    }
    error.rate[trial, sp.idx - 1L] <- testError[[cl]][2L, 3L, trial, (sp.idx - 1L)*length.d + best.oob.error]
  }
}

df <- data.frame(classifier = rep("RerF", length(ind)), error.mean = apply(error.rate, 2, mean), error.sem = apply(error.rate, 2, function(x) sd(x)/sqrt(length(x))), lambda = 2:5)

cl <- "frc"
ind <- 1:length(params[[cl]][[8L]]$sparsity)
length.d <- length(params[[cl]][[8L]]$d)
error.rate <- matrix(0, num.trials, length(ind))
for (trial in 1:num.trials) {
  for (sp.idx in ind) {
    best.oob.error <- which(OOBError[[cl]][2L, 3L, trial, ((sp.idx - 1L)*length.d + 1L):(sp.idx*length.d)] == min(OOBError[[cl]][2L, 3L, trial, ((sp.idx - 1L)*length.d + 1L):(sp.idx*length.d)]))
    if (length(best.oob.error) > 1L) {
      best.oob.auc <- which(OOBAUC[[cl]][2L, 3L, trial, (sp.idx - 1L)*length.d + best.oob.error] == max(OOBAUC[[cl]][2L, 3L, trial, (sp.idx - 1L)*length.d + best.oob.error]))
      if (length(best.oob.auc) > 1L) {
        best.oob.auc <- sample(best.oob.auc, 1L)
      }
      best.oob.error <- best.oob.error[best.oob.auc]
    }
    error.rate[trial, sp.idx] <- testError[[cl]][2L, 3L, trial, (sp.idx - 1L)*length.d + best.oob.error]
  }
}

df2 <- data.frame(classifier = rep("F-RC", length(ind)), error.mean = apply(error.rate, 2, mean), error.sem = apply(error.rate, 2, function(x) sd(x)/sqrt(length(x))), lambda = 2:5)

df <- rbind(df, df2)
df$classifier <- factor(df$classifier, levels = rev(levels(df$classifier)))


panel.label <- textGrob("(A)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- 2
xmax <- 5
x.range <- xmax - xmin
xmin <- xmin - 0.1*x.range
xmax <- xmax + 0.1*x.range
x.range <- xmax - xmin
label.x <- xmin - 0.45*x.range
ymin <- min(df$error.mean)
ymax <- max(df$error.mean)
y.range <- ymax - ymin
ymin <- ymin - 0.1*y.range
ymax <- ymax + 0.1*y.range
label.y <- ymax + 0.45*y.range

p1 <- ggplot(df, aes(x = lambda, y = error.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  geom_errorbar(aes(ymin=error.mean-error.sem, ymax=error.mean+error.sem, colour = classifier), size = line.width, width = 0.1) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = 2:5, labels = paste0(as.character(2:5), "/p")) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = c(0, 0.04, 0.08)) +
  xlab(expression(lambda)) +
  ylab("Error Rate") +
  ggtitle("Sparse Parity", subtitle = bquote(n == .(n))) +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width),
        text = element_text(size = 14),
        aspect.ratio = 1,
        title = element_text(size = 14),
        legend.position = "bottom") +
  scale_color_manual(values = color.map[c(3, 4)]) +
  guides(color = guide_legend(title = NULL))

gt1 <- ggplot_gtable(ggplot_build(p1))
gt.legend <- gt1$grobs[[which(sapply(gt1$grobs, function(x) x$name == "guide-box"))]]
p1 <- p1 + guides(color = F)
gt1 <- ggplot_gtable(ggplot_build(p1))
gt1$layout$clip[] <- "off"

# Make plots for Orthant

load("~/RandomerForest/R/Results/2017.10.01/Orthant_2017_10_01.RData")

p <- 6L
n <- 400L
num.trials <- dim(testError[[1]])[3]
classifiers <- c("rf", "rerf", "frc")
class.map <- list(rerf = "RerF", rf = "RF", frc = "F-RC")


# color.map <- c("#41ab5d", "#4292c6", "#f768a1")

cl <- "rerf"
ind <- 2:length(params[[cl]][[7L]]$sparsity)
length.d <- length(params[[cl]][[7L]]$d)
error.rate <- matrix(0, num.trials, length(ind))
for (trial in 1:num.trials) {
  for (sp.idx in ind) {
    best.oob.error <- which(OOBError[[cl]][1L, 3L, trial, ((sp.idx - 1L)*length.d + 1L):(sp.idx*length.d)] == min(OOBError[[cl]][1L, 3L, trial, ((sp.idx - 1L)*length.d + 1L):(sp.idx*length.d)]))
    if (length(best.oob.error) > 1L) {
      best.oob.auc <- which(OOBAUC[[cl]][1L, 3L, trial, (sp.idx - 1L)*length.d + best.oob.error] == max(OOBAUC[[cl]][1L, 3L, trial, (sp.idx - 1L)*length.d + best.oob.error]))
      if (length(best.oob.auc) > 1L) {
        best.oob.auc <- sample(best.oob.auc, 1L)
      }
      best.oob.error <- best.oob.error[best.oob.auc]
    }
    error.rate[trial, sp.idx - 1L] <- testError[[cl]][3L, 3L, trial, (sp.idx - 1L)*length.d + best.oob.error]
  }
}

df <- data.frame(classifier = rep("RerF", length(ind)), error.mean = apply(error.rate, 2, mean), error.sem = apply(error.rate, 2, function(x) sd(x)/sqrt(length(x))), lambda = 2:5)

cl <- "frc"
ind <- 1:length(params[[cl]][[7L]]$sparsity)
length.d <- length(params[[cl]][[7L]]$d)
error.rate <- matrix(0, num.trials, length(ind))
for (trial in 1:num.trials) {
  for (sp.idx in ind) {
    best.oob.error <- which(OOBError[[cl]][1L, 3L, trial, ((sp.idx - 1L)*length.d + 1L):(sp.idx*length.d)] == min(OOBError[[cl]][1L, 3L, trial, ((sp.idx - 1L)*length.d + 1L):(sp.idx*length.d)]))
    if (length(best.oob.error) > 1L) {
      best.oob.auc <- which(OOBAUC[[cl]][1L, 3L, trial, (sp.idx - 1L)*length.d + best.oob.error] == max(OOBAUC[[cl]][1L, 3L, trial, (sp.idx - 1L)*length.d + best.oob.error]))
      if (length(best.oob.auc) > 1L) {
        best.oob.auc <- sample(best.oob.auc, 1L)
      }
      best.oob.error <- best.oob.error[best.oob.auc]
    }
    error.rate[trial, sp.idx] <- testError[[cl]][1L, 3L, trial, (sp.idx - 1L)*length.d + best.oob.error]
  }
}

df2 <- data.frame(classifier = rep("F-RC", length(ind)), error.mean = apply(error.rate, 2, mean), error.sem = apply(error.rate, 2, function(x) sd(x)/sqrt(length(x))), lambda = 2:5)

df <- rbind(df, df2)
df$classifier <- factor(df$classifier, levels = rev(levels(df$classifier)))

panel.label <- textGrob("(B)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- 2
xmax <- 5
x.range <- xmax - xmin
xmin <- xmin - 0.1*x.range
xmax <- xmax + 0.1*x.range
x.range <- xmax - xmin
label.x <- xmin - 0.45*x.range
ymin <- min(df$error.mean)
ymax <- max(df$error.mean)
y.range <- ymax - ymin
ymin <- ymin - 0.1*y.range
ymax <- ymax + 0.1*y.range
label.y <- ymax + 0.44*y.range

p2 <- ggplot(df, aes(x = lambda, y = error.mean)) +
  geom_line(aes(colour = classifier), size = line.width) +
  geom_errorbar(aes(ymin=error.mean-error.sem, ymax=error.mean+error.sem, colour = classifier), size = line.width, width = 0.1) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = 2:5, labels = paste0(as.character(2:5), "/p")) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = c(0, 0.2, 0.4), labels = c("0.00", "0.20", "0.40")) +
  xlab(expression(lambda)) +
  ylab("Error Rate") +
  ggtitle("Orthant", subtitle = bquote(n == .(n))) + 
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  guides(color=F) +
  theme_classic() +
  theme(plot.margin = unit(c(0.1, 0, 0.075, 0), "in"), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, title = element_text(size = 14)) +
  scale_color_manual(values = color.map[c(3, 4)])

gt2 <- ggplot_gtable(ggplot_build(p2))
gt2$layout$clip[] <- "off"

psub <- ggarrange(gt1, gt2, ncol = 2L)
pall <- ggarrange(psub, gt.legend, nrow = 2L, heights = c(1, 1/6))
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/sparsity_sensitivity_synthetic_data.pdf", plot = pall, width = 4, height = 2.5, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/sparsity_sensitivity_synthetic_data.png", plot = pall, width = 4, height = 2.5, units = "in")
