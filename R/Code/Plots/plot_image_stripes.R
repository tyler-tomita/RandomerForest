library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)

ntrain <- 10L
trial <- 1L
dataPath <- "~/R/Data/Image_stripes/"
D <- as.matrix(read.table(file = paste0(dataPath, "Train/Image_stripes_train_set_n", ntrain, "_trial", trial, ".csv"), header = F, sep = ",", quote = "", row.names = NULL))
p <- ncol(D) - 1L
Xtrain <- D[, 1:p]
Ytrain <- as.integer(D[, p + 1L])

Xi <- melt(matrix(Xtrain[1, ], 20L, 20L))

p1 <- ggplot(Xi, aes(x = Var2, y = Var1, fill = factor(value))) +
  geom_raster() +
  scale_fill_manual(breaks = levels(factor(Xi$value)),
                    values = c("#A9A9A9", "black")) +
  theme(plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        aspect.ratio = 1) +
  guides(fill = F) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

Xi <- melt(matrix(Xtrain[2, ], 20L, 20L))

p2 <- ggplot(Xi, aes(x = Var2, y = Var1, fill = factor(value))) +
  geom_raster() +
  scale_fill_manual(breaks = levels(factor(Xi$value)),
                    values = c("#A9A9A9", "black")) +
  theme(plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        aspect.ratio = 1) +
  guides(fill = F) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

Xi <- melt(matrix(Xtrain[3, ], 20L, 20L))

p3 <- ggplot(Xi, aes(x = Var2, y = Var1, fill = factor(value))) +
  geom_raster() +
  scale_fill_manual(breaks = levels(factor(Xi$value)),
                    values = c("#A9A9A9", "black")) +
  theme(plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        aspect.ratio = 1) +
  guides(fill = F) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

Xi <- melt(matrix(Xtrain[6, ], 20L, 20L))

p4 <- ggplot(Xi, aes(x = Var2, y = Var1, fill = factor(value))) +
  geom_raster() +
  scale_fill_manual(breaks = levels(factor(Xi$value)),
                    values = c("#A9A9A9", "black")) +
  theme(plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        aspect.ratio = 1) +
  guides(fill = F) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

Xi <- melt(matrix(Xtrain[7, ], 20L, 20L))

p5 <- ggplot(Xi, aes(x = Var2, y = Var1, fill = factor(value))) +
  geom_raster() +
  scale_fill_manual(breaks = levels(factor(Xi$value)),
                    values = c("#A9A9A9", "black")) +
  theme(plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        aspect.ratio = 1) +
  guides(fill = F) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

Xi <- melt(matrix(Xtrain[8, ], 20L, 20L))

p6 <- ggplot(Xi, aes(x = Var2, y = Var1, fill = factor(value))) +
  geom_raster() +
  scale_fill_manual(breaks = levels(factor(Xi$value)),
                    values = c("#A9A9A9", "black")) +
  theme(plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        aspect.ratio = 1) +
  guides(fill = F) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

# gt = ggplotGrob(p1)

# Select plot panel only
#   gt = gt[6,4]    # Using index notation; OR
# gt = gtable::gtable_filter(gt, "panel")  

# Draw it
# grid.newpage()
# grid.draw(gt)

pall1 <- grid.arrange(p1, p2, p3, ncol = 3, top = textGrob("Class 0", gp = gpar(fontsize = 16)))
pall2 <- grid.arrange(p4, p5, p6, ncol = 3, top = textGrob("Class 1", gp = gpar(fontsize = 16)))
pall3 <- ggarrange(pall1, pall2, ncol = 1L, nrow = 2L)


load("~/RandomerForest/R/Results/2018.01.02/Image_stripes_2018_01_02.RData")

length.n <- dim(testError[[1]])[2]
num.trials <- dim(testError[[1]])[1]
classifiers <- c("rf", "rerf", "strerf", "control")
class.map <- list(rerf = "RerF", rf = "RF", strerf = "S-RerF", control = "Control")

ns <- c(10L, 20L, 50L)

# color.map <- c("#41ab5d", "#4292c6", "#f768a1")
# color.map <- brewer.pal(5L, "Set2")
color.map <- c("CCF" = "#6a51a3",
               "RF" = "#4a484c",
               "F-RC" = "#318dde",
               "RerF" = "#F41711",
               "XGBoost" = "#E78AC3")
color.map2 <- brewer.pal(8L, "Dark2")
# color.map2 <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")[c(2L, 3L, 1L, 4L)]
line.width <- 1
marker.size <- 1

classifier <- vector("character", length(classifiers)*length.n*num.trials)
num.obs <- vector("integer", length(classifiers)*length.n*num.trials)
test.error <- vector("double", length(classifiers)*length.n*num.trials)
train.time <- vector("double", length(classifiers)*length.n*num.trials)
tree.strength <- vector("double", length(classifiers)*length.n*num.trials)
tree.correlation <- vector("double", length(classifiers)*length.n*num.trials)
for (i in 1:length(classifiers)) {
  cl <- classifiers[i]
  classifier[((i - 1L)*length.n*num.trials + 1L):(i*length.n*num.trials)] <- class.map[[cl]]
  for (k in 1:length.n) {
    num.obs[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- ns[k]
    test.error[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- testError[[cl]][, k]
    train.time[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- trainTime[[cl]][, k]
    tree.strength[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- treeStrength[[cl]][, k]
    tree.correlation[(((i - 1L)*length.n + (k - 1L))*num.trials + 1L):(((i - 1L)*length.n + k)*num.trials)] <- treeCorr[[cl]][, k]
  }
}

df <- data.frame(classifier = factor(classifier, levels = sapply(names(class.map), function(x) class.map[[x]])))
df$n <- num.obs
df$error <- test.error
df$train.time <- train.time
df$tree.strength <- tree.strength
df$tree.correlation <- tree.correlation

df.mean <- aggregate(cbind(error,train.time,tree.strength,tree.correlation)~classifier+n, df, mean)
names(df.mean)[3:ncol(df.mean)] <- paste0(names(df)[3:ncol(df)], ".mean")
df.sem <- aggregate(cbind(error,train.time,tree.strength,tree.correlation)~classifier+n, df, function(x) sd(x)/sqrt(length(x)))
names(df.sem)[3:ncol(df.sem)] <- paste0(names(df)[3:ncol(df)], ".sem")
df.mean <- cbind(df.mean, df.sem[3:ncol(df.sem)])
df.mean$classifier <- factor(df.mean$classifier, levels = c("Control", "RF", "RerF", "S-RerF"))

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

panel.label <- textGrob("(B)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- ns[1L]
xmax <- ns[3L]
ymin <- 0
ymax <- 0.5
x.range <- xmax - xmin
xmin <- xmin - 0.1*x.range
xmax <- xmax + 0.1*x.range
x.range <- xmax - xmin
label.x <- xmin - 0.31*x.range
label.y <- ymax

p7 <- ggplot(df.mean, aes(x = n, y = error.mean)) +
  geom_line(aes(color = classifier, linetype = classifier), size = line.width) +
  geom_errorbar(aes(ymin=error.mean-error.sem, ymax=error.mean+error.sem, color = classifier, linetype = classifier), size = line.width, width = line.width*2) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = ns, labels = ns) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = c(0, 0.25, 0.5), labels = c("0.00", "0.25", "0.50")) +
  xlab("n") +
  ylab("Error Rate") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  # guides(color=guide_legend(title="Algorithm")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 16), legend.title = element_blank(), legend.position = "right") +
  scale_color_manual(values = c("Control" = "#6a51a3", "RF" = "#4a484c", "RerF" = "#F41711", "S-RerF" = "#F41711"), labels = c("Control", "RF", "RerF", "S-RerF")) +
  scale_linetype_manual(values = c("Control" = "solid", "RF" = "solid", "RerF" = "dotted", "S-RerF" = "solid"), labels = c("Control", "RF", "RerF", "S-RerF"))

gt7 <- ggplot_gtable(ggplot_build(p7))
gt7$layout$clip[] <- "off"

pall4 <- ggarrange(pall3, gt7, nrow = 2L, heights = c(1, 0.6), labels = c("(A)", ""), common.legend = F)

ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/image_stripes.pdf", plot = pall4, width = 5, height = 5.5, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/image_stripes.png", plot = pall4, width = 5, height = 5.5, units = "in")
