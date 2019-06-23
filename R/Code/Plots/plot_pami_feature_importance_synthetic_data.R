rm(list = ls())
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

load("~/RandomerForest/R/Results/2017.10.15/Trunk_feature_importance_2017_10_15.RData")
color.map <- c("#af8dc3", "#f7f7f7", "#7fbf7b")
# color.map2 <- brewer.pal(5L, "Set2")
color.map2 <- c("CCF" = "#6a51a3",
               "RF" = "#4a484c",
               "F-RC" = "#318dde",
               "RerF" = "#F41711",
               "XGBoost" = "#E78AC3")
line.width <- 1
marker.size <- 2

p <- 10L
top <- 10L

mu0 <- -1/sqrt(1:p)
mu1 <- 1/sqrt(1:p)
sigma <- rep(1, p)
L.bayes <- rep(0, length(top))

# randomer forest
for (i in 1:top) {
  proj <- matrix(c(rep(i, p), 1:p, rep(0, p)), nrow = p)
  nnz <- length(feature.imp$rerf$proj[[i]])/2
  nz.idx <- feature.imp$rerf$proj[[i]][(1:nnz)*2L - 1L]
  nz.coef <- feature.imp$rerf$proj[[i]][(1:nnz)*2L]
  proj[nz.idx, 3L] <- nz.coef
  L.bayes[i] <- 1 - pnorm(abs((proj[, 3L]/sqrt(sum(proj[, 3L]^2)))%*%(mu1 - mu0))*sqrt(sigma%*%(proj[, 3L]/sqrt(sum(proj[, 3L]^2)))^2)/2)
  if (i == 1L) {
    df <- as.data.frame(proj)
    names(df) <- c("projection", "feature.idx", "coefficient")
  } else {
    df <- rbind(df, data.frame(projection = proj[, 1L], feature.idx = proj[, 2L], coefficient = proj[, 3L]))
  }
}
df$projection <- as.factor(df$projection)
df$feature.idx <- as.factor(df$feature.idx)
df$coefficient <- as.factor(df$coefficient)
df$feature.idx <- with(df, factor(feature.idx, levels = rev(levels(feature.idx))))

panel.label <- textGrob("(A)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- 0.5
xmax <- 10.5
x.range <- xmax - xmin
label.x <- xmin - 0.15*x.range
ymin <- 0.5
ymax <- 10.5
y.range <- ymax - ymin
label.y <- ymax + 0.25*y.range

p1 <- ggplot(df, aes(projection, feature.idx)) +
  geom_tile(aes(fill = coefficient)) +
  geom_vline(aes(xintercept = as.numeric(projection) + 0.5)) +
  xlab("Projection") +
  ylab("Dimension") +
  ggtitle("RerF", subtitle = "Top 10 Projections") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, legend.position = "bottom", legend.background = element_rect(fill = "grey90")) +
  scale_fill_manual(values = color.map, labels = c("-1", "0", "+1"))

gt1 <- ggplot_gtable(ggplot_build(p1))
gt.legend1 <- gt1$grobs[[which(sapply(gt1$grobs, function(x) x$name == "guide-box"))]]
p1 <- p1 + guides(fill = F)
gt1 <- ggplot_gtable(ggplot_build(p1))
gt1$layout$clip[] <- "off"

df2 <- data.frame(projection = 1:top, importance = (feature.imp$rerf$imp/feature.imp$rerf$imp[1L])[1:top], bayes.error = L.bayes, classifier = rep("RerF", top))

# pall <- ggarrange(p2, p1, ncol = 1L, nrow = 2L, heights = c(0.52, 1))
# 
# pall <- ggarrange(pall, ncol = 1L, nrow = 1L, labels = "B")
# 
# ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/feature_importance_synthetic_data_rerf.png", plot = pall, width = 6, height = 4.6, units = "in")
# ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/feature_importance_synthetic_data_rerf.pdf", plot = pall, width = 6, height = 4.6, units = "in")

# random forest
L.bayes <- rep(0, length(top))
for (i in 1:top) {
  proj <- matrix(c(rep(i, p), 1:p, rep(0, p)), nrow = p)
  nnz <- length(feature.imp$rf$proj[[i]])/2
  nz.idx <- feature.imp$rf$proj[[i]][(1:nnz)*2L - 1L]
  nz.coef <- feature.imp$rf$proj[[i]][(1:nnz)*2L]
  proj[nz.idx, 3L] <- nz.coef
  L.bayes[i] <- 1 - pnorm(abs(proj[, 3L]%*%(mu1 - mu0))*sqrt(sigma%*%proj[, 3L]^2)/2)
  if (i == 1L) {
    df <- as.data.frame(proj)
    names(df) <- c("projection", "feature.idx", "coefficient")
  } else {
    df <- rbind(df, data.frame(projection = proj[, 1L], feature.idx = proj[, 2L], coefficient = proj[, 3L]))
  }
}
df$projection <- as.factor(df$projection)
df$feature.idx <- as.factor(df$feature.idx)
df$coefficient <- as.factor(df$coefficient)

df$feature.idx = with(df, factor(feature.idx, levels = rev(levels(feature.idx))))

panel.label <- textGrob("(B)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- 0.5
xmax <- 10.5
x.range <- xmax - xmin
label.x <- xmin - 0.15*x.range
ymin <- 0.5
ymax <- 10.5
y.range <- ymax - ymin
label.y <- ymax + 0.25*y.range

p2 <- ggplot(df, aes(projection, feature.idx)) +
  geom_tile(aes(fill = coefficient)) +
  geom_vline(aes(xintercept = as.numeric(projection) + 0.5)) +
  xlab("Projection") +
  ylab("Dimension") +
  ggtitle("RF", subtitle = "Top 10 Projections") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, legend.background = element_rect(fill = "grey90")) +
  scale_fill_manual(values = color.map[2:3], labels = c("0", "+1")) +
  guides(fill = F)

gt2 <- ggplot_gtable(ggplot_build(p2))
gt2$layout$clip[] <- "off"

df2 <- rbind(df2, data.frame(projection = 1:top, importance = (feature.imp$rf$imp/feature.imp$rf$imp[1L])[1:top], bayes.error = L.bayes, classifier = rep("RF", top)))

panel.label <- textGrob("(C)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- 0.5
xmax <- 10.5
x.range <- xmax - xmin
label.y <- xmax + 0.15*x.range
ymin <- min(df2$importance)
ymax <- max(df2$importance)
y.range <- ymax - ymin
label.y <- ymax + 0.32*y.range

p3 <- ggplot(df2, aes(x = projection + 0.02, y = importance, color = classifier)) +
  geom_line(size = line.width) +
  geom_point(size = marker.size) +
  scale_x_continuous(limits = c(0.85, 10.15), breaks = 1:10) +
  scale_y_continuous(limits = c(ymin, ymax)) +
  xlab("Projection") +
  ylab("Gini Importance") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  theme_classic() +
  # theme(axis.line = element_line(size = line.width), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = 14), aspect.ratio = 1) +
  theme(plot.margin = unit(c(0.575, 0, 0.075, 0), "in"), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, legend.position = "bottom") +
  scale_color_manual(values = color.map2[c(4, 2)])

gt3 <- ggplot_gtable(ggplot_build(p3))
gt.legend2 <- gt3$grobs[[which(sapply(gt3$grobs, function(x) x$name == "guide-box"))]]
p3 <- p3 + guides(color = F)
gt3 <- ggplot_gtable(ggplot_build(p3))
gt3$layout$clip[] <- "off"

panel.label <- textGrob("(D)", gp=gpar(fontsize=14, fontface="bold"))
xmin <- 0.5
xmax <- 10.5
x.range <- xmax - xmin
label.y <- xmax + 0.15*x.range
ymin <- min(df2$bayes.error)
ymax <- max(df2$bayes.error)
y.range <- ymax - ymin
label.y <- ymax + 0.32*y.range

p4 <- ggplot(df2, aes(x = projection + 0.02, y = bayes.error, color = classifier)) +
  geom_line(size = line.width) +
  geom_point(size = marker.size) +
  scale_x_continuous(limits = c(0.85, 10.15), breaks = 1:10) +
  scale_y_continuous(limits = c(ymin, ymax)) +
  xlab("Projection") +
  ylab("Bayes Error") +
  annotation_custom(panel.label, xmin = label.x, xmax = label.x, ymin = label.y, ymax = label.y) +
  theme_classic() +
  # theme(axis.line = element_line(size = line.width), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = 14), aspect.ratio = 1) +
  theme(plot.margin = unit(c(0.575, 0, 0.075, 0), "in"), axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1) +
  scale_color_manual(values = color.map2[c(4, 2)]) +
  guides(color = F)

gt4 <- ggplot_gtable(ggplot_build(p4))
gt4$layout$clip[] <- "off"

psub1 <- ggarrange(gt1, gt2, gt.legend1, nrow = 3L, heights = c(1, 1, 1/8))
psub2 <- ggarrange(gt3, gt4, gt.legend2, nrow = 3L, heights = c(1, 1, 1/8))
# psub3 <- ggarrange(gt.legend1, gt.legend2, nrow = 2L)
pall <- ggarrange(psub1, psub2, ncol = 2L, widths = c(1, 1))

ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/feature_importance_synthetic_data.png", plot = pall, width = 5.5, height = 6, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/feature_importance_synthetic_data.pdf", plot = pall, width = 5.5, height = 6, units = "in")
