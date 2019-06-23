rm(list = ls())
library(rerf)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

line.width <- 0.75
marker.size <- 0.75
color.map <- brewer.pal(5L, "Set2")

load("~/RandomerForest/R/Results/2017.10.28/vector_angle_probability.RData")
df$p <- log2(df$p)
xmin <- min(df$p)
xmax <- max(df$p)

theta <- 1
tv <- "sparse"
df.plot <- df[(df$target.vector == tv) & (df$theta == theta), ]
ps <- unique(df.plot$p)
p1 <- ggplot(df.plot, aes(x = p, y = prob.success)) +
  geom_line(aes(colour = random.vector), size = line.width) +
  geom_errorbar(aes(ymin=prob.success-sem.success, ymax=prob.success+sem.success, colour = random.vector), size = line.width, width = 0.1) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = ps, labels =2^ps) +
  xlab("p") +
  ylab(bquote("P("*theta <= .(theta)*degree*")")) +
  ggtitle(expression(lambda^"*" == 1/p)) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.title = element_text(size = 14)) +
  scale_color_manual(values = color.map[c(2, 1)], labels = c("RR-RF", "RerF"))

theta <- 10
tv <- "sparse"
df.plot <- df[(df$target.vector == tv) & (df$theta == theta), ]
ps <- unique(df.plot$p)
p2 <- ggplot(df.plot, aes(x = p, y = prob.success)) +
  geom_line(aes(colour = random.vector), size = line.width) +
  geom_errorbar(aes(ymin=prob.success-sem.success, ymax=prob.success+sem.success, colour = random.vector), size = line.width, width = 0.1) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = ps, labels =2^ps) +
  xlab("p") +
  ylab(bquote("P("*theta <= .(theta)*degree*")")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.title = element_text(size = 14)) +
  scale_color_manual(values = color.map[c(2, 1)], labels = c("RR-RF", "RerF"))

theta <- 22.5
tv <- "sparse"
df.plot <- df[(df$target.vector == tv) & (df$theta == theta), ]
ps <- unique(df.plot$p)
p3 <- ggplot(df.plot, aes(x = p, y = prob.success)) +
  geom_line(aes(colour = random.vector), size = line.width) +
  geom_errorbar(aes(ymin=prob.success-sem.success, ymax=prob.success+sem.success, colour = random.vector), size = line.width, width = 0.1) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = ps, labels =2^ps) +
  xlab("p") +
  ylab(bquote("P("*theta <= .(theta)*degree*")")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.title = element_text(size = 14)) +
  scale_color_manual(values = color.map[c(2, 1)], labels = c("RR-RF", "RerF"))

theta <- 45
tv <- "sparse"
df.plot <- df[(df$target.vector == tv) & (df$theta == theta), ]
ps <- unique(df.plot$p)
p4 <- ggplot(df.plot, aes(x = p, y = prob.success)) +
  geom_line(aes(colour = random.vector), size = line.width) +
  geom_errorbar(aes(ymin=prob.success-sem.success, ymax=prob.success+sem.success, colour = random.vector), size = line.width, width = 0.1) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = ps, labels =2^ps) +
  xlab("p") +
  ylab(bquote("P("*theta <= .(theta)*degree*")")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.title = element_text(size = 14)) +
  scale_color_manual(values = color.map[c(2, 1)], labels = c("RR-RF", "RerF"))

theta <- 1
tv <- "dense"
df.plot <- df[(df$target.vector == tv) & (df$theta == theta), ]
ps <- unique(df.plot$p)
p5 <- ggplot(df.plot, aes(x = p, y = prob.success)) +
  geom_line(aes(colour = random.vector), size = line.width) +
  geom_errorbar(aes(ymin=prob.success-sem.success, ymax=prob.success+sem.success, colour = random.vector), size = line.width, width = 0.1) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = ps, labels =2^ps) +
  xlab("p") +
  ylab(bquote("P("*theta <= .(theta)*degree*")")) +
  ggtitle(expression(lambda^"*" == 1)) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.title = element_text(size = 14)) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = ps, labels =2^ps) +
  scale_color_manual(values = color.map[c(2, 1)], labels = c("RR-RF", "RerF"))

theta <- 10
tv <- "dense"
df.plot <- df[(df$target.vector == tv) & (df$theta == theta), ]
ps <- unique(df.plot$p)
p6 <- ggplot(df.plot, aes(x = p, y = prob.success)) +
  geom_line(aes(colour = random.vector), size = line.width) +
  geom_errorbar(aes(ymin=prob.success-sem.success, ymax=prob.success+sem.success, colour = random.vector), size = line.width, width = 0.1) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = ps, labels =2^ps) +
  xlab("p") +
  ylab(bquote("P("*theta <= .(theta)*degree*")")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.title = element_text(size = 14)) +
  scale_color_manual(values = color.map[c(2, 1)], labels = c("RR-RF", "RerF"))

theta <- 22.5
tv <- "dense"
df.plot <- df[(df$target.vector == tv) & (df$theta == theta), ]
ps <- unique(df.plot$p)
p7 <- ggplot(df.plot, aes(x = p, y = prob.success)) +
  geom_line(aes(colour = random.vector), size = line.width) +
  geom_errorbar(aes(ymin=prob.success-sem.success, ymax=prob.success+sem.success, colour = random.vector), size = line.width, width = 0.1) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = ps, labels =2^ps) +
  xlab("p") +
  ylab(bquote("P("*theta <= .(theta)*degree*")")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.title = element_text(size = 14)) +
  scale_color_manual(values = color.map[c(2, 1)], labels = c("RR-RF", "RerF"))

theta <- 45
tv <- "dense"
df.plot <- df[(df$target.vector == tv) & (df$theta == theta), ]
ps <- unique(df.plot$p)
p8 <- ggplot(df.plot, aes(x = p, y = prob.success)) +
  geom_line(aes(colour = random.vector), size = line.width) +
  geom_errorbar(aes(ymin=prob.success-sem.success, ymax=prob.success+sem.success, colour = random.vector), size = line.width, width = 0.1) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = ps, labels =2^ps) +
  xlab("p") +
  ylab(bquote("P("*theta <= .(theta)*degree*")")) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 1, plot.title = element_text(size = 14)) +
  scale_color_manual(values = color.map[c(2, 1)], labels = c("RR-RF", "RerF"))

pall <- ggarrange(p1, p5, p2, p6, p3, p7, p4, p8, nrow = 4L, ncol = 2L, common.legend = T, heights = c(1, 0.9, 0.9, 0.9), legend = "right", labels = paste0("(", LETTERS[c(1,5,2,6,3,7,4,8)], ")"))
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/vector_angle_probability.pdf", plot = pall, width = 6.5, height = 9, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/vector_angle_probability.png", plot = pall, width = 6.5, height = 9, units = "in")
