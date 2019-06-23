rm(list = ls())
library(ggplot2)
library(ggpubr)
library(scales)
library(RColorBrewer)

# color.map <- brewer.pal(5L, "Set2")
color.map <- c("CCF" = "#6a51a3",
               "RF" = "#4a484c",
               "F-RC" = "#318dde",
               "RerF" = "#F41711",
               "XGBoost" = "#E78AC3")
line.width <- 0.75
marker.size <- 0.75

load(file="~/RandomerForest/R/Results/2018.03.11/exp0027.Rdata")
ds <- cbind(dataset = "MNIST (60000x784)", ress1)

load(file="~/RandomerForest/R/Results/2018.03.11/exp0036.Rdata")
ds <- rbind(ds, cbind(dataset = "Higgs(250000x31)", ress1))

load(file="~/RandomerForest/R/Results/2018.03.11/exp0037.Rdata")
ds <- rbind(ds, cbind(dataset = "p53(31159x5409)", ress1))

lab <- labs(title="Forest Algorithm Strong Scaling", x="Number of Threads", y="Relative Performance")

cols <- c("Ideal"="#000000", "RerF"="#F41711", "XGBoost"="#E78AC3", "Ranger"="#6a51a3")

data.text <- data.frame(x = -Inf, y = Inf, lab = paste0("(", LETTERS[4:6], ")"), dataset = unique(ds$dataset))

p1 <- ggplot(ds,aes(x=Cores_Used, y=Speed_Up)) +
  geom_line(aes(color = Line_Type), size = line.width) +
  facet_grid(. ~ dataset) +
  geom_text(aes(x, y, label = lab), data = data.text, hjust = 0, vjust = 1) +
  scale_x_continuous(breaks = c(0,10,20,30,40)) +
  xlab("Number of Threads") +
  ylab("Relative Performance") +
  ggtitle("Strong Scaling") +
  scale_color_manual(values=cols, labels=c("Ideal", "Ranger", "RerF", "XGBoost")) +
  guides(color = F) +
  theme(aspect.ratio = 1.2,
      legend.position = "bottom") +
  guides(color = guide_legend(title = NULL))

gt1 <- ggplot_gtable(ggplot_build(p1))
gt.legend <- gt1$grobs[[which(sapply(gt1$grobs, function(x) x$name == "guide-box"))]]
p1 <- p1 + guides(color = F)


load(file="~/RandomerForest/R/Results/2018.03.11/exp0033.Rdata")
ds <- cbind(dataset = "MNIST (60000x784)", ress1)

load(file="~/RandomerForest/R/Results/2018.03.11/exp0034.Rdata")
ds <- rbind(ds, cbind(dataset = "Higgs(250000x31)", ress1))

load(file="~/RandomerForest/R/Results/2018.03.11/exp0035.Rdata")
ds <- rbind(ds, cbind(dataset = "p53(31159x5409)", ress1))

lab <- labs(title=paste("Median Per Tree Training Time vs\nNumber of Threads"), x="Number of Threads", y="Time (s)")

# cols <- c("Ideal"="#000000", "RerF"="#009E73", "XGBoost"="#E69F00", "Ranger"="#0072B2", "RF"="#CC79A7")

data.text <- data.frame(x = -Inf, y = Inf, lab = paste0("(", LETTERS[1:3], ")"), dataset = unique(ds$dataset))

p2 <- ggplot(ds,aes(x=Cores_Used, y=Time_Sec)) +
  geom_line(aes(color = Line_Type), size = line.width) +
  facet_grid(. ~ dataset) +
  geom_text(aes(x, y, label = lab), data = data.text, hjust = 0, vjust = 1) +
  scale_x_continuous(breaks = c(0,10,20,30,40)) +
  scale_y_continuous(limits = c(0, 15)) +
  xlab("Number of Threads") +
  ylab("Time (s)") +
  ggtitle("Parallel Execution Time") +
  scale_color_manual(values=cols, labels=c("Ideal", "Ranger", "RerF", "XGBoost")) +
  theme(aspect.ratio = 1.2) +
  guides(color = F)

pall <- ggarrange(p2, p1, gt.legend, nrow = 3L, heights = c(1, 1, 1/8))

ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/speed_scaling.pdf", plot = pall, width = 4, height = 5, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/speed_scaling.png", plot = pall, width = 4, height = 5, units = "in")
