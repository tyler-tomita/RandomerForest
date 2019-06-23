library(ggplot2)
library(RColorBrewer)

fileName <- "~/RandomerForest/R/Results/2018.08.14/predTimes.csv"

mydata <- read.csv(file=fileName, header=FALSE)

leg <- theme(legend.text = element_text(size = 12), legend.title=element_blank(), plot.title = element_text(size = 16,  face="bold"), plot.subtitle = element_text(size = 12),axis.title.x = element_text(size=12), axis.text.x = element_text(size=12), axis.title.y = element_text(size=12), axis.text.y = element_text(size=12))

cols <- c("FP"="#F41711", "RerF"="#F41711", "XGBoost"="#E78AC3", "Ranger"="#6a51a3")
shapes <- c("FP"=19, "RerF"=1, "XGBoost"=19, "Ranger"=19)

p <- ggplot(mydata, aes(x=V2, y=V3, group=V1, color=V1, shape=V1)) + geom_jitter(width=0.3, size=1)

p <- p + guides(fill=guide_legend(title="System"))
p <- p + labs(x = "Dataset", y = "Test Prediction Time (s)", size=14)
p <- p + scale_y_log10(breaks=c(.01,.1,1,10,100,1000))
p <- p + scale_color_manual(values=cols, breaks=c("Ranger","RerF","XGBoost","FP"), labels=c("Ranger","RerF","XGBoost","RerF-Pack"))
p <- p + scale_shape_manual(values=shapes, breaks=c("Ranger","RerF","XGBoost","FP"), labels=c("Ranger","RerF","XGBoost","RerF-Pack"))
p <- p + leg
p <- p + theme(aspect.ratio=0.5)

ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/forest_packing.png", plot = p, width = 6, height = 3, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/forest_packing.pdf", plot = p, width = 6, height = 3, units = "in")
