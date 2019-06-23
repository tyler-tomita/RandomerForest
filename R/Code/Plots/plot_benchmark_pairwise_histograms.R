rm(list = ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
# library(xtable)
library(gridExtra)
source("~/RandomerForest/R/Code/Utils/GetCatMap.R")
source("~/RandomerForest/R/Code/Utils/GetFolds.R")


filePath <- "~/RandomerForest/R/Results/2018.02.07/"
filePathCCF <- "~/RandomerForest/Results/2018.06.12/"
contents <- list.files(filePath)
contents <- contents[!grepl("frc", contents) & !grepl("rerfc", contents) & !grepl("rr-rf", contents) & !grepl("rr_rf", contents) & !grepl("xgb", contents)]
catfiles <- list.files("~/tmp/uci/processed/categorical_map/")
load(paste0(filePath, contents[[1L]]))
fieldNames <- ls()
fieldNames <- fieldNames[(fieldNames != "filePath") & (fieldNames != "contents") & (fieldNames != "multiplot") &
                         (fieldNames != "catfiles") & (fieldNames != "GetCatMap") & (fieldNames != "filePathCCF") &
                         (fieldNames != "GetFolds")]

res <- vector("list", length(fieldNames))
names(res) <- fieldNames

line.width <- 0.75
marker.size <- 1
# color.map <- brewer.pal(5L, "Set2")
color.map <- c("RerF" = "#F41711",
               "RF" = "#4a484c",
               "XGBoost" = "#E78AC3",
               "RR-RF" = "#318dde",
               "CCF" = "#6a51a3")
color.map <- color.map[names(color.map) != "RerF"]

for (f in contents) {
  ds <- strsplit(f, "_2018_02_07.RData")[[1L]][1L]
  load(paste0(filePath, f))
  for (fn in fieldNames) {
    res[[fn]][[strsplit(f, "_2018_02_07.RData")[[1L]]]] <- get(fn)[[1L]]
  }
  load(paste0(filePath, paste0(ds, "_xgb_2018_02_07.RData")))
  res[["testError"]][[ds]]$xgb <- get("testError")[[1L]]
  res[["bestIdx"]][[ds]]$xgb <- get("colSample")[[1L]]
  
  if (file.exists(paste0(filePathCCF, paste0(ds, "_ccf_2018_06_12.csv")))) {
    res[["testError"]][[ds]]$ccf <- c(as.matrix(read.table(paste0(filePathCCF, ds, "_ccf_2018_06_12.csv"), header = F, sep = ',', quote = '', row.names = NULL)))
    res[["bestIdx"]][[ds]]$ccf <- NULL
  } else {
    res[["testError"]][[ds]]$ccf <- NULL
    res[["bestIdx"]][[ds]]$ccf <- NULL
    print(ds)
  }
  if (file.exists(paste0(filePathCCF, paste0(ds, "_ccf_naive_cat_2018_06_12.csv")))) {
    res[["testError"]][[ds]]$ccfn <- c(as.matrix(read.table(paste0(filePathCCF, ds, "_ccf_naive_cat_2018_06_12.csv"), header = F, sep = ',', quote = '', row.names = NULL)))
    res[["bestIdx"]][[ds]]$ccfn <- NULL
  } else {
    res[["testError"]][[ds]]$ccfn <- NULL
    res[["bestIdx"]][[ds]]$ccfn <- NULL
    print(ds)
  }
  
  if (file.exists(paste0(filePath, paste0(ds, "_rr_rf_2018_02_07.RData")))) {
    load(paste0(filePath, paste0(ds, "_rr_rf_2018_02_07.RData")))
    res[["testError"]][[ds]]$`rr-rf` <- get("testError")[[1L]][[1L]]
    res[["bestIdx"]][[ds]]$`rr-rf` <- get("bestIdx")[[1L]][[1L]]
    if (length(testError[[ds]]$`rr-rf`) == 0) {
      print(ds)
    }
  } else {
    print(ds)
    res[["testError"]][[ds]]$`rr-rf` <- NULL
    res[["bestIdx"]][[ds]]$`rr-rf` <- NULL
  }
}

classifiers <- c("rf", "rerf", "xgb", "rr-rf", "ccfn")
dataSets <- names(res$testError)
error.matrix <- matrix(NA, nrow=length(res$testError)*5L, ncol=length(classifiers))
mean.diff <- matrix(NA, nrow = length(res$testError), ncol=length(classifiers)-1L)
mean.error <- matrix(NA, nrow = length(res$testError), ncol = length(classifiers))
sem.error <-  matrix(NA, nrow = length(res$testError), ncol = length(classifiers))
chance.error <- matrix(NA, nrow=length(res$testError)*5L, ncol=1L)
rownames(error.matrix) <- rep(dataSets, each=5L)
colnames(error.matrix) <- classifiers
rownames(mean.diff) <- dataSets
colnames(mean.diff) <- classifiers[classifiers!="rerf"]
rownames(mean.error) <- dataSets
colnames(mean.error) <- classifiers
rownames(sem.error) <- dataSets
colnames(sem.error) <- classifiers
rownames(chance.error) <- rep(dataSets, each=5L)

n <- integer(length(dataSets))
p <- integer(length(dataSets))
pcat <- integer(length(dataSets))
for (i in seq.int(length(dataSets))) {
ds <- dataSets[i]
D <- as.matrix(read.table(paste0("~/tmp/uci/processed/data/", ds, ".csv"), header = F, sep = ",", quote = "", row.names = NULL))
Y <- as.integer(D[, ncol(D)])
n[i] <- length(Y)
if (paste0(ds, "_catmap.txt") %in% catfiles) {
  cat.map <- GetCatMap(paste0("~/tmp/uci/processed/categorical_map/", ds, "_catmap.txt"))
  pcat[i] <- length(cat.map)
  p[i] <- cat.map[[1L]][1L] - 1L + pcat[i]
} else {
  p[i] <- ncol(D) - 1L
}
fold <- GetFolds(paste0("~/tmp/uci/processed/cv_partitions/", ds, "_partitions.txt"))
for (k in 1:5) {
  chance.error[(i-1L)*5L + k, ] <- 1 - max(tabulate(Y[fold[[k]]]))/length(Y[fold[[k]]])
}

for (j in seq.int(length(classifiers))) {
  cl <- classifiers[j]
  if (is.null(res$testError[[ds]][[cl]]) | (length(res$testError[[ds]][[cl]]) == 0) | (any(res$bestIdx[[ds]][[cl]] == 0)) | any(is.na(res$bestIdx[[ds]][[cl]]))) {
    mean.error[i, j] <- NA
    sem.error[i, j] <- NA
    next
  }
  if (cl == "xgb" || cl == "ccf" || cl == "ccfn") {
    error.matrix[((i-1L)*5L + 1L):(i*5), j] <- res$testError[[ds]][[cl]]
    mean.error[i, j] <- mean(res$testError[[ds]][[cl]])
    sem.error[i, j] <- sd(res$testError[[ds]][[cl]])/sqrt(length(res$testError[[ds]][[cl]]))
  } else {
    error.matrix[((i-1L)*5L + 1L):(i*5), j] <- sapply(1:length(res$bestIdx[[ds]][[cl]]),
                                                    function (x) res$testError[[ds]][[cl]][x, res$bestIdx[[ds]][[cl]][x]])
    mean.error[i, j] <- mean(sapply(1:length(res$bestIdx[[ds]][[cl]]),
                                    function (x) res$testError[[ds]][[cl]][x, res$bestIdx[[ds]][[cl]][x]]))
    sem.error[i, j] <- sd(sapply(1:length(res$bestIdx[[ds]][[cl]]),
                                 function (x) res$testError[[ds]][[cl]][x, res$bestIdx[[ds]][[cl]][x]]))/sqrt(length(res$bestIdx[[ds]][[cl]]))
  }
}
for (j in colnames(mean.error)) {
  mean.error[i, j] <- mean(error.matrix[((i-1L)*5L + 1L):(i*5), j]/chance.error[((i-1L)*5L + 1L):(i*5), ])
  sem.error[i, j] <- sd(error.matrix[((i-1L)*5L + 1L):(i*5), j]/chance.error[((i-1L)*5L + 1L):(i*5), ])
  if (j != "rerf") {
    mean.diff[i, j] <- mean((error.matrix[((i-1L)*5L + 1L):(i*5), j] - error.matrix[((i-1L)*5L + 1L):(i*5), "rerf"])/chance.error[((i-1L)*5L + 1L):(i*5), ])
  }
}
}

no.na <- apply(mean.error, 1, function(x) !any(is.na(x)))
chance.error <- chance.error[no.na, , drop=FALSE]
mean.error <- mean.error[no.na, ]
sem.error <- sem.error[no.na, ]
mean.diff <- mean.diff[no.na, ]
colnames(sem.error) <- colnames(mean.error)
n <- n[no.na]
p <- p[no.na]
pcat <- pcat[no.na]
dataSets <- dataSets[no.na]
pcat2 <- rep(pcat, each = 5L)
error.matrix <- error.matrix[apply(error.matrix, 1, function(x) !any(is.na(x))), ]

### make table of error rates for each benchmark dataset ###

# df.error <- as.data.frame(cbind(as.character(n), as.character(p - pcat), as.character(pcat), matrix(paste0("$", as.character(t(sapply(seq_along(n), function(x) round(mean.error[x, c("rerf", "rf", "xgb", "rr-rf", "ccfn")], max(1, floor(log10(n[x]))))))), " \\pm ", as.character(sapply(seq_along(n), function(x) round(sem.error[x, c("rerf", "rf", "xgb", "rr-rf", "ccf")], max(1, floor(log10(n[x])))))), "$"), nrow(mean.error), ncol(mean.error))))
# rownames(df.error) <- dataSets
# colnames(df.error) <- c("n", "p_{num}", "p_{cat}", "RerF", "RF", "XGBoost", "RR-RF", "CCF")
# print(xtable(df.error, type = "latex"), file = "~/benchmark_error_table.tex")

### ###

colnames(mean.error) <- c("RF", "RerF", "XGBoost", "RR-RF", "CCF")
colnames(error.matrix) <- c("RF", "RerF", "XGBoost", "RR-RF", "CCF")
colnames(mean.diff) <- c("RF", "XGBoost", "RR-RF", "CCF")
df <- data.frame(error.diff = mean.diff[, "RF"],
                 error.mean = apply(mean.error[, c("RerF", "RF")], 1, mean),
                 pair = rep("RF", nrow(mean.error)),
                 type = ifelse(pcat > 0, "Categorical", "Numeric"),
                 type2 = rep("All", nrow(mean.error)))
for (cl in colnames(mean.error)[!(colnames(mean.error) %in% c("RF", "RerF"))]) {
  df <- rbind(df, data.frame(error.diff = mean.diff[, cl],
                             error.mean = apply(mean.error[, c("RerF", cl)], 1, mean),
                             pair = rep(cl, nrow(mean.error)),
                             type = ifelse(pcat > 0, "Categorical", "Numeric"),
                             type2 = rep("All", nrow(mean.error))))
}

df$error.diff[df$error.diff > 0] <- sqrt(df$error.diff[df$error.diff > 0])
df$error.diff[df$error.diff < 0] <- -sqrt(-df$error.diff[df$error.diff < 0])

data.text <- data.frame(x = c(Inf, Inf), y = c(-0.55, 0.55), lab = c("RerF\nworse", "RerF\nbetter"))
data.text2 <- data.frame(x = 0, y = -1/2, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 == 0, c("RerF", "RF")], 1, diff)/chance.error[pcat2 == 0, ], alternative = "greater")$p.value, 2), digits = 1, format = "e")))
data.text3 <- data.frame(x = 0, y = -2/3, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 == 0, c("RerF", "XGBoost")], 1, diff)/chance.error[pcat2 == 0, ], alternative = "greater")$p.value, 2), digits = 1, format = "e")))
data.text4 <- data.frame(x = 0, y = -5/6, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 == 0, c("RerF", "RR-RF")], 1, diff)/chance.error[pcat2 == 0, ], alternative = "greater")$p.value, 2), digits = 1, format = "e")))
data.text5 <- data.frame(x = 0, y = -1, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 == 0, c("RerF", "CCF")], 1, diff)/chance.error[pcat2 == 0, ], alternative = "greater")$p.value, 2), digits = 1, format = "e")))

mgn1 <- 0.1
mgn2 <- 0.5
mgn3 <- 0
ar <- 3

p1.scatter <- ggplot(df[df$type == "Numeric", ], aes(sqrt(error.mean), error.diff)) +
  geom_point(aes(color = pair), size = marker.size) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  geom_hline(yintercept = 0, size = line.width) +
  geom_text(aes(x, y, label = lab), data = data.text, hjust = 1, vjust = 0.5, alpha = 1) +
  geom_text(aes(x, y, label = lab), data = data.text2, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[1L]) +
  geom_text(aes(x, y, label = lab), data = data.text3, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[2L]) +
  geom_text(aes(x, y, label = lab), data = data.text4, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[3L]) +
  geom_text(aes(x, y, label = lab), data = data.text5, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[4L]) +
  xlab(expression(sqrt("Mean Error"))) +
  ylab(expression(sqrt("Difference in Error"))) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width),
        text = element_text(size = 14),
        aspect.ratio = 1,
        legend.title = element_blank(),
        plot.margin = unit(c(0.1, 0, 0.1, 0.4), "in")) +
  scale_color_manual(values = color.map) +
  guides(color = F)

p1.hist <- ggplot(df[df$type == "Numeric", ]) +
  geom_density(aes(x = error.diff, color = pair), size = line.width) +
  # geom_density(aes(x = error.diff, color = type2), size = line.width) +
  geom_vline(xintercept = 0, size = line.width) +
  scale_x_continuous(limits = c(-1, 1)) +
  coord_flip() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = ar,
        plot.margin = unit(c(mgn1, 0, mgn2, mgn3), "in"),
        panel.background = element_blank()) +
  scale_color_manual(values = color.map) +
  guides(color = F)


gt1 <- grid.arrange(p1.scatter, p1.hist, ncol = 2L, widths = c(1, 1/3), top = "Numeric")
gt1$layout$clip[] <- "off"

data.text2 <- data.frame(x = 0, y = -1/2, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 > 0, c("RerF", "RF")], 1, diff)/chance.error[pcat2 > 0, ], alternative = "greater")$p.value, 2), digits = 1, format = "e")))
data.text3 <- data.frame(x = 0, y = -2/3, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 > 0, c("RerF", "XGBoost")], 1, diff)/chance.error[pcat2 > 0, ], alternative = "greater")$p.value, 2), digits = 1, format = "e")))
data.text4 <- data.frame(x = 0, y = -5/6, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 > 0, c("RerF", "RR-RF")], 1, diff)/chance.error[pcat2 > 0, ], alternative = "greater")$p.value, 2), digits = 1, format = "e")))
data.text5 <- data.frame(x = 0, y = -1, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 > 0, c("RerF", "CCF")], 1, diff)/chance.error[pcat2 > 0, ], alternative = "greater")$p.value, 2), digits = 1, format = "e")))

p2.scatter <- ggplot(df[df$type == "Categorical", ], aes(sqrt(error.mean), error.diff)) +
  geom_point(aes(color = pair), size = marker.size) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  geom_hline(yintercept = 0, size = line.width) +
  geom_text(aes(x, y, label = lab), data = data.text2, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[1L]) +
  geom_text(aes(x, y, label = lab), data = data.text3, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[2L]) +
  geom_text(aes(x, y, label = lab), data = data.text4, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[3L]) +
  geom_text(aes(x, y, label = lab), data = data.text5, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[4L]) +
  xlab(expression(sqrt("Mean Error"))) +
  ylab(expression(sqrt("Difference in Error"))) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width),
        text = element_text(size = 14),
        aspect.ratio = 1,
        legend.title = element_blank(),
        plot.margin = unit(c(0.1, 0, 0.1, 0.4), "in")) +
  scale_color_manual(values = color.map) +
  guides(color = F)

p2.hist <- ggplot(df[df$type == "Categorical", ]) +
  geom_density(aes(x = error.diff, color = pair), size = line.width) +
  geom_vline(xintercept = 0, size = line.width) +
  scale_x_continuous(limits = c(-1, 1)) +
  coord_flip() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = ar,
        plot.margin = unit(c(mgn1, 0, mgn2, mgn3), "in"),
        panel.background = element_blank()) +
  scale_color_manual(values = color.map) +
  guides(color = F)

gt2 <- grid.arrange(p2.scatter, p2.hist, ncol = 2L, widths = c(1, 1/3), top = "Categorical")
gt2$layout$clip[] <- "off"

data.text2 <- data.frame(x = 0, y = -1/2, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[, c("RerF", "RF")], 1, diff)/chance.error, alternative = "greater")$p.value, 2), digits = 1, format = "e")))
data.text3 <- data.frame(x = 0, y = -2/3, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[, c("RerF", "XGBoost")], 1, diff)/chance.error, alternative = "greater")$p.value, 2), digits = 1, format = "e")))
data.text4 <- data.frame(x = 0, y = -5/6, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[, c("RerF", "RR-RF")], 1, diff)/chance.error, alternative = "greater")$p.value, 2), digits = 1, format = "e")))
data.text5 <- data.frame(x = 0, y = -1, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[, c("RerF", "CCF")], 1, diff)/chance.error, alternative = "greater")$p.value, 2), digits = 1, format = "e")))

p3.scatter <- ggplot(df, aes(sqrt(error.mean), error.diff)) +
  geom_point(aes(color = pair), size = marker.size) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  geom_hline(yintercept = 0, size = line.width) +
  geom_text(aes(x, y, label = lab), data = data.text2, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[1L]) +
  geom_text(aes(x, y, label = lab), data = data.text3, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[2L]) +
  geom_text(aes(x, y, label = lab), data = data.text4, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[3L]) +
  geom_text(aes(x, y, label = lab), data = data.text5, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[4L]) +
  xlab(expression(sqrt("Mean Error"))) +
  ylab(expression(sqrt("Difference in Error"))) +
  theme_classic() +
  theme(axis.line = element_line(size = line.width),
        text = element_text(size = 14),
        aspect.ratio = 1,
        legend.title = element_blank(),
        plot.margin = unit(c(0.1, 0, 0.1, 0.4), "in")) +
  scale_color_manual(values = color.map) +
  guides(color = F)

p3.hist <- ggplot(df) +
  geom_density(aes(x = error.diff, color = pair), size = line.width) +
  geom_vline(xintercept = 0, size = line.width) +
  scale_x_continuous(limits = c(-1, 1)) +
  coord_flip() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = ar,
        plot.margin = unit(c(mgn1, 0, mgn2, mgn3), "in"),
        panel.background = element_blank(),
        legend.position = "bottom") +
  scale_color_manual(values = color.map) +
  guides(color = guide_legend(title.position = "top", title = "Alg Compared to RerF"))

gt3 <- ggplot_gtable(ggplot_build(p3.hist))
# gt.legend <- gt3$grobs[[which(sapply(gt3$grobs, function(x) x$name == "guide-box"))]]
gt.legend <- cowplot::get_legend(p3.hist)
# gt.legend$grobs[[1]]$grobs[[2]]$label <- "Alg Compared to RerF"
# gt.legend$grobs[[1]]$grobs[[3]]$gp$fill <- NA
# gt.legend$grobs[[1]]$grobs[[6]]$gp$fill <- NA
# gt.legend$grobs[[1]]$grobs[[9]]$gp$fill <- NA
p3.hist <- p3.hist + guides(color = F)

gt3 <- grid.arrange(p3.scatter, p3.hist, ncol = 2L, widths = c(1, 1/3), top = "All")
gt3$layout$clip[] <- "off"

pall1 <- ggarrange(gt1, gt2, gt3, ncol = 3L, labels = paste0("(", LETTERS[1:3], ")"))

pall2 <- ggarrange(pall1, gt.legend, nrow = 2L, heights = c(1, 1/4))

ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/error_histogram_benchmarks.png", plot = pall2, width = 12, height = 3.5, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/error_histogram_benchmarks.pdf", plot = pall2, width = 12, height = 3.5, units = "in")
