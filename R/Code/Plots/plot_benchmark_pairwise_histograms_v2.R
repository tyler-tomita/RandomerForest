# Plot histograms of error relative to RF for each algorithm and plot heatmap of each algorithm's frequency of performance rank

rm(list = ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(gridExtra)
library(grid)
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
color.map <- color.map[names(color.map) != "RF"]

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
colnames(mean.diff) <- classifiers[classifiers!="rf"]
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
    if (j != "rf") {
      mean.diff[i, j] <- mean((error.matrix[((i-1L)*5L + 1L):(i*5), j] - error.matrix[((i-1L)*5L + 1L):(i*5), "rf"])/chance.error[((i-1L)*5L + 1L):(i*5), ])
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

colnames(mean.error) <- c("RF", "RerF", "XGBoost", "RR-RF", "CCF")
colnames(error.matrix) <- c("RF", "RerF", "XGBoost", "RR-RF", "CCF")
colnames(mean.diff) <- c("RerF", "XGBoost", "RR-RF", "CCF")
df <- data.frame(error.diff = mean.diff[, "RerF"],
                 pair = rep("RerF", nrow(mean.error)),
                 type = ifelse(pcat > 0, "Categorical", "Numeric"),
                 type2 = rep("All", nrow(mean.error)))
for (cl in colnames(mean.error)[!(colnames(mean.error) %in% c("RF", "RerF"))]) {
  df <- rbind(df, data.frame(error.diff = mean.diff[, cl],
                             pair = rep(cl, nrow(mean.error)),
                             type = ifelse(pcat > 0, "Categorical", "Numeric"),
                             type2 = rep("All", nrow(mean.error))))
}

df$error.diff[df$error.diff > 0] <- sqrt(df$error.diff[df$error.diff > 0])
df$error.diff[df$error.diff < 0] <- -sqrt(-df$error.diff[df$error.diff < 0])

panel.label <- data.frame(x = -1.25, y = 4, lab = "(A)")
data.text <- data.frame(x = c(-0.5, 0.5), y = c(Inf, Inf), lab = c("RF worse", "RF better"))
data.text2 <- data.frame(x = -1, y = 3, lab = paste0("MRES = ", formatC(signif(mean(apply(error.matrix[pcat2 == 0, c("RF", "RerF")], 1, diff)/chance.error[pcat2 == 0, ])), digits = 1, format = "e")))
data.text3 <- data.frame(x = -1, y = 5/2, lab = paste0("MRES = ", formatC(signif(mean(apply(error.matrix[pcat2 == 0, c("RF", "XGBoost")], 1, diff)/chance.error[pcat2 == 0, ])), digits = 1, format = "e")))
data.text4 <- data.frame(x = -1, y = 2, lab = paste0("MRES = ", formatC(signif(mean(apply(error.matrix[pcat2 == 0, c("RF", "RR-RF")], 1, diff)/chance.error[pcat2 == 0, ])), digits = 1, format = "e")))
data.text5 <- data.frame(x = -1, y = 3/2, lab = paste0("MRES = ", formatC(signif(mean(apply(error.matrix[pcat2 == 0, c("RF", "CCF")], 1, diff)/chance.error[pcat2 == 0, ])), digits = 1, format = "e")))
# data.text2 <- data.frame(x = -1, y = 3, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 == 0, c("RF", "RerF")], 1, diff)/chance.error[pcat2 == 0, ], alternative = "less")$p.value, 2), digits = 1, format = "e")))
# data.text3 <- data.frame(x = -1, y = 5/2, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 == 0, c("RF", "XGBoost")], 1, diff)/chance.error[pcat2 == 0, ], alternative = "less")$p.value, 2), digits = 1, format = "e")))
# data.text4 <- data.frame(x = -1, y = 2, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 == 0, c("RF", "RR-RF")], 1, diff)/chance.error[pcat2 == 0, ], alternative = "less")$p.value, 2), digits = 1, format = "e")))
# data.text5 <- data.frame(x = -1, y = 3/2, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 == 0, c("RF", "CCF")], 1, diff)/chance.error[pcat2 == 0, ], alternative = "less")$p.value, 2), digits = 1, format = "e")))

mgn1 <- 0.1
mgn2 <- 0.5
mgn3 <- 0
ar <- 3

p1.hist <- ggplot(df[df$type == "Numeric", ]) +
  geom_vline(xintercept = 0, size = line.width, linetype = "dashed") +
  geom_density(aes(x = error.diff, color = pair), size = line.width) +
  # geom_density(aes(x = error.diff, color = type2), size = line.width) +
  # geom_text(aes(x, y, label = lab), data = panel.label, hjust = 1, vjust = 0.5, fontface = "bold") +
  geom_text(aes(x, y, label = lab), data = data.text, hjust = 0.5, vjust = 1) +
  geom_text(aes(x, y, label = lab), data = data.text2, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[1L]) +
  geom_text(aes(x, y, label = lab), data = data.text3, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[2L]) +
  geom_text(aes(x, y, label = lab), data = data.text4, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[3L]) +
  geom_text(aes(x, y, label = lab), data = data.text5, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[4L]) +
  scale_x_continuous(limits = c(-1, 1)) +
  ylim(0, 3.5) +
  xlab(expression("sign(RES)"~sqrt("|"~RES~"|"))) +
  ylab("Density Estimate") +
  ggtitle("Numeric") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = color.map) +
  guides(color = guide_legend(title.position = "top", title = "Alg Compared to RerF"))

  # gt <- ggplot_gtable(ggplot_build(p1.hist))
  # gt$layout$clip[gt$layout$name == "panel"] <- "off"
  # grid.draw(gt)

  data.text2 <- data.frame(x = -1, y = 3, lab = paste0("MRES = ", formatC(signif(mean(apply(error.matrix[pcat2 > 0, c("RF", "RerF")], 1, diff)/chance.error[pcat2 > 0, ])), digits = 1, format = "e")))
  data.text3 <- data.frame(x = -1, y = 5/2, lab = paste0("MRES = ", formatC(signif(mean(apply(error.matrix[pcat2 > 0, c("RF", "XGBoost")], 1, diff)/chance.error[pcat2 > 0, ])), digits = 1, format = "e")))
  data.text4 <- data.frame(x = -1, y = 2, lab = paste0("MRES = ", formatC(signif(mean(apply(error.matrix[pcat2 > 0, c("RF", "RR-RF")], 1, diff)/chance.error[pcat2 > 0, ])), digits = 1, format = "e")))
  data.text5 <- data.frame(x = -1, y = 3/2, lab = paste0("MRES = ", formatC(signif(mean(apply(error.matrix[pcat2 > 0, c("RF", "CCF")], 1, diff)/chance.error[pcat2 > 0, ])), digits = 1, format = "e")))
# data.text2 <- data.frame(x = -1, y = 3, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 > 0, c("RF", "RerF")], 1, diff)/chance.error[pcat2 > 0, ], alternative = "less")$p.value, 2), digits = 1, format = "e")))
# data.text3 <- data.frame(x = -1, y = 5/2, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 > 0, c("RF", "XGBoost")], 1, diff)/chance.error[pcat2 > 0, ], alternative = "less")$p.value, 2), digits = 1, format = "e")))
# data.text4 <- data.frame(x = -1, y = 2, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 > 0, c("RF", "RR-RF")], 1, diff)/chance.error[pcat2 > 0, ], alternative = "less")$p.value, 2), digits = 1, format = "e")))
# data.text5 <- data.frame(x = -1, y = 3/2, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 > 0, c("RF", "CCF")], 1, diff)/chance.error[pcat2 > 0, ], alternative = "less")$p.value, 2), digits = 1, format = "e")))

p2.hist <- ggplot(df[df$type == "Categorical", ]) +
  geom_vline(xintercept = 0, size = line.width, linetype = "dashed") +
  geom_density(aes(x = error.diff, color = pair), size = line.width) +
  # geom_density(aes(x = error.diff, color = type2), size = line.width) +
  geom_text(aes(x, y, label = lab), data = data.text2, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[1L]) +
  geom_text(aes(x, y, label = lab), data = data.text3, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[2L]) +
  geom_text(aes(x, y, label = lab), data = data.text4, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[3L]) +
  geom_text(aes(x, y, label = lab), data = data.text5, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[4L]) +
  scale_x_continuous(limits = c(-1, 1)) +
  ylim(0, 3.5) +
  xlab(expression("sign(RES)"~sqrt("|"~RES~"|"))) +
  ylab("Density Estimate") +
  ggtitle("Categorical") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = color.map) +
  guides(color = guide_legend(title.position = "top", title = "Alg Compared to RerF"))

  data.text2 <- data.frame(x = -1, y = 3, lab = paste0("MRES = ", formatC(signif(mean(apply(error.matrix[, c("RF", "RerF")], 1, diff)/chance.error[, ])), digits = 1, format = "e")))
  data.text3 <- data.frame(x = -1, y = 5/2, lab = paste0("MRES = ", formatC(signif(mean(apply(error.matrix[, c("RF", "XGBoost")], 1, diff)/chance.error[, ])), digits = 1, format = "e")))
  data.text4 <- data.frame(x = -1, y = 2, lab = paste0("MRES = ", formatC(signif(mean(apply(error.matrix[, c("RF", "RR-RF")], 1, diff)/chance.error[, ])), digits = 1, format = "e")))
  data.text5 <- data.frame(x = -1, y = 3/2, lab = paste0("MRES = ", formatC(signif(mean(apply(error.matrix[, c("RF", "CCF")], 1, diff)/chance.error[, ])), digits = 1, format = "e")))
# data.text2 <- data.frame(x = -1, y = 3, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[, c("RF", "RerF")], 1, diff)/chance.error, alternative = "less")$p.value, 2), digits = 1, format = "e")))
# data.text3 <- data.frame(x = -1, y = 5/2, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[, c("RF", "XGBoost")], 1, diff)/chance.error, alternative = "less")$p.value, 2), digits = 1, format = "e")))
# data.text4 <- data.frame(x = -1, y = 2, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[, c("RF", "RR-RF")], 1, diff)/chance.error, alternative = "less")$p.value, 2), digits = 1, format = "e")))
# data.text5 <- data.frame(x = -1, y = 3/2, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[, c("RF", "CCF")], 1, diff)/chance.error, alternative = "less")$p.value, 2), digits = 1, format = "e")))


p3.hist <- ggplot(df) +
  geom_vline(xintercept = 0, size = line.width, linetype = "dashed") +
  geom_density(aes(x = error.diff, color = pair), size = line.width) +
  geom_text(aes(x, y, label = lab), data = data.text2, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[1L]) +
  geom_text(aes(x, y, label = lab), data = data.text3, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[2L]) +
  geom_text(aes(x, y, label = lab), data = data.text4, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[3L]) +
  geom_text(aes(x, y, label = lab), data = data.text5, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[4L]) +
  scale_x_continuous(limits = c(-1, 1)) +
  ylim(0, 3.5) +
  xlab(expression("sign(RES)"~sqrt("|"~RES~"|"))) +
  ylab("Density Estimate") +
  ggtitle("All") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  # theme(axis.text = element_blank(),
  #       # axis.title = element_blank(),
  #       # legend.title = element_blank(),
  #       # axis.ticks = element_blank(),
  #       aspect.ratio = ar,
  #       # plot.margin = unit(c(mgn1, 0, mgn2, mgn3), "in"),
  #       panel.background = element_blank()) +
  scale_color_manual(values = color.map) +
  guides(color = guide_legend(title.position = "top", title = "Alg Compared to RF"))

pall <- ggarrange(p1.hist, p2.hist, p3.hist, ncol=3L, common.legend = TRUE, legend = "bottom")

# ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/error_histogram_benchmarks_v2.png", plot = pall, width = 12, height = 3.5, units = "in")



# alg rank heatmap
# alg.rank <- t(apply(mean.error, 1, order))
# colnames(alg.rank) <- colnames(mean.error)
alg.rank <- t(apply(mean.error, 1, function(rw) rank(rw, ties.method="min")))
alg.rank <- alg.rank[, ncol(alg.rank):1]
alg.rank <- cbind(alg.rank[, colnames(alg.rank) != "RerF"], alg.rank[, "RerF", drop = FALSE])

# Numeric Datasets
rank.df <- melt(apply(alg.rank[pcat == 0L, ], 2, tabulate), value.name = "Frequency", varnames = c("Rank", "Algorithm"))

mgn_top <- 0.2
mgn_right <- 0
mgn_bottom <- 0
mgn_left <- 0
ar <- 1/2

# data.text1 <- data.frame(x = 5.5, y = 4, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 == 0, c("RerF", "RF")], 1, diff)/chance.error[pcat2 == 0, ], alternative = "greater")$p.value, 2), digits = 1, format = "e")))
# data.text2 <- data.frame(x = 5.5, y = 3, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 == 0, c("RerF", "XGBoost")], 1, diff)/chance.error[pcat2 == 0, ], alternative = "greater")$p.value, 2), digits = 1, format = "e")))
# data.text3 <- data.frame(x = 5.5, y = 2, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 == 0, c("RerF", "RR-RF")], 1, diff)/chance.error[pcat2 == 0, ], alternative = "greater")$p.value, 2), digits = 1, format = "e")))
# data.text4 <- data.frame(x = 5.5, y = 1, lab = paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 == 0, c("RerF", "CCF")], 1, diff)/chance.error[pcat2 == 0, ], alternative = "greater")$p.value, 2), digits = 1, format = "e")))

labs1 <- levels(rank.df$Algorithm)
labs2 <- rep("", length(labs1))
for (i in 1:(length(labs2) - 1L)) {
  labs2[i] <- paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 == 0L, c("RerF", labs1[i])], 1, diff)/chance.error[pcat2 == 0L, ], alternative = "greater")$p.value, 2), digits = 1, format = "e"))
}
labs2[i + 1L] <- ""
rank.df$Algorithm <- as.numeric(rank.df$Algorithm)
p.rank_num <- ggplot(rank.df, aes(Rank, Algorithm)) +
  geom_tile(aes(fill = Frequency)) +
  geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5), size = 1) +
  scale_x_discrete(limits = as.character(1:5), expand=c(0,0)) +
  # scale_y_discrete(expand=c(0,0)) +
  scale_y_continuous(breaks = 1:length(labs1),
                     labels = labs1,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(labs2),
                                         labels = labs2),
                     expand=c(0,0)) +
  ylab("") +
  ggtitle("Numeric") +
  # geom_text(aes(x, y, label = lab), data = data.text1, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[1L]) +
  # geom_text(aes(x, y, label = lab), data = data.text2, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[2L]) +
  # geom_text(aes(x, y, label = lab), data = data.text3, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[3L]) +
  # geom_text(aes(x, y, label = lab), data = data.text4, hjust = 0, vjust = 0.5, fontface = "bold", color = color.map[4L]) +
  theme(plot.margin = unit(c(mgn_top, mgn_right, mgn_bottom, mgn_left), "in"),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient2(limits = c(min(rank.df$Frequency), max(rank.df$Frequency)),
                       midpoint = (max(rank.df$Frequency) - min(rank.df$Frequency))/2,
                       low="#40004b", mid="#f7f7f7", high="#00441b")
p.rank_num
# Categorical Datasets
rank.df <- melt(apply(alg.rank[pcat > 0L, ], 2, tabulate), value.name = "Frequency", varnames = c("Rank", "Algorithm"))

mgn_top <- 0.2
mgn_right <- 0
mgn_bottom <- 0
mgn_left <- 0
ar <- 1/2

labs2 <- rep("", length(labs1))
for (i in 1:(length(labs2) - 1L)) {
  labs2[i] <- paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[pcat2 > 0L, c("RerF", labs1[i])], 1, diff)/chance.error[pcat2 > 0L, ], alternative = "greater")$p.value, 2), digits = 1, format = "e"))
}
labs2[i + 1L] <- ""
rank.df$Algorithm <- as.numeric(rank.df$Algorithm)
p.rank_cat <- ggplot(rank.df, aes(Rank, Algorithm)) +
  geom_tile(aes(fill = Frequency)) +
  geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5), size = 1) +
  scale_x_discrete(limits = as.character(1:5), expand=c(0,0)) +
  # scale_y_discrete(expand=c(0,0)) +
  scale_y_continuous(breaks = 1:length(labs1),
                     labels = labs1,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(labs2),
                                         labels = labs2),
                     expand=c(0,0)) +
  ylab("") +
  ggtitle("Categorical") +
  theme(plot.margin = unit(c(mgn_top, mgn_right, mgn_bottom, mgn_left), "in"),
       legend.position = "bottom",
       plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient2(limits = c(min(rank.df$Frequency), max(rank.df$Frequency)),
                      midpoint = (max(rank.df$Frequency) - min(rank.df$Frequency))/2,
                      low="#40004b", mid="#f7f7f7", high="#00441b")
p.rank_cat

# All Datasets
rank.df <- melt(apply(alg.rank, 2, tabulate), value.name = "Frequency", varnames = c("Rank", "Algorithm"))

mgn_top <- 0.2
mgn_right <- 0
mgn_bottom <- 0
mgn_left <- 0
ar <- 1/2

labs2 <- rep("", length(labs1))
for (i in 1:(length(labs2) - 1L)) {
  labs2[i] <- paste0("p = ", formatC(signif(wilcox.test(apply(error.matrix[, c("RerF", labs1[i])], 1, diff)/chance.error[, ], alternative = "greater")$p.value, 2), digits = 1, format = "e"))
}
labs2[i + 1L] <- ""
rank.df$Algorithm <- as.numeric(rank.df$Algorithm)
p.rank_all <- ggplot(rank.df, aes(Rank, Algorithm)) +
  geom_tile(aes(fill = Frequency)) +
  geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5), size = 1) +
  scale_x_discrete(limits = as.character(1:5), expand=c(0,0)) +
  # scale_y_discrete(expand=c(0,0)) +
  scale_y_continuous(breaks = 1:length(labs1),
                     labels = labs1,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(labs2),
                                         labels = labs2),
                     expand=c(0,0)) +
  ylab("") +
  ggtitle("All") +
  theme(plot.margin = unit(c(mgn_top, mgn_right, mgn_bottom, mgn_left), "in"),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient2(limits = c(min(rank.df$Frequency), max(rank.df$Frequency)),
                       midpoint = (max(rank.df$Frequency) - min(rank.df$Frequency))/2,
                       low="#40004b", mid="#f7f7f7", high="#00441b")
p.rank_all

# gt2 <- ggplot_gtable(ggplot_build(p.rank))
#   gt2$layout$clip[gt2$layout$name == "panel"] <- "off"
#   grid.draw(gt2)

pall2 <- ggarrange(p.rank_num, p.rank_cat, p.rank_all, ncol=3L, vjust = 4, common.legend = F)

p.plot <- ggarrange(pall, pall2, nrow=2L, heights = c(1, 1), labels = c("(A)", "(B)"))

ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/error_histogram_benchmarks_v2.png", plot = p.plot, width = 16, height = 6.5, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/error_histogram_benchmarks_v2.pdf", plot = p.plot, width = 16, height = 6.6, units = "in")
