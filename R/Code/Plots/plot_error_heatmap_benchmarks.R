rm(list = ls())
library(reshape2)
library(ggplot2)
library(ks)
library(meda)
library(ggridges)
library(ggpubr)
source("~/RandomerForest/R/Code/Utils/GetCatMap.R")

filePath <- "~/RandomerForest/R/Results/2017.10.03/"

contents <- list.files(filePath)
contents <- contents[!grepl("frc", contents) & !grepl("rerfc", contents) & !grepl("rr-rf", contents)]
catfiles <- list.files("~/tmp/uci/processed/categorical_map/")
load(paste0(filePath, contents[[1L]]))
fieldNames <- ls()
fieldNames <- fieldNames[(fieldNames != "filePath") & (fieldNames != "contents") & (fieldNames != "multiplot") &
                           (fieldNames != "catfiles") & (fieldNames != "GetCatMap")]
res <- vector("list", length(fieldNames))
names(res) <- fieldNames
classifiers <- vector("list", length(contents))

for (f in contents) {
  load(paste0(filePath, f))
  classifiers[[strsplit(f, "_2017_10_03.RData")[[1L]]]] <- names(get(fieldNames[[1L]])[[1L]])
  for (fn in fieldNames) {
    res[[fn]][[strsplit(f, "_2017_10_03.RData")[[1L]]]] <- get(fn)[[1L]]
  }
}

# classifiers <- unique(unlist(classifiers))
# classifiers <- c("rerf", "rerfr", "rf", "frc", "frank")
classifiers <- rev(c("rf", "rerfr", "frc"))
class.map <- list(rerf = "RerF", rerfr = "RerF(r)", rf = "RF", frc = "F-RC", frank = "Frank")
class.map2 <- class.map[names(class.map) != "rf"]

# color.map <- c("#41ab5d", "#4292c6", "#f768a1")
color.map <- brewer.pal(5L, "Set2")
line.width <- 0.75

dataSets <- names(res$testError)
mean.error <- matrix(NA, nrow = length(res$testError), ncol = length(classifiers))
rownames(mean.error) <- dataSets
colnames(mean.error) <- classifiers
train.time <- matrix(NA, nrow = length(res$testError), ncol = length(classifiers))
row.names(train.time) <- dataSets
colnames(train.time) <- classifiers
tree.strength <- matrix(NA, nrow = length(res$testError), ncol = length(classifiers))
row.names(tree.strength) <- dataSets
colnames(tree.strength) <- classifiers
tree.correlation <- matrix(NA, nrow = length(res$testError), ncol = length(classifiers))
row.names(tree.correlation) <- dataSets
colnames(tree.correlation) <- classifiers

chance.error <- double(length(dataSets))
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
  chance.error[i] <- 1 - max(tabulate(Y))/length(Y)
  # clNames <- names(res$testError[[ds]])
  for (j in seq.int(length(classifiers))) {
    cl <- classifiers[j]
    mean.error[i, j] <- mean(sapply(1:length(res$bestIdx[[ds]][[cl]]),
                                    function (x) res$testError[[ds]][[cl]][x, res$bestIdx[[ds]][[cl]][x]]))
    train.time[i, j] <- mean(sapply(1:length(res$bestIdx[[ds]][[cl]]),
                                    function (x) res$trainTime[[ds]][[cl]][x, res$bestIdx[[ds]][[cl]][x]]))
    tree.strength[i, j] <- mean(sapply(1:length(res$bestIdx[[ds]][[cl]]),
                                       function (x) res$treeStrength[[ds]][[cl]][x, res$bestIdx[[ds]][[cl]][x]]))
    tree.correlation[i, j] <- mean(sapply(1:length(res$bestIdx[[ds]][[cl]]),
                                          function (x) res$treeCorr[[ds]][[cl]][x, res$bestIdx[[ds]][[cl]][x]]))
  }
}

no.na <- apply(mean.error, 1, function(x) !any(is.na(x)))
chance.error <- chance.error[no.na]
mean.error <- mean.error[no.na, ]
train.time <- train.time[no.na, ]
tree.strength <- tree.strength[no.na, ]
tree.correlation <- tree.correlation[no.na, ]
n <- n[no.na]
p <- p[no.na]
pcat <- pcat[no.na]
dataSets <- dataSets[no.na]

norm.rel.error <- -sapply(colnames(mean.error), function(x) apply(mean.error[, colnames(mean.error) != x], 1, diff))/chance.error
colnames(norm.rel.error) <- c("RerF(r) - RF", "F-RC - RF", "F-RC - RerF(r)")
df <- melt(norm.rel.error, varnames = c("dataset", "classifier"), value.name = "error")
df$type <- rep(paste0("All ", length(dataSets)), nrow(df))
df <- rbind(df, cbind(melt(norm.rel.error[pcat > 0L, ], varnames = c("dataset", "classifier"), value.name = "error"), data.frame(type = rep(paste0(sum(pcat > 0L), " Categorical"), sum(pcat > 0L)))))
df <- rbind(df, cbind(melt(norm.rel.error[pcat == 0L, ], varnames = c("dataset", "classifier"), value.name = "error"), data.frame(type = rep(paste0(sum(pcat == 0L), " Numeric"), sum(pcat == 0L)))))
df$type <- factor(df$type)
df$type <- factor(df$type, levels = rev(levels(df$type)))

breaks <- c(-1, -0.75, seq(-0.5, 0.5, 0.1), 0.75, 1)

p <- ggplot(df, aes(error)) +
  geom_histogram(position = "identity", size = 0.5, breaks = breaks, alpha = 0.5, color = "#8da0cb", fill = "#8da0cb") +
  facet_grid(classifier ~ type) +
  xlab("Relative Error") +
  ylab("Count") +
  ggtitle(paste0(nrow(norm.rel.error), " UCI Datasets")) +
  # theme_bw() +
  theme(text = element_text(size = 14), aspect.ratio = 0.5)
  # scale_fill_manual(values = color.map[c(3, 1)]) +
  # scale_color_manual(values = color.map[c(3, 1)])
  
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/error_heatmap_benchmarks.png", plot = pall, width = 10, height = 2, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/error_heatmap_benchmarks.pdf", plot = pall, width = 10, height = 2, units = "in")

# norm.error <- mean.error/chance.error

# norm.rel.error <- (mean.error[, colnames(mean.error) != "rf"] - mean.error[, "rf"])/chance.error
# colnames(norm.rel.error) <- classifiers[classifiers != "rf"]

# all datasets
x <- rep(0, (length(breaks) - 1L)*3L*length(classifiers[classifiers != "rf"]))
y <- rep(0, (length(breaks) - 1L)*3L*length(classifiers[classifiers != "rf"]))
classifier <- rep(0, (length(breaks) - 1L)*3L*length(classifiers[classifiers != "rf"]))
for (j in seq_along(classifiers[classifiers != "rf"])) {
  cl <- classifiers[classifiers != "rf"][j]
  h <- hist(norm.rel.error[, cl], breaks = breaks)
  x[((j - 1L)*(length(breaks) - 1L)*3L + 1L):(j*(length(breaks) - 1L)*3L)] <- c(sapply(1:(length(h$breaks) - 1L), function(x) c(h$breaks[x], rep(h$breaks[x + 1L], 2L))))
  y[((j - 1L)*(length(breaks) - 1L)*3L + 1L):(j*(length(breaks) - 1L)*3L)] <- c(sapply(h$counts, function(x) c(rep(x, 2L), NA)))
  classifier[((j - 1L)*(length(breaks) - 1L)*3L + 1L):(j*(length(breaks) - 1L)*3L)] <- rep(class.map[[cl]], (length(breaks) - 1L)*3L)
}

df <- data.frame(x = x, y = y, classifier = factor(classifier, levels = c("RerF", "RerF(r)", "F-RC", "Frank")))

p1 <- ggplot(df, aes(x = x, y = y, color = classifier)) +
  geom_line(size = line.width) +
  xlab("Error Relative to RF") +
  ylab("Frequency") +
  ggtitle(paste0(nrow(norm.rel.error), " UCI Datasets"), subtitle = "all datasets") +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 0.5) +
  scale_color_manual(values = color.map[c(1, 5, 3, 4)])

# categorical datasets
x <- rep(0, (length(breaks) - 1L)*3L*length(classifiers[classifiers != "rf"]))
y <- rep(0, (length(breaks) - 1L)*3L*length(classifiers[classifiers != "rf"]))
classifier <- rep(0, (length(breaks) - 1L)*3L*length(classifiers[classifiers != "rf"]))
for (j in seq_along(classifiers[classifiers != "rf"])) {
  cl <- classifiers[classifiers != "rf"][j]
  h <- hist(norm.rel.error[pcat > 0L, cl], breaks = breaks)
  x[((j - 1L)*(length(breaks) - 1L)*3L + 1L):(j*(length(breaks) - 1L)*3L)] <- c(sapply(1:(length(h$breaks) - 1L), function(x) c(h$breaks[x], rep(h$breaks[x + 1L], 2L))))
  y[((j - 1L)*(length(breaks) - 1L)*3L + 1L):(j*(length(breaks) - 1L)*3L)] <- c(sapply(h$counts, function(x) c(rep(x, 2L), NA)))
  classifier[((j - 1L)*(length(breaks) - 1L)*3L + 1L):(j*(length(breaks) - 1L)*3L)] <- rep(class.map[[cl]], (length(breaks) - 1L)*3L)
}

df <- data.frame(x = x, y = y, classifier = classifier)

p2 <- ggplot(df, aes(x = x, y = y, color = classifier)) +
  geom_line(size = line.width) +
  xlab("Error Relative to RF") +
  ylab("Frequency") +
  ggtitle(paste0(sum(pcat > 0L), " UCI Datasets"), subtitle = "categorical datasets") +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 0.5) +
  scale_color_manual(values = color.map[c(1, 5, 3, 4)])

# numeric datasets
x <- rep(0, (length(breaks) - 1L)*3L*length(classifiers[classifiers != "rf"]))
y <- rep(0, (length(breaks) - 1L)*3L*length(classifiers[classifiers != "rf"]))
classifier <- rep(0, (length(breaks) - 1L)*3L*length(classifiers[classifiers != "rf"]))
for (j in seq_along(classifiers[classifiers != "rf"])) {
  cl <- classifiers[classifiers != "rf"][j]
  h <- hist(norm.rel.error[pcat == 0L, cl], breaks = breaks)
  x[((j - 1L)*(length(breaks) - 1L)*3L + 1L):(j*(length(breaks) - 1L)*3L)] <- c(sapply(1:(length(h$breaks) - 1L), function(x) c(h$breaks[x], rep(h$breaks[x + 1L], 2L))))
  y[((j - 1L)*(length(breaks) - 1L)*3L + 1L):(j*(length(breaks) - 1L)*3L)] <- c(sapply(h$counts, function(x) c(rep(x, 2L), NA)))
  classifier[((j - 1L)*(length(breaks) - 1L)*3L + 1L):(j*(length(breaks) - 1L)*3L)] <- rep(class.map[[cl]], (length(breaks) - 1L)*3L)
}

df <- data.frame(x = x, y = y, classifier = classifier)

p3 <- ggplot(df, aes(x = x, y = y, color = classifier)) +
  geom_line(size = line.width) +
  xlab("Error Relative to RF") +
  ylab("Frequency") +
  ggtitle(paste0(sum(pcat == 0L), " UCI Datasets"), subtitle = "non-categorical datasets") +
  theme_classic() +
  theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 0.5) +
  scale_color_manual(values = color.map[c(1, 5, 3, 4)])


# df <- data.frame(classifier = c(sapply(class.map[class.map != "RF"], function(x) rep(x, length(dataSets)))),
#                  error = c(norm.rel.error),
#                  categorical = rep(pcat > 0L, length(class.map) - 1L))
# 
# p1 <- ggplot(df, aes(x = error, fill = classifier, color = classifier)) +
#   geom_histogram(position = "identity", size = 0.5, breaks = c(-1, -0.75, -0.5, seq(-0.4, 0.4, 0.05), 0.5, 0.75, 1), alpha = 0.5) +
#   xlab("Error Relative to RF") +
#   ylab("Frequency") +
#   ggtitle(paste0(nrow(norm.rel.error), " UCI Datasets")) +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 0.5, plot.title = element_text(face = "bold")) +
#   scale_fill_manual(values = color.map[c(3, 1)]) +
#   scale_color_manual(values = color.map[c(3, 1)])
# 
# p2 <- ggplot(df[df$categorical, ], aes(x = error, fill = classifier, color = classifier)) +
#   geom_histogram(position = "identity", size = 0.5, breaks = c(-1, -0.75, -0.5, seq(-0.4, 0.4, 0.05), 0.5, 0.75, 1), alpha = 0.5) +
#   xlab("Error Relative to RF") +
#   ylab("Frequency") +
#   ggtitle(paste0(sum(pcat > 0L), " UCI Datasets"), subtitle = "categorical datasets") +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 0.5, plot.title = element_text(face = "bold")) +
#   scale_fill_manual(values = color.map[c(3, 1)]) +
#   scale_color_manual(values = color.map[c(3, 1)])
# 
# p3 <- ggplot(df[!df$categorical, ], aes(x = error, fill = classifier, color = classifier)) +
#   geom_histogram(position = "identity", size = 0.5, breaks = c(-1, -0.75, -0.5, seq(-0.4, 0.4, 0.05), 0.5, 0.75, 1), alpha = 0.5) +
#   xlab("Error Relative to RF") +
#   ylab("Frequency") +
#   ggtitle(paste0(sum(pcat == 0L), " UCI Datasets"), subtitle = "non-categorical datasets") +
#   theme_classic() +
#   theme(axis.line = element_line(size = line.width), text = element_text(size = 14), aspect.ratio = 0.5, plot.title = element_text(face = "bold")) +
#   scale_fill_manual(values = color.map[c(3, 1)]) +
#   scale_color_manual(values = color.map[c(3, 1)])
#   
# p1 <- ggplot(df, aes(x = error, y = classifier, height = ..density..)) +
#   geom_density_ridges(stat = "binline", bins = 8, aes(fill = classifier)) +
#   xlab("Error Relative to RF") +
#   ylab("") +
#   ggtitle(paste0(nrow(norm.rel.error), " UCI Datasets")) +
#   scale_fill_manual(values = color.map)
# 
# p2 <- ggplot(df[df$categorical, ], aes(x = error, y = classifier, height = ..density..)) +
#   geom_density_ridges(stat = "binline", bins = 8, aes(fill = classifier)) +
#   xlab("Error Relative to RF") +
#   ylab("") +
#   ggtitle(paste0(nrow(norm.rel.error), " UCI Datasets")) +
#   scale_fill_manual(values = color.map)
# 
# p3 <- ggplot(df[!df$categorical, ], aes(x = error, y = classifier, height = ..density..)) +
#   geom_density_ridges(stat = "binline", bins = 8, aes(fill = classifier)) +
#   xlab("Error Relative to RF") +
#   ylab("") +
#   ggtitle(paste0(nrow(norm.rel.error), " UCI Datasets")) +
#   scale_fill_manual(values = color.map)

# hmap <- d1heat(norm.rel.error, trunc = NA, breaks = seq(-1, 1, 0.25))
# count.min <- min(hmap$dat$Count)
# count.max <- max(hmap$dat$Count)
# 
# p1 <- plot(hmap, bincount = T, limits = c(count.min, count.max)) +
#   geom_hline(yintercept = 1.5, size = line.width) +
#   xlab("Error Relative to RF") +
#   ylab("") +
#   ggtitle(paste0(nrow(norm.rel.error), " UCI Datasets")) +
#   scale_y_discrete(labels = unlist(class.map)) +
#   theme(plot.margin = unit(c(0.095, 0, 0.02, 0), units = "npc"))
# 
# hmap <- d1heat(norm.rel.error[pcat > 0L, ], trunc = NA, breaks = seq(-1, 1, 0.25))
# p2 <- plot(hmap, bincount = T, limits = c(count.min, count.max)) +
#   geom_hline(yintercept = 1.5, size = line.width) +
#   xlab("Error Relative to RF") +
#   ylab("") +
#   ggtitle(paste0(sum(pcat > 0L), " UCI Datasets"), subtitle = "categorical datasets") +
#   scale_y_discrete(labels = unlist(class.map))
# 
# hmap <- d1heat(norm.rel.error[pcat == 0L, ], trunc = NA, breaks = seq(-1, 1, 0.25))
# p3 <- plot(hmap, bincount = T, limits = c(count.min, count.max)) +
#   geom_hline(yintercept = 1.5, size = line.width) +
#   xlab("Error Relative to RF") +
#   ylab("") +
#   ggtitle(paste0(sum(pcat == 0L), " UCI Datasets"), subtitle = "non-categorical datasets") +
#   scale_y_discrete(labels = unlist(class.map))
# 
# 
pall <- ggarrange(p1, p2, p3, ncol = 3L, widths = c(1, 1, 1), legend = "right", common.legend = T, labels = "AUTO")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/error_heatmap_benchmarks.png", plot = pall, width = 10, height = 2, units = "in")
ggsave(filename = "~/RandomerForest/Draft/PAMI/Figures/error_heatmap_benchmarks.pdf", plot = pall, width = 10, height = 2, units = "in")
