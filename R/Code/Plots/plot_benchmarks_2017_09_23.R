rm(list = ls())
library(reshape2)
library(ggplot2)
library(ks)

filePath <- "~/RerF/Results/2017.09.23/"

contents <- list.files(filePath)
load(paste0(filePath, contents[[1L]]))
fieldNames <- ls()
fieldNames <- fieldNames[fieldNames != "filePath" || fieldNames != "contents"]
res <- vector("list", length(fieldNames))
names(res) <- fieldNames
classifiers <- vector("list", length(contents))

for (f in contents) {
  load(paste0(filePath, f))
  classifiers[[strsplit(f, "_2017_09_23.RData")[[1L]]]] <- names(get(fieldNames[[1L]])[[1L]])
  for (fn in fieldNames) {
    res[[fn]][[strsplit(f, "_2017_09_23.RData")[[1L]]]] <- get(fn)[[1L]]
  }
}

classifiers <- unique(unlist(classifiers))

mean.error <- matrix(NA, nrow = length(res$testError), ncol = length(classifiers))

dataSets <- names(res$testError)
chance.error <- double(length(dataSets))
for (i in seq.int(length(res$testError))) {
  ds <- dataSets[i]
  D <- as.matrix(read.table(paste0("~/tmp/uci/processed/data/", ds, ".csv"), header = F, sep = ",", quote = "", row.names = NULL))
  Y <- as.integer(D[, ncol(D)])
  chance.error[i] <- 1 - max(tabulate(Y))/length(Y)
  clNames <- names(res$testError[[ds]])
  for (j in seq.int(length(clNames))) {
    cl <- clNames[j]
    mean.error[i, j] <- mean(sapply(1:length(res$bestIdx[[ds]][[cl]]), function (x) res$testError[[ds]][[cl]][x, res$bestIdx[[ds]][[cl]][x]]))
  }
}

no.na <- apply(mean.error, 1, function(x) !any(is.na(x)))
chance.error <- chance.error[no.na]
mean.error <- mean.error[no.na, ]

norm.rel.error <- (mean.error[, 2:ncol(mean.error)] - mean.error[, 1L])/chance.error
colnames(norm.rel.error) <- classifiers[2:length(classifiers)]

# mn <- min(norm.rel.error)
# mx <- max(norm.rel.error)
# eval.points <- seq(mn, mx, (mx - mn)/100)
# kscdf <- sapply(lapply(colnames(norm.rel.error), FUN = function(cname) kcde(x = norm.rel.error[, cname], eval.points = eval.points)), FUN = function(x) x$estimate)
# for (j in 1:ncol(kscdf)) {
#   kscdf[kscdf[, j] == 0 & cumsum(kscdf[, j]) > 0, j] <- 1
# }
# colnames(kscdf) <- colnames(norm.rel.error)
# # kscdf <- kcde(x = norm.rel.error[, 1], eval.points = eval.points)
# ccol <- rainbow(length(classifiers) - 1L)
# cdf.plot <- melt(kscdf)
# names(cdf.plot) <- c("relative.error", "classifier", "kcde")
# cdf.plot$relative.error <- eval.points[cdf.plot$relative.error]
# ggplot(cdf.plot, aes(x = relative.error, y = kcde)) +
#   scale_color_manual(values = ccol) +
#   geom_line(aes(group = classifier, colour = classifier))
# 
# df <- melt(as.matrix(norm.rel.error))
# names(df) <- c("ind","classifier","value")
# ggplot(df, aes(x=value)) + 
#   scale_color_manual(values=ccol) +
#   geom_density(aes(group=classifier, colour=classifier))

