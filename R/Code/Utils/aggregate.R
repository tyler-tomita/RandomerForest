rm(list = ls())

# files.frc <- list.files(pattern = "frc")
files.rerfc <- list.files(path = "~/RandomerForest/R/Results/2017.10.03/", pattern = "2017_10_03.RData")
# files.rrrfr <- list.files(pattern = "rr-rfr")
for (f in files.rerfc) {
  load(paste0("~/RandomerForest/R/Results/2017.10.03/", f))
  if (f == files.rerfc[1L]) {
    fieldNames <- ls()
    fieldNames <- fieldNames[(fieldNames != "files.rerfc") & (fieldNames != "f")]
  }
  res <- vector("list", length(fieldNames))
  names(res) <- fieldNames
  for (fn in fieldNames) {
    res[[fn]] <- get(fn)
  }
  # if (!("rf" %in% names(res[[fieldNames[1]]][[1]]))) {
    ds <- names(get(fn))
    load(paste0(ds, "_2017_09_23.RData"))
    res2 <- vector("list", length(fieldNames))
    for (fn in fieldNames) {
      res2[[fn]] <- get(fn)
      clNames <- names(res[[fn]][[ds]])
      # clNames <- c("rerfc", "rerfcr")
      for (cl in clNames) {
        res2[[fn]][[ds]][[cl]] <- res[[fn]][[ds]][[cl]]
      }
      assign(fn, res2[[fn]])
    }
    save(testError, OOBError, OOBAUC, trainTime, OOBTime, testTime, treeStrength, treeCorr, numNodes, bestIdx, params,
         file = paste0("~/RandomerForest/R/Results/2017.10.03/", ds, "_2017_10_03_aggregated.RData")) 
  # }
}
