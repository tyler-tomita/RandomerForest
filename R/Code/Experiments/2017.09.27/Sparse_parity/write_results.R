rm(list=ls())
options(scipen = 999)

load("~/R-RerF/Results/2017.06.13/Sparse_parity/Sparse_parity_rerf_2017_06_13_partial.RData")

ps <- c(3, 10, 20)
ns <- list(c(10, 100, 1000), c(100, 1000, 10000), c(1000, 5000, 10000))

for (j in 1:3) {
  p <- ps[j]
  for (i in 1:3) {
    n <- ns[[j]][i]
    write.table(testError[i,j,], file = paste("~/R-RerF/Results/2017.06.13/Sparse_parity/Sparse_parity_testError_n",n,"_p",p,".dat", sep = ""), row.names = FALSE, col.names = FALSE)
    write.table(trainTime[i,j,], file = paste("~/R-RerF/Results/2017.06.13/Sparse_parity/Sparse_parity_trainTime_n",n,"_p",p,".dat", sep = ""), row.names = FALSE, col.names = FALSE)
  }
}