rm(list=ls())
options(scipen = 999)

load("~/R-RerF/Results/2017.06.13/Orthant/Orthant_rerf_2017_06_13.RData")

ps <- c(2,4,6)
ns <- list(c(20,200,400), c(80,400,4000), c(400,2000,4000))

for (j in 1:3) {
  p <- ps[j]
  for (i in 1:3) {
    n <- ns[[j]][i]
    write.table(testError[i,j,], file = paste("~/R-RerF/Results/2017.06.13/Orthant/Orthant_testError_n",n,"_p",p,".dat", sep = ""), row.names = FALSE, col.names = FALSE)
    write.table(trainTime[i,j,], file = paste("~/R-RerF/Results/2017.06.13/Orthant/Orthant_trainTime_n",n,"_p",p,".dat", sep = ""), row.names = FALSE, col.names = FALSE)
  }
}