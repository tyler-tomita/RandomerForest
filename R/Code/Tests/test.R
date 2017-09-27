
sink("candidate.log")
trees <- 500
numT <- 50
numC <- 25

pvals <- c(.5,2,.5,2, .5,2,.5,1.35)

#The 8th test must be manually changed.  This is the long MNIST.
numTests <- c(rep(numT, 7), 50)
cores <- c(rep(numC, 7), 25)
numTrees <- c(rep(trees, 7), 500)
#This skips the 8th test.  If you don't want to skip the test #then make skipMNIST FALSE
skipMNIST <- TRUE

trainSets <- c("Orthant_train.csv", "Sparse_parity_train.csv", "Trunk_train.csv", "mnist_train.csv")
testSets <- c("Orthant_test.csv", "Sparse_parity_test.csv", "Trunk_test.csv", "mnist_test.csv")

pass <- TRUE
#RerF_baseline <- "baseline_RFR.R"
#RerF_candidate <- "candidate_RFR.R"

RerF_baseline <- "candidate_RFR.R"
RerF_candidate <- "working_rfr.R"

testDS <- function(RerF_baseline, RerF_candidate, trainSet, testSet, numTests, numTrees, cores, p1){
    X <- read.csv(file(trainSet), header=FALSE)
    Y <- as.numeric(X[,ncol(X)])
    X <- as.matrix(X[,1:(ncol(X)-1)])

    Ymod <- 1 - as.numeric(min(levels(as.factor(Y))))
    Y <- Y+Ymod

    Xtest <- read.csv(file(testSet), header=FALSE)
    Ytest <- as.numeric(Xtest[,ncol(Xtest)])+Ymod
    Xtest <- as.matrix(Xtest[,1:(ncol(Xtest)-1)])

    pass <- TRUE

    if(exists("comp_err")){rm("comp_err", envir=.GlobalEnv)}
    if(exists("comp_errOOB")){rm("comp_errOOB", envir=.GlobalEnv)}
    if(exists("comp_rfr")){rm("comp_rfr", envir=.GlobalEnv)}
    source(RerF_baseline)
    #    source("~/dropbox/gitRepos/R-RerF/code/rfr_function.R")
    ##################  Baseline, mtry = p^.5

    ptmtrain_baseline <- NA
    ptmtest_baseline <- NA
    error_baseline <- NA
    OOBerror_baseline <- NA
    NodeSize_baseline <- NA
    memSize_baseline <- NA

    for (i in 1:numTests){
        initMem <- sum(gc(reset=TRUE)[9:10])
        ptm <- proc.time()
        forest<-rfr(X,Y,trees=numTrees, MinParent=2L, MaxDepth="inf", stratify=TRUE, COOB=TRUE, NumCores=cores,  options=c(ncol(X), ceiling(ncol(X)^p1),1L, 1/ncol(X)), seed = sample(1:10000,1)) 
        ptmtrain_baseline[i]<- (proc.time() - ptm)[3]
        if (length(forest)!=numTrees){
            print(paste("rfr failed to make desired number of trees.  ", length(forest) , " trees made.  ", numTrees , " trees requested."))
            pass <- FALSE
        }
        temp_size <- 0
        for(z in 1:length(forest)){
            temp_size <- temp_size + length(forest[[z]]$ClassProb[,1])
            for(q in 1:length(forest[[z]]$ClassProb[,1])){
                if (sum(forest[[z]]$ClassProb[q,]) == 0){
                    print("A node with 0 samples exists")
                    pass <- FALSE
                }
            }
        }
        NodeSize_baseline[i] <- temp_size/length(forest)
        capture.output(OOBmat_temp <- OOBpredict(X,Y,forest, cores))
        OOBmat_cols <- length(OOBmat_temp[[length(OOBmat_temp)]][1,])
        numWrong<- 0L
        numTotal<- 0L

        for(k in 1:nrow(X)){ 
            if(any(OOBmat_temp[[length(OOBmat_temp)]][k,3:OOBmat_cols]!=0)){
                if(order(OOBmat_temp[[length(OOBmat_temp)]][k,3:OOBmat_cols],decreasing=T)[1L]!=Y[k]){
                    numWrong <- numWrong+1L
                }
                numTotal<-numTotal+1L
            }
        }
        OOBerror_baseline[i] <- 100*numWrong/numTotal
        ptm <- proc.time()
        error_baseline[i] <- error_rate(Xtest,Ytest,forest, NumCores = cores)
        ptmtest_baseline[i]<- (proc.time() - ptm)[3]
        memSize_baseline[i] <- sum(gc()[9:10]) - initMem
    }
    #####################################################################################################################
    #####################################################################################################################
    if(exists("comp_err")){rm("comp_err", envir=.GlobalEnv)}
    if(exists("comp_errOOB")){rm("comp_errOOB", envir=.GlobalEnv)}
    if(exists("comp_rfr")){rm("comp_rfr", envir=.GlobalEnv)}
    source(RerF_candidate)
    #source("~/dropbox/gitRepos/R-RerF/rfr_function.R")
    #source("~/dropbox/gitRepos/R-RerF/code/rfr_function.R")

    ##################  Candidate, mtry = p^.5

    ptmtrain_candidate <- NA
    ptmtest_candidate <- NA
    error_candidate <- NA
    OOBerror_candidate <- NA
    NodeSize_candidate <- NA
    memSize_candidate <- NA

    for (i in 1:numTests){
        initMem <- sum(gc(reset=TRUE)[9:10])
        ptm <- proc.time()
        forest<-rfr(X,Y,trees=numTrees, MinParent=2L, MaxDepth="inf", stratify=TRUE, COOB=TRUE, NumCores=cores,  options=c(ncol(X), ceiling(ncol(X)^p1),1L, 1/ncol(X)), seed = sample(1:10000,1)) 
        ptmtrain_candidate[i]<- (proc.time() - ptm)[3]
        if (length(forest)!=numTrees){
            print(paste("rfr failed to make desired number of trees.  ", length(forest) , " trees made.  ", numTrees , " trees requested."))
            pass <- FALSE
        }
        temp_size <- 0
        for(z in 1:length(forest)){
            temp_size <- temp_size + length(forest[[z]]$ClassProb[,1])
            for(q in 1:length(forest[[z]]$ClassProb[,1])){
                if (sum(forest[[z]]$ClassProb[q,]) == 0){
                    print("A node with 0 samples exists")
                    pass <- FALSE
                }
            }
        }
        NodeSize_candidate[i] <- temp_size/length(forest)
        capture.output(OOBmat_temp <- OOBpredict(X,Y,forest, cores))
        OOBmat_cols <- length(OOBmat_temp[[length(OOBmat_temp)]][1,])
        numWrong<- 0L
        numTotal<- 0L

        for(k in 1:nrow(X)){ 
            if(any(OOBmat_temp[[length(OOBmat_temp)]][k,3:OOBmat_cols]!=0)){
                if(order(OOBmat_temp[[length(OOBmat_temp)]][k,3:OOBmat_cols],decreasing=T)[1L]!=Y[k]){
                    numWrong <- numWrong+1L
                }
                numTotal<-numTotal+1L
            }
        }
        OOBerror_candidate[i] <- 100*numWrong/numTotal
        ptm <- proc.time()
        error_candidate[i] <- error_rate(Xtest,Ytest,forest,NumCores = cores)
        ptmtest_candidate[i]<- (proc.time() - ptm)[3]
        memSize_candidate[i] <- sum(gc()[9:10]) - initMem
    }
    NodeSizeRatio = median(NodeSize_candidate)/median(NodeSize_baseline)
    memSizeRatio = median(memSize_baseline)/median(memSize_candidate)
    OOBerrorRatio = median(OOBerror_baseline)/median(OOBerror_candidate)
    errorRatio = median(error_baseline)/median(error_candidate)
    ptmtrainRatio = median(ptmtrain_baseline)/median(ptmtrain_candidate)
    ptmtestRatio = median(ptmtest_baseline)/median(ptmtest_candidate)

    print("")
    print(paste(trainSet, "pmod= ", p1))
    wt <- wilcox.test(NodeSize_baseline, NodeSize_candidate)$p.value
    print(paste("Node Size, pval(2sided)= ", wt, ", ratio(C/B)= ", NodeSizeRatio, ", pass=", wt>.001))
    if(!is.nan(wt)){
        if(wt<.001){
            pass <- FALSE
        }
    }

    wt <- wilcox.test(OOBerror_baseline, OOBerror_candidate)$p.value
    print(paste("OOBerror, pval(2sided)= ", wt, ", ratio(B/C)= ", OOBerrorRatio, ", pass=", wt>.001))
    if(!is.nan(wt)){
        if(wt<.001){
            pass <- FALSE
        }
    }

    wt <- wilcox.test(error_baseline, error_candidate)$p.value
    print(paste("Error, pval(2sided)= ", wt, ", ratio(B/C)= ", errorRatio, ", pass=", wt>.001))
    if(!is.nan(wt)){
        if(wt<.001){
            pass <- FALSE
        }
    }

    wt <- wilcox.test(ptmtrain_baseline, ptmtrain_candidate, alternative="l")$p.value
    print(paste("Train Time, pval(l)= ", wt, ", ratio(B/C)= ", ptmtrainRatio, ", pass=", wt>.001))
    if(!is.nan(wt)){
        if(wt<.001){
            pass <- FALSE
        }
    }

    wt <- wilcox.test(ptmtest_baseline, ptmtest_candidate, alternative="l")$p.value
    print(paste("Test Time, pval(l)= ", wt, ", ratio(B/C)= ", ptmtestRatio, ", pass=", wt>.001))
    if(!is.nan(wt)){
        if(wt<.001){
            pass <- FALSE
        }
    }

    print(paste("mem ratio 1 (B/C): ", memSizeRatio))
    flush.console()
    return(pass)
}


lh <- NA
for (z in 1:(2*length(trainSets)-skipMNIST)){
    gc()
    lh[z] <- testDS(RerF_baseline,RerF_candidate, trainSets[ceiling(z/2)], testSets[ceiling(z/2)], numTests[z], numTrees[z], cores[z], pvals[z])
    if (lh[z]==FALSE){
        pass <- FALSE
        print(paste(trainSets[ceiling(z/2)], " failed"))
    }
}

if(pass){
    print("All tests passed.")
}else{
    print("Candidate failed.")
}
sink()

