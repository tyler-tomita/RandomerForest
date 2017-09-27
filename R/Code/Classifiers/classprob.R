classprob <- function(forest, Xtest) {
  if(!is.null(Forest$forest)){
    Forest<-Forest$forest
  }
  n <- nrow(Xtest)
  nTrees <- length(Forest)
  nClasses <- length(forest[[1]]$ClassProb[1,])
  classProb <- matrix(0, nrow = n, ncol = nClasses)
  for(i in 1:n){
    cp <- vector(mode = "numeric", length = nClasses)
    for(j in 1:nTrees){
      Tree <- Forest[[j]]
      currentNode <- 1L
      while(Tree$Children[currentNode]!=0L){
        rotX <- 0
        for(s in 1:(length(Tree$matA[[currentNode]])/2)){
          rotX<-rotX+Tree$matA[[currentNode]][2*s]*Xtest[i,Tree$matA[[currentNode]][2*s-1]]
        }
        if(rotX<=Tree$CutPoint[currentNode]){
          currentNode <- Tree$Children[currentNode,1L]
        }else{
          currentNode <- Tree$Children[currentNode,2L]
        }
      }
      cp <- cp + Tree$ClassProb[currentNode,]
    }
    classProb[i, ] <- cp
  }
  return(classProb/nTrees)
}