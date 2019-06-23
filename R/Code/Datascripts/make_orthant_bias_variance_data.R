# make orthant bias variance data
rm(list = ls())

p <- 6L
ns <- c(200L, 400L, 2000L, 4000L)
n.test <- 10000L
num.trials <- 100L

# distribution parameters
num.orthants = 2L^p;
num.nodes = num.orthants*2L - 1L;
num.levels = p + 1L;
num.classes = num.orthants;
tree.node = 1:num.nodes
tree.cutpoint = rep(as.double(NA), num.nodes)
tree.cutvar = rep(as.integer(NA), num.nodes)
tree.children = matrix(rep(as.integer(NA), num.nodes*2L), num.nodes, 2L)
tree.parent = matrix(rep(as.integer(NA), num.nodes*2L), num.nodes, 2L)
for (level in 1:num.levels) {
  nodes = (2L^(level - 1L)):(2L^level - 1L)
  if (level != num.levels) {
    for (k in 1:length(nodes)) {
      node <- nodes[k]
      tree.cutpoint[node] <- 0
      tree.cutvar[node] <- level
      tree.children[node, ] <- c(node*2L,node*2L + 1L)
      tree.parent[tree.children[node, ]] <- node
    }
  }
}


# test set
X <- matrix(runif(n.test*p, min = -1, max = 1), n.test, p)
Y <- rep(0L, n.test)
posteriors <- matrix(0, n.test, num.classes)
node.observations <- vector("list", num.nodes)
node.observations[[1L]] <- 1:n.test
for (level in 1:num.levels) {
  nodes <- (2L^(level - 1L)):(2L^level - 1L)
  if (level != num.levels) {
    for (k in 1:length(nodes)) {
      node <- nodes[k]
      move.left <- X[node.observations[[node]], tree.cutvar[node]] <= tree.cutpoint[node]
      node.observations[[tree.children[node, 1L]]] <- node.observations[[node]][move.left]
      node.observations[[tree.children[node, 2L]]] <- node.observations[[node]][!move.left]
    }
  } else {
    for (k in 1:length(nodes)) {
      node <- nodes[k]
      Y[node.observations[[node]]] <- k
    }
  }
}
for (k in 1:num.classes) {
  posteriors[, k] <- as.double(Y == k)
}

write.table(cbind(X, Y),
            file = paste0("~/R/Data/Orthant_bias_variance/Test/Orthant_bias_variance_test.csv"),
            sep = ",", row.names = F, col.names = F)

write.table(posteriors,
            file = paste0("~/R/Data/Orthant_bias_variance/Test/Orthant_bias_variance_test_posteriors.csv"),
            sep = ",", row.names = F, col.names = F)


for (i in 1:length(ns)) {
  n.train <- ns[i]
  print(paste0("n = ", n.train))
  for (trial in 1:num.trials) {
    # training set
    X <- matrix(runif(n.train*p, min = -1, max = 1), n.train, p)
    Y <- rep(0L, n.train)
    node.observations <- vector("list", num.nodes)
    node.observations[[1L]] <- 1:n.train
    for (level in 1:num.levels) {
      nodes <- (2L^(level - 1L)):(2L^level - 1L)
      if (level != num.levels) {
        for (k in 1:length(nodes)) {
          node <- nodes[k]
          move.left <- X[node.observations[[node]], tree.cutvar[node]] <= tree.cutpoint[node]
          node.observations[[tree.children[node, 1L]]] <- node.observations[[node]][move.left]
          node.observations[[tree.children[node, 2L]]] <- node.observations[[node]][!move.left]
        }
      } else {
        for (k in 1:length(nodes)) {
          node <- nodes[k]
          Y[node.observations[[node]]] <- k
        }
      }
    }
    
    write.table(cbind(X, Y),
                file = paste0("~/R/Data/Orthant_bias_variance/Train/Orthant_bias_variance_train_n", n.train, "_trial", trial, ".csv"),
                sep = ",", row.names = F, col.names = F)
  }
}