source("../src/phySigs.R")

library(graph)

sim_path  <- "../data/simulations/"
instances <- list.files(sim_path, pattern="*.P_eps.txt")

parseCountMatrix <- function(filename) {
  feat_mat <- t(read.csv(filename, header = TRUE, row.names=1, check.names=FALSE))
  return(as.data.frame(feat_mat))
}

parseTree <- function(filename, V) {
  tree <- new("graphNEL", nodes=V, edgemode="directed")
  t <- read.table(filename)
  for (i in 1:nrow(t)) {
    tree <- addEdge(as.character(t[i,1]), as.character(t[i,2]), tree, 1)
  }
  return(tree)
}

exp_list_sims <- list()
for (instance in instances) {
  exp_list_sims[[instance]] <- list()
  print(instance)
  s <- unlist(strsplit(instance, "_"))
  treeFilename <- paste0(sim_path, s[[1]], "_", s[[2]], ".tree")
  
  feat_mat2 <- parseCountMatrix(paste0(sim_path, instance))
  tree <- parseTree(treeFilename, rownames(feat_mat2))
  
  nrEdges <- ncol(edgeMatrix(tree))
  exp_list_sims[[instance]] <- allTreeExposures(tree, feat_mat2, row.names(signatures.cosmic))
}

# infer errors for all trees
error_list_sims <- list()
BIC_list_sims <- list()
for (instance in instances) {
  error_list_sims[[instance]] <- list()
  BIC_list_sims[[instance]] <- list()
  
  print(instance)
  s <- unlist(strsplit(instance, "_"))
  treeFilename <- paste0(sim_path, s[[1]], "_", s[[2]], ".tree")
  
  feat_mat2 <- parseCountMatrix(paste0(sim_path, instance))
  tree <- parseTree(treeFilename, rownames(feat_mat2))

  nrEdges  <- length(nodes(tree)) - 1
  for (k in 0:nrEdges) {
    print(paste(instance, nrEdges, k))
    error <- getError(feat_mat2, 
                      exp_list_sim[[instance]][[as.character(k)]], 
                      row.names(signatures.cosmic))
    error_list_sims[[instance]][[as.character(k)]] <- error
    BIC_list_sims[[instance]][[as.character(k)]] <- getBIC(feat_mat2, 
                                                           exp_list_sims[[instance]][[as.character(k)]], 
                                                           row.names(signatures.cosmic))
  }
}