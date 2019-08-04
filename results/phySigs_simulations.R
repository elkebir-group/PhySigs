source("../src/phySigs.R")
load("sims.RData")

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

#exp_list_sims <- list()
for (instance in instances) {
  if (!(instance %in% names(exp_list_sims))) {
    exp_list_sims[[instance]] <- list()
    print(instance)
    s <- unlist(strsplit(instance, "_"))
    treeFilename <- paste0(sim_path, s[[1]], "_", s[[2]], ".tree")
    
    feat_mat2 <- parseCountMatrix(paste0(sim_path, instance))
    tree <- parseTree(treeFilename, rownames(feat_mat2))
    
    nrEdges <- ncol(edgeMatrix(tree))
    exp_list_sims[[instance]] <- allTreeExposures(tree, feat_mat2, row.names(signatures.cosmic))
  }
}

# l <- foreach (i=1:length(instances), combine=c) %dopar% {
#   instance <- instances[i]
#   #if (!(instance %in% names(exp_list_sims))) {
#     result <- list()
#     print(instance)
#     s <- unlist(strsplit(instance, "_"))
#     treeFilename <- paste0(sim_path, s[[1]], "_", s[[2]], ".tree")
#     
#     feat_mat2 <- parseCountMatrix(paste0(sim_path, instance))
#     tree <- parseTree(treeFilename, rownames(feat_mat2))
#     
#     nrEdges <- ncol(edgeMatrix(tree))
#     list(c(instance, list(allTreeExposures(tree, feat_mat2, row.names(signatures.cosmic)))))
#   #}
# }


#for (i in 61:180) { instance <- l[[i]][[1]][[1]][[1]]; exp_list_sims[[instance]] <- l[[i]][[1]][[2]]}

# infer errors for all trees
error_list_sims <- list()
BIC_list_sims <- list()
for (instance in names(exp_list_sims)) {
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
                      exp_list_sims[[instance]][[as.character(k)]], 
                      row.names(signatures.cosmic))
    error_list_sims[[instance]][[as.character(k)]] <- error
    BIC_list_sims[[instance]][[as.character(k)]] <- getBIC(feat_mat2, 
                                                           exp_list_sims[[instance]][[as.character(k)]], 
                                                           row.names(signatures.cosmic))
  }
}

sink("sims.tsv")
cat("instance", "\t", "n", "\t", "k", "\t", "min_mut", "\t", "max_mut", "\t", "RSS", "\t", "BIC", "\n", sep="")
for (instance in names(exp_list_sims)) {
  s <- unlist(strsplit(instance, "_"))
  treeFilename <- paste0(sim_path, s[[1]], "_", s[[2]], ".tree")
  
  feat_mat2 <- parseCountMatrix(paste0(sim_path, instance))
  tree <- parseTree(treeFilename, rownames(feat_mat2))

  nrEdges  <- length(nodes(tree)) - 1
  for (k in 0:nrEdges) {
    exp_mat <- exp_list[[instance]][[as.character(k)]]
    mut_list <- list()
    for (col in names(exp_mat)) {
      s <- unlist(strsplit(col, ";"))
      count <- 0
      for (node in s) {
        count <- count + sum(feat_mat2[node,])
      }
      mut_list[[col]] <- count
    }
    
    cat(instance, "\t", 
        nrEdges + 1, "\t",
        k, "\t",
        min(as.numeric(mut_list)), "\t",
        max(as.numeric(mut_list)), "\t",
        error_list_sims[[instance]][[as.character(k)]], "\t",
        BIC_list_sims[[instance]][[as.character(k)]], "\n", sep="")
  }
}
sink()

# Best BIC
sink("sims_opt.tsv")
cat("instance", "\t", "n", "\t", "subtree", "\t", "mutation_count", "\t", "RSS", sep="")
for (sig in row.names(signatures.cosmic)) {
  cat("\t", sig, sep="")
}
cat("\n", sep="")
for (instance in names(exp_list_sims)) {
  s <- unlist(strsplit(instance, "_"))
  treeFilename <- paste0(sim_path, s[[1]], "_", s[[2]], ".tree")
  
  feat_mat <- parseCountMatrix(paste0(sim_path, instance))
  tree <- parseTree(treeFilename, rownames(feat_mat))

  nrEdges  <- length(nodes(tree)) - 1
    
  min_bic <- Inf
  k_star <- -1
  for (k in 0:nrEdges) {
    bic <- BIC_list_sims[[instance]][as.character(k)][[1]]
    if (bic < min_bic) {
      min_bic <- bic
      k_star <- k
    }
  }
    
  exp_mat <- exp_list_sims[[instance]][[as.character(k_star)]]
    
  mut_list <- list()
  for (col in names(exp_mat)) {
    s <- unlist(strsplit(col, ";"))
    count <- 0
    for (node in s) {
      count <- count + sum(feat_mat[node,])
    }
    mut_list[[col]] <- count
    
    cat(instance, "\t", 
        nrEdges + 1, "\t",
        col, "\t",
        mut_list[[col]], "\t",
        getError(feat_mat[s,], exp_mat[col], row.names(signatures.cosmic)), sep="")
    
    for (sig in row.names(signatures.cosmic)) {
      cat("\t", exp_mat[sig,col])
    }
    cat("\n")
  }
}
sink()
