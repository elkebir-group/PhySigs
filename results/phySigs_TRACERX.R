source("../src/phySigs.R")
load("tracerx.RData")

library(deconstructSigs)
library(readxl)
library(graph)

getFeatures <- function(input_mat) {
  # Use deconstructSigs to convert SNVs to 96 Features
  mut.to.sigs.input(mut.ref = input_mat, 
                    sample.id = "Sample", 
                    chr = "chr", 
                    pos = "pos", 
                    ref = "ref", 
                    alt = "alt")
}

getPatientFeatures <- function(sigs.input, patient, tree) {
  pf <- sigs.input[paste0(patient, ":", nodes(tree)),]
  rownames(pf) <- nodes(tree)
  return(pf)  
}

getPhenotypeMatrix <- function(file_name) {
  # Get information on individuals
  return(read_excel(file_name, sheet = "TableS2", skip=1))
}

getTreeMatrix <- function(file_name) {
  # Get information on trees
  return(read_excel(file_name, sheet = "TableS7", skip=1))
}

getSNVMatrix <- function(file_name) {
  # Get raw input data on SNVs
  input_mat = read_excel("appendix_2.xlsx", sheet = "cleaned_variants", col_types = c("text", "text", "text", "text", "numeric", "text", "text", "numeric"))
  input_mat <- as.data.frame(subset(input_mat, include ==1))
  input_mat$include <- NULL
  input_mat$chr <- paste("chr", input_mat$chr, sep="") 
  return(input_mat)
}

getLungSigs <- function() {
  # Create vector of signatures for each type of lung cancer
  # See https://cancer.sanger.ac.uk/signatures_v2/matrix.png
  sigs_adeno <- c("Signature.1", "Signature.2", "Signature.4", "Signature.5", "Signature.6", "Signature.13", "Signature.17")
  sigs_squamous <- c("Signature.1", "Signature.2", "Signature.4", "Signature.5", "Signature.13")#, "Signature.unknown")
  sigs_all <- union(sigs_adeno, sigs_squamous)
  return(sigs_all)
}

getTrees <- function(patient, tree_mat) {
  pt <- subset(tree_mat, SampleID==patient) # patient tree
  
  # We make a tree for each option in the Tracerx data
  trees <- append(vector(), as.character(pt['PrimaryTreeStructure']))
  if (!is.na(pt['AdditionalTreeStructures'])){
    extra_trees <- unlist(strsplit(as.character(pt['AdditionalTreeStructures']), ","))
    for (item in extra_trees){
      trees <- append(trees, item)
    }
  }
  node_ids <- strsplit(pt$TreeClusters, ",")
  
  i <- 1
  res_trees <- list()
  for (adj in as.list(trees)) {
    # Initialize a tree 
    tree <- new("graphNEL", nodes=as.character(node_ids[[1]]), edgemode="directed")
    
    # Add edges to the tree
    edges <- unlist(strsplit(adj, ";"))
    eAttrs <- list()
    el <- list()
    for (e in edges){
      endpoints <- unlist(strsplit(e, "->"))
      tree <- addEdge(endpoints[1], endpoints[2], tree, 1)
    }
    
    res_trees[[i]] <- tree
    i <- i + 1
  }
  return(res_trees)
}

sigs_all  <- getLungSigs()
ind_mat   <- getPhenotypeMatrix("appendix_2.xlsx")
tree_mat  <- getTreeMatrix("appendix_2.xlsx")
input_mat <- getSNVMatrix("appendix_2.xlsx")

# infer exposures for all trees
sigs.input <- getFeatures(input_mat)
exp_list <- list()
for (patient in tree_mat$SampleID) {
  exp_list[[patient]] <- list()
  trees <- getTrees(patient, tree_mat)
  i <- 1
  for (t in trees) {
    feat_mat <- getPatientFeatures(sigs.input, patient, t)
    nrEdges  <- length(nodes(t)) - 1
    exp_list[[patient]][[i]] <- allTreeExposures(t, feat_mat, sigs_all)
    i <- i + 1
  }
}

# infer errors for all trees
error_list <- list()
BIC_list <- list()
for (patient in tree_mat$SampleID) {
  error_list[[patient]] <- list()
  BIC_list[[patient]] <- list()
  trees <- getTrees(patient, tree_mat)
  i <- 1
  for (t in trees) {
    feat_mat <- getPatientFeatures(sigs.input, patient, t)
    nrEdges  <- length(nodes(t)) - 1
    error_list[[patient]][[i]] <- list()
    BIC_list[[patient]][[i]] <- list()
    for (k in 0:nrEdges) {
      print(paste(patient, nrEdges, i, k))
      error <- getError(feat_mat, exp_list[[patient]][[i]][[as.character(k)]], sigs_all)
      error_list[[patient]][[i]][[as.character(k)]] <- error
      BIC_list[[patient]][[i]][[as.character(k)]] <- getBIC(feat_mat, exp_list[[patient]][[i]][[as.character(k)]], sigs_all)
    }
    i <- i + 1
  }
}

# plot results
pdf("tracerx_results.pdf", width=10, height=10)
for (patient in tree_mat$SampleID) {
  trees <- getTrees(patient, tree_mat)
  pd <- subset(ind_mat, TRACERxID==patient) # patient data
  i <- 1
  for (t in trees) {
    feat_mat <- getPatientFeatures(sigs.input, patient, t)
    nrEdges  <- length(nodes(t)) - 1
    x <- 0:nrEdges
    y <- as.numeric(error_list[[patient]][[i]])
    y2 <- as.numeric(BIC_list[[patient]][[i]])
    plot(x, y, type="b", xlab = "#removed edges", ylab = "") # first plot
    par(new = TRUE)
    plot(x, y2, type = "b", axes = FALSE, bty = "n", xlab = "", ylab = "", col="red")
    axis(side=1, at = pretty(range(x)))
    axis(side=2, at = pretty(range(y)))
    axis(side=4, at = pretty(range(y2)), col="red")
    mtext("RSS", side=2, line=3)
    mtext("BIC", side=4, line=3, col="red")
    
    for (k in 0:nrEdges) {
      print(paste(patient, nrEdges, i, k))
      error <- error_list[[patient]][[i]][[as.character(k)]]
      exp_mat <- exp_list[[patient]][[i]][[as.character(k)]]

      tree_idx <- 0
      if (length(trees) > 1) {
        tree_idx <- i
      }
      
      title <- paste("k =", k, ";", round(error, digits=4), "--", "Age:", pd['Age'],"--","Status:", 
                     pd['Smoking status'], paste("(",pd['Pack years'],")",sep="" ), 
                     sep=' ')
      
      plotTree(patient, title, t, feat_mat, exp_mat, tree_idx=tree_idx)
    }
    i <- i + 1
  }
}
dev.off()

# tracerx Best BIC
sink("tracerx_opt.tsv")
cat("patient", "\t", "n", "\t", "tree", "\t", "subtree", "\t", "mutation_count", "\t", "RSS", sep="")
for (sig in sigs_all) {
  cat("\t", sig, sep="")
}
cat("\n", sep="")
for (patient in tree_mat$SampleID) {
  trees <- getTrees(patient, tree_mat)
  i <- 1
  
  for (t in trees) {
    feat_mat <- getPatientFeatures(sigs.input, patient, t)
    nrEdges  <- length(nodes(t)) - 1
    
    min_bic <- Inf
    k_star <- -1
    for (k in 0:nrEdges) {
      bic <- BIC_list[[patient]][[i]][as.character(k)][[1]]
      if (bic < min_bic) {
        min_bic <- bic
        k_star <- k
      }
    }
    
    exp_mat <- exp_list[[patient]][[i]][[as.character(k_star)]]
    
    mut_list <- list()
    for (col in names(exp_mat)) {
      s <- unlist(strsplit(col, ";"))
      count <- 0
      for (node in s) {
        count <- count + sum(feat_mat[node,])
      }
      mut_list[[col]] <- count
      
      cat(patient, "\t", 
          nrEdges + 1, "\t",
          i, "\t", 
          col, "\t",
          mut_list[[col]], "\t",
          getError(feat_mat[s,], exp_mat[col], sigs_all), sep="")
      
      for (sig in sigs_all) {
        cat("\t", exp_mat[sig,col])
      }
      cat("\n")
    }

    i <- i + 1
  }
}
sink()

# tracerx summary
sink("tracerx.tsv")
cat("patient", "\t", "n", "\t", "tree", "\t", "k", "\t", "min_mut", "\t", "max_mut", "\t", "RSS", "\t", "BIC", "\n", sep="")
for (patient in tree_mat$SampleID) {
  trees <- getTrees(patient, tree_mat)
  i <- 1
  for (t in trees) {
    feat_mat <- getPatientFeatures(sigs.input, patient, t)
    nrEdges  <- length(nodes(t)) - 1
    for (k in 0:nrEdges) {
      exp_mat <- exp_list[[patient]][[i]][[as.character(k)]]
      mut_list <- list()
      for (col in names(exp_mat)) {
        s <- unlist(strsplit(col, ";"))
        count <- 0
        for (node in s) {
          count <- count + sum(feat_mat[node,])
        }
        mut_list[[col]] <- count
      }
      
      cat(patient, "\t", 
          nrEdges + 1, "\t",
          i, "\t", 
          k, "\t",
          min(as.numeric(mut_list)), "\t",
          max(as.numeric(mut_list)), "\t",
          error_list[[patient]][[i]][[as.character(k)]], "\t",
          BIC_list[[patient]][[i]][[as.character(k)]], "\n", sep="")
    }
    i <- i + 1
  }
}
sink()
