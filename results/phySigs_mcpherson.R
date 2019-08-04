setwd("/Users/sarahchristensen/Documents/UIUC/Bioinformatics/Signatures/phySigs/results/")
source("../src/phySigs.R")
#load("mcpherson.RData")

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

getPatients <- function(samples) {
  patients <- strsplit(as.character(samples), "[:]")
  patients <- sapply(patients, "[", 1 )
  patients <- unique(patients)
  
  return(patients)
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

getTreeMatrix <- function(patients) {
  
  # Import trees
  trees <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(trees) <- c("SampleID", "V1", "V2")
  
  for (i in 1:length(patients)){
    p_tree <- read.table(paste("../data/mcpherson/trees/patient", patients[i], ".tree", sep=""), 
                         header = FALSE, sep = " ")
    p_tree$V2 <- sapply(strsplit(as.character(p_tree$V2), ""), "[", 1 )
    p_tree$SampleID <- rep(patients[i],nrow(p_tree))
    p_tree <- subset(p_tree, V1!=V2)
    trees <- rbind(trees, p_tree)
  }
  return(trees)
}

getSNVMatrix <- function(file_name) {
  # Get raw input data on SNVs
  input_mat = as.data.frame(read.csv(file_name, colClasses = c("character", "character", "numeric", "character", "character")))
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

getOVSigs <- function() {
  # Create vector of signatures for OV cancer
  # See https://cancer.sanger.ac.uk/signatures_v2/matrix.png
  sigs_ovary <- c("Signature.1", "Signature.3", "Signature.5")
  #sigs_additional <- c("Signature.2", "Signature.13")
  sigs_additional <- c("Signature.2", "Signature.5", "Signature.6", "Signature.8", "Signature.10", "Signature.13", "Signature.17", "Signature.18", "Signature.20", "Signature.26", "Signature.30")
  sigs_all <- union(sigs_ovary, sigs_additional)
  return(sigs_all)
}

getTrees <- function(patient, tree_mat) {
  
  pt <- subset(tree_mat, SampleID==patient) # patient tree entries
  node_ids <- as.character(unique(union(pt$V1, pt$V2)))

  # Initialize a tree 
  tree <- new("graphNEL", nodes=node_ids, edgemode="directed")
    
  # Add edges to the tree
  for (i in 1:nrow(pt)){
    tree <- addEdge(as.character(pt$V1[i]), as.character(pt$V2[i]), tree, 1)
  }

  return(tree)
}

###############################################################################
###############################################################################
sigs_all  <- getOVSigs()
input_mat <- getSNVMatrix("../data/mcpherson/snv.csv")
patients <- getPatients(unique(input_mat$Sample))
tree_mat  <- getTreeMatrix(patients)

# infer exposures for all trees
sigs.input <- getFeatures(input_mat)
exp_list <- list()
for (patient in patients) {
  exp_list[[patient]] <- list()
  t <- getTrees(patient, tree_mat)
  feat_mat <- getPatientFeatures(sigs.input, patient, t)
  nrEdges  <- length(nodes(t)) - 1
  exp_list[[patient]][[1]] <- allTreeExposures(t, feat_mat, sigs_all, "default")
}

# infer errors for all trees
error_list <- list()
BIC_list <- list()
for (patient in patients) {
  i <- 1
  error_list[[patient]] <- list()
  BIC_list[[patient]] <- list()
  t <- getTrees(patient, tree_mat)
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
}

# plot results
pdf("mcpherson_results.pdf", width=10, height=10)
for (patient in patients) {
  t <- getTrees(patient, tree_mat)
  # pd <- subset(ind_mat, TRACERxID==patient) # patient data

  i <- 1
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
      
    title <- paste("k =", k, ";", round(error, digits=4), sep=' ')
      
    plotTree(patient, title, t, feat_mat, exp_mat, tree_idx=tree_idx)
  }
}
dev.off()

# mcpherson Best BIC
sink("mcpherson_opt.tsv")
cat("patient", "\t", "n", "\t", "tree", "\t", "subtree", "\t", "mutation_count", "\t", "RSS", sep="")
for (sig in sigs_all) {
  cat("\t", sig, sep="")
}
cat("\n", sep="")
for (patient in patients) {
  t <- getTrees(patient, tree_mat)
  i <- 1

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
}
sink()

# mcpherson summary
sink("mcpherson.tsv")
cat("patient", "\t", "n", "\t", "tree", "\t", "k", "\t", "min_mut", "\t", "max_mut", "\t", "RSS", "\t", "BIC", "\n", sep="")
for (patient in patients) {
  t <- getTrees(patient, tree_mat)
  i <- 1
    
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
}
sink()

save.image("mcpherson.RData")
