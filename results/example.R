library(deconstructSigs)
library(graph)
setwd("/Users/sarahchristensen/Documents/UIUC/Bioinformatics/Signatures/phySigs/results/")
source("../src/phySigs.R")

# Set signatures of interest
S <- c("Signature.1", "Signature.3")

# Input CSV file 
tree_matrix <- read.csv(file="../data/example/tree.csv", header=TRUE, sep=",")

# V1/V2 contains clone label for first/second endpoint of edge
V1 <- as.character(tree_matrix$V1)
V2 <- as.character(tree_matrix$V2)

# Get nodes in the tree
node_ids <- unique(union(V1, V2))

# Initialize a tree 
T <- new("graphNEL", nodes=node_ids, edgemode="directed")

# Add edges to the tree
for (i in 1:nrow(tree_matrix)){
  T <- addEdge(V1[i], V2[i], T, 1)
}

# Input CSV file
input_mat <- as.data.frame(read.csv(file="../data/example/snv.csv", 
                           colClasses = c("character", "character", "numeric", "character", "character")))

# Use deconstructSigs to convert SNVs to 96 Features
P <- mut.to.sigs.input(mut.ref = input_mat, 
                    sample.id = "Sample", 
                    chr = "chr", 
                    pos = "pos", 
                    ref = "ref", 
                    alt = "alt")

# Normalize feature matrix
P_Norm <- normalizeFeatureMatrix(P, "genome")

# Recover lowest error exposure matrix for every possible number of clusters. 
E_list <- allTreeExposures(T, P_Norm, S)

# Get error for exposure matrix from best set of k clusters.
k <- 1
error <- getError(P_Norm, E_list[[k]], S)

# Get BIC for exposure matrix from best set of k clusters. 
bic <- getBIC(P_Norm, E_list[[k]], S)

# Get tree figure with pie chart nodes showing exposures.
pdf("example_results.pdf", width=10, height=10)
tumorID <- "Tumor1"
title <- "Example Plot for PhySigs"
plotTree(tumorID, title, T, P_Norm, E_list[[k]], tree_idx=1)
dev.off()


