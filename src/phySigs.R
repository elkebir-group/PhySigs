library(deconstructSigs)
library(graph)
library(Rgraphviz)
library(RColorBrewer)

#sigs.input['CRUK0001:4',]+sigs.input['CRUK0001:5',]

treeExposures <- function(tree, feat_mat, k, sigs_filter) {
  V <- nodes(tree)
  E <- edgeMatrix(tree)
  E_list <- list()
  nrEdges <- dim(E)[[2]]
  for (i in 1:nrEdges) {
    E_list[[i]] <- cbind(E[1,i], E[2,i])
  }
  
  C <- combn(E_list, k)
  
  best_CC <- NULL
  best_error <- Inf
  best_sample_exp <- NULL
  for (idx in 1:dim(C)[2]) {
    cpy_tree <- tree
    if (k > 0) {
      for (i in 1:k) {
        edge <- C[i,idx][[1]]
        cpy_tree <- removeEdge(V[edge[1]], V[edge[2]], cpy_tree)
      }
    }
    
    sample_exp <- data.frame(matrix(0L, nrow = length(sigs_filter), ncol = 0))
    row.names(sample_exp)<- sigs_filter
    
    CC <- connComp(cpy_tree)
    error <- 0
    for (CCC in CC) {
      feat_mat_CCC <- feat_mat[CCC[1],]
      if (length(CCC) >= 2) {
        for (i in 2:length(CCC)) {
          feat_mat_CCC <- feat_mat_CCC + feat_mat[CCC[i],]
        }  
      }
      row.names(feat_mat_CCC) <- paste(CCC,collapse=";")
      # print(dim(feat_mat_CCC))
      
      # Get exposure for sample
      sample_exp_CCC = whichSignatures(tumor.ref = feat_mat_CCC, 
                                       signatures.ref = signatures.cosmic, 
                                       associated = sigs_filter,
                                       contexts.needed = TRUE,
                                       signature.cutoff = 0.0001,
                                       tri.counts.method = 'default')
      
      # Add any unknown signatures
      sample_exp_CCC$weights$Signature.unknown <- sample_exp_CCC[["unknown"]]
      active <- sample_exp_CCC$weights[sigs_filter]
      sample_exp <- cbind(sample_exp, t(active))
    }
    
    error <- getError(feat_mat, sample_exp, sigs_filter)
    
    if (error < best_error) {
      best_error <- error
      best_CC <- CC
      best_sample_exp <- sample_exp
    }
  }
  
  print(paste("k:", k, "; error:", best_error))
  return(best_sample_exp)
}

getError <- function(feat_mat, exp_mat, sigs_filter) {
  expanded_exp_mat <- expand(exp_mat)
  expanded_exp_mat <- expanded_exp_mat[, row.names(feat_mat)]
  for (node in row.names(feat_mat)) {
    expanded_exp_mat[node] <- expanded_exp_mat[node] * sum(feat_mat[node,])
  }
  
  diff_mat <- feat_mat - (t(as.matrix(expanded_exp_mat)) %*% as.matrix(signatures.cosmic[sigs_filter,]))
  
  error <- sum(abs(diff_mat)^2)
  
  return(error)
}

getBIC <- function(feat_mat, exp_mat, sigs_filter) {
  error <- getError(feat_mat, exp_mat, sigs_filter)
  n <- prod(dim(feat_mat))
  k <- prod(dim(exp_mat))
  
  return(n * log(error/n) + k * log(n))
}

expand <- function(exp_mat) {
  res_exp_mat <- data.frame(matrix(0L, nrow = length(row.names(exp_mat)), ncol = 0))

  for (col in names(exp_mat)) {
    s <- unlist(strsplit(col, ";"))
    for (node in s) {
      res_exp_mat[as.character(node)] <- exp_mat[col]
    }
  }
  
  return(res_exp_mat)
}

allTreeExposures <- function(tree, feat_mat, sigs_filter) {
  exp_list <- list()
  nrEdges  <- length(nodes(tree)) - 1
  for (k in 0:nrEdges) { exp_list[[as.character(k)]] <- treeExposures(tree, feat_mat, k, sigs_filter) }
  return(exp_list)
}

plotTree <- function(patient, title, tree, feat_mat, exp_mat, tree_idx=0) {
  expanded_exp_mat <- expand(exp_mat)
  expanded_exp_mat <- expanded_exp_mat[, nodes(tree)]
  
  # Sub graphs -- disabled for now
  # subTList <- vector(mode="list", length=length(names(exp_mat)))
  # i <- 1
  # for (C in names(exp_mat)) {
  #   print(nodes(tree))
  #   print(strsplit(C, "[;]")[[1]])
  #   subTList[[i]] <- list(graph=subGraph(strsplit(C, "[;]")[[1]], tree))
  #   i <- i + 1
  # }
  
  eAttrs <- list()
  eAttrs$color <- list()
  V <- nodes(tree)
  E <- edgeMatrix(tree)
  nrEdges <- dim(E)[[2]]
  for (i in 1:nrEdges) {
    source <- V[[E[1,i]]]
    target <- V[[E[2,i]]]
    
    ok <- FALSE
    for (C in names(exp_mat)) {
      CC <- strsplit(C, "[;]")[[1]]

      if ((source %in% CC) && (target %in% CC)) {
        ok <- TRUE
      }
    }
    if (ok) {
      eAttrs$color[[paste(source, "~",  target, sep="")]] <- "black"
    }
    else {
      eAttrs$color[[paste(source, "~",  target, sep="")]] <- "red"
    }
  }
  
  # Set color palette
  fill <- brewer.pal(8, "Set1")
  names(fill) <- row.names(exp_mat)
  
  # Format pie chart nodes
  g1layout <- agopen(tree, name="foo")
  makeNodeDrawFunction <- function(x, fill, patient, feat_mat) {
    force(x)
    function(node, ur, attrs, radConv) {
      
      # Remove labels
      names(x) <- vector(mode="character", length=length(names(x)))
      
      # Get node locations
      nc <- getNodeCenter(node)
      
      # Get consistent color scheme
      pal <- fill[which(x > 0)]
      
      # Make plots
      pieGlyph(x[which(x > 0)], xpos=getX(nc), ypos=getY(nc), radius=getNodeRW(node), col=pal)
      text(getX(nc), getY(nc), paste(name(node), ": ", sum(feat_mat[name(node), ]), 
                                     #"\\n",
                                     #frob_mat_p[[paste0(patient, ":", name(node))]], 
                                     sep=""), cex=.3, col="black", font=2)
    }
  }
  drawFunc <- apply(expanded_exp_mat, 2, makeNodeDrawFunction, fill=fill, patient=patient, feat_mat=feat_mat) #mut_counts=mut_counts, frob_mat_p=frob_mat_p)
  
  # Make plot
  title_id <- paste(patient, sep="")
  if (tree_idx > 0) {
    title_id <- paste(patient, letters[tree_idx], sep="")
  }
  
  # Truncate signature names
  sig_ids <- strsplit(as.character(row.names(exp_mat)), "[.]")
  sig_ids <- as.character(sapply(sig_ids, "[", 2 ))
  
  plot(tree, #subGList=subTList,
       drawNode=drawFunc, edgeAttrs = eAttrs, 
       attrs=list(node=list(height=.75, width=.75, fontsize=2), 
                  edge=list(fontsize=5)), mai=c(0.15, 0.15, 0.15, 1), 
       main=paste("ID:", title_id, "--", title, sep=' '))
  legend("topright", inset=c(-0.09,.03), title="Signatures", sig_ids, fill=fill, xpd = TRUE)
}
