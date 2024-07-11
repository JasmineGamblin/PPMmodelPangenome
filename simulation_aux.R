# Code used to simulate the PPM model developped by Gamblin, Lambert and Blanquart.
#
# Code organization:
# 1. dependencies
# 2. functions used for tree simulation
# 3. functions used for genes simulation
# 4. post-treatment
# 
# This code is distributed under the GNU GPL license.
# 
# Author: Jasmine Gamblin





############################################ DEPENDENCIES ######################################################

library(ape)
library(castor)
library(data.table) # for matrix compression
library(parallel) # for mclapply





############################################ TREE SIMULATION ###################################################

# store in Newick format
newick <- function(tree) {
  if (length(tree$subtrees) == 1) {
    return(paste(tree$subtrees, ":", toString(tree$branch.length), sep = ""))
  }
  else {
    return(paste("(", newick(tree$subtrees[[1]]), ",", newick(tree$subtrees[[2]]),"):", toString(tree$branch.length), sep = ""))
  }
}

# make tree recursion
make.tree <- function(t, times, leaves, d, b) {
  n <- length(leaves)
  if (n == 1) {
    return(list(branch.length = d-t, subtrees = toString(leaves[1])))
  }
  else {
    coal <- times[1]
    if (n == 2) {
      ind.left <- sample(1:2, 1)
      times.left <- c()
    }
    else {
      p <- rbeta(1, b+1, b+1)
      n.left <- sample(1:(n-1), 1, prob = sapply(1:(n-1), function(k) choose(n, k)*p**k*(1-p)**(n-k)/(1-p**n-(1-p)**n)))
      ind.left <- sample(1:n, n.left, replace = F)
      times.left <- sample(2:(n-1), n.left-1, replace = F)
      if (n==3 && n.left==2) {times.left <- 2}
    }
    return(list(branch.length = coal-t, subtrees = list(make.tree(coal, sort(times[times.left]), leaves[ind.left], d, b),
                                                        make.tree(coal, sort(times[-c(1,times.left)]), leaves[-ind.left], d, b))))
  }
}

# simulate beta-splitting ultrametric tree, with uniformly distributed coalescence times
sim.tree <- function(nb.leaves, height, b=-0.5) {
  coal_times <- c(0, sort(runif(nb.leaves-2, 0, height)))
  tree <- make.tree(0, coal_times, 1:nb.leaves, height, b)
  tree <- read.tree(text = paste(newick(tree),";", sep=""))
  return(tree)
}





############################################ GENES SIMULATION ##################################################

# simulate N0 genes
sim.N0 <- function(N0, e, err = "cst") {
  rep <- matrix(rep(1,N0*G), nrow = N0)
  if (err=="cst") {
    rep <- sapply(seq(ncol(rep)), function(j) sapply(seq(nrow(rep)), function(i) rbinom(1,1,1-e)))
  } else if (err=="exp") {
    par_ber <- rexp(nrow(rep), rate = 1/e)
    rep <- sapply(seq(ncol(rep)), function(j) sapply(seq(nrow(rep)), function(i) rbinom(1,1,1-par_ber[i])))
  }
  selec <- which(rowSums(rep) == 0)
  if (length(selec) > 0) {
    return(rep[-selec,])
  } else {
    return(rep)
  }
}

# simulate N1 genes
compute.loss <- function(length, l1) {
  loss <- rexp(1, rate = l1)
  return(if (loss < length) 0 else 1)
}
rep.N1 <- function(l1) {
  nodes <- c()
  nodes[root] <- 1 # gene is present at root
  for (i in seq(length(tree$edge.length))) {
    mother <- tree$edge[i,1]
    child <- tree$edge[i,2]
    if (nodes[mother]) {
      branch.length <- tree$edge.length[i]
      presence <- compute.loss(branch.length, l1)
    }
    else {
      presence <- 0
    }
    nodes[child] <- presence
  }
  return(nodes[1:G])
}
sim.N <- function(N, l, e, err="cst") {
  rep <- t(sapply(seq(N), function(i) rep.N1(l)))
  if (err=="cst") {
    rep <- sapply(seq(ncol(rep)), function(j) sapply(seq(nrow(rep)), function(i) if (rep[i,j] == 1) rbinom(1,1,1-e) else 0))
  } else if (err=="exp") {
    par_ber <- rexp(nrow(rep), rate = 1/e)
    rep <- sapply(seq(ncol(rep)), function(j) sapply(seq(nrow(rep)), function(i) if (rep[i,j] == 1) rbinom(1,1,1-par_ber[i]) else 0))
  }
  if (is.null(dim(rep))) {dim(rep) <- c(1, G)}
  selec <- which(rowSums(rep) == 0)
  if (length(selec) > 0) {
    rep <- rep[-selec,]
    if (is.null(dim(rep))) {dim(rep) <- c(1, G)}
  }
  return(rep)
}

# choose I1 immigration events
compute.im1 <- function(length, i1) {
  t <- 0
  im <- c()
  dt <- rexp(1, rate = i1)
  while (t+dt < length) {
    t <- t + dt
    im <- c(im, t)
    dt <- rexp(1, rate = i1)
  }
  return(im)
}

# simulate I1 genes
rep.I1 <- function(edge, time, l1) {
  nodes <- rep(0, 2*G-1)
  first_node <- tree$edge[edge,2]
  presence_first_node <- compute.loss(tree$edge.length[edge]-time, l1)
  if (presence_first_node) {
    nodes[first_node] <- 1 # gene is present at first node
    for (i in seq(length(tree$edge.length))) {
      mother <- tree$edge[i,1]
      if (nodes[mother]) {
        child <- tree$edge[i,2]
        branch.length <- tree$edge.length[i]
        presence <- compute.loss(branch.length, l1)
        nodes[child] <- presence
      }
    }
  }
  return(nodes[1:G])
}
sim.I1 <- function(i1, l1, e, err="cst") {
  I1_im <- lapply(tree$edge.length, function(l) compute.im1(l, i1))
  I1_edges <- unlist(sapply(seq(length(I1_im)), function(i) rep(i, length(I1_im[[i]]))))
  I1_times <- unlist(I1_im)
  rep <- t(matrix(unlist(mclapply(seq(length(I1_edges)), function(i)
    rep.I1(I1_edges[i], I1_times[i], l1), mc.cores=1)), ncol=length(I1_edges)))
  if (err=="cst") {
    rep <- sapply(seq(ncol(rep)), function(j) sapply(seq(nrow(rep)), function(i) if (rep[i,j] == 1) rbinom(1,1,1-e) else 0))
  } else if (err=="exp") {
    par_ber <- rexp(nrow(rep), rate = 1/e)
    rep <- sapply(seq(ncol(rep)), function(j) sapply(seq(nrow(rep)), function(i) if (rep[i,j] == 1) rbinom(1,1,1-par_ber[i]) else 0))
  }
  selec <- which(rowSums(rep) == 0)
  if (length(selec) > 0) {
    return(rep[-selec,])
  } else {
    return(rep)
  }
}

# simulate N2 genes
compute.gain_loss <- function(presence, length, g2, l2) {
  t <- 0
  p <- presence
  dt <- rexp(1, rate = if (p) l2 else g2)
  while (t+dt < length) {
    t <- t + dt
    p <- 1-p
    dt <- rexp(1, rate = if (p) l2 else g2)
  }
  return(p)
}

# choose I2 immigration events
compute.im2 <- function(length, i2) {
  t <- 0
  im <- c()
  dt <- rexp(1, rate = i2)
  while (t+dt < length) {
    t <- t + dt
    im <- c(im, t)
    dt <- rexp(1, rate = i2)
  }
  return(im)
}

# simulate I2 genes
cut.tree <- function(time) {
  cut <- split_tree_at_height(tree, height = time)
  return(sapply(cut$subtrees, function(s) {
    n <- s$new2old_clade[s$tree$Nnode + 2]
    if (is.na(n)) s$new2old_clade[1] else n
  }))
}
rep.I2 <- function(time, g2, l2) {
  nodes <- rep(-1,2*G-1)
  first_nodes <- cut.tree(time)
  nodes[first_nodes] <- sapply(first_nodes, function(n)
    compute.gain_loss(0, depths[n]-time, g2, l2)) # gene is present or absent at each first nodes
  for (i in seq(length(tree$edge.length))) {
    mother <- tree$edge[i,1]
    if (nodes[mother] != -1) {
      child <- tree$edge[i,2]
      branch.length <- tree$edge.length[i]
      presence <- compute.gain_loss(nodes[mother], branch.length, g2, l2)
      nodes[child] <- presence
    }
  }
  nodes[nodes == -1] <- 0
  return(nodes[1:G])
}
sim.I2 <- function(i2, g2, l2, e1, e2, err="cst") {
  I2_im <- compute.im2(H, i2)
  rep <- t(matrix(unlist(mclapply(I2_im, function(t) rep.I2(t, g2, l2), mc.cores=1)), ncol=length(I2_im)))
  if (err=="cst") {
    rep <- sapply(seq(ncol(rep)), function(j) sapply(seq(nrow(rep)), function(i)
      if (rep[i,j] == 0) rbinom(1,1,e1) else rbinom(1,1,1-e2)))
  } else if (err=="exp") {
    par_ber1 <- rexp(nrow(rep), rate = 1/e1)
    par_ber2 <- rexp(nrow(rep), rate = 1/e2)
    rep <- sapply(seq(ncol(rep)), function(j) sapply(seq(nrow(rep)), function(i)
      if (rep[i,j] == 0) rbinom(1,1,par_ber1[i]) else rbinom(1,1,1-par_ber2[i])))
  }
  selec <- which(rowSums(rep) == 0)
  if (length(selec) > 0) {
    return(rep[-selec,])
  } else {
    return(rep)
  }
}





############################################ POST-TREATMENT ####################################################

# compression
compress.rep <- function(rep) {
  rep <- as.data.table(rep)
  rep <- rep[, Count := .N, by = names(rep)]
  rep <- data.frame(unique(rep), check.names = F) # merge identical gene patterns
  return(as.matrix(rep))
}






