# Code used to simulate the PPM model developped by Gamblin, Lambert and Blanquart.
#
# Script organization:
# 1. dependencies
# 2. functions used for tree simulation
# 3. functions used for genes simulation
# 4. post-treatment
# 5. example of use
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

# simulate I1 genes while recording their branch and time of arrival
sim.I1_time <- function(i1, l1, e, err="cst") {
  I1_im <- lapply(tree$edge.length, function(l) compute.im1(l, i1))
  I1_edges <- unlist(sapply(seq(length(I1_im)), function(i) rep(i, length(I1_im[[i]]))))
  I1_times <- unlist(I1_im)
  rep <- cbind(t(matrix(unlist(mclapply(seq(length(I1_edges)), function(i)
    rep.I1(I1_edges[i], I1_times[i], l1), mc.cores=1)), ncol=length(I1_edges))),
    I1_edges, tree$edge.length[I1_edges]-I1_times)
  if (err=="cst") {
    rep[,1:G] <- sapply(seq(ncol(rep[,1:G])), function(j) sapply(seq(nrow(rep)), function(i)
      if (rep[i,j] == 1)rbinom(1,1,1-e) else 0))
  } else if (err=="exp") {
    par_ber <- rexp(nrow(rep), rate = 1/e)
    rep[,1:G] <- sapply(seq(ncol(rep[,1:G])), function(j) sapply(seq(nrow(rep)), function(i)
      if (rep[i,j] == 1) rbinom(1,1,1-par_ber[i]) else 0))
  }
  selec <- which(rowSums(rep[,1:G]) == 0)
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

# simulate I2 genes while counting the number of introductions/gene
# genes present in the pool but not gained yet in a lineage are encoded as -2
compute.gain_loss_intro <- function(presence, length, g2, l2) {
  intro <- 0
  t <- 0
  p <- presence
  dt <- rexp(1, rate = if (p==1) l2 else g2)
  while (t+dt < length) {
    t <- t + dt
    if (p==-2) {
      p <- 1
      intro <- 1
    }
    else {p <- 1-p}
    dt <- rexp(1, rate = if (p==1) l2 else g2)
  }
  return(c(p, intro))
}
rep.I2_intro <- function(time, g2, l2) {
  intro <- 0
  nodes <- rep(-1,2*G-1)
  first_nodes <- cut.tree(time)
  first_gain_loss <- sapply(first_nodes, function(n)
    compute.gain_loss_intro(-2, depths[n]-time, g2, l2)) # gene is present or absent at each first node
  nodes[first_nodes] <- first_gain_loss[1,]
  intro <- sum(first_gain_loss[2,])
  for (i in seq(length(tree$edge.length))) {
    mother <- tree$edge[i,1]
    if (nodes[mother] != -1) {
      child <- tree$edge[i,2]
      branch.length <- tree$edge.length[i]
      presence <- compute.gain_loss_intro(nodes[mother], branch.length, g2, l2)
      nodes[child] <- presence[1]
      intro <- intro + presence[2]
    }
  }
  nodes[nodes == -1] <- 0
  return(c(nodes[1:G], intro))
}
sim.I2_intro <- function(i2, g2, l2) {
  I2_im <- compute.im2(H, i2)
  rep <- t(matrix(unlist(mclapply(I2_im, function(t) rep.I2_intro(t, g2, l2), mc.cores=1)), ncol=length(I2_im)))
  return(rep)
}

# simulate I2 genes while recording their importation time
sim.I2_time <- function(i2, g2, l2, e1, e2, err="cst") {
  I2_im <- compute.im2(H, i2)
  rep <- cbind(t(matrix(unlist(mclapply(I2_im, function(t) rep.I2(t, g2, l2), mc.cores=1)), ncol=length(I2_im))), H-I2_im)
  if (err=="cst") {
    rep[,1:G] <- sapply(seq(ncol(rep[,1:G])), function(j) sapply(seq(nrow(rep)), function(i)
      if (rep[i,j] == 0) rbinom(1,1,e1) else rbinom(1,1,1-e2)))
  } else if (err=="exp") {
    par_ber1 <- rexp(nrow(rep), rate = 1/e1)
    par_ber2 <- rexp(nrow(rep), rate = 1/e2)
    rep[,1:G] <- sapply(seq(ncol(rep[,1:G])), function(j) sapply(seq(nrow(rep)), function(i)
      if (rep[i,j] == 0) rbinom(1,1,par_ber1[i]) else rbinom(1,1,1-par_ber2[i])))
  }
  selec <- which(rowSums(rep[,1:G]) == 0)
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





############################################ USE EXAMPLE #######################################################

# simulate an ultrametric random tree with 20 leaves and height 1
tree <- sim.tree(20, 1)

# or load a reconstructed species tree (in Newick format)
# tree <- read.tree(file = "tree.nwk")


# compute some tree statistics
depths <- get_all_distances_to_root(tree)
H <- max(depths)
G <- tree$Nnode + 1
root <- G + 1


# set parameters
N0 <- 200
l0 <- 0.1
i1 <- 50
l1 <- 2
i2 <- 150.
g2 <- 5.
l2 <- 10.
eps0 <- 0.015
eps1 <- 0.015
eps2 <- 0.015

# or draw them randomly
# N0 <- sample(500:1500, 1)
# l0 <- rlunif(1, 0.1/L, 10/L)
# i1 <- runif(1, 500/L, 1500/L)
# l1 <- rlunif(1, 5/L, 500/L)
# i2 <- runif(1, 500/H, 1500/H)
# g2 <- rlunif(1, 5/L, 500/L)
# l2 <- rlunif(1, 250/L, 25000/L)
# eps0 <- rexp(1, rate = 100)
# eps1 <- rexp(1, rate = 100)
# eps2 <- rexp(1, rate = 100)


# simulate gene patterns along the tree
rep_N0 <- sim.N(N0, l0, eps0, err="cst")
rep_N1 <- sim.N(round(i1/l1), l1, eps1, err="cst")
rep_I1 <- sim.I1(i1, l1, eps1, err="cst")
rep_I2 <- sim.I2(i2, g2, l2, eps2, eps2, err="cst")


# build presence/absence matrix
rep <- rbind(rep_N0, rep_N1, rep_I1, rep_I2)
rep <- cbind(rep, c(rep(0, nrow(rep_N0)), rep(1, nrow(rep_N1)+nrow(rep_I1)), rep(2, nrow(rep_I2))))
colnames(rep) <- c(tree$tip.label, "cat")
print(table(rep[,"cat"])) # print pangenome composition


# store presence/absence matrix in csv format, compatible with C++ inference tool
# write.table(t(rep[,1:G]), file = "pa_matrix.txt", sep = ",", quote = F)






