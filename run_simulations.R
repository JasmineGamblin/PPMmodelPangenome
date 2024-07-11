# Script allowing to simulate the PPM model developped by Gamblin, Lambert and Blanquart.
#
# Rely on functions implemented in simulation_aux.R
# 
# This code is distributed under the GNU GPL license.
# 
# Author: Jasmine Gamblin





############################################ USE EXAMPLE #######################################################

library(KScorrect) # for rlunif (log-uniform distribution)

# simulate an ultrametric random tree with 20 leaves and height 1
tree <- sim.tree(20, 1)
write.tree(tree, "tree.nwk")

# or load a reconstructed species tree (in Newick format)
# tree <- read.tree(file = "tree.nwk")


# compute some tree statistics
depths <- get_all_distances_to_root(tree)
H <- max(depths)
G <- tree$Nnode + 1
root <- G + 1
L <- sum(tree$edge.length)


# set parameters
# N0 <- 200
# l0 <- 0.1
# i1 <- 50
# l1 <- 2
# i2 <- 150.
# g2 <- 5.
# l2 <- 10.
# eps0 <- 0.015
# eps1 <- 0.015
# eps2 <- 0.015

# or draw them randomly
set.seed(1)
N0 <- sample(200:800, 1)
l0 <- rlunif(1, 0.1/L, 10/L)
i1 <- runif(1, 500/L, 1500/L)
l1 <- rlunif(1, 5/L, 50/L)
i2 <- runif(1, 200/H, 800/H)
g2 <- rlunif(1, 5/L, 500/L)
l2 <- rlunif(1, 250/L, 25000/L)
eps0 <- rexp(1, rate = 100)
eps1 <- rexp(1, rate = 100)
eps2 <- rexp(1, rate = 100)


# simulate gene patterns along the tree
source(file = "simulation_aux.R")
set.seed(1)
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
write.table(t(rep[,1:G]), file = "pa_matrix.txt",
            sep = ",", quote = F)


