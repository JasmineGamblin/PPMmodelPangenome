# Inference tool for the PPM model

## Description
This repository contains the source code used in the manuscript 'Persistent, Private and Mobile genes: a model for gene dynamics in evolving pangenomes' by Jasmine Gamblin, Amaury Lambert and François Blanquart.

## Files
- `inference.cpp`: main file to perform inference.
- `fonctions.cpp` and `fonctions.hpp`: auxiliary functions for inference procedure.
- `objets.cpp` and `objets.hpp`: objects definitions for inference procedure.
- `nelder_mead_dfoptim.h`: Nelder-Mead optimzation algorithm, as implemented in the R package `dfoptim`
- `simulation.R`: R code used to simulate the PPM model

## Installation and Compilation

### Requirements
- GCC (or other C++ compiler)
- `libtbb` (Intel Threading Building Blocks)

### Installation

1. **Clone repository**:
```sh
git clone https://github.com/JasmineGamblin/PPMmodelPangenome
cd PPMmodelPangenome
```

2. **Install TBB**:

Sur Ubuntu/Debian :
```sh
sudo apt-get install libtbb-dev
```
Sur Fedora :
```sh
sudo dnf install tbb-devel
```
Sur macOS :
```sh
brew install tbb
```

### Compilation
```sh
make
```

## Run inference
To run an inference with a random starting point, use the following command:
```sh
./inference seed "tree.nwk" "pa_matrix.txt" par_nb "mle_param.txt" "inf_cat.txt"
```
where:
- `seed` is the random seed
- `tree.nwk` is the file containing the species tree
- `pa_matrix.txt` is the file containing the presence/absence matrix
- `par_nb` is the number of free parameters (7 or 9, depending if the error parameter is the same for all categories or not)
- `mle_param.txt` is the file where the paramter estimates will be stored
- `inf_cat.txt` is the file where the inferred category for each gene will be stored


To run an inference with a chosen starting point, use instead (example with 9 free parameters):
```sh
./inference 0 "tree.nwk" "pa_matrix.txt" 9 "mle_param.txt" "inf_cat.txt" N0 l0 i1 l1 g2 l2 eps0 eps1 eps2
```
where `N0`, `l0`, `i1`, `l1`, `g2`, `l2`, `eps0`, `eps1` and `eps2` are replaced by the chosen initial values.
