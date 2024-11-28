# Inference tool for the PPM model

## Description
This repository contains the source code used in the manuscript 'Persistent, Private and Mobile genes: a model for gene dynamics in evolving pangenomes' by Jasmine Gamblin, Amaury Lambert and François Blanquart.

## Files
- folder `source`: source code for inference tool
- folder `test`: test files containing simulated tree and matrix to run a test inference
- `Makefile`: compilation file
- `simulation_aux.R`, `run_simulations.R`: R scripts used to simulate the PPM model

## Installation and Compilation

### 1. Requirements
- GCC (or other C++ compiler)
- `libtbb` (Intel Threading Building Blocks)

Install TBB on Ubuntu/Debian:
```sh
sudo apt-get install libtbb-dev
```
on Fedora :
```sh
sudo dnf install tbb-devel
```
on macOS :
```sh
brew install tbb
```

### 2. Clone repository
```sh
git clone https://github.com/JasmineGamblin/PPMmodelPangenome
cd PPMmodelPangenome
```

### 3. Compilation
```sh
make
```

## Running an inference

### 1. Input file format
- Species tree `tree.nwk`: Newick format with leaf labels (handles node labels by ignoring it). Tree must be ultrametric. Example:
```txt
(genome2:0.35,genome1:0.35):0;
```

- Presence/absence matrix `pa_matrix.txt`: Matrix in csv format (comma-separated values) where rows are genomes and columns are genes. Matrix must contain row names (genome IDs matching the tree leaf labels) and column names (which are ignored). Values must be either 0 or 1. Example:
```txt
gene1,gene2,gene3
genome1,1,0,1
genome2,0,1,1
```

### 2. Execute code
To run an inference with a random starting point, use the following command:
```sh
./inference seed "test/tree.nwk" "test/pa_matrix.txt" "test/mle_param.txt" "test/inf_cat.txt"
```
where:
- `seed` is the random seed (must be an integer different from 0, as 0 indicates that the user is chosing the starting point)
- `tree.nwk` and `pa_matrix.txt` are the input files
- `mle_param.txt` and `inf_cat.txt` are the output files


To run an inference with a chosen starting point, use instead:
```sh
./inference 0 "test/tree.nwk" "test/pa_matrix.txt" "test/mle_param.txt" "test/inf_cat.txt" N0 l0 i1 l1 g2 l2 s_10 s_01
```
where `N0`, `l0`, `i1`, `l1`, `g2`, `l2`, `s_10`, and `s_01` are replaced by the chosen initial values.

Inference should take around 5 minutes on a laptop with the provided test data (20 genomes, 949 genes), but using a cluster is recommended for bigger datasets.

### 3. Output file format
- Parameter estimates `mle_param.txt`: values are stored in the following order: `seed`, `N0`, `l0`, `i1`, `l1`, `i2`, `g2`, `l2`, `s_10`, `s_01`, and the maximum log-likelihood value reached

- Inferred gene categories `inf_cat.txt`: file containing inferred category number for each gene (0 for Persistent, 1 for Private and 2 for Mobile), separated by blank spaces
