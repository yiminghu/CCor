# CCor
A network-based similarity measure for clustering genes into co-expression modules

## Introduction
CCor is a whole genome network-based similarity measure, that makes use of information of all the genes in the network. Application on microarray data and theoretical properties can be found at http://onlinelibrary.wiley.com/doi/10.1111/biom.12508/abstract. 

## Prerequisites

CCor_mCCor.R computes the CCor and mCCor. Before use, please open command line and change directory to where the ccor.C and mccor.C are saved. Then please run ‘R CMD SHLIB ccor.c’ and ‘R CMD SHLIB mccor.c’ in command line, which will generate the following files: ccor.o, ccor.so, mccor.o and mccor.so. We suggest using command line instead of RStudio to run CCor_mCCor.R too (change directory to the same directory as above). Since when the number of genes is large, e.g. larger than 1000, an error will happen with RStudio.

Concentration.R computes the derived concentrating value in the concentration inequality.

Illustration.R computes the within and between-module concentrating values shown in Inference section of the paper.

simulation.R and simulation_tom.R are used to perform simulation and generate Table 2 and 3.

