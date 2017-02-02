# CCor
A network-based similarity measure for clustering genes into co-expression modules

## Introduction
CCor is a whole genome network-based similarity measure, that makes use of information of all the genes in the network. Application on microarray data and theoretical properties can be found at http://onlinelibrary.wiley.com/doi/10.1111/biom.12508/abstract. 

## Prerequisites
This package was developed under R 3.1.2 and requires the following packages
* miscTools
* ROCR
* GeneNet
* space
* WGCNA
To install these packages, you can conveniently use the **install.packages** function in R terminal. For example,
```
install.packages("miscTools")
```

## Input Data
* A matrix with each row being the expression of a gene, e.g. microarray or RNAseq data

## Setup and Usage Example
1) Clone this repository
```
git clone https://github.com/yiminghu/AnnoPred.git
```
2) Compile the C codes for efficient calculation of CCor and mCCor
```
R CMD SHLIB ccor.c
R CMD SHLIB mccor.c
```
This will generate the following files: ccor.o, ccor.so, mccor.o and mccor.so. We suggest using command line instead of RStudio when calculating CCor and mCCor. Since when the number of genes is large, e.g. larger than 1000, a memory error will happen with RStudio.

3) CCor_mCCor.R computes the CCor and mCCor. 

4) Concentration.R computes the derived concentrating value in the concentration inequality.

5) Illustration.R computes the within and between-module concentrating values shown in Inference section of the paper.

6) simulation.R and simulation_tom.R are used to perform simulation and generate Table 2 and 3.

## Some tips on processing Microarray data
Since most of the microarray data are deposited on GEO <a href="https://www.ncbi.nlm.nih.gov/geo/">GEO</a>
