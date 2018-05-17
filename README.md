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
git clone https://github.com/yiminghu/CCor.git
```
2) Compile the C codes for efficient calculation of CCor and mCCor
```
R CMD SHLIB ccor.c
R CMD SHLIB mccor.c
```
This will generate the following files: ccor.o, ccor.so, mccor.o and mccor.so. We suggest using command line instead of RStudio when calculating CCor and mCCor. Since when the number of genes is large, e.g. larger than 1000, a memory error will happen with RStudio.

3) CCor_mCCor.R computes the CCor and mCCor. 

4) Concentration.R computes the derived concentrating value in the concentration inequality.

5) Illustration.R computes the within and between-module concentrating values shown in Theoretical Results section of the <a href="https://www.ncbi.nlm.nih.gov/pubmed/26953524">paper</a>.

6) simulation.R and simulation_tom.R are used to perform simulation and generate Table 2 and 3 of the <a href="https://www.ncbi.nlm.nih.gov/pubmed/26953524">paper</a>.

## Some tips on processing Microarray data
Since most of the microarray data are deposited on <a href="https://www.ncbi.nlm.nih.gov/geo/">GEO</a>, we therefore provide some codes for conveniently querying, processing and transforming data from GEO to the ideal format.

1) install necessary pacakges
* GEOquery
* affy
* hgu133plus2.db
* annotate
Since some of the packages are deposited on <a href="https://www.bioconductor.org/">bioconductor</a> instead of <a href="https://cran.r-project.org/">CRAN</a>, we will need to use the following commands to install them.
```
source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
```
2) After install and load necessary packages, we can start query gene expression data on GEO with accession number, take <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37418">GSE37418</a> as an example.
```
getGEOSuppFiles("GSE37418")
```
3) Then we unzip the raw data files into a directory named "data"
```
untar("GSE37418/GSE37418_RAW.tar", exdir="data")
cels <- list.files("data/", pattern = "[gz]")
sapply(paste("data", cels, sep="/"), gunzip)
setwd('data/')
```
4) Read in all CEL files and normalize data with RMA method (Highly recommended!!!)
```
AB <- ReadAffy()
ESet <- rma(AB,normalize=T) 
ESet.filtered <- nsFilter(ESet, require.entrez=FALSE, remove.dupEntrez=FALSE)
ESet.filtered$filter.log
```
5) Map probes to gene symbols (convenient for future analysis)
```
dat <- exprs(ESet.filtered$eset)
prb.list <- row.names(dat)
enID <- select(hgu133plus2.db, prb.list, c("SYMBOL", "ENTREZID", "GENENAME"))
enID1 <- enID[!duplicated(enID[,1]), 3]
gs_map <- enID1[!is.na(enID1)]
dat_map <- dat[!is.na(enID1),]
```
6) Remove duplicated observations
```
dupgene <- names(table(gs_map))[table(gs_map)>1]
indupgene <- names(table(gs_map))[table(gs_map)==1]
expfile.indup <- matrix(0,length(indupgene),ncol(dat_map))
for(i in 1:length(indupgene)){
  expfile.indup[i,] <- as.numeric(dat_map[which(gs_map==indupgene[i],arr.ind=T),])
}
expfile.dup <- matrix(0,length(dupgene),ncol(dat_map))
for(i in 1:length(dupgene)){
  expfile.dup[i,] <- apply(as.matrix(dat_map[which(gs_map==dupgene[i],arr.ind=T),]),2,mean)
}
exp.output <- rbind(expfile.indup,expfile.dup)
row.names(exp.output) <- c(indupgene,dupgene)
```
7) Write normalized expression data and corresponding gene symbols to files
```
write.table(exp.output,'tumor.expr.txt',quote=F,col.names=F,row.names=F)
write.table(row.names(exp.output),'tumor.gene.txt',quote=F,col.names=F,row.names=F)
```




