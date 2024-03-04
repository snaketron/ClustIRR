# ClustIRR: Clustering of Immune Receptor Repertoires

## Overview 
ClustIRR is a quantitative method for *clustering* of immune receptor 
repertoires (IRRs). The algorithm identifies groups of T or B cell 
receptors (TCRs or BCRs) that likely have the same antigen specificity. 
This is achieved by comparing *global* and *local* sequence features of the 
complementarity determining regions (CDRs) of TCRs and BCRs. Once the 
specificity groups are identified in a set of IRRs, the algorithm identifies
clusters and performs a quantitative comparison of the cluster abundances 
between biological conditions. 

## How to use ClustIRR
ClustIRR is an R-package available from Bioconductor: 

https://bioconductor.org/packages/ClustIRR/

To install this package, start R and enter:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ClustIRR")
```

To use the development version of ClustIRR you can install it from Github.
To do this, start R and enter:

```r
library(devtools)
install_github("snaketron/ClustIRR")
```


Case studies are provided in the directory /vignettes

## Workflow & output 

![clustirr workflow](/inst/extdata/logo.png)

