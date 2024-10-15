# ClustIRR: Clustering of Immune Receptor Repertoires

## Overview 
ClustIRR is a quantitative method designed for the clustering of immune 
receptor repertoires (IRRs). It identifies groups of T cell or B cell 
receptors (TCRs or BCRs) that share antigen specificity across different 
biological conditions or longitudinal samples. Using a Bayesian hierarchical 
model, ClustIRR performs differential cluster occupancy (DCO) analysis to 
detect specificity clusters that are expanding or contracting in response 
to each biological condition.

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

