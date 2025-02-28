# ClustIRR: Clustering of Immune Receptor Repertoires

![Build Status](https://github.com/<OWNER>/<REPO>/actions/workflows/r.yml/badge.svg)

## Overview 
ClustIRR analyzes repertoires of B- and T-cell receptors. It begins the 
analysis by identifying communities (i.e., specificity groups) of immune 
receptors with similar specificities, based on the sequences of their 
complementarity determining regions (CDRs). Next, it employs a Bayesian 
probabilistic models to quantify differential community occupancy (DCO) 
between repertoires, allowing the identification of expanding or contracting 
communities in response to factors such as infection or cancer treatment.

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

