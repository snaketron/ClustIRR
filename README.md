# ClustIRR: Clustering of Immune Receptor Repertoires

![Build Status](https://github.com/snaketron/ClustIRR/actions/workflows/r.yml/badge.svg)

## Overview

ClustIRR is a computational framework for uncertainty-aware quantitative analysis of immune receptor repertoires (IRRs). It maps multiple repertoires across replicates, time points, or conditions onto a joint graph induced by immune receptor sequence similarity, detects communities on the joint graph (CJs), and uses these communities as reference structures for Bayesian analyses of differential community occupancy. It can also integrate single-cell gene expression data to relate community expansion to transcriptional activation signatures.

This framework is described in the eLife paper: Kitanovski et al., “Uncertainty-aware quantitative analysis of the structure and dynamics of T cell receptor repertoires” (eLife, 2026, DOI: 10.7554/eLife.112537). See `inst/CITATION` for the citation.

ClustIRR is designed for applications such as quantitative changes after antigen challenge, longitudinal tracking during cancer immunotherapy, and comparison of V(D)J biases across human and murine repertoires.

## How to use ClustIRR

ClustIRR is an R package available from Bioconductor:

https://bioconductor.org/packages/ClustIRR/

To install from Bioconductor, start R and enter:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ClustIRR")
```

To install the development version from GitHub:

```r
library(devtools)
install_github("snaketron/ClustIRR")
```

Case studies are provided in the `vignettes/` directory.

## Workflow & output

![clustirr workflow](/inst/extdata/logo.png)
