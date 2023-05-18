# gliphR

Clustering of T-cell receptors (TCRs) to find groups of TCRs with similar 
specificity. Quantitative comparison of TCR cluster structure between 
biological conditions. Includes low memory mode to reduce memory usage 
for large input files.

A completely revised version of the 
[turboGliph](https://github.com/HetzDra/turboGliph) package focused on the 
clustering part of the Gliph algorithm. Input parameters are reduced to
a minimum for user convenience.

Includes R implementations of the 
[Gliph](https://github.com/immunoengineer/gliph) and 
[Gliph2](http://50.255.35.37:8080) clustering algorithms.
Introduces a new version of the Gliph algorithm called `gliph3`,
designed for single-cell sequencing input data.

See the package vignette for further information.

### References

**1**: Glanville, Jacob, et al. "Identifying specificity groups in the 
T cell receptor repertoire." Nature 547.7661 (2017): 94.<br>
**2**: Huang, Huang, et al. "Analyzing the Mycobacterium tuberculosis immune 
response by T-cell receptor clustering with GLIPH2 and genome-wide antigen 
screening." Nature Biotechnology 38.10 (2020): 1194-1202.<br>

