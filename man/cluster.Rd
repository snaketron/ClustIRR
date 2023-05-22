\name{cluster_irr}
\alias{cluster_irr}
\title{Clustering of immune receptor repertoires}
\usage{
cluster_irr(
  data_sample,
  data_ref,
  version = 2,
  ks = c(2, 3, 4),
  cores = 1,
  control = list(B = 1000,
                 global_max_dist = 1,
                 local_max_fdr = 0.05,
                 local_min_ove = 2,
                 local_min_o = 3,
                 trim_flank_aa = 0,
                 global_pairs = NULL,
                 low_mem = FALSE)
)
}
\arguments{
\item{data_sample}{data.frame, Complementarity determining regions 3 (CDR3s)
observed in an immune receptor repertoire (IRR). The data.frame can have either
one column or two columns:
\itemize{
 \item 1 column: data_sample contains CDR3s from a single chain: 'CDR3b' or
 'CDR3a'
 \item 2 columns: data_sample contains CDR3s from both chains (paired):
    'CDR3b' and 'CDR3a'
}
'CDR3b' stands for CDR3 sequences derived from i) beta chain T cell receptors
or ii) from gamma chain T cell receptors; or iii) from heavy chain B cell
receptors

'CDR3a' stands for CDR3 sequences derived from i) alpha chain T cell receptors
or ii) from delta chain T cell receptors; or iii) from light chain B cell
receptors
}
\item{data_ref}{data.frame, reference database of CDR3 sequences. This
data.frame must have the same structure (number of columns and column names)
as data_sample}
\item{version}{integer, version of the algorithm: version=1, 2 or 3 (default)}
\item{ks}{integer or integer vector, motif lengths ks=2:4 (default)}
\item{cores}{integer, number of CPU cores, cores=1 (default).}
\item{control}{list, a named list of auxiliary parameters
to control algorithm's behavior. See the details below:
\itemize{
  \item {B}, integer, number of bootstrap samples. B is used only by
  version=1. B=1000 (default)
  \item {global_max_dist}, number, Hamming distance (HD) threshold to consider
  two CDR3s as globally clustered. CDR3s \eqn{a} and \eqn{b} are globally
  clustered if HD(\eqn{a}, \eqn{b}) \eqn{\leq} global_max_dist.
  global_max_dist=1 (default)
  \item {local_max_fdr}, numeric, maximum False Discovery Rate (FDR) for the
  detection of enriched motifs. local_max_fdr=0.05 (default)
  \item {local_min_ove}, numeric, minimum fold change between observed and
  expected relative abundances for the detection of enriched motifs.
  local_min_ove=2 (default)
  \item {local_min_o}, numeric, minimum absolute frequency of a motif in the
  data_sample in order for the motif to be used in the enrichment analysis.
  local_min_o=1 (default)
  \item {trim_flank_aa}, integer, how many amino acids should be trimmed
  from the flanks of all CDR3 sequences (only used for local motif
  clustering. trim_flank_aa=0 (default))
  \item {low_mem}, logical, allows low memory mode for global clustering. This
  will lead to increase in the CPU time but lead to a lower memory footprint.
  low_mem=FALSE (default)
  \item {global_pairs}, matrix, pre-computed global pairs. If global_pairs
  is provided by the user, then global clustering is not performed. Instead
  the CDR3 pairs from global_pairs are used as global clustering pairs.
  global_pairs is a character matrix with 3 columns. First two columns contain
  pairs of CDR3 sequences. These are considered globally clustered. The third
  column contains information about the TCR chain of each pair of CDR3s: TRA or
  TRB. global_pairs=NULL (default)
}}
}
\value{
The output is a list of the following elements:
\itemize{
   \item clust: clustering results for each TCR chain
   \item edges: local + clustering connections between rows in data_sample
   \item data_sample: provided input data_sample
   \item version: provided input version
   \item ks: provided input ks
   \item cores: provided input cores
   \item control: provided input control
}
}

\description{
This algorithm finds groups of TCRs that likely have similar peptide:MHC
specificity. This is done with two strategies: local motif clustering and
global CDR3 sequence clustering.

SK: here we need two paragraphs of how each algorithm (local vs. global)
works, as well as *key* differences between v=1, 2 and 3
}

\examples{
# load package input data
data("hs_CD8_ref")
data_sample <- data.frame(CDR3b = hs_CD8_ref[1:(10^3), "CDR3b"])
data_ref <- data.frame(CDR3b = hs_CD8_ref[1:(10^4), "CDR3b"])

# detect cores
# SK: we use the R-package future for multicore execution. parallel
# is not available on Win OS + we do not declare it in the DESCRIPTION
# file. Also, we do not use cores later on in the main cluster call, so
# why declare it here?
# cores <- future::availableCores()

# run analysis
out <- cluster_irr(
data_sample = data_sample,
data_ref = data_ref
)

}
