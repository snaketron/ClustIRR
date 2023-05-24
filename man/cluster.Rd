\name{cluster_irr}
\alias{cluster_irr}
\title{Clustering of immune receptor repertoires}
\usage{
cluster_irr(
  data_sample,
  data_ref,
  version = 3,
  ks = c(2, 3, 4),
  cores = 1,
  control = list(B = 1000,
                 global_max_dist = 1,
                 local_max_fdr = 0.05,
                 local_min_ove = 2,
                 local_min_o = 1,
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
}
}
}
\value{
The output is an object with 'clust_irr' structure. This object is formed by
two sublists ('clust'= clustering results; and 'inputs'= original user inputs):
\itemize{
   \item clust, list, clustering results
   \itemize{
      \item CDR3b, list, contains global+local clustering results for CDR3s
      from chain TRB
      \itemize{
        \item local, list, local clustering results
          \itemize{
            \item m, data.frame, with motif enrichment results

            For version=2 and version=3, data.frame m has the following set of
            columns: 'motif': motif sequence; 'f_sample' & 'f_ref': observed
            motif counts in sample and reference data; n_sample & n_ref: number of
            all observed motifs in sample and reference set; k: motif length;
            ove, ove_ci_l95, ove_ci_h95: mean observed/expected relative motif
            frequency and 95\% confidence intervals of ove; p_value: p-value
            from Fisher's exact test; fdr: false discovery rate, i.e. adjusted
            p-value by Benjamini & Hochberg correction; pass: logical value
            indicating whether a motifs are enriched (pass=TRUE) given the
            user-defined thresholds in control.

            For version = 1, data.frame m has the following columns:
            motif = motif sequence; k: motif length, f_sample: motif counts in
            sample; mean_f_ref, min_f_ref, max_f_ref: mean, minimum and maximum
            motif counts in reference set computed based on B simulations; ove:
            observed (f_sample)/expected (mean_f_ref) motif relative frequency;
            p: probability of observing a motif count in the reference set that
            is at least as high or higher as in the sample set; fdr: false
            discovery rate, i.e. adjusted p-value by Benjamini & Hochberg
            correction; pass: logical value indicating whether motifs are
            enriched (pass=TRUE) given the user-defined thresholds in control.

            \item lp, data.frame, enriched motifs are linked to their original
            CDR3 sequences and shown as rows in the data.frame

            The data.frame lp contains the following columns: cdr3 & cdr3_core:
            CDR3 sequence and the corresponding core of the CDR3 sequence
            obtained by trimming 'trim_flank_aa' amino acids (user-defined
            parameter in control). if trim_flank_aa=0, then cdr3=cdr3_core;
            motif: enriched motif from cdr3_core.
          }
        \item global, matrix, global clustering results. Pairs of globally
        similar CDR3s are shown in each row of the matrix (analogous to lp)
      }
      \item CDR3a, list, contains global+local clustering results for CDR3s
      from chain TRA (has the same structure as in CDR3b)
    }
   \item inputs, list, user provided inputs (see input arguments: data_sample,
   version, ks, cores, control)
}
}

\description{
This algorithm finds groups of TCRs that likely have similar peptide:MHC
specificity. This is done with two strategies: local motif clustering and
global CDR3 sequence clustering.

SK: here we need two paragraphs of how each algorithm (local vs. global)
works, as well as *key* differences between version=1, 2 and 3
}

\examples{
# load package input data
data("CD8")
data_sample <- data.frame(CDR3b = CD8[1:5000, "CDR3b"])
data_ref <- CD8

# run analysis
out <- cluster_irr(data_sample = data_sample,
                   data_ref = data_ref,
                   version = 3,
                   ks = 4,
                   cores = 1,
                   control = list(
                      B = 1000,
                      global_max_dist = 1,
                      local_max_fdr = 0.05,
                      local_min_ove = 2,
                      local_min_o = 1,
                      trim_flank_aa = 3,
                      global_pairs = NULL,
                      low_mem = FALSE))

# output class
class(out)

# output structure
str(out)

# inspect local motif enrichment results
knitr::kable(head(out$clust$CDR3b$local$m))

# inspect which CDR3bs are globally similar
knitr::kable(head(out$clust$CDR3b$global))
}
