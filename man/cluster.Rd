\name{cluster_irr}
\alias{cluster_irr}
\title{Clustering of immune receptor repertoires}
\usage{
cluster_irr(
  s,
  r,
  version = 3,
  ks = 4,
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
\item{s}{data.frame, Complementarity determining regions 3 (CDR3s)
observed in an immune receptor repertoire (IRR). The data.frame can have either
one column or two columns:
\itemize{
 \item One column: s contains CDR3s from a single chain: \emph{CDR3b}
 or \emph{CDR3a}
 \item Two columns: s contains CDR3s from both chains (paired):
    \emph{CDR3b} and \emph{CDR3a}
}
\emph{CDR3b} stands for CDR3 sequences derived from one of the following:
\itemize{ \item \eqn{\beta} chain T cell receptors
\item \eqn{\gamma} chain T cell receptors
\item Heavy chain B cell receptors }

\emph{CDR3a} stands for CDR3 sequences derived from one of the following:
\itemize{
\item \eqn{\alpha} chain T cell receptors
\item \eqn{\delta} chain T cell receptors
\item Light chain B cell receptors
}

}
\item{r}{data.frame, reference (or control) repertoire of CDR3 sequences. Must
have the same structure (number of columns and column names) as \code{s}}
\item{version}{integer, version of the algorithm: \code{version} = 1, 2 or
3 (default)}
\item{ks}{integer or integer vector, motif lengths. \code{ks} = 4 (default)}
\item{cores}{integer, number of CPU cores, \code{cores} = 1 (default).}
\item{control}{list, a named list of auxiliary parameters
to control algorithm's behavior. See the details below:
\itemize{\item \code{B} - integer, number of bootstrap samples. Only used by
  version = 1. \code{B} = 1000 (default)
  \item \code{global_max_dist} - number, Hamming distance (HD) threshold to
  consider two CDR3s as globally clustered. CDR3s are globally clustered if
  HD(\eqn{a}, \eqn{b}) \eqn{\leq} \code{global_max_dist}.
  \code{global_max_dist} = 1 (default)
  \item \code{local_max_fdr} - numeric, maximum False Discovery Rate (FDR) for
  the detection of enriched motifs. \code{local_max_fdr} = 0.05 (default)
  \item \code{local_min_ove} - numeric, minimum fold change between observed
  and expected relative abundances for the detection of enriched motifs.
  \code{local_min_ove} = 2 (default)
  \item \code{local_min_o} - numeric, minimum absolute frequency of a motif in
  the s in order for the motif to be used in the enrichment analysis.
  \code{local_min_o} = 1 (default)
  \item \code{trim_flank_aa} - integer, how many amino acids should be trimmed
  from the flanks of all CDR3 sequences (only used for local clustering.
  \code{trim_flank_aa} = 0 (default))
  \item \code{low_mem} - logical, allows low memory mode for global clustering.
  This will lead to increase in the CPU time but lead to a lower memory
  footprint. \code{low_mem} = \code{FALSE} (default)
  \item \code{global_pairs} - matrix, pre-computed global pairs. If
  \code{global_pairs} is provided by the user, then global clustering is not
  performed. Instead the CDR3 pairs from global_pairs are used as global
  clustering pairs. \code{global_pairs} is a character matrix with 3 columns.
  The first two columns contain pairs of CDR3 sequences. These are considered
  globally clustered. The third column contains information about the TCR chain
  of each pair of CDR3s: \emph{TRA} or \emph{TRB}. \code{global_pair} =
  \code{NULL} (default)
}
}
}
\value{
The output is an object with \code{clust_irr} structure. This object is formed
by two sublists:
\item{clust}{list, contains clustering results
\itemize{
      \item \code{CDR3b} - list, contains global+local clustering results for
      CDR3s from chain \emph{TRB}
      \itemize{
        \item \code{local} - list, local clustering results
          \itemize{
            \item \code{m} - data.frame, motif enrichment results

            For version = 2 and version = 3, \code{m} has the following set of
            columns:
            \itemize{
            \item \code{motif} - motif sequence
            \item \code{f_s} - observed motif counts in \code{s}
            \item \code{f_r} - observed motif counts in \code{r}
            \item \code{n_s} - number of all observed motifs in \code{s}
            \item \code{n_r} - number of all observed motifs in \code{r}
            \item \code{k} - motif length
            \item \code{ove} - mean observed/expected relative motif frequency
            \item \code{ove_ci_l95} - 95\% confidence intervals of ove (lower
            boundary)
            \item \code{ove_ci_h95} - 95\% confidence intervals of ove (upper
            boundary)
            \item \code{p_value} - p-value from Fisher's exact test
            \item \code{fdr} - false discovery rate, i.e. adjusted p-value by
            Benjamini & Hochberg correction
            \item \code{pass} - logical value indicating whether a motifs are
            enriched (\code{pass=TRUE}) given the user-defined thresholds in
            control}

            For version = 1, \code{m} has the following columns:
            \itemize{
            \item \code{motif} - motif sequence
            \item \code{k} - motif length
            \item \code{f_s} - motif counts in \code{s}
            \item \code{mean_f_r} - mean motif counts in \code{r}
            computed based on \code{B} simulations
            \item \code{min_f_r} - maximum motif counts in \code{r}
            computed based on \code{B} simulations
            \item \code{max_f_r} - maximum motif counts in \code{r}
            computed based on \code{B} simulations
            \item \code{ove} - observed (\code{f_s)} / expected
            (\code{mean_f_r)} motif relative frequency
            \item \code{p_value} - probability of observing a motif count in
            the \code{r} that is at least as high or higher as in the \code{s}
            \item \code{fdr} - false discovery rate, i.e. adjusted p-value by
            Benjamini & Hochberg correction
            \item \code{pass} - logical value indicating whether motifs are
            enriched (pass=TRUE) given the user-defined thresholds in control
            }

            \item \code{lp} - data.frame, enriched motifs are linked to their
            original CDR3 sequences and shown as rows in the data.frame

            Contains the following columns:
            \itemize{
            \item \code{cdr3} - CDR3 sequence
            \item \code{cdr3_core} - corresponding core of the CDR3 sequence,
            obtained by trimming \code{trim_flank_aa} amino acids (user-defined
            parameter). If \code{trim_flank_aa} = 0, then
            \code{cdr3} = \code{cdr3_core}
            \item \code{motif} - enriched motif from \code{cdr3_core}
            }
          }
        \item \code{global} - matrix, global clustering results. Pairs of
        globally similar CDR3s are shown in each row of the matrix
        (analogous to \code{lp})
      }
      \item \code{CDR3a} - list, contains global+local clustering results for
      CDR3s from chain \emph{TRA} (same structure as in \code{CDR3b})
    }
}

\item{inputs}{list, contains all user provided inputs (see \bold{Arguments})}
}

\description{
This algorithm finds groups of TCRs that likely have similar peptide:MHC
specificity.

Two strategies are used:
\enumerate{
    \item Local clustering
    \item Global clustering
}

\bold{Local clustering procedure}

\enumerate{
\item CDR3 processing steps
    \itemize{
        \item Get set of non-redundant CDR3s (only in version = 1 or 2)
        \item Trim CDR3 ends
    }
\item Motif processing steps
    \itemize{
          \item motif frequencies in data set \code{s} (\eqn{f_s}) and \code{r}
          (\eqn{f_r})
          \item total number of motifs in data set \code{s} (\eqn{n_s}) and
          \code{r} (\eqn{n_r})
          \item ratio of observed vs. expected motif counts using the following
          formula: OvE=\eqn{(f_s/n_s)/(f_r/n_r})
          \item In version 2 and 3, probability \eqn{p_i} of finding the
          observed or a larger OvE for motif \eqn{i} given that the
          null hypothesis is true is computed with the Fisher's exact test
          \item In version 1, \eqn{p_i} is computed by bootstrapping
          \item classify motif \eqn{i} as \code{pass=TRUE} if the motif passes
          all filters specified in the user-provided control list, otherwise
          as \code{pass=FALSE}
    }

}


\bold{Global clustering}

ClustIRR employs two approaches for global clustering:

\enumerate{
\item For equal-length CDR3 sequences \eqn{i} and \eqn{j} we compute the
Hamming distance \eqn{d_{ij}}. If \eqn{d_{ij}\leq} \code{global_max_dist}
(user-defined input), then \eqn{i} and \eqn{j} are globally similar.
\item The user can provide a matrix of globally similar CDR3 sequence pairs,
computed by a complementary approach (e.g. tcrdist)
}
}

\examples{
# load package input data
data("CDR3ab")
s <- data.frame(CDR3b = CDR3ab[1:1000, "CDR3b"])
r <- data.frame(CDR3b = CDR3ab[1:5000, "CDR3b"])

# run analysis
out <- cluster_irr(s = s,
                   r = r,
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

# inspect motif enrichment results
knitr::kable(head(out$clust$CDR3b$local$m))

# inspect which CDR3bs are globally similar
knitr::kable(head(out$clust$CDR3b$global))
}
