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
 \item One column: data_sample contains CDR3s from a single chain: \emph{CDR3b} 
 or \emph{CDR3a}
 \item Two columns: data_sample contains CDR3s from both chains (paired):
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
\item{data_ref}{data.frame, reference database of CDR3 sequences. Must have the 
same structure (number of columns and column names) as \code{data_sample}}
\item{version}{integer, version of the algorithm: \code{version} = 1, 2 or 
3 (default)}
\item{ks}{integer or integer vector, motif lengths. \code{ks} = 2:4 (default)}
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
  the data_sample in order for the motif to be used in the enrichment analysis.
  \code{local_min_o} = 1 (default)
  \item \code{trim_flank_aa} - integer, how many amino acids should be trimmed
  from the flanks of all CDR3 sequences (only used for local motif
  clustering. \code{trim_flank_aa} = 0 (default))
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
The output is an object with \code{clust_irr} structure. This object is formed by
two sublists:
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
            \item \code{f_sample} - observed motif counts in sample data
            \item \code{f_ref} - observed motif counts in reference data
            \item \code{n_sample} - number of all observed motifs in sample set
            \item \code{n_ref} - number of all observed motifs in reference set
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
            enriched (\code{pass=TRUE}) given the user-defined thresholds in control
            }
            
            For version = 1, \code{m} has the following columns:
            \itemize{
            \item \code{motif} - motif sequence
            \item \code{k} - motif length
            \item \code{f_sample} - motif counts in sample
            \item \code{mean_f_ref} - mean motif counts in reference set 
            computed based on \code{B} simulations
            \item \code{min_f_ref} - maximum motif counts in reference set 
            computed based on \code{B} simulations
            \item \code{max_f_ref} - maximum motif counts in reference set 
            computed based on \code{B} simulations
            \item \code{ove} - observed (\code{f_sample)} / expected 
            (\code{mean_f_ref)} motif relative frequency
            \item \code{p_value} - probability of observing a motif count in 
            the reference set that is at least as high or higher as in the 
            sample set
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
    
\item{inputs}{list, contains user provided inputs: 
\itemize{
\item \code{data_sample} - input CDR3 sequences
\item \code{version} - \emph{ClustIRR} version that was used
\item \code{ks} - kmer-count to be searched for 
\item \code{cores} - used CPU cores
\item \code{control} - other control parameters, details see \bold{Arguments})
}
}

}

\description{
This algorithm finds groups of TCRs that likely have similar peptide:MHC
specificity.

Two strategies are used: 
\enumerate{
    \item Local motif clustering
    \item Global CDR3 sequence clustering
}

\bold{Local motif clustering procedure}

\enumerate{
\item CDR3 processing steps 
    \itemize{
        \item Get set of non-redundant CDR3s (only in version = 1 or 2)
        \item Trim CDR3 ends
    }
\item Motif processing steps
    \itemize{
          \item motif frequencies in sample (\eqn{f_{sample}}) and reference 
          (\eqn{f_{ref}})
          \item total number of motifs in sample (\eqn{n_{sample}}) and 
          reference (\eqn{n_{ref}})
          \item ratio of observed vs. expected motif counts using the following
          formula:
          
          \eqn{OvE} = \eqn{(f_{sample}/n_{sample})/(f_{ref}/n_{ref}})
          \item In version 2 and 3, probability \eqn{p_i} of finding the 
          observed or a larger \eqn{OvE} for motif \eqn{i} given that the 
          null hypothesis is true is computed with the Fisher's exact test
          \item In version 1, \eqn{p_i} is computed by bootstrapping
          \item classify motif \eqn{i} as \code{pass=TRUE} if the motif passes 
          all filters specified in the user-provided control list, otherwise 
          as \code{pass=FALSE}
    }

}


\bold{Global CDR3 sequence clustering}

ClustIRR provides two possible approaches for global clustering:

\enumerate{
\item For equal-length CDR3 sequences \eqn{i} and \eqn{j} compute their Hamming
distance \eqn{d_{ij}}

Sequences with \eqn{d_{ij}\leq} \code{global_max_dist}
(user-defined input) are considered globally similar
\item  Alternatively, the user can provide a matrix of globally similar CDR3 
sequences computed by a complementary approach (e.g. \emph{tcrdist})
}

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
