\name{cluster_irr}
\alias{cluster_irr}
\title{Clustering of immune receptor repertoires}
\usage{
cluster_irr(
    s,
    r,
    version = 2,
    ks = 4,
    cores = 1,
    control = list(global_max_dist = 1,
                   local_max_fdr = 0.05,
                   local_min_ove = 2,
                   local_min_o = 1,
                   trim_flank_aa = 0,
                   global_pairs = NULL,
                   low_mem = FALSE))
}
\arguments{
\item{s}{data.frame, complementarity determining region 3 (CDR3) amino acid
sequences observed in an immune receptor repertoire (IRR). The data.frame can 
have either one column or two columns:
\itemize{
\item One column: s contains CDR3s from a single chain: \emph{CDR3b}, 
\emph{CDR3a}, \emph{CDR3g}, \emph{CDR3d}, \emph{CDR3h} or \emph{CDR3l}
\item Two columns: s contains CDR3s from both chains (paired), for instance:
    \itemize{
        \item \emph{CDR3b} and \emph{CDR3a} [for \eqn{\alpha\beta} TCRs]
        \item \emph{CDR3g} and \emph{CDR3d} [for \eqn{\gamma\delta} TCRs]
        \item \emph{CDR3h} and \emph{CDR3l} [for heavy/ligh chain BCRs]
    }
}
}
\item{r}{data.frame, reference (or control) repertoire of CDR3 sequences. Must
have the same structure (number of columns and column names) as \code{s}}
\item{version}{integer, version of the algorithm: \code{version} = 1 or 2 
(default)}
\item{ks}{integer or integer vector, motif lengths. \code{ks} = 4 (default)}
\item{cores}{integer, number of CPU cores, \code{cores} = 1 (default).}
\item{control}{list, a named list of auxiliary parameters
to control algorithm's behavior. See the details below:
\itemize{
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
The output is an S4 object of class \code{clust_irr}. This object contains 
two sublists:
\item{clust}{list, contains clustering results for each TCR/BCR chain. The
results are stored in separate sub-list named appropriately (e.g. CDR3a, CDR3b, 
CDR3g, etc.). In the following we who the typical structure of these lists:
\itemize{
    \item \code{local} - list, local clustering results
    \itemize{
        \item \code{m} - data.frame, motif enrichment results with columns:
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
            \item \code{lp} - data.frame, enriched motifs are linked to their
            original CDR3 sequences and shown as rows in the data.frame with
            the following columns:
            \itemize{
                \item \code{cdr3} - CDR3 amino acid sequence
                \item \code{cdr3_core} - core portion of the CDR3 sequence,
                obtained by trimming \code{trim_flank_aa} amino acids (user-
                defined parameter). If \code{trim_flank_aa} = 0, then 
                \code{cdr3} = \code{cdr3_core}
            \item \code{motif} - enriched motif from \code{cdr3_core}
            }
        }
        \item \code{global} - matrix, global clustering results. Pairs of
        globally similar CDR3s are shown in each row of the matrix
        (analogous to \code{lp})
    }
}
\item{inputs}{list, contains all user provided inputs (see \bold{Arguments})}
}

\description{
This algorithm finds groups of TCRs or BCRs with similar specificity. Two 
clustering strategies are employed:
\enumerate{
    \item Local clustering
    \item Global clustering
}

\bold{Local clustering}

\enumerate{
\item CDR3 processing steps
    \itemize{
        \item each row of s and r is considered as a CDR3 sequence from an 
        individual T- or B-cell (version = 2, default). If version=1 is
        specified, then we compute the set of non-redundant CDR3s from s and r
        and use them for clustering.
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
        \item probability \eqn{p_i} of finding the observed or a larger OvE for 
        motif \eqn{i} given that the null hypothesis is true is computed with 
        the Fisher's exact test
        \item classify motif \eqn{i} as \code{pass=TRUE} if the motif passes
        all filters specified in the user-provided control list, otherwise
        as \code{pass=FALSE}
    }
}


\bold{Global clustering}

The default ClustIRR algorithm for global clustering is simple. For each pair 
of equal-length CDR3 sequences \eqn{i} and \eqn{j} we compute the Hamming 
distance \eqn{d_{ij}}. If \eqn{d_{ij}\leq} \code{global_max_dist} (user-defined 
input), then \eqn{i} and \eqn{j} are globally similar.

Alternatively, the user can provide a matrix of globally similar CDR3 sequence 
pairs, computed by a complementary approachs such as TCRdist.
}

\examples{
# load package input data
data("CDR3ab")
s <- data.frame(CDR3b = CDR3ab[1:1000, "CDR3b"])
r <- data.frame(CDR3b = CDR3ab[1:5000, "CDR3b"])

# artificially enrich motif 'RQWW' inside sample dataset
base::substr(x = s$CDR3b[1:20], start = 6, stop = 9) <- "RQWW"

# add an artificial clonal expansion of two sequences to the sample dataset
s <- base::rbind(s, base::data.frame(CDR3b = rep(x = c("CATSRAAKPDGLRALETQYF",
                                    "CATSRAAKPDRQWWLSTQYF"),
                                    times = 15)))

# run analysis
out <- cluster_irr(s = s,
                   r = r,
                   version = 2,
                   ks = 4,
                   cores = 1,
                   control = list(
                        global_max_dist = 1,
                        local_max_fdr = 0.05,
                        local_min_ove = 2,
                        local_min_o = 1,
                        trim_flank_aa = 3,
                        global_pairs = NULL,
                        low_mem = FALSE))

# output class
base::class(out)

# output structure
utils::str(out)

# inspect motif enrichment results
knitr::kable(utils::head(slot(out, "clust")$CDR3b$local$m))

# inspect which CDR3bs are globally similar
knitr::kable(utils::head(slot(out, "clust")$CDR3b$global))

# plot graph
plot_graph(out)
}
