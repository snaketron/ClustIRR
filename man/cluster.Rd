\name{cluster_irr}
\alias{cluster_irr}
\title{Clustering of immune receptor repertoires}
\usage{
cluster_irr(
    s,
    r,
    ks = 4,
    cores = 1,
    control = list(global_smart = FALSE,
                   global_max_hdist = 1,
                   global_min_identity = 0.7,
                   local_max_fdr = 0.05,
                   local_min_o = 1,
                   trim_flank_aa = 3,
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
        \item \emph{CDR3h} and \emph{CDR3l} [for heavy/light chain BCRs]
    }
}
}
\item{r}{data.frame, reference (or control) repertoire of CDR3 sequences. Must
have the same structure (number of columns and column names) as \code{s}. If
this is not specified or set to NULL, then ClustIRR performs only global 
clustering using sample \code{s}}
\item{ks}{integer or integer vector, motif lengths. \code{ks} = 4 (default)}
\item{cores}{integer, number of CPU cores, \code{cores} = 1 (default).}
\item{control}{list, a named list of auxiliary parameters
to control algorithm's behavior. See the details below:
\itemize{
\item \code{global_smart} - logical, should we use smart global clustering 
based of BLOSUM62 scores (slower but more accurate; default) or less smart 
global clustering based on Hamming distances (faster but less accurate)
\item \code{global_min_identity} - probability, what is the minimum sequence
identity between a pair of CDR3 sequences for them to even be considered for
global similarity inspection (default = 0.7; 70 percent identity) This input
is only used if \code{global_smart} = TRUE.
\item \code{global_max_hdist} - integer, if \code{global_smart}=FALSE, then 
\code{global_max_hdist} defines a Hamming distance (HD) threshold, i.e. two 
CDR3s as globally similar if their Hamming distance is smaller or equal to 
\code{global_max_hdist} HD(\eqn{a}, \eqn{b}) \eqn{\leq} \code{global_max_hdist}.
\code{global_max_hdist} = 1 (default)
\item \code{local_max_fdr} - numeric, maximum False Discovery Rate (FDR) for
the detection of enriched motifs. \code{local_max_fdr} = 0.05 (default)
\item \code{local_min_o} - numeric, minimum absolute frequency of a motif in
the s in order for the motif to be used in the enrichment analysis.
\code{local_min_o} = 1 (default)
\item \code{trim_flank_aa} - integer, how many amino acids should be trimmed
from the flanks of all CDR3 sequences (only used for local clustering.
\code{trim_flank_aa} = 3 (default))
\item \code{low_mem} - logical, allows low memory mode for global clustering.
This will lead to increase in the CPU time but lead to a lower memory
footprint. \code{low_mem} = \code{FALSE} (default)
\item \code{global_pairs} - data.frame, pre-computed global pairs. If
\code{global_pairs} is provided by the user, then global clustering is not
performed. Instead the CDR3 pairs from \code{global_pairs} are used as global
clustering pairs. \code{global_pairs} is a data.frame matrix with 4 columns.
The first two columns, named, from_cdr3 and to_cdr3 contain pairs of CDR3 
sequences that are considered globally similar. The third column, called
weight, contains a similarity weight. If weights are not available they should 
be set to 1. The fourth column, called chain, contains the chain immune 
receptor in which the CDR3s are found: CDR3b or CDR3a [for \eqn{\alpha\beta} 
TCRs]; CDR3g or CDR3d [for \eqn{\gamma\delta} TCRs]; or CDR3h or CDR3l 
[for heavy/ligh chain BCRs].
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
        \item \code{global} - data.frame, global clustering results. Pairs of
        globally similar CDR3s are shown in each row (analogous to \code{lp}).
        If \code{global_smart}=FALSE in the control, then global clustering is 
        done based on Hamming distances and the remaining columns of this 
        data.frame are not important. Else, if \code{global_smart}=FALSE, then
        the remainig columns are relevant, i.e. global similarity scores are 
        shown for the complete CDR3 sequence pairs (column \code{weight}) or 
        their core (trimmed) CDR3 sequence part (column \code{cweight}). The 
        column \code{max_len} stores the the maximum length in each pair of 
        CDR3 sequences, and is used to normalize the scores \code{weight} and 
        \code{cweight}: the normalized scores are shown in the columns 
        \code{nweight} and \code{ncweight}.
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
The aim of local clustering is to find motifs (contiguous k-mers of the CDR 
sequence) that are overrepresented in repertoire \code{s} compared to 
repertoire \code{r}. This is an outline of the local clustering procedure:
\enumerate{
\item Trim CDR3 flanks based on \code{control$trim_flank_aa}
\item For each motif found in \code{s} compute the following:
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
        \item if a motif passes the criteria defined in control list, set flag 
        \eqn{\text{pass}_i=\text{T}}, else \eqn{\text{pass}_i=\text{F}}
    }
}

\bold{Global clustering} The aim of global clustering is to find similar
CDR3 sequences in repertoire \code{s}. This is an outline of the global 
clustering approaches implemented in ClustIRR:

If \code{global_smart}=FALSE

For each pair of equal-length CDR3 sequences \eqn{i} and \eqn{j} we compute 
the Hamming distance \eqn{d_{ij}}. If \eqn{d_{ij}\leq} \code{global_max_hdist} 
(user-defined input), then \eqn{i} and \eqn{j} are globally similar.

If \code{global_smart}=TRUE

We find pairs of CDR3 sequences that satisfy a minimum sequence identity 
defined by \code{global_min_identity}. Each pair of sequences are aligned.
Then, we compute the length of the longest of the two CDR3 sequences (column 
\code{max_len}). Then, we compute four alignment scores: a) \code{weight} -
BLOSUM62 score of the aligned CDR3 sequences; b) \code{nweight} - BLOSUM62 
score of the aligned CDR3 sequences normalized by \code{max_len}; c) 
\code{cweight} - BLOSUM62 score of the aligned CDR3 cores (untrimmed part of 
the CDR3 sequences); d) \code{ncweight} - BLOSUM62 score of the aligned CDR3 
cores (untrimmed part of the CDR3 sequences) normalized by the length of the 
longest core sequence.

Alternatively, the user can provide a matrix of globally similar CDR3 sequence 
pairs, computed by a complementary approachs such as TCRdist.
}

\examples{
# load package input data
data("CDR3ab")
s <- data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], clone_size = 1)
r <- data.frame(CDR3b = CDR3ab[1:500, "CDR3b"], clone_size = 1)

# artificially enrich motif 'RQWW' inside sample dataset
substr(x = s$CDR3b[1:20], start = 6, stop = 9) <- "RQWW"

# add an artificial clonal expansion of two sequences to the sample dataset
s <- rbind(s, data.frame(CDR3b = c("CATSRAAKPDGLRALETQYF",
                                   "CATSRAAKPDRQWWLSTQYF"),
                         clone_size = 10))

# run analysis
out <- cluster_irr(s = s,
                   r = r,
                   ks = 4,
                   cores = 1,
                   control = list(
                       global_smart = TRUE,
                       global_max_hdist = 1,
                       local_max_fdr = 0.05,
                       local_min_o = 1,
                       trim_flank_aa = 3,
                       global_pairs = NULL,
                       low_mem = FALSE))

# output class
class(out)

# output structure
str(out)

# inspect motif enrichment results
knitr::kable(head(slot(out, "clust")$CDR3b$local$m))

# inspect which CDR3bs are globally similar
knitr::kable(head(slot(out, "clust")$CDR3b$global))

# get graph
g <- get_graph(out)

# plot graph
plot_graph(g)

# plot graph as visgraph
plot_graph(g, as_visnet = TRUE)
}
