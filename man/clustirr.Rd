\name{cluster_irr}
\alias{cluster_irr}
\title{Clustering of immune receptor repertoires (IRRs)}
\description{
This algorithm performs clustering of CDR3 sequences to find groups of immune
receptors (IRs = T-cell and B-cell receptors) with similar specificity. The 
output are pairs of similar CDR3 sequences and their similarity scores. 
}
\usage{
cluster_irr(s,
            control = list(gmi = 0.7,
                           trim_flank_aa = 3,
                           db_dist = 0,
                           db_custom = NULL))
}
\arguments{
\item{s}{a data.frame with complementarity determining region 3 (CDR3) amino 
acid sequences observed in IRR clones. Each rows in the data.frame represents 
a clone. The data.frame has the following columns:
\itemize{
\item{\code{sample}: name of the IRR (e.g. 'A')}
\item{\code{clone_size}: frequency of cells in the clone (=clonal expansion)}
\item{\code{CDR3?}: amino acid CDR3 sequence. Replace '?' with the appropriate 
name of the chain (e.g. CDR3a for CDR3s from TCR\eqn{\alpha} chain; or CDR3d 
for CDR3s from TCR\eqn{\delta} chain}. Meanwhile, if paired CDR3s from both 
chains are available, then you can provide both in separate columns. In the 
following see a list of possible chains:.
\itemize{
    \item \emph{CDR3b} and \emph{CDR3a} [for \eqn{\alpha\beta} TCRs]
    \item \emph{CDR3g} and \emph{CDR3d} [for \eqn{\gamma\delta} TCRs]
    \item \emph{CDR3h} and \emph{CDR3l} [for heavy/light chain BCRs]
}
}
}
\item{control}{a named list of auxiliary parameters to control the 
algorithm's behavior. See the details below:
\itemize{
\item \code{gmi}: probability, what is the global minimum sequence identity 
(gmi) between a pair of CDR3 sequences for them to even be considered for
alignment and scoring (default = 0.7; 70 percent identity).
\item \code{trim_flank_aa}: integer, how many amino acids should be trimmed
from the flanks of all CDR3 sequences. \code{trim_flank_aa} = 3 (default).
\item \code{db_custom}: additional database (data.frame) that allows us to 
map CDR3 sequences from the input (\code{s}) to their cognate antigens. The 
structure of \code{db_custom} must be identical to that in data(vdjdb). 
ClustIRR will use the internal databases if \code{db_custom}=NULL (default).
Three databases (\bold{data only from human CDR3}) are integrated in ClustIRR: 
\code{VDJdb}, \code{tcr3d} and \code{mcpas-tcr}. 
\item \code{db_dist}: CDR3 sequences having edit distance \eqm{\leq} 
\code{edit_dist} compared to database CDR3 sequences are annotated with the 
known specificities of the database CDR3 sequences.
}
}
}
\value{
The output is an \code{S4} object of class \code{clust_irr}. This object 
contains two sublists:
\itemize{
\item{\code{clust}, list, contains clustering results for each IR chain. 
The results are stored as data.frame in separate sub-list named appropriately 
(e.g. CDR3a, CDR3b, CDR3g, etc.). Each row in the data.frames contains a pair 
of CDR3s.

The remaining columns contain similarity scores for the complete CDR3 sequences 
(column \code{weight}) or their cores (column \code{cweight}). The columns 
\code{max_len} and \code{max_clen} store the length of the longer CDR3 sequence
and core in the pair, and these used to normalize the scores \code{weight} and 
\code{cweight}: the normalized scores are shown in the columns \code{nweight} 
and \code{ncweight}}

\item{\code{inputs}, list, contains all user provided inputs (see Arguments)}
}
}
\details{
IRRs, such as T-cell receptor repertoires, are made up of T-cells which
are distributed over T-cell clones. TCR clones with \bold{identical} pairs
of e.g. CDR3\eqn{\alpha} and CDR3$\eqn{\beta} sequences most likely recognize 
the same sets of antigens. Meanwhile, TCR clones with \bold{similar} pairs 
of CDR3\eqn{\alpha} and CDR3$\eqn{\beta} sequences may also share common
specificity. ClustIRR aims to quantify the similarity between pairs of TCR 
clones based on the similarities of their CDR3s sequences.

TODO: copy from vignette

}

\examples{
# load package input data
data("CDR3ab", package = "ClustIRR")
s <- data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], sample = "A", clone_size = 1)

# run analysis
c <- cluster_irr(s = s)

# output class
class(c)

# output structure
str(c)

# inspect which CDR3bs are similar
knitr::kable(head(slot(c, "clust")$CDR3b))
}
