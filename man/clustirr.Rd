\name{cluster_irr}
\alias{cluster_irr}
\title{Clustering of immune receptor repertoires (IRRs)}
\description{
\code{cluster_irr} computes similarities between immune receptors (IRs = 
T-cell and B-cell receptors) based on their CDR3 sequences.
}
\usage{
cluster_irr(s,
            meta = NULL,
            control = list(gmi = 0.7,
                           trim_flank_aa = 3,
                           db_dist = 0,
                           db_custom = NULL))
}
\arguments{
\item{s}{a data.frame with complementarity determining region 3 (CDR3) amino 
acid sequences observed in IRR clones (data.frame rows). The data.frame has 
the following columns (IR clone features):
\itemize{
\item{\code{sample}: name of the IRR (e.g. 'A')}
\item{\code{clone_size}: cell count in the clone (=clonal expansion)}
\item{\code{CDR3?}: amino acid CDR3 sequence. Replace '?' with the appropriate 
name of the immmune receptor chain (e.g. CDR3a for CDR3s from TCR\eqn{\alpha} 
chain; or CDR3d for CDR3s from TCR\eqn{\delta} chain}. Meanwhile, if paired 
CDR3s from both chains are available, then you can provide both in separate 
columns e.g.:
\itemize{
    \item \emph{CDR3b} and \emph{CDR3a} [for \eqn{\alpha\beta} TCRs]
    \item \emph{CDR3g} and \emph{CDR3d} [for \eqn{\gamma\delta} TCRs]
    \item \emph{CDR3h} and \emph{CDR3l} [for heavy/light chain BCRs]
}
}
}
\item{meta}{data.frame with meta-data for each clone, which may contain
clone-specific data, such as, V/J genes, cell-type (e.g. CD8+, CD4+), nut
also repertoire-specific data, such as, biological condition, HLA type, 
age, etc. This data will be used to annotate the graph nodes and help 
downstream analyses.}
\item{control}{auxiliary parameters to control the algorithm's behavior. See 
the details below:
\itemize{
\item \code{gmi}: the minimum sequence identity between a pair of CDR3 
sequences for them to even be considered for alignment and scoring 
(default = 0.7; 70 percent identity).
\item \code{trim_flank_aa}: how many amino acids should be trimmed from 
the flanks of all CDR3 sequences to isolate the \bold{CDR3 cores}. 
\code{trim_flank_aa} = 3 (default).
\item \code{db_custom}: additional database (data.frame) which allows us 
to annotate CDR3 sequences from the input (\code{s}) with their cognate 
antigens. The structure of \code{db_custom} must be identical to that in 
\code{data(vdjdb, package = "ClustIRR")}. ClustIRR will use the internal 
databases if \code{db_custom}=NULL (default). Three databases (\bold{data 
only from human CDR3}) are integrated in ClustIRR: \code{VDJdb}, \code{TCR3d}
and \code{McPAS-TCR}. 
\item \code{db_dist}: we compute edit distances between CDR3 sequences from 
\code{s} and from a database (e.g. VDJdb). If a particular distance is 
smaller than or equal to \code{db_dist} (default = 0), then we annotate 
the CDR3 from \code{s} with the specificity of the database CDR3 sequence.
}
}
}
\value{
The output is an \code{S4} object of class \code{clust_irr}. This object 
contains two sublists:
\itemize{
\item{\code{clust}, list, contains clustering results for each receptor chain. 
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
ClustIRR performs the following steps.
\enumerate{ 
\item Compute similarities between clones within each repertoire 
\eqn{\rightarrow} \strong{the function \code{cluster_irr} performs this step}
\item Construct a graph from each TCR repertoire
\item Construct a joint similarity graph (\eqn{J})
\item Detect communities in \eqn{J}
\item Analyze Differential Community Occupancy (DCO)
\itemize{
    \item Between individual TCR repertoires with model \eqn{M}
    \item Between groups of TCR repertoires from biological conditions
    with model \eqn{M_h}
}
\item Inspect results
}
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
