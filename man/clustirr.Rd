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
name of the IR chain (e.g. CDR3a for CDR3s from TCR\eqn{\alpha} chain; or 
CDR3d for CDR3s from TCR\eqn{\delta} chain}. Meanwhile, if paired CDR3s from 
both chains are available, then you can provide both in separate columns e.g.:
\itemize{
    \item \emph{CDR3b} and \emph{CDR3a} [for \eqn{\alpha\beta} TCRs]
    \item \emph{CDR3g} and \emph{CDR3d} [for \eqn{\gamma\delta} TCRs]
    \item \emph{CDR3h} and \emph{CDR3l} [for heavy/light chain BCRs]
}
}
}
\item{meta}{data.frame with meta-data for each IR clone, which may contain
data such as, V/J genes, biological condition, age, etc. This data will be
used to annotate the graph nodes and help downstream analyses.}
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
are distributed over T-cell clones. TCR clones with \bold{identical} 
pairs of CDR3\eqn{\alpha} and CDR3\eqn{\beta} sequences most likely 
recognize the same sets of antigens. Meanwhile, TCR clones with 
\bold{similar} pairs of CDR3\eqn{\alpha} and CDR3\eqn{\beta} sequences 
may also share common specificity. ClustIRR aims to quantify the 
similarity between pairs of TCR clones based on the similarities of their 
CDR3s sequences.

How to compute a similarity score between a pair of CDR3 sequences?

1. Align pairs of sequences

2. Score alignment with BLOSUM62 substitution matrix and gap open/exten costs

3. Normalize alignment scores by alignment length

4. Compute the normalized alignment score of the CDR3 \bold{cores}.

CDR3 \bold{cores} are the central parts of the CDR3 loop, which tend to have 
high probability of making a contact with the antigen. ClustIRR constructs the 
CDR3 cores by trimming few residues (defined by \code{control$trim_flanks}) 
from both ends of each CDR3 sequence.

For large IRRs with many clones, step 1 requires significant computational 
resources. To mitigate this challenge, we employ a screening step in which 
similar sequence pairs selected. In short, each CDR3 is used as a query in 
a \bold{fast} protein-BLAST search as implemented in the R-package blaster, 
while the remaining CDR3s are considered as a database of amino acid sequences 
against which the query is compared. CDR3 sequences which share at least 70\% 
sequence identity (user parameter \code{control$gmi}) with the query are 
selected, and only these are scored according to steps 2-4.
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
