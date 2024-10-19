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
\item \code{db_dist}: we compute edit distances between CDR3 sequences from 
the sample and from a database (e.g. VDJdb). If a particular distance is 
smaller than or equal to \code{edit_dist}, then we annotate the sample CDR3
with the specificity of the database CDR3 sequence.
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

Pair of sequences, \eqn{a} and \eqn{b}, are aligned with the Needleman-
Wunsch algorithm. The output is an alignment score (\eqn{\omega}). 
Identical or similar CDR3 sequence pairs get a large positive \eqn{\omega}, 
and dissimilar CDR3 sequence pairs get a low (or even negative) \eqn{\omega}.

To make sure that \eqn{\omega} is comparable across pairs of CDR3s with
different lengths, ClustIRR divides (normalizes) \eqn{\omega}
by the length of the longest CDR3 sequences in each pair: 

\eqn{\bar{\omega} = \dfrac{\omega}{\max(|a|, |b|)}}

where \eqn{|a|} and \eqn{|b|} are the lengths of CDR3 sequences \eqn{a}
and \eqn{b}; and \eqn{\bar{\omega}} is the normalized alignment score.

The CDR3 \bold{cores}, which represent the central parts of the CDR3 loop 
and tend to have high probability of making a contact with the antigen, 
are also compared. ClustIRR constructs the CDR3 cores by trimming few 
residues (defined by \code{control$trim_flanks}) from either end of each 
CDR3 sequences. These are then aligned and scored based on the same 
algorithm, yielding for each pair of CDR3 cores a normalized alignment 
scores \eqn{\bar{\omega}_c}.

To compute the similarity score between TCR clones, \eqn{i} and \eqn{j}, we
compute the average alignment score from their CDR3\eqn{\alpha} and
CDR3\eqn{\beta} alignment scores (\bold{notice:} in the next, I will use
TCR\eqn{\alpha\beta} as an example. However, this approach can also be used
to compare TCR\eqn{\gamma\delta} or BCR\emph{IgH-IgL} clones): 


\eqn{\omega(i,j)=\dfrac{{\bar{\omega}}^\alpha + {\bar{\omega}}^\beta}{2}
\qquad\text{or}\qquad
\omega(i,j)=\dfrac{{\bar{\omega}}^\alpha_c + {\bar{\omega}}^\beta_c}{2}},


where \eqn{\bar{\omega}^\alpha} and \eqn{\bar{\omega}^\beta} are the
alignment scores for the CDR3\eqn{\alpha} and CDR3\eqn{\beta} sequences,
respectively; and \eqn{\bar{\omega}^\alpha_c} and \eqn{\bar{\omega}^\beta_c} 
are the alignment scores for the CDR3\eqn{\alpha} and CDR3\eqn{\beta} cores,
respectively. Based on this metric, the contributions of CDR3\eqn{\alpha}
and CDR3\eqn{\beta} towards the overall similarity of the TCR clones are
assigned equal weights.

\bold{Important note}: If the data contains CDR3 sequences from only one
chain, such as CDR3\eqn{\beta}, then \eqn{\omega(i,j)} is defined as

\eqn{\omega(i,j)={\bar{\omega}}^\beta \qquad\text{or}\qquad 
\omega(i,j)={\bar{\omega}}^\beta_c}

ClustIRR provides two additional metrics for computing
similarity scores between TCR clones, including a \emph{strict metric}, 
which assigns high similarity score to a pair of TCR clones only if both of
their CDR3\eqn{\alpha} and CDR3\eqn{\beta} sequence pairs are similar

\eqn{\omega^s(i,j) = \min({\bar{\omega}}^\alpha, {\bar{\omega}}^\beta)
    \qquad\text{or}\qquad
    \omega^s(i,j) = \min({\bar{\omega}}^\alpha_c, {\bar{\omega}}^\beta_c)},

and a \emph{loose metric}, which assigns high similarity score to
a pair of TCR clones if either of their CDR3\eqn{\alpha} and CDR3\eqn{\beta}
sequence pairs are similar 

\eqn{\omega^l(i,j) = \max({\bar{\omega}}^\alpha, {\bar{\omega}}^\beta)
    \qquad\text{or}\qquad
    \omega^l(i,j) = \max({\bar{\omega}}^\alpha_c, {\bar{\omega}}^\beta_c)},

\bold{This strategy is computationally very expensive!}

For large IRRs with \eqn{n>10^6} this algorithm requires significant
computational resources. To mitigate this challenge, we employ a
screening step in which dissimilar sequences pairs are flagged. In
short, each CDR3 is used as a query in a \bold{fast} protein-BLAST search
as implemented in the R-package blaster, while the remaining CDR3s are
considered as a database of amino acid sequences against which the query
is compared. CDR3 sequences which share at least 70\% sequence identity
(user parameter \code{control$gmi}) with the query are selected, and only
these are aligned with query CDR3. For the remaining CDR3 pairs we
assume \eqn{\omega(i,j)=0}.

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
