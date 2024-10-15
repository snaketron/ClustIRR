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
            control = list(global_min_identity = 0.7,
                           trim_flank_aa = 3,
                           global_hamming = FALSE,
                           global_max_hdist = 1,
                           low_mem = FALSE))
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
\item \code{global_min_identity} - probability, what is the minimum sequence
identity between a pair of CDR3 sequences for them to even be considered for
global similarity inspection (default = 0.7; 70 percent identity) This input
is only used if \code{global_hamming} = TRUE.
\item \code{trim_flank_aa} - integer, how many amino acids should be trimmed
from the flanks of all CDR3 sequences. \code{trim_flank_aa} = 3 (default)
\item \code{global_hamming} - logical (default = FALSE), should we use global 
clustering basedbased of Hamming distances (faster but less reliable).
\item \code{global_max_hdist} - integer, if \code{global_hamming}=TRUE, then 
\code{global_max_hdist} defines a Hamming distance (HD) threshold, i.e. two 
CDR3s as globally similar if their Hamming distance is smaller or equal to 
\code{global_max_hdist} HD(\eqn{a}, \eqn{b}) \eqn{\leq} \code{global_max_hdist}.
\code{global_max_hdist} = 1 (default)
\item \code{low_mem} - logical, allows low memory mode for global clustering.
This will lead to increase in the CPU time but lead to a lower memory
footprint. \code{low_mem} = \code{FALSE} (default)
}
}
}
\value{
The output is an \code{S4} object of class \code{clust_irr}. This object 
contains two sublists:
\itemize{
\item{\code{clust}, list, contains clustering results for each IR chain. 
The results are stored as data.frame in separate sub-list named appropriately 
(e.g. CDR3a, CDR3b, CDR3g, etc.). 

Each row in the data.frames contains a pair of globally similar CDR3s.
If \code{global_hamming}=TRUE in the control, then global clustering is 
done based on Hamming distances, then the remaining columns of this 
data.frame are not important. 

Otherwise, the remainig columns are relevant, i.e. global similarity scores 
are shown for the complete CDR3 sequences (column \code{weight}) or their 
cores (column \code{cweight}). The column \code{max_len} stores the the maximum 
length in each pair of CDR3 sequences, and is used to normalize the scores 
\code{weight} and \code{cweight}: the normalized scores are shown in the 
columns \code{nweight} and \code{ncweight}}

\item{\code{inputs}, list, contains all user provided inputs (see \bold{Arguments})}
}
}
\details{
We assume that two CDR3 having identical sequences will recognize the same 
sets of antigens. Meanwhile, we also know that similar CDR3 sequences (e.g., 
CDR3 sequences which differ by only one amino acid) may also have common 
specificity. `r Biocpkg("ClustIRR")` employs clustering to find groups of 
CDR3s which are **globally similar** (based on the complete CDR3 sequences). 

For this we employs one of two strategies.

\bold{Strategy A: if control$global_hamming = FALSE (default)}

Assume we have \eqn{n} CDR3 sequences in a given IRR, \code{cluster_irr} 
performs \eqn{n^2} pairwise Needleman-Wunsch alignments. For a pair \eqn{i} 
of CDR3 sequences, it computes a BLOSUM62 similarity score \eqn{b_i}. 

To ensure that the similarity scores are comparable between pairs of CDR3s 
with different lengths, we divide \eqn{b_i} with \eqn{l_i}, the maximum 
length of the longer CDR3 sequence in \eqn{i}: 

\eqn{b'_i = \dfrac{b_i}{l_i} = \dfrac{b_i}{\max(m_i, k_i)}},

where \eqn{m_i} and \eqn{k_i} are the lengths of the two CDR3 sequences in 
pair \eqn{i}; and \eqn{b'_i} is the normalized similarity score. We compute 
\eqn{b_i} and \eqn{b'_i} for two parts of the CDR3 sequence:
\itemize{
\item{\bold{complete CDR3}: \eqn{b_i}=\code{weight} and 
\eqn{b'_i}=\code{nweight}}
\item{\bold{CDR3 core}: \eqn{b_i}=\code{cweight} and 
\eqn{b'_i}=\code{ncweight}}
}

\bold{This strategy is computationally very expensive!}

The Needleman-Wunsch algorithm has a time complexity \eqn{O(mk)}. Empirically 
we know that \eqn{m \approx k}, and the typical values of \eqn{m} and \eqn{k} 
range between 6 and 20 (mean = 13). This operation has to be performed 
\eqn{n^2} times. Therefore, the time complexity of this procedure is 
approximately \eqn{O(m^2n^2)}.

For large IRRs with \eqn{n>10^6}, which are routinely generated nowadays, 
this algorithm requires significant computational resources. To mitigate 
this challenge, we only consider pairs of CDR3 sequences with sequence 
identity \eqn{\geq} 70\%. This threshold may be modified by the user (by 
modifying the control parameter \code{global_min_identity}). For the 
remaining pairs we consider \eqn{b'_i=0}.

\bold{Strategy B: if control$global_hamming = TRUE}

Alternatively, \code{clust_irr} quantifies the dissimilarity between 
pairs of CDR3 sequences using Hamming distances. Two CDR3 sequences with 
Hamming distance \eqn{\leq x} are considered globally similar, where \eqn{x} 
is an user-defined threshold (input control$global_max_hdist, default = 1). 

This is a computationally and conceptually simple strategy, which also has 
drawbacks. First, CDR3 sequences with different lengths are by definition 
dissimilar. Second, the properties of the replaced amino acids are not 
considered by the Hamming distance.

In summary, for each pair of equal-length CDR3 sequences \eqn{i} and \eqn{j} 
we compute the Hamming distance \eqn{d_{ij}}. If \eqn{d_{ij}\leq} 
\code{global_max_hdist} (user-defined input), then \eqn{i} and \eqn{j} are 
globally similar.
}

\examples{
# load package input data
data("CDR3ab")
s <- data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], sample = "A", clone_size = 1)

# run analysis
c <- cluster_irr(s = s)

# output class
class(c)

# output structure
str(c)

# inspect which CDR3bs are globally similar
knitr::kable(head(slot(c, "clust")$CDR3b))
}
