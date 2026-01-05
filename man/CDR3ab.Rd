\docType{data}
\name{Datasets}
\alias{CDR3ab}
\alias{D1}
\alias{D2}
\title{Datasets \code{CDR3ab}, \code{D1} and \code{D2} with TCR\eqn{\alpha\beta} 
mock repertoires}

\format{
CDR3ab is a \code{data.frame} with rows as TCR clones and 6 columns 
\itemize{
\item \code{CDR3a}: CDR3\eqn{\alpha} amino acid sequence
\item \code{TRAV}: variable (V) gene of TCR\eqn{\alpha}
\item \code{TRAV}: joining (J) gene of TCR\eqn{\alpha}
\item \code{CDR3b}: CDR3\eqn{\beta} amino acid sequence
\item \code{TRBV}: variable (V) gene of TCR\eqn{\beta}
\item \code{TRBV}: joining (J) gene of TCR\eqn{\beta}
}
}
\source{
\href{http://50.255.35.37:8080/}{GLIPH version 2}
}
\usage{
# For the raw data with 10,000 TCR clones
data(CDR3ab)

# For dataset D1
data(D1)

# For dataset D2
data(D2)
}
\description{
TCR\eqn{\alpha\beta} repertoire with 10,000 T-cells (rows). Each T-cell 
has the following features: amino acid sequences of their complementarity 
determining region 3 (CDR3); and variable (V) and joining (J) gene names 
for TCR chains \eqn{\alpha} and \eqn{\beta}. 

Important remark: this is a mock dataset, all CDR3 sequences and the genes were 
sampled from a larger set of CDR3\eqn{\beta} sequences and genes of naive 
CD8+ T cells in humans.

About dataset \code{D1}:

We used this data to create dataset D1: three TCR\eqn{\alpha\beta} 
repertoires \code{a}, \code{b}, and \code{c}, each with 500 TCR 
clones. We simulated clonal expansion with increasing degree in TCR 
repertoires \code{b} and \code{c}. The TCR repertoires are stored as 
element of a list. For each TCR repertoires we have a metadata: 
\code{ma}, \code{mb}, and \code{mc}.

About dataset \code{D2}:

We used \code{b} and \code{c} from D1 to create dataset D2. Three 
TCR\eqn{\alpha\beta} repertoires (replicates) were generated from 
two biological conditions: A (A1, A2, A3) and C (C1, C2, C3), each 
with 2,000 T-cells.
}
\examples{
data("CDR3ab")
data("D1")
data("D2")
}
\value{
data(CDR3ab) loads the object CDR3ab, which is a data.frame with six columns
(3 for TCR\eqn{\alpha} and 3 for TCR\eqn{\beta}) and rows for each TCR clone
(see details).
}
\keyword{datasets}
