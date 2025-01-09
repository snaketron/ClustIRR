\docType{data}
\name{CDR3ab}
\alias{CDR3ab}
\title{Mock data set of complementarity determining region 3 (CDR3)
sequences and variable (V) and joining (J) genes from the \eqn{\alpha} 
and \eqn{\beta} chains of 10,000 T cell receptors.}

\format{
\code{data.frame} with 10,000 rows and 6 columns 
\itemize{
\item \code{CDR3a}: CDR\eqn{\alpha} amino acid sequence
\item \code{TRAV}: variable (V) gene of TCR\eqn{\alpha}
\item \code{TRAV}: joining (J) gene of TCR\eqn{\alpha}
\item \code{CDR3b}: CDR\eqn{\beta} amino acid sequence
\item \code{TRBV}: variable (V) gene of TCR\eqn{\beta}
\item \code{TRBV}: joining (J) gene of TCR\eqn{\beta}
}
}
\source{
\href{http://50.255.35.37:8080/}{GLIPH version 2}
}
\usage{
data(CDR3ab)
}
\description{
T-cell receptor (TCR) repertoire with 10,000 T-cells (rows). Each T-cell 
has the following features: amino acid sequences of their complementarity 
determining region 3 (CDR3); and variable (V) and joining (J) gene names 
for TCR chains \eqn{\alpha} and \eqn{\beta}. 

Important remark: this is a mock dataset, all CDR3 sequences and the genes were 
sampled from a larger set of CDR3\eqn{\beta} sequences and genes of naive 
CD8+ T cells in humans.
}
\examples{
data("CDR3ab")
}
\value{
data(CDR3ab) loads the object CDR3ab, which is a data.frame with six columns
(3 for TCR\eqn{\alpha} and 3 for TCR\eqn{\beta}) and 10,000 rows (see details).
}
\keyword{datasets}
