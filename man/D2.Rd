\docType{data}
\name{D2}
\alias{D2}
\title{Synthetic dataset \code{D2} with six TCR\eqn{\alpha\beta} repertoires}

\description{
\code{D2} contains six synthetic TCR\eqn{\alpha\beta} repertoires
generated from two biological conditions, \code{A} and \code{C}. Each
condition comprises three replicate repertoires (\code{A1}--\code{A3}
and \code{C1}--\code{C3}), each containing 2,000 T cells distributed across 
approximately 450 clonotypes.

The dataset is stored as a single data.frame, which serves as the
primary input for \code{clustirr}. Each row represents a clonotype and
contains the paired amino acid sequences of the
complementarity-determining region 3 (CDR3\eqn{\alpha\beta}), the
repertoire identifier (column \code{sample}), and the clonotype
abundance (column \code{clone_size}).
}

\format{
A data.frame with one row per clonotype and the following columns:

\itemize{
  \item \code{CDR3a}: CDR3\eqn{\alpha} amino acid sequence.
  \item \code{CDR3b}: CDR3\eqn{\beta} amino acid sequence.
  \item \code{sample}: repertoire identifier (\code{A1}--\code{A3} or
  \code{C1}--\code{C3}).
  \item \code{clone_size}: clonotype abundance, represented as the
  number of cells assigned to the clonotype.
}
}

\source{
Simulated dataset (no source).
}

\usage{
data(D2)
str(D2)
}

\examples{
data(D2)
str(D2)
aggregate(clone_size~sample, data = D2, FUN = sum)
table(D2$sample)
}

\value{
Loads the dataset \code{D2}.
}

\keyword{datasets}