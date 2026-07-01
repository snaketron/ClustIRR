\docType{data}
\name{D1}
\alias{D1}
\title{Synthetic dataset \code{D1} with three TCR\eqn{\alpha\beta} repertoires}

\description{
\code{D1} contains three synthetic TCR\eqn{\alpha\beta} repertoires
(\code{a}, \code{b}, and \code{c}), each comprising 500 clonotypes.
The dataset is stored as a list with two data.frames, \code{s} and
\code{meta}.

The data.frame \code{s} contains a clonotype table with 1,500 clonotypes
(500 from each repertoire). Each clonotype is represented by the paired
amino acid sequences of the complementarity-determining region 3
(CDR3\eqn{\alpha\beta}), the repertoire identifier (column
\code{sample}), and the clonotype expansion (column
\code{clone_size}).

The data.frame \code{meta} contains additional metadata for each
clonotype, including the clonotype identifier, variable gene usage,
cell subset (e.g., CD4+ or CD8+), patient HLA type, and age.
}

\format{
A list with two data.frames:

\itemize{
  \item \code{s}: the primary \code{clustirr} input, containing one row
  per clonotype with the following columns:
  \itemize{
    \item \code{CDR3a}: CDR3\eqn{\alpha} amino acid sequence.
    \item \code{CDR3b}: CDR3\eqn{\beta} amino acid sequence.
    \item \code{sample}: repertoire identifier.
    \item \code{clone_size}: clonotype abundance, represented as the
    number of cells assigned to the clonotype.
  }

  \item \code{meta}: metadata for each clonotype, including clonotype
  identifiers, V gene usage, cell subset, HLA type, and age.
}
}

\source{
https://zenodo.org/records/18925758
}

\usage{
data(D1)
str(D1)
}

\examples{
data(D1)
}

\value{
Loads the dataset \code{D1}.
}

\keyword{datasets}