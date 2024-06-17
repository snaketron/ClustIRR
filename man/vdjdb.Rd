\docType{data}
\name{vdjdb}
\alias{vdjdb}
\title{CDR3 sequences and their matching epitopes obtained from VDJdb}

\format{
\code{data.frame} with columns: 
\enumerate{
\item \code{CDR3a}: CDR3a amino acid sequence
\item \code{CDR3b}: CDR3b amino acid sequence
\item \code{CDR3g}: CDR3g amino acid sequence -> NA
\item \code{CDR3d}: CDR3d amino acid sequence -> NA
\item \code{CDR3h}: CDR3h amino acid sequence -> NA
\item \code{CDR3l}: CDR3l amino acid sequence -> NA
\item \code{CDR3_species}: CDR3 species (e.g. human, mouse, ...)
\item \code{Antigen_species}: antigen species
\item \code{Antigen_gene}: antigen gene
\item \code{Reference}: Reference (Pubmed ID)
}
}

\source{
\href{https://vdjdb.cdr3.net/}{VDJdb, June 2024}
}
\usage{
data(vdjdb)
}
\description{
data.frame with unpaired \code{CDR3a} or \code{CDR3b} sequences and their
matching epitopes obtained from VDJdb. The remaining CDR3 columns are set 
to NA. For data processing details see the script inst/script/get_vdjdb.R}
\examples{
data("vdjdb")
}
\value{
data(vdjdb) loads the object vdjdb
}
\keyword{datasets}
