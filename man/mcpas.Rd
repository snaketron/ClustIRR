\docType{data}
\name{McPAS-TCR}
\alias{McPAS-TCR}
\title{CDR3 sequences and their matching epitopes obtained from McPAS-TCR}

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
\href{http://friedmanlab.weizmann.ac.il/McPAS-TCR/}{McPAS-TCR, June 2024}
}
\usage{
data(mcpas)
}
\description{
data.frame with \code{CDR3a} and/or \code{CDR3b} sequences and their 
matching antigenic epitopes obtained from McPAS-TCR. The remaining CDR3 
columns are set to NA. For data processing details see the script 
inst/script/get_mcpastcr.R
}
\examples{
data(mcpas)
}
\value{
data(mcpas) loads the object McPAS-TCR
}
\keyword{datasets}
