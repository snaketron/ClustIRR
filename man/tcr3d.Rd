\docType{data}
\name{TCR3d}
\alias{TCR3d}
\title{CDR3 sequences and their matching epitopes obtained from TCR3d}

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
\item \code{Reference}: Reference ID
}
}

\source{
\href{https://tcr3d.ibbr.umd.edu/}{TCR3d, June 2024}
}
\usage{
data(tcr3d)
}
\description{
data.frame with paired \code{CDR3a} and \code{CDR3b} CDR3 sequences and their 
matching epitopes obtained from TCR3d. The remaining CDR3 columns are set to 
NA. The antigenic epitopes come from cancer antigens and from viral antigens. 
For data processing details see the script inst/script/get_tcr3d.R
}
\examples{
data("tcr3d")
}
\value{
data(tcr3d) loads the object tcr3d
}
\keyword{datasets}
