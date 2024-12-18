\name{BLOSUM62}
\alias{BLOSUM62}
\docType{data}
\title{
BLOSUM62 matrix
}
\description{
Predefined scoring matrix for amino acid or nucleoitide alignments.
}
\usage{data("BLOSUM62")}
\format{
BLOSUM62 is a square symmetric matrix. Rows and columns are identical 
single letters, representing nucleotide or amino acid. Elements are 
integer coefficients (substitution scores).
}
\details{
BLOSUM62 was obtained from NCBI (the same matrix used by the stand-
alone BLAST software). 
}
\source{
See https://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62
}
\references{
See https://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62
}
\examples{
data(BLOSUM62, package = "ClustIRR")
BLOSUM62
}
\keyword{datasets}
