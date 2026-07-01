\docType{data}
\name{D_CMV_FLU}
\alias{D_CMV_FLU}
\title{Processed single-cell TCR and gene expression dataset of antigen 
stimulated T cells}

\description{
\code{D_CMV_FLU} is a processed single-cell TCR and gene expression
dataset derived from a publicly available Parse Biosciences Evercode TCR
v4 dataset. It combines two antigen-stimulated human T-cell repertoires:

\itemize{
\item \code{CMV_A0201}: HLA-A0201+ donor T cells stimulated with
cytomegalovirus (CMV) antigens.
\item \code{Flu}: donor T cells stimulated with influenza M1 (Flu)
antigens.
}

The \code{CMV_A0201} and \code{Flu} repertoires contain 154 and 1,540
clonotypes, respectively. The dataset is stored as a list with three
data frames: \code{s}, \code{meta}, and \code{gex}.

The data frame \code{s} is the primary input to \code{clustirr} and
contains one row per clonotype with paired CDR3\eqn{\alpha\beta} amino
acid sequences, repertoire identifier (sample name), and clonotype
abundance.

The data frame \code{meta} contains clonotype metadata, including the
Parse Biosciences clonotype identifier and a comma-separated list of
cell barcodes assigned to each clonotype.

The data frame \code{gex} contains one row per cell with log-normalized
expression values for 12 T-cell activation genes and a
\code{cell_id} column linking each cell to the corresponding clonotype
through the cell barcodes stored in \code{meta}.
}

\format{
A list with three data.frames:

\itemize{

\item \code{s}: the primary \code{clustirr} input, containing one row
per clonotype with the following columns:
\itemize{
\item \code{CDR3a}: CDR3\eqn{\alpha} amino acid sequence.
\item \code{CDR3b}: CDR3\eqn{\beta} amino acid sequence.
\item \code{sample}: antigen-stimulation condition (\code{CMV_A0201} or 
\code{Flu}).
\item \code{clone_size}: clonotype abundance, represented by the number
of cells assigned to the clonotype.
}

\item \code{meta}: processed clonotype metadata, including Parse
Biosciences clonotype identifiers and a comma-separated list of cell
barcodes for each clonotype.

\item \code{gex}: one row per cell containing log-normalized expression
values for 12 T-cell activation genes together with the column
\code{cell_id}, which links each cell to its clonotype.
}
}

\source{
Performance of Evercode™ TCR v4 in Antigen-Stimulated Human T Cells, 
https://www.parsebiosciences.com/datasets/performance-of-evercode-tcr-v4-in-
antigen-stimulated-human-t-cells/; Parse Biosciences; Accessed 2026 June 30
}



\usage{
data(D_CMV_FLU)
str(D_CMV_FLU)
}

\examples{
data(D_CMV_FLU)
str(D_CMV_FLU)
}

\value{
Dataset \code{D_CMV_FLU}.
}

\keyword{datasets}