\name{get_cdr3_cosine_similarity}
\alias{get_cdr3_cosine_similarity}
\title{Compute and Visualize Cosine Similarity (CS)}
\description{
This function calculates pairwise CSs between repertoires on CDR3 sequence level
and generates a labeled heatmap visualization. If two chains are compared,
only perfect matches between both sequences will be considered as overlap.
Cell counts of identical sequences within a repertoire will be aggregated before
comparison. If only one chain is supplied, but two chains exist in the clust_irr
input, cell counts of identical sequences of this chain within the repertoire 
will be aggregated. Here CS measures the cosine of the angle between vectors of
cell counts (positive integers or zero), providing a value between 0 (perfect 
dissimilarity) and 1 (perfect similarity).
}
\details{
The function performs three main operations:
\enumerate{
  \item Computes CS between CDR3 sequences between repertoires using: 
  \deqn{\text{CS}(x, y) = \frac{x \cdot y}{||x|| \cdot ||y||}}
  \item Creates a heatmap with CS values annotated in each cell
}
Red and blue tiles in the heatmap indicates high and low CS, respectively.
}
\usage{
get_cdr3_cosine_similarity(clust_irr,
                           chains)
}
\arguments{
\item{clust_irr}{An S4 object of class clust_irr.}
\item{chains}{Which chains to compare. Can be one or both of the existing pairs
in the clust_irr object. If not set, all detected chains will be compared.}
}
\value{
A list containing two elements:
\itemize{
  \item \code{g}: A \code{ggplot2} heatmap displaying the CS
  \item \code{cs}: A data frame in long format with columns:
    \itemize{
      \item \code{i}: First repertoire
      \item \code{j}: Second repertoire
      \item \code{CS}: CS (rounded to 2 decimals in plot labels)
    }
}
}
\examples{
# load package input data
data("CDR3ab", package = "ClustIRR")
s <- data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], 
                CDR3a = CDR3ab[1:100, "CDR3a"],
                sample = rep(c("A", "B"), 50), 
                clone_size = rpois(100, 1)+1)

# run analysis
c <- clustirr(s = s)

# Compute similarities and plot
result <- get_cdr3_cosine_similarity(c)

# Display heatmap
print(result$g)

# Inspect similarity values
head(result$cs)

}