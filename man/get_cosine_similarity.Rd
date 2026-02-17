\name{get_cosine_similarity}
\alias{get_cosine_similarity}
\title{Compute and Visualize Cosine Similarity (CS)}
\description{
This function calculates pairwise CSs between columns (repertoires) 
of the community occupancy matrix and generates a labeled heatmap 
visualization. Here CS  measures the cosine of the angle between 
vectors of cell counts (positive integers or zero), providing a 
value between 0 (perfect dissimilarity) and 1 (perfect similarity).
}
\details{
The function performs three main operations:
\enumerate{
  \item Computes CS between all column pairs using: 
  \deqn{\text{CS}(x, y) = \frac{x \cdot y}{||x|| \cdot ||y||}}
  \item Creates a heatmap with CS values annotated in each cell
}
Red and blue tiles in the heatmap indicates high and low CS, respectively.
}
\usage{
get_cosine_similarity(com)
}
\arguments{
\item{com}{Community occupancy matrix where columns represent repertoires for 
CS computation. Rows should are communities. Entries are the number of
cells in each repertoire and community.}
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
# Create a sample matrix
mat <- matrix(rpois(n=15, lambda = 4), nrow = 5, 
              dimnames = list(NULL, c("A", "B", "C")))

# Compute similarities and plot
result <- get_cosine_similarity(mat)

# Display heatmap
print(result$g)

# Inspect similarity values
head(result$cs)

}