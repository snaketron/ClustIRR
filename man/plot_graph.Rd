\name{plot_graph}
\alias{plot_graph}
\title{Plot ClustIRR graph}
\usage{plot_graph(clust_irr, as_visnet=FALSE)}
\arguments{
\item{clust_irr}{S4 object of type \code{clust_irr}, result of \code{clust_irr} 
function}
\item{as_visnet}{logical, if as_visnet=TRUE we plot an interactive graph with 
visNetwork. If as_visnet=FALSE, we plot a static graph with igraph.}}
\value{
The output is an igraph plot.

The size of the vertices increases linearly as the logarithm of the degree of 
the clonal expansion (number of cells per clone) in the corresponding clones.
}

\description{
This this function visualizes a graph. The input is \code{clust_irr} object
created by the function \code{cluster_irr}.
}

\examples{
# load package input data
data("CDR3ab")
s <- data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], clone_size = 1)
r <- data.frame(CDR3b = CDR3ab[1:500, "CDR3b"], clone_size = 1)

# artificially enrich motif 'RWGW' inside sample dataset
substr(x = s$CDR3b[1:20], start = 6, stop = 9) <- "RWGW"

# add an artificial clonal expansion of two sequences to the sample dataset
s <- rbind(s, data.frame(CDR3b = c("CATSRADKPDGLDALETQYF", 
                                   "CATSRAAKPDGLAALSTQYF"),
                         clone_size = 5))

# run analysis
out <- cluster_irr(s = s,
                   r = r,
                   ks = 4,
                   cores = 1,
                   control = list(trim_flank_aa = 3))

# plot graph with vertices as clones
p1 <- plot_graph(out, as_visnet=FALSE)
p1
}
