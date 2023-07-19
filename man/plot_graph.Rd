\name{plot_graph}
\alias{plot_graph}
\title{Get graph object}
\usage{plot_graph(clust_irr, expand_clones = FALSE)}
\arguments{
\item{clust_irr}{object of type \code{clust_irr}, result of \code{clust_irr} 
function, formed by two sublists \itemize{
\item \code{clust} - list,  clustering results
\item \code{inputs} - list, original user inputs
}
}
\item{expand_clones}{logical, show clonally expanded cells as distinct nodes
in the graph, including connections between them. 
\code{expand_clones} = \code{FALSE} (default)
}
}
\value{
The output is an interactive \emph{visNetwork} graph.
}


\description{
This is a helper function for visualization of an \code{clustIRR} object.
}

\examples{
# load package input data
data("CDR3ab")
s <- data.frame(CDR3b = CDR3ab[1:1000, "CDR3b"])
r <- data.frame(CDR3b = CDR3ab[1:5000, "CDR3b"])

# artificially enrich motif 'RQWW' inside sample dataset
substr(x = s$CDR3b[1:20], start = 6, stop = 9) <- "RQWW"

# add the artificial clonal expansion of two sequences to the sample dataset
clones <- data.frame(CDR3b = rep(x = c("CATSRAAKPDGLAALETQYF",
                                       "CATSRAAKPDGLAALSTQYF"),
                                 times = 15))
s <- rbind(s, clones)

# run analysis
out <- cluster_irr(s = s,
                   r = r,
                   version = 3,
                   ks = 4,
                   cores = 1,
                   control = list(
                      B = 1000,
                      global_max_dist = 1,
                      local_max_fdr = 0.05,
                      local_min_ove = 2,
                      local_min_o = 1,
                      trim_flank_aa = 3,
                      global_pairs = NULL,
                      low_mem = FALSE))


# plot graph
plot_graph(out)

}
