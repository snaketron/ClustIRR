\name{get_graph}
\alias{get_graph}
\title{Get graph object}
\usage{get_graph(clust_irr)}
\arguments{
\item{clust_irr}{object of type \code{clust_irr}, result of \code{clust_irr} 
function, formed by two sublists \itemize{
\item \code{clust} - list,  clustering results
\item \code{inputs} - list, original user inputs
}
}
}
\value{
The output is an \code{igraph} object, that can be used to visualize the graph,
for example using \emph{visNetwork} or \emph{cytoscape}.

The graph consists of the following edge attributes:
\itemize{
\item \code{from} - ID of the source CDR3 sequence in the sample dataset
\item \code{to} - ID of the target CDR3 sequence in the sample dataset
\item \code{from_cdr3} - AA of the source CDR3 sequence in the sample dataset
\item \code{to_cdr3} - AA of the target CDR3 sequence in the sample dataset
\item \code{motif} - If type \emph{local}, shared motif between sequences. 
\code{NA}, if type \emph{global}
\item \code{type} - Type of connection, \emph{global} for globally clustered
or \emph{local} for locally clustered
\item \code{chain} - Chain the sequence belongs to, e.g., \emph{CDR3b}
}
}


\description{
This is a helper function for visualization preparation of an \code{clustIRR}
object.
}

\examples{
# load package input data
data("CDR3ab")
s <- base::data.frame(CDR3b = CDR3ab[1:1000, "CDR3b"])
r <- base::data.frame(CDR3b = CDR3ab[1:5000, "CDR3b"])

# artificially enrich motif 'RQWW' inside sample dataset
base::substr(x = s$CDR3b[1:20], start = 6, stop = 9) <- "RQWW"

# add an artificial clonal expansion of two sequences to the sample dataset
clones <- base::data.frame(CDR3b = rep(x = c("CATSRAAKPDGLAALETQYF",
                                       "CATSRAAKPDGLAALSTQYF"),
                                       times = 15))
s <- base::rbind(s, clones)

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

# create igraph object
out_ig <- get_graph(out)

# output class
base::class(out_ig)

# extract edges
utils::head(igraph::as_data_frame(out_ig, what = "edges"))

# extract nodes
utils::head(igraph::as_data_frame(out_ig, what = "vertices"))


}
