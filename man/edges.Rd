\name{get_edges}
\alias{get_edges}
\title{Get cluster edges data.frame}
\usage{get_edges(clust_irr)}
\arguments{
\item{clust_irr}{object of type \code{clust_irr}, result of \code{clust_irr} 
function, formed by two sublists \itemize{
\item \code{clust} - list,  clustering results
\item \code{inputs} - list, original user inputs
}
}
}
\value{
The output is an data.frame of the local and global edges, that can be
used as input for visualization, for example with \emph{igraph}.
\itemize{
\item \code{from_ID} - ID of the source CDR3 sequence in the sample dataset
\item \code{to_ID} - ID of the target CDR3 sequence in the sample dataset
\item \code{from} - AA of the source CDR3 sequence in the sample dataset
\item \code{to} - AA of the target CDR3 sequence in the sample dataset
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
data("CD8")
data_sample <- data.frame(CDR3b = CD8[1:5000, "CDR3b"])
data_ref <- CD8

# run analysis
out <- cluster_irr(data_sample = data_sample,
                   data_ref = data_ref,
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

# output class
class(out)

# output structure
str(out)

# inspect local motif enrichment results
knitr::kable(head(out$clust$CDR3b$local$m))

# inspect which CDR3bs are globally similar
knitr::kable(head(out$clust$CDR3b$global))

# create data.frame
out_df <- get_edges(out)

# TODO visualize data with igraph/visNetwork

}
