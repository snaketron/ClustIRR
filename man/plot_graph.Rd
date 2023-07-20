\name{plot_graph}
\alias{plot_graph}
\title{Plot clustIRR graph}
\usage{plot_graph(clust_irr, expand_clones = FALSE)}
\arguments{
\item{clust_irr}{object of type \code{clust_irr}, result of \code{clust_irr} 
function, formed by two sublists \itemize{
\item \code{clust} - list,  clustering results
\item \code{inputs} - list, original user inputs
}
}
\item{expand_clones}{logical, show clonally expanded cells as distinct nodes
in the graph, including connections between them. Only relevant for version 3
of the algorithm. \code{expand_clones} = \code{FALSE} (default)
}
}
\value{
The output is an interactive \emph{visNetwork} graph. 

Nodes represent cells, edges represent clustered local or global similarities.

Informations about shared motifs and/or global sequences can be accessed via 
hovering over nodes and edges (copy&paste possible). Nodes can be selected 
directly, via sequence and by shared local or global connections. 

Colors and symbols of each graph are explained in a legend. 
Connections between chains are represented by different colors. 
Clone counts are represented by different symbols and sizes, cells 
sharing the exact same sequences (clonal expansions) are represented by dots 
expanded in size. Edge color represents local connections, global connections 
or both between cells. Edge thickness represents number of shared local and/or 
global connections between cells.

Clonally expanded cells can be represented by single dots, containing a list of
in between shared motifs and sequences, or multiple dots with edges between 
them by setting \code{expand_clones} to \code{TRUE}.

\itemize{
  \item{\code{scroll}} - click & drag background or edges
  \item{\code{zoom}} - use the mousewheel
  \item{\code{show information}} - hover-over node or edge
  \item{\code{select node directly}} - click on node
  \item{\code{select node by id}} - dropdown menu to select node via CDR3 
  sequence, sorted alphabetically
  \code{select node by group} - dropdown menu to select node via shared motif 
  or global sequence, sorted alphabetically
  \item{\code{move node}} - click & drag node
  \item{\code{export graph}} - use console, viewer and \code{export} button
  
}

}


\description{
This is a helper function for visualization of an \code{clustIRR} object.
}

\examples{
# load package input data
data("CDR3ab")
s <- base::data.frame(CDR3b = CDR3ab[1:1000, "CDR3b"])
r <- base::data.frame(CDR3b = CDR3ab[1:5000, "CDR3b"])

# artificially enrich motif 'RWGW' inside sample dataset
base::substr(x = s$CDR3b[1:20], start = 6, stop = 9) <- "RWGW"

# add an artificial clonal expansion of two sequences to the sample dataset
s <- rbind(s, base::data.frame(CDR3b = rep(x = c("CATSRADKPDGLDALETQYF",
                                       "CATSRAAKPDGLAALSTQYF"),
                                       times = 5)))

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

# plot graph with expanded clones as individual cells/nodes
plot_graph(out, expand_clones = TRUE)

}
