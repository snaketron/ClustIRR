\name{plot_graph}
\alias{plot_graph}
\title{Plot ClustIRR graph}
\usage{plot_graph(g,
           select_by = "Ag_species",
           as_visnet = FALSE,
           show_singletons = TRUE, 
           node_opacity = 1)
}
\arguments{
\item{g}{Object returned by the functions \code{get_graph} or 
\code{get_joint_graph}}
\item{as_visnet}{logical, if \code{as_visnet}=TRUE we plot an interactive graph 
with visNetwork. If \code{as_visnet}=FALSE, we plot a static graph with igraph.}
\item{select_by}{character string, two values are possible: "Ag_species" or
"Ag_gene". This only has an effect if \code{as_visnet} = TRUE, i.e. if the graph
is interactive. It will allow the user to highligh clones (nodes) in the graph
that are associated with a specific antigenic specie or gene. The mapping 
between CDR3 and antigens is extracted from databases, such as, VDJdb, 
McPAS-TCR and TCR3d. This mapping is done by the function \code{get_graph}. If
none of the clones in the graph are matched to a CDR3, then the user will
have no options to select/highlight.}
\item{show_singletons}{logical, if \code{show_singletons}=TRUE we plot all 
vertices. If \code{show_singletons}=FALSE, we plot only vertices connected by 
edges.}
\item{node_opacity}{probability, controls the opacity of node colors. Lower 
values corresponding to more transparent colors.}
}
\value{
The output is an \code{igraph} or \code{visNetwork} plot.

The size of the vertices increases linearly as the logarithm of the degree of 
the clonal expansion (number of cells per clone) in the corresponding clones.
}

\description{
This function visualizes a graph. The main input is \code{g} object
created by the function \code{get_graph}.
}

\examples{
# load package input data
data("CDR3ab", package = "ClustIRR")
s <- data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], sample = "A", clone_size = 1)

# run ClustIRR analysis
out <- cluster_irr(s = s,
                   control = list(gmi = 0.7,
                                  trim_flank_aa = 3,
                                  global_hamming = FALSE,
                                  global_max_hdist = 1,
                                  low_mem = FALSE))

# get graph
g <- get_graph(clust_irr = out)

# plot graph with vertices as clones
plot_graph(g, as_visnet=FALSE, show_singletons=TRUE, node_opacity = 0.8)
}
