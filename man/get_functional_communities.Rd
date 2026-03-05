\name{get_functional_communities}
\alias{get_functional_communities}
\title{Detect functionally similar communities}
\description{
Detect functionally similar nodes and communities in a supplied 
detect_communities output. Adds new binary columns to node_summary and 
communities and adds index vector of functional communities to the communities
object.
}
\details{
Searches for antigen gene and/or species matches based on the supplied vectors
of antigens and species. Adds binary functional_node column to node_summary and 
functional_community column if any of the antigens or species are found.
}
\usage{
get_functional_communities(com, ag_gene, ag_species)
}
\arguments{
\item{com}{Community occupancy matrix where columns represent repertoires for 
CS computation. Rows are communities. Entries are the number of cells in each 
repertoire and community.}
\item{ag_gene}{Vector of genes}
\item{ag_species}{Vector of species}
}
\value{
The output is a list with the folling elements:
\item{community_occupancy_matrix}{matrix, rows are communities, columns 
are repertoires, matrix entries are numbers of cells in each community and 
repertoire.}
\item{community_summary}{data.frame, rows are communities and their properties
are provided as columns with added column for functionality.}
\item{node_summary}{data.frame, rows are nodes (clones) and their properties 
are provided as columns. Contains additional columns for functionality and 
ag_gene as well as ag_species hits.}
\item{graph}{igraph object, processed graph object.}
\item{graph_structure_quality}{graph modularity and quality (only for Leiden) 
measure of the strength of division of the graph into communities.}
\item{input_config}{list, inputs provided as list.}
\item{functional_communities}{Numeric vector of functional community indices.}
}
\examples{
# test data
data("CDR3ab", package = "ClustIRR")
s <- data.frame(CDR3b = CDR3ab[, "CDR3b"],
                CDR3a = CDR3ab[, "CDR3a"],
                clone_size = 1,
                sample = rep(c("a", "b", "c", "d"), 125))[1:500,]

# run ClustIRR analysis
c <- clustirr(s = s)

com <- detect_communities(graph = c$graph,
                          algorithm = "leiden",
                          resolution = 1, 
                          weight = "ncweight",
                          metric = "average",
                          chains = c("CDR3b", "CDR3a"))

ag_gene <- c("MLANA", "MART")
ag_species <- c("HomoSapiens")

com <- get_functional_communities(com, ag_gene, ag_species)

# communities found to be specific
com$functional_communities

}