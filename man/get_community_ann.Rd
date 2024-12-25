\name{get_community_ann}
\alias{get_community_ann}
\title{Get antigen annotations for communities}
\description{
Annotate the \code{node_summary} and \code{community_summary} data 
created by the function \code{detect_communities} with antigen annotations 
from databases such as VDJdb.}
\details{
The user has to provide a vector of antigen species (e.g. 
\code{ag_species} = c("EBV", "CMV")) and/or a vector of antigen 
genes (e.g. \code{ag_genes} = "MLANA"). Furthermore, the user has 
to provide nodes (\code{node_summary} data.frame created by the 
function \code{detect_communities}).}
\usage{
get_community_ann(node_summary, community_summary, ag_species, ag_genes)
}
\arguments{
\item{node_summary}{\code{node_summary} data.frame}
\item{community_summary}{\code{community_summary} data.frame}
\item{ag_species}{antigen species, character vector, e.g. "EBV".}
\item{ag_genes}{antigen genes, character vector, e.g. "MLANA".}
}
\value{
The output is a data.frame with communities and their antigen annotations,
provided for each chain and annotations database.
}
\examples{
# load package input data
data("CDR3ab", package = "ClustIRR")
a <- data.frame(CDR3a = CDR3ab[1:300, "CDR3a"],
                CDR3b = CDR3ab[1:300, "CDR3b"],
                clone_size = 1,
                sample = "a")

b <- data.frame(CDR3a = CDR3ab[201:400, "CDR3a"],
                  CDR3b = CDR3ab[201:400, "CDR3b"],
                  clone_size = 1,
                  sample = "b")
b$clone_size[1] <- 20

# run ClustIRR analysis
c <- c(cluster_irr(s = a), cluster_irr(s = b))

# get joint graph
jg <- get_joint_graph(clust_irrs = c)

# detect communities
gcd <- detect_communities(graph = jg$graph, 
                          algorithm = "leiden",
                          resolution = 1,
                          weight = "ncweight",
                          metric = "average",
                          chains = c("CDR3a", "CDR3b"))

# look at outputs
names(gcd)

# annotate communities. Which communities contain clones (nodes) which
# contain CDR3s specific for antigens from the EBV
an <- get_community_ann(node_summary = gcd$node_summary, 
                        community_summary = gcd$community_summary, 
                        ag_species = "EBV",
                        ag_genes = "MLANA")
kable(an)
}
