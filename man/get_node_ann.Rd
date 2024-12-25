\name{get_node_ann}
\alias{get_node_ann}
\title{Get antigen annotations for nodes}
\description{
Use \code{node_summary} data created by the function \code{detect_communities} 
to perform antigen annotation based databases such as VDJdb.}
\details{
The user has to provide antigen name (e.g. \code{ag_name} = "EBV") and 
antigen source (\code{ag_source} = "species" or \code{ag_source} = "gene"). 
Furthermore, the user has to provide nodes (\code{node_summary} data.frame 
created by the function \code{detect_communities}).}
\usage{
get_node_ann(node_summary, ag_name, ag_source)
}
\arguments{
\item{node_summary}{\code{node_summary} data.frame}
\item{ag_name}{antigen name, character, e.g. "EBV".}
\item{ag_source}{antigen source, character, "species" or "gene".}
}
\value{
The output is a data.frame with nodes and their antigen annotations,
provided for each immune receptor chain and annotation database.
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

# annotate nodes. Which nodes contain CDR3s specific for antigens 
# from the EBV
an <- get_node_ann(node_summary = gcd$node_summary, 
                   ag_name = "EBV",
                   ag_source = "species")
kable(an)
}
