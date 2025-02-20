\name{decode_communities}
\alias{decode_communities}
\title{Decode graph communities}
\description{
Given a graph and community ID, the function explores the heterogeneity 
within communities and summarizes its features.
}
\details{

}
\usage{
decode_communities(community_id,
                   graph,
                   edge_at,
                   node_at)
}
\arguments{
\item{graph}{\code{igraph} object that has been analyzed by graph-based 
community detection methods as implemented in \code{detect_communities}}
\item{community_id}{which community should be decoded?}
\item{edge_at}{which edge attributes should be used to decode the community}
\item{node_at}{which node attributes should be used to decode the community}
}
\value{

}
\examples{

}
