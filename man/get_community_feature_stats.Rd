\name{get_community_feature_stats}
\alias{get_community_feature_stats}
\title{Compute descriptive statistics of a community node feature}
\description{
Given a \code{node_summary} object containing community assignments and node-
level attributes, this function computes descriptive statistics of a specified 
node feature within a single community.
}
\details{
The function requires a \code{node_summary} data.frame containing a
\code{community} column and a valid \code{node_feature} column.

The type of statistics returned depends on whether the node feature is numeric
or categorical:

\itemize{
\item For numeric features, the function computes mean, median, sum, and the
number of nodes.
\item For categorical (character, factor, or logical) features, the function
computes counts and proportions for each feature level.
}

The community identifier is inferred from the \code{community} column and is
assumed to be constant within \code{node_summary}.
}
\usage{
get_community_feature_stats(node_summary, node_feature)
}
\arguments{
\item{node_summary}{A data.frame containing node-level attributes and a
\code{community} column.}
\item{node_feature}{Character scalar specifying the node attribute for which
statistics should be computed.}
}
\value{
A \code{data.frame} summarizing feature statistics for the given community:

\itemize{
\item For numeric features: columns \code{community}, \code{feature},
\code{feature_mean}, \code{feature_median}, \code{feature_sum},
\code{feature_type}, and \code{n}.
\item For categorical features: columns \code{community}, \code{feature},
\code{feature_count}, \code{feature_prop}, \code{feature_type}, and \code{n}.
}
}
\examples{
# load package input data
data("D1", package = "ClustIRR")

# run ClustIRR analysis
c <- clustirr(s = rbind(D1$a, D1$b), meta = rbind(D1$ma, D1$mb))

# detect communities
gcd <- detect_communities(graph = c$graph, 
                          algorithm = "leiden",
                          metric = "average",
                          resolution = 1,
                          weight = "ncweight",
                          iterations = 100,
                          chains = c("CDR3a", "CDR3b"))

# numeric feature
get_community_feature_stats(node_summary = gcd$node_summary[gcd$node_summary$community == 1, ],
                            node_feature = "age")

# categorical feature
get_community_feature_stats(node_summary = gcd$node_summary[gcd$node_summary$community == 1, ],
                            node_feature = "HLA_A")
}
