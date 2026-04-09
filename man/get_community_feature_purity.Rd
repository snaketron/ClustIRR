\name{get_community_feature_purity}
\alias{get_community_feature_purity}
\title{Compute community purity with respect to a node feature}
\description{
Compute the purity of graph communities with respect to a node feature.

Given a \code{node_summary} object returned by \code{detect_communities} 
containing a \code{community} column, this function evaluates how homogeneous 
each community is with respect to a specified node feature.

If the feature is numeric (e.g., gene expression, clone size), purity is 
measured using the coefficient of variation (CV). If the feature is categorical 
(e.g., tissue type, sample, or cell type), purity is measured using 
Gini impurity (GI) and Shannon's entropy (H). For binary categorical features, 
the signed dominance score (D) is additionally reported to indicate which 
class dominates a community.
}
\details{
The function verifies that \code{node_summary} contains a \code{community}
column and that the requested \code{node_feature} exists. The feature type is 
determined automatically.

For categorical features, purity metrics are computed per community:

\itemize{
\item \strong{Gini impurity (GI)} measures the probability that two randomly 
selected nodes from the community have different feature labels.
\item \strong{Shannon's entropy (H)} measures how evenly feature labels are 
distributed within the community.
\item \strong{Signed dominance (D)} is computed when the feature has exactly 
two classes and quantifies which class dominates the community.
}

Lower values of GI and H indicate higher purity (i.e., communities dominated 
by a single feature value). The dominance score ranges from -1 to 1, where 
positive and negative values indicate dominance of one or the other class.

For numeric features, the coefficient of variation (CV) is computed:

\itemize{
\item \strong{Mean} and \strong{standard deviation (sd)} are computed per community.
\item \strong{CV = sd / mean} measures the relative dispersion of the feature 
within the community.
}

Lower CV values indicate more homogeneous communities with respect to the 
numeric feature.
}
\usage{
get_community_feature_purity(node_summary, node_feature)
}
\arguments{
\item{node_summary}{
A \code{data.frame} returned by \code{detect_communities} containing node 
attributes and a column named \code{community}.
}

\item{node_feature}{
Character scalar specifying the node attribute for which community purity 
should be evaluated. The feature must be either numeric or categorical.
}
}
\value{
A \code{data.frame} summarizing purity statistics for each community.

For categorical features the returned columns are:

\itemize{
\item \code{community} – community identifier
\item \code{GI} – Gini impurity
\item \code{H} – Shannon's entropy
\item \code{D} – signed dominance score (only meaningful for binary features)
\item \code{n} – number of nodes in the community
}

For numeric features the returned columns are:

\itemize{
\item \code{community} – community identifier
\item \code{mean} – mean feature value
\item \code{sd} – standard deviation
\item \code{n} – number of nodes in the community
\item \code{cv} – coefficient of variation
}
}
\examples{
# only run if blast is installed
if(has_blast()) {
    # load package input data
    data("D1", package = "ClustIRR")
    
    # run ClustIRR analysis
    c <- clustirr(s = rbind(D1$a, D1$b), meta = rbind(D1$ma, D1$mb))
    
    # detect communities
    gcd <- detect_communities(graph = c$graph, 
                              algorithm = "leiden",
                              metric = "average",
                              resolution = 1,
                              iterations = 100,
                              chains = c("CDR3a", "CDR3b"))
    
    # categorical feature
    get_community_feature_purity(node_summary = gcd$node_summary, 
                                 node_feature = "HLA_A")
    
    # another categorical feature
    get_community_feature_purity(node_summary = gcd$node_summary, 
                                 node_feature = "TRAV")
    
    # numeric feature
    get_community_feature_purity(node_summary = gcd$node_summary, 
                                 node_feature = "age")
}
}
