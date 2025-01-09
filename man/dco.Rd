\name{dco}
\alias{dco}
\title{Model-based differential community occupancy (DCO) analysis}
\description{
This algorithm takes as input a community matrix, and quantifies the relative
enrichment/depletion of individual communities in each sample using a Bayesian
hierarchical model.
}
\usage{
dco(community_occupancy_matrix, mcmc_control, compute_delta=TRUE, groups = NA)
}
\arguments{
\item{community_occupancy_matrix}{matrix, rows are communities, columns 
are repertoires, matrix entries are numbers of cells in each community and 
repertoire.}
\item{mcmc_control}{list, configurations for the Markov Chain Monte Carlo 
(MCMC) simulation.}
\itemize{
\item mcmc_warmup = 750; number of MCMC warmups
\item mcmc_iter = 1500; number of MCMC iterations
\item mcmc_chains = 4; number of MCMC chains
\item mcmc_cores = 1; number of computer cores 
\item mcmc_algorithm = "NUTS"; which MCMC algorithm to use
\item adapt_delta = 0.95; MCMC step size
\item max_treedepth = 12; the max value, in exponents of 2, of what the binary 
tree size in NUTS should have.
}
\item{compute_delta}{should delta be computed by the Stan model? This may be 
take up extra memory.}
\item{groups}{vector with integers \eqn{\geq1}, one for each repertoire (column
in community_occupancy_matrix). This specifies the biological group of each 
repertoire (e.g. for cancer repertoire we may specify the index 1, and for 
normal repertoires the index 2). If this vector is specified, ClustIRR will 
employ a hierarchical model, modeling the dependence between the repertoires 
within each group. Else (which is the default setting in ClustIRR), ClustIRR 
will treat the repertoires as independent samples by employing a simpler model.}
}
\value{
The output is a list with the folling elements:
\item{fit}{model fit (stan object)}
\item{posterior_summary}{nested list with data.frames, summary of model 
parameters, including their means, medians, 95\% credible intervals, etc.
Predicted observations (y_hat), which are useful for posterior predictive
checks are also provided.}
\item{community_occupancy_matrix}{matrix, rows are communities, columns are 
repertoires, matrix entries are numbers of cells in each community and 
repertoire.}
\item{mcmc_control}{mcmc configuration inputs provided as list.}
\item{compute_delta}{the input compute_delta.}
\item{groups}{the input groups.}
}

\examples{
# load package input data
data("CDR3ab", package = "ClustIRR")
a <- data.frame(CDR3a = CDR3ab[1:500, "CDR3a"],
                  CDR3b = CDR3ab[1:500, "CDR3b"],
                  clone_size = 1,
                  sample = "a")

b <- data.frame(CDR3a = CDR3ab[401:900, "CDR3a"],
                  CDR3b = CDR3ab[401:900, "CDR3b"],
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

# look at the community matrix
head(gcd$community_occupancy_matrix)

# look at the community summary
head(gcd$community_summary$wide)

# look at the node summary
head(gcd$node_summary)

# differential community occupancy analysis
dco <- dco(community_occupancy_matrix = gcd$community_occupancy_matrix)

names(dco)
}
