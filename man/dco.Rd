\name{dco}
\alias{dco}
\title{Model-based differential community occupancy (DCO) analysis}
\usage{
detect_communities(community_matrix, mcmc_control)
}
\arguments{
\item{community_matrix}{matrix, rows are communities, columns are repertoires,
matrix entries are numbers of cells in each community and repertoire.}
\item{mcmc_control}{list, configurations for the Markov Chain Monte Carlo 
(MCMC) simulation.}
\itemize{
\item mcmc_warmup = 750; number of MCMC warmups
\item mcmc_iter = 1500; number of MCMC iterations
\item mcmc_chains = 4; number of MCMC chains
\item mcmc_chains = 1; number of computer cores 
\item mcmc_algorithm = "NUTS"; which MCMC algorithm to use
\item adapt_delta = 0.95; MCMC step size
\item max_treedepth = 12; the max value, in exponents of 2, of what the binary 
tree size in NUTS should have.
}
}
\value{
The output is a list with the folling elements:
\item{fit}{stan object, model fit}
\item{posterior_summary}{nested list with data.frames, summary of model 
parameter fits}
\item{community_matrix}{matrix, rows are communities, columns are repertoires,
matrix entries are numbers of cells in each community and repertoire}
\item{mcmc_control}{list, mcmc configuration inputs provided as list.}
}
\description{
This algorithm takes as input a community matrix, and quantifies the relative
enrichment/depletion of individual communities in each sample using a Bayesian
hierarchical model.
}

\examples{
# load package input data
data("CDR3ab")
s_1 <- data.frame(CDR3a = CDR3ab[1:500, "CDR3a"],
                  CDR3b = CDR3ab[1:500, "CDR3b"],
                  clone_size = 1)

s_2 <- data.frame(CDR3a = CDR3ab[401:900, "CDR3a"],
                  CDR3b = CDR3ab[401:900, "CDR3b"],
                  clone_size = 1)

s_2$clone_size[1] <- 20

# run ClustIRR analysis
c_1 <- cluster_irr(s = s_1, control = list(global_smart=TRUE))
c_2 <- cluster_irr(s = s_2, control = list(global_smart=TRUE))

# get joint graph
cs <- c(c_1, c_2)
names(cs) <- c("A", "B")
jg <- get_joint_graph(clust_irrs = cs)

# detect communities
coms <- detect_communities(graph = jg$graph, 
                           algorithm = "leiden",
                           resolution = 1,
                           weight_type = "ncweight",
                           chains = c("CDR3a", "CDR3b"))

# look at outputs
names(coms)

# look at the community matrix
head(coms$community_matrix)

# look at the community summary
head(coms$community_summary)

# look at the node summary
head(coms$node_summary)

# differential community occupancy analysis
dco <- dco(community_matrix = coms$community_matrix)

names(dco)

}
