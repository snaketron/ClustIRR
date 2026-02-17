\name{save_graph}
\alias{save_graph}
\title{Save ClustIRR graph}
\usage{save_graph(graph,
            file_name,
            output_folder)
}
\arguments{
\item{graph}{Object returned by the function \code{plot_graph}.}
\item{file_name}{Name of the exported graph.}
\item{output_folder}{Path, where to store the graph.}
}
\value{
This function saves an interactive \code{visNetwork} plot as self-contained
html file. If \code{output_folder} is not existent, it will be saved in the 
working directory.
}

\description{
This function saves an interactive \code{visNetwork} object as html file
of a \code{clustirr} graph.
}


\examples{
# load package input data
data("CDR3ab", package = "ClustIRR")
s <- data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], sample = "A", clone_size = 1)

# run ClustIRR analysis
c <- clustirr(s = s)

# plot graph with vertices as clones
g <- plot_graph(c, as_visnet=FALSE, show_singletons=TRUE, node_opacity = 0.8)

\dontrun{
# Save the graph to a temporary directory
my_temp_dir <- tempdir()
save_graph(g, file_name = "test_graph", output_folder = my_temp_dir)
}
}
