\name{get_ag_species_hits}
\alias{get_ag_species_hits}
\title{Annotate antigen species hits in node summary}
\description{
Annotates clones (nodes) in the node summary output from the function 
\code{detect_communities} with specificity to given antigen species. 
Adds new binary columns indicating hits and computes antigen-specific 
cellular statistics per community and sample.
}
\details{
Searches for antigen species matches (e.g., "EBV" or "CMV") in annotation 
columns from specified databases (e.g., "vdjdb"). For each match, adds a 
new column to \code{node_summary} (1 = hit, 0 = no hit) and computes:
\itemize{
  \item Repertoire-level: Total cells/clones per sample
  \item Community-level: Cells/clones per community
  \item Antigen-specific: Cells/clones per antigen-community
}
Requires \code{node_summary} from \code{detect_communities} with database 
annotation columns).
}
\usage{
get_ag_species_hits(node_summary, 
                    db = "vdjdb", 
                    ag_species)
}
\arguments{
\item{node_summary}{Node summary data.frame from \code{detect_communities}}
\item{db}{Annotation database (e.g., "vdjdb", "mcpas", "tcr3d")}
\item{ag_species}{Antigen species (e.g., "EBV", "CMV", "SARS-CoV-2")}
}
\value{
List containing:
\item{node_summary}{Input with added antigen hit columns. See \code{new_columns}
for the names of the added columns}
\item{new_columns}{Names of added columns}
\item{table_summary}{Aggregated data.frame with columns:
  \itemize{
    \item \code{sample}, \code{community}, \code{ag_key}
    \item \code{rep_cells}, \code{rep_clones} (repertoire totals)
    \item \code{com_cells}, \code{com_clones} (community totals)
    \item \code{ag_cells}, \code{ag_clones} (antigen-specific counts)
  }
}
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
c <- clustirr(s = rbind(a, b))

# detect communities
gcd <- detect_communities(graph = c$graph, 
                          algorithm = "leiden",
                          resolution = 1,
                          weight = "ncweight",
                          chains = c("CDR3a", "CDR3b"))

ag <- get_ag_species_hits(node_summary = gcd$node_summary,
                          db = "vdjdb",
                          ag_species = "EBV")

head(ag)
}