\name{attach_clonotypes}
\alias{attach_clonotypes}
\title{Attach Clonotypes to Existing Communities}
\description{
Assigns clonotypes to existing communities by matching their amino
acid sequences against the sequences of nodes assigned to communities
using BLASTP.
}
\usage{
attach_clonotypes(s, node_summary, control, keep_best = TRUE)
}
\arguments{
  \item{s}{
    A data frame containing clonotypes to attach. This data.frame must 
    have an identical as the main input provided to the function 
    \code{clustirr()}. Each row represents a clonotype.
  }
  \item{node_summary}{
    A data frame of nodes and their properties, obtained from the
    \code{node_summary} element of the output of
    \code{detect_communities()}. Rows represent nodes (clones)
    assigned to communities. The data frame must include the relevant
    amino acid sequence columns and a \code{community} column.
  }
  \item{control}{
    A list of auxiliary parameters controlling sequence matching.
    The following elements are used:
    \itemize{
      \item \code{blast_gmi}: minimum sequence identity between a
      query sequence and a node sequence for the match to be retained
      (default = 0.9, corresponding to 90\% identity).
      \item \code{blast_cores}: number of CPU cores to use with BLAST
      (default = 1).
      \item \code{trim_flank_aa}: number of amino acids trimmed from
      the flanks of CDR3 sequences when calculating alignment scores
      for the CDR3 cores (default = 3).
    }
    These parameters can be supplied using the same control list as
    in \code{clustirr()}.
  }
  \item{keep_best}{
    Logical. If \code{TRUE}, retains only the highest-ranked match
    for each query clonotype within each receptor chain, ranked by
    normalized BLAST score in decreasing order. Defaults to
    \code{TRUE}.
  }
}
\details{
  This function can be used to attach additional clonotypes to
  communities identified by \code{detect_communities()}.
  The \code{node_summary} argument should be obtained from the
  \code{node_summary} component of the output of that function.

  For each receptor chain returned by \code{get_chains()}, the function
  constructs a protein BLAST database from the sequences in
  \code{node_summary} and queries it using the corresponding
  sequences in \code{s}.

  Matches with percent identity below
  \code{100 * control$blast_gmi} are discarded. Alignments that do
  not cover the full length of both query and subject sequences are
  also discarded.

  Retained matches are ranked by normalized BLAST score, calculated
  as the raw alignment score divided by alignment length. If
  \code{keep_best = TRUE}, only the highest-ranked match per query
  clonotype is retained within each chain.
}
\value{
  Match data frame (one for each chain) contains the following columns:
  \describe{
    \item{\code{from}}{
      Row index of the query clonotype in \code{s}.
    }
    \item{\code{to}}{
      Row index of the matched node in \code{node_summary}.
    }
    \item{\code{len}}{
      Length of the BLAST alignment.
    }
    \item{\code{weight}}{
      Raw BLAST alignment score.
    }
    \item{\code{nweight}}{
      Raw BLAST alignment score normalized by alignment length.
    }
    \item{\code{cweight}}{
      BLOSUM62-based alignment score calculated by
      \code{get_bscore}.
    }
    \item{\code{clen}}{
      Alignment length used to normalize the BLOSUM62-based score.
    }
    \item{\code{ncweight}}{
      BLOSUM62-based alignment score normalized by \code{clen}.
    }
    \item{\code{chain}}{
      Receptor chain associated with the match.
    }
    \item{\code{community}}{
      Community identifier of the matched node, obtained from
      \code{node_summary$community}.
    }
  }
}
\seealso{
  \code{clustirr}, \code{detect_communities}
}
