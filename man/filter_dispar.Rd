\name{filter_dispar}
\alias{filter_dispar}
\docType{package}
\title{
Build Disease Co-morbidity Network Based on Disparity Filter
}
\description{
This is a filtering function that removes insignificant co-morbidity links/edges from the raw co-morbidity network generated from individual disease progression data. The observed link/edge weights are compared with the network null model to determine each link/edge's p value, as in (Serrano, et al. 2009).
}
\usage{
filter_dispar(A, sig_lv)	
}
\arguments{
  \item{A}{adjacency matrix of raw progression network}
  \item{sig_lv}{the significance level}
}
\examples{
\dontrun{
  res <- filter_dispar(as.matrix(raw_net),sig_lv = 0.05)
  dispar_net <- res$A_filtered
  sub_net <- dispar_net[1:30,1:30]
  ind <- which(rowSums(sub_net)>0 | colSums(sub_net)>0)
  sub_net <- sub_net[ind,ind]
  g1 <- graph_from_adjacency_matrix(sub_net)
  g1 <- set_vertex_attr(g1, "label", value = dzLevel[V(g1)])
  g1$layout <- layout_with_fr
  plot(g1)
}
}

\references{
Serrano, M. A., Boguna, M., & Vespignani, A. (2009). Extracting the multiscale backbone of complex weighted networks. Proceedings of the national academy of sciences, 106(16), 6483-6488.
}