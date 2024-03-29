\name{filter_phi}
\alias{filter_phi}
\docType{package}
\title{
Build Disease Co-morbidity Network Based on Phi-correlation Measure
}
\description{
This is a filtering function that removes insignificant co-morbidity links/edges from the raw co-morbidity network generated from individual disease progression data. Links/edges  with significant Phi-correlation coefficient greater than 0 are kept as the resulting network.
}
\usage{
filter_phi(A, sig_lv)	
}
\arguments{
  \item{A}{adjacency matrix of raw progression network}
  \item{sig_lv}{the significance level}
}
\examples{
\dontrun{
  res <- filter_phi(as.matrix(raw_net),sig_lv = 0.05)
  phi_net <- res$A_filtered
  sub_net <- phi_net[1:30,1:30]
  ind <- which(rowSums(sub_net)>0 | colSums(sub_net)>0)
  sub_net <- sub_net[ind,ind]
  g1 <- graph_from_adjacency_matrix(sub_net)
  g1 <- set_vertex_attr(g1, "label", value = dzLevel[V(g1)])
  g1$layout <- layout_with_fr
  plot(g1)
}
}

\references{
Hidalgo, C. A., Blumm, N., Barabasi, A. L., & Christakis, N. A. (2009). A dynamic network approach for the study of human phenotypes. PLoS computational biology, 5(4), e1000353.
}