\name{comoNet_from_prog_sp}
\alias{comoNet_from_prog_sp}
\docType{package}
\title{
Create Raw Progression Comorbidity Network from Individual-level Progression Sequence
}
\description{
The function constructs the directed, weighted disease comorbidity network by counting the potential progression links in the individual pregression sequence.
}
\usage{
comoNet_from_prog_sp(ipseq, dz.lv, bidir = F, conditional = T)	
}

\arguments{
  \item{ipseq}{ list of individual progression sequences }
  \item{dz.lv}{ the level of the diagnostic code the defines the diseases }
  \item{bidir}{ True (default); whether to include bidirectional links for diseases didnosed at the same date/time }
  \item{conditional}{ False (default); whether to calculate conditional, directed, weighted disease comorbidity network. If True, the function will output an adjacency array instead of an adjacency matrix. }
  
}

\examples{
\dontrun{
    # Convert to individual progression sequence from demo dataset
    require(comoNet)
    require(Matrix)
    head(diagnostic_data)
    ipseq <- get_ipseq(df = diagnostic_data, var_id = "ID", var_t = "ageday", 
                       var_dz = c("DZ1","DZ2","DZ3","DZ4","DZ5"))
    
    dzLevel <- sort(unique(
      c(diagnostic_data$DZ1,diagnostic_data$DZ2, diagnostic_data$DZ3, diagnostic_data$DZ4, diagnostic_data$DZ5)
    ))
    res <- comoNet_from_prog_sp(ipseq = ipseq, dz.lv = dzLevel, bidir = F, conditional = F)
    raw_net <- res[[1]]
    image(raw_net)
}
}


\references{
Fotouhi, B., Momeni, N., Riolo, M. A., & Buckeridge, D. L. (2018). Statistical methods for constructing disease comorbidity networks from longitudinal inpatient data. Applied network science, 3(1), 46.
}
