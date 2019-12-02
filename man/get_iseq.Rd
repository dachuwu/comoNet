\name{get_iseq}
\alias{get_iseq}
\docType{package}
\title{
Get Individual Disease Progression Sequence
}
\description{

}
\usage{
get_iseq(x, df, 
           var_id = "ID", var_t = "date", 
           var_diag = c("ICD_1","ICD_2","ICD_3","ICD_4","ICD_5"))	
}
\arguments{
  \item{x}{An id}
  \item{df}{A dataframe}
  \item{var_id}{The variable name of personal identification}
  \item{var_t}{The variable name of date/time of diagnosis}
  \item{var_diag}{The variable names of diagnosed diseases (a character vector)}
  
}
\value{
A dataframe of individual disease progression sequence with two variable: on is the type of the disease, and the other is the date/time of the diagnosis.
}
\examples{
\dontrun{

}
}
