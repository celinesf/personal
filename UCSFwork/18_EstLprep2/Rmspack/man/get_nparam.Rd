\name{get_nparam}
\alias{get_nparam}

\title{
  Function to obtain the position in a list of a particular
  keyword/argument/parameter. 
}

\description{
  The function \code{\link{get_nparam}} outputs the position in the list
  \code{listparam} of the keyword specified in the argument
  \code{"keyword"}. \cr 
  This function is called by the \R function \code{\link{order_param}}. 
}

\usage{
get_nparam(keyword, listparam)
}

\arguments{
  \item{keyword}{
    The name of a keyword/argument/parameter for which one requires the
    position in the list \code{listparam}. 
  }
  \item{listparam}{
    The array of possible keywords/arguments/parameters in the required
    order. \cr
     For function   \code{\link{simulate_data}}:\cr
    \code{c("}\emph{\code{theta_1}}\code{", "}\emph{\code{M_present}}\code{", "}\emph{\code{theta_2}}\code{", "}\emph{\code{theta_A}}\code{", "}\emph{\code{T_split}}\code{", "}\emph{\code{T_change}}\code{", "}\emph{\code{M_change}}\code{", "}\emph{\code{rho}}\code{", "}\emph{\code{nregions}}\code{")}.\cr
    For function  \code{\link{estimate_IMc}}:\cr
    \code{c("}\emph{\code{theta_1}}\code{", "}\emph{\code{M_present}}\code{", "}\emph{\code{theta_2}}\code{", "}\emph{\code{theta_A}}\code{", "}\emph{\code{T_split}}\code{", "}\emph{\code{T_change}}\code{", "}\emph{\code{M_change}}\code{", "}\emph{\code{rho}}\code{", "}\emph{\code{nregions}}\code{", "}\emph{\code{howmany}}\code{", "}\emph{\code{parallel}}\code{")}..
  }
}

\value{
  The function \code{\link{get_nparam}} outputs an integer between
  \code{[1,length(listparam)]}. \cr 
  If the \code{keyword} is not in \code{listparam}, then
  \code{\link{get_nparam}} returns \code{"0"}. 
}

\author{
  Celine Becquet - \email{celine.becquet@gmail.com}. 
}

\seealso{
  The function \code{\link{get_nparam}} is called by the function
  \code{\link{order_param}}, which in turn is called by the functions
  \code{\link{simulate_data}} and \code{\link{estimate_IMc}}. 
}

\examples{
### Create the input for the function.
listparam=c("theta_1","M_present","theta_2","theta_A","T_split","T_change","M_change","rho", "nregions")

get_nparam("theta_2", listparam) ## ->"3"
get_nparam("r", listparam)     ## ->"0"
}

\keyword{misc}