\name{order_param}
\alias{order_param}

\title{
  Function to obtain an ordered matrix of values and information for
  specific keywords/arguments/parameters. 
}

\description{
  The function \code{\link{order_param}} extracts the numerical
  information from an unordered vector of strings of characters,
  \code{param}, into a matrix following the order of
  keywords/arguments/parameters specified by \code{listparam}. \cr 
  This function is called by the \R functions \code{\link{simulate_data}}
  and \code{\link{estimate_IMc}}. 
}
\usage{
order_param(param, listparam)
}

\arguments{
  \item{param}{
    The unsorted vector of strings of characters with the values specified
    for various keywords/arguments/parameters. \cr 
    The keyword/argument/parameter names is the first non-blank "word"
    the string of characters. \cr
    See \code{\link{param_sim}} and \code{\link{param_est}} for more
    information on those values and keywords/arguments/parameters. 
  }
  \item{listparam}{
    The array of possible keywords/arguments/parameters in the required
    order. \cr
     For function   \code{\link{simulate_data}}:\cr
    \code{c("}\emph{\code{theta_1}}\code{", "}\emph{\code{M_present}}\code{", "}\emph{\code{theta_2}}\code{", "}\emph{\code{theta_A}}\code{", "}\emph{\code{T_split}}\code{", "}\emph{\code{T_change}}\code{", "}\emph{\code{M_change}}\code{", "}\emph{\code{rho}}\code{", "}\emph{\code{nregions}}\code{")}.\cr
    For function  \code{\link{estimate_IMc}}:\cr
    \code{c("}\emph{\code{theta_1}}\code{", "}\emph{\code{M_present}}\code{", "}\emph{\code{theta_2}}\code{", "}\emph{\code{theta_A}}\code{", "}\emph{\code{T_split}}\code{", "}\emph{\code{T_change}}\code{", "}\emph{\code{M_change}}\code{", "}\emph{\code{rho}}\code{", "}\emph{\code{nregions}}\code{", "}\emph{\code{howmany}}\code{", "}\emph{\code{parallel}}\code{")}.
  }
}

\value{
  The function \code{\link{order_param}} outputs a matrix of three columns
  with the numerical values specified after the
  keywords/arguments/parameters in \code{param}. \cr 
  The rows follow the order specified by the list \code{listparam}. 
}

\author{
  Celine Becquet - \email{celine.becquet@gmail.com}. 
}

\seealso{
  The functions \code{\link{simulate_data}} and
  \code{\link{estimate_IMc}} call the function \code{\link{order_param}},
  which in turn calls the function \code{\link{get_nparam}}. 
}

\examples{
### Write the file of parameter values in the local directory.
data(simulation_files) # download the data

write.table(file="sim", x=param_sim, row.name=FALSE, col.names=FALSE, quote=FALSE)
# Creates the file "sim" containing (with comments):
read.csv("sim",header =FALSE,sep="", comment="#")

### Create the inputs for the function.
param=scan("sim", comment.char="#", what=c("c", "n", "n", "n"), allowEscapes=FALSE, fill=TRUE, sep="\n", strip.white=TRUE, quiet=TRUE)
param

## Matrix of ordered parameter values.
listparam=c("theta_1", "M_present", "theta_2", "theta_A", "T_split", "T_change", "M_change", "rho", "nregions")
order_param(param, listparam) 

## A new possible order.
listparam=c("theta_1", "theta_2", "theta_A", "T_split", "T_change", "M_present", "M_change", "rho", "nregions")
order_param(param, listparam) 

# Clean up the directory.
unlink("sim")
}

\keyword{misc}