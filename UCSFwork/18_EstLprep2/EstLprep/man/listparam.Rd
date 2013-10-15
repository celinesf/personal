\name{listparam}
\alias{listparam}
\docType{data}
\title{
  List of all possible keywords used in "paramfile"
}
\description{
  Lists the 13 possible keywords one can use in the file \code{\link{paramfile}}.
}
\usage{data(listparam)}
\format{
  The 13 keywords are as follows and need to start a line in file  \code{\link{paramfile}}:
  \itemize{
    \item \emph{\code{theta_1}}\cr
    \item \emph{\code{M_present}}\cr
    \item \emph{\code{theta_2}}\cr 
    \item \emph{\code{theta_A}}\cr
    \item \emph{\code{T_split}}\cr
    \item \emph{\code{T_change}}\cr
    \item \emph{\code{M_change}}\cr
    \item \emph{\code{rho}}\cr
    \item \emph{\code{nregion}}s\cr
    \item  \emph{\code{type}}\cr
    \item  \emph{\code{stats}}\cr
    \item \emph{\code{howmany}}\cr
    \item \emph{\code{parallel}}
  }
}
\seealso{
  The list \code{\link{listparam}} is used as input by the functions
 \code{\link{order_param}} and
  \code{\link{check_param}} called from \code{\link{simulate_data}} and
  \code{\link{input_EstL}}.\cr  
  \code{\link{liststats}}, \code{\link{msout}}, \code{\link{regfile}},
  \code{\link{paramfile}} are other data and file examples used
  by the functions \code{\link{simulate_data}} and/or
  \code{\link{input_EstL}}.\cr
}
\references{
  See \code{References} section in \code{\link{EstLprep-package}}.
}
\examples{
data(listparam)
## The list of possible keywords in paramfile
listparam
}
\keyword{datasets}
