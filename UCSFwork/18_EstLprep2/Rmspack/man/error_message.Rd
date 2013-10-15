\name{error_message}
\alias{error_message}

\title{
  Function to output error messages.
}

\description{
  The function \code{\link{error_message}} pastes together and prints
  strings of characters, value of parameters and other relevant
  information into an error message. \cr 
  This function is called by the \R functions checking the format of data
  files. 
}

\usage{
error_message(error_message, separator)
}

\arguments{
  \item{error_message}{
    The vector of strings of characters and values of \R objects to
    output into an error message. 
  }
  \item{separator}{
    The vector of characters used as argument \code{sep} of the
    function \code{\link{paste}}. 
  }
  ATTENTION: The function \code{\link{error_message}} requires
  \code{length(error_message)<=length(separator)-1}. 
} 

\value{
  The function \code{\link{error_message}} outputs a string of
  characters with the following information:\cr
  \code{"}\code{error_message[1]separator[1]error_message[2]separator[2] \dots separator[length(separator)]error_message[length(error_message)]}\code{"}.
}

\note{
  \item{ ATTENTION: It is the user's responsibility to mind the following
    restrictions:}{
    The function \code{\link{error_message}} requires
    \code{length(error_message)<=length(separator)-1}.
  }
}

\author{
  Celine Becquet - \email{celine.becquet@gmail.com}. 
}

\seealso{
  The function \code{\link{error_message}} is called by the
  functions:\cr 
  \code{\link{get_rho}}, \code{\link{check_param}}, 
  \code{\link{check_reg}}, and \code{\link{get_ms_output}}. 
}

\examples{
### Create the inputs for the function.
theta=5
paramfile="simulation.par"
error_message=c("PROBLEM: The required argument theta=",theta," was not found in the file of argument ", paramfile, "");
separator=c("'", "'", "'", "'");

error_message(error_message=error_message, separator=separator);
}

\keyword{error}
\keyword{print}
