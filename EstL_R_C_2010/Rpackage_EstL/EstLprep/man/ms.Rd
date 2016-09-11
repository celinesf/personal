\name{ms}
\alias{ms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate ms-like random samples
}
\description{
The function \code{\link{ms}} calls Hudson's C program \code{ms} to simulate
independent data sets for a genomic region of interest.\cr  
\code{\link{ms}} uses as input an \code{ms}-like command 
line defined by the user and output the data simulated by \code{ms}.
}
\usage{
ms(cmdline)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cmdline}{
The string of character with all the tags and values required
by the C program \code{ms} (cf.   \url{http://pps-spud.uchicago.edu/~cbecquet/msdoc.pdf}). \cr
This string   is of the same format as output by the function \code{\link{cmd_line}}.\cr
For usage "\code{ms(".\ms")}".\cr
ATTENTION:  \code{cmdline} needs to have a minimum
of FIVE "words": \cr
"\code{./ms} \eqn{n} \eqn{H} \code{-t} \eqn{\theta_1}".
\tabular{ll}{
\eqn{n} \tab : The total number of chromosomes sampled for the genomic region considered  (\eqn{n>0}). \cr 
\eqn{H} \tab : The number of independent data sets to simulate with
different gene genealogy samples. \cr
\eqn{\theta_1} \tab : The region-specific population mutation rate per
generation for the population considered (\eqn{\theta>0}).  
}
}
}

\value{
The function \code{\link{ms}} outputs the following \code{msframe}:\cr
an \R list  with the information usually output by Hudson's C program \code{ms}.
\item{$cmdline_seeds}{ The vector with two entries: the \code{ms} command
line and the three random seeds used to run the C program \code{ms}.  
}
\item{$segsites}{ The vector of \eqn{H} integers \eqn{S_h}, the number of
segregating sites found in the simulated data set  \eqn{h \in [1, H]}.
} 
\item{$positions}{ The vector of \eqn{H} lists of
positions of the segregating sites found in the \eqn{H} independent data
sets.  
\tabular{ll}{ 
\code{$positions[[}\eqn{h}\code{]]} \tab : The list
of positions of the \eqn{S_h} segregating sites found in the simulated
data set \eqn{h}. \cr
\code{$positions[[}\eqn{h}\code{]][}\eqn{s}\code{]} \tab :
The retalive position (i.e., \eqn{p \in [0,1]})  of the segregating site \eqn{s \in [1, S_h]}.  }
}  
\item{$haplotypes}{ The matrix \eqn{H \times n} of haplotypes simulated for the \eqn{H}
independent data sets. \cr 
\tabular{ll}{
\code{$haplotypes[}\eqn{h,}\code{]} \tab : The \eqn{n} haplotypes
simulated for the data set \eqn{h}. \cr 
\code{$haplotypes[}\eqn{h,i}\code{]} \tab : The
string of \eqn{S_h} characters (\code{"0"} and \code{"1"} for ancestral and derived, respectively) corresponding to the alleles found at the
\eqn{S_h} segregating sites in the data set \eqn{h} (\eqn{i \in [1,n]}). }
}
}
\references{
  See \code{References} section in \code{\link{EstLprep-package}}.
}
\author{
  Celine Becquet - \email{celine.becquet@gmail.com}. 
}
\note{
ATTENTION: It is the user's responsibility to mind the following
restrictions: \cr 
-> The user needs to be comfortable with the use of the
C program \code{ms} and needs to be familiar with the documentation at
\url{http://pps-spud.uchicago.edu/~cbecquet/msdoc.pdf}. \cr 
-> In general, mind the restrictions listed in
\code{\link{cmd_line}}.}

\seealso{
See the documentation of Dick Hudson's program
\code{ms}:
\url{http://pps-spud.uchicago.edu/~cbecquet/msdoc.pdf}. \cr 
The  function \code{\link{cmd_line}} generates inputs for the function
\code{\link{ms}}. \cr
The function  \code{\link{simulate_data}} calls the function \code{\link{ms}}. \cr 
The following functions input, output
and/or manipulate the object output by the function \code{\link{ms}}: \cr 
\code{\link{get_loci}}, \code{\link{get_ms_output}},
\code{\link{out_ms_loci}}, \code{\link{output_ms}},
\code{\link{stats_pop}} and \code{\link{write_ms_output}}. \cr 
The lists of definitions of the symbols and parameters
mentioned in this file are found in \code{\link{EstLprep-package}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
msline="./ms 8 2 -t 1 -r 0.1 9999 -I 2 4 4 10 -ej 1 1 2 -eN 1 3"
ms(msline)

# ms("./ms") # for ms usage
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{programming}