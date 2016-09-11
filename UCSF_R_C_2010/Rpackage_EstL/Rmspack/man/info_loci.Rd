\name{info_loci}
\alias{info_loci}
\docType{data}

\title{
  Example of a file with information on the loci from multi-locus genomic
  regions. 
}

\description{
  The file \code{\link{info_loci}} provides an example of the file
  \code{locifile}, which contains information on the loci for the
  multi-locus genomic regions (as specified with \eqn{Y>1} in the file
  \code{regfile}, see \code{\link{info_region}}) required to simulate
  data with the \R function \code{\link{simulate_data}} or estimate
  parameters of extensions of the Isolation-with-migration model with the
  \R function \code{\link{estimate_IMc}}. 
}

\usage{
data(info_loci)
}

\format{
  The file \code{\link{info_loci}} provides an example of the file
  \code{locifile}, which contains information on the loci for the
  multi-locus genomic regions (as specified with \eqn{Y>1} in the file
  \code{regfile}, see \code{\link{info_region}}) required to simulate
  data with the \R function \code{\link{simulate_data}} or estimate
  parameters of extensions of the Isolation-with-migration model with the
  \R function \code{\link{estimate_IMc}}. \cr 
  The first line of the file is ignored. The information on the loci
  starts on the SECOND line. \cr 
  Each locus is described by six values: 
  \tabular{ll}{
    \eqn{r} \tab 
    : The multi-locus genomic region number, \eqn{r} \eqn{\in}
    \eqn{[1,R]}, that this locus is part of. \cr 
    \tab \eqn{R} is the number of independent genomic regions considered
    as specified with the keyword \emph{\code{nregions}} in the file
    \code{paramfile} (see \code{\link{param_sim}} and
    \code{\link{param_est}}). \cr 
    \tab ATTENTION: The multi-locus genomic region numbers, \eqn{r}, need
    to be in order, i.e., the information for the loci for the first
    multi-locus genomic region (\eqn{a}) are from line 2 to \eqn{Y_a+1},
    for the second multi-locus genomic region (\eqn{b}), the loci
    information starts on line \eqn{Y_a+2} to \eqn{Y_a+Y_b+1} \dots \eqn{a}
    and \eqn{b} \eqn{\in} \eqn{[1,R]}. \cr 
    \eqn{y} \tab 
    : The locus number, \eqn{y} \eqn{\in} \eqn{[1,Y]}. \cr 
    \tab \eqn{Y} is the total number of loci spanning the multi-locus
    genomic region \eqn{r} as specified in the file \code{regfile} (see
    \code{\link{info_region}}). \cr 
    \tab ATTENTION: The loci numbers, \eqn{y}, need to be in order, i.e.,
    information for locus 1 of the first multi-locus genomic region
    \eqn{r} is on line 2, locus 2 on line 3 \dots locus \eqn{Y} on line
    \eqn{Y+1}. \cr 
    \eqn{n_1y} \tab 
    : The sample size from population 1 for the locus \eqn{y}.\cr
    \tab \eqn{n_1y}\eqn{\le}\eqn{n_1}, where \eqn{n_1} for the multi-locus
    genomic region \eqn{r} is specified in the file \code{regfile} (see
    \code{\link{info_region}}). \cr 
    \eqn{n_2y} \tab 
    : The sample size from population 2 for the locus \eqn{y}.\cr
    \tab \eqn{n_2y}\eqn{\le}\eqn{n_2}, where \eqn{n_2} for the multi-locus
    genomic region \eqn{r} is specified in the file \code{regfile} (see
    \code{\link{info_region}}). \cr 
    \eqn{z_sy} \tab 
    : The start position of the locus \eqn{y} in bp.\cr
    \tab \eqn{z_s1}\eqn{\ge}\eqn{z_s}, where \eqn{z_s} is the start position
    of the multi-locus genomic region \eqn{r} specified in the file
    \code{regfile} (see \code{\link{info_region}}). \cr 
    \eqn{z_ey} \tab 
    : The end position of the locus \eqn{y} in bp.\cr
    \tab \eqn{z_eY}\eqn{\le}\eqn{z_e} where \eqn{z_e} is the end position
    of the multi-locus genomic region \eqn{r} specified in the file
    \code{regfile} (see \code{\link{info_region}}). 
  }
}

\note{
  \item{ ATTENTION: It is the user's responsibility to mind the following
    restrictions:}{
    -> The information on the first locus of the first multi-locus genomic
    region starts on the SECOND line of the file. \cr 
    -> The loci numbers, \eqn{y}, need to be in order, i.e., information
    for locus 1 of the first multi-locus genomic region is on line 2, locus
    2 on line 3 \dots locus \eqn{Y} on line \eqn{Y+1}. \cr 
    -> The multi-locus genomic region numbers, \eqn{r}, need to be in
    order, i.e., the information for the loci for the first multi-locus
    genomic region (\eqn{a}) are from line 2 to \eqn{Y_a+1}, for the second
    multi-locus genomic region (\eqn{b}), the loci information starts on
    line \eqn{Y_a+2} to \eqn{Y_a+Y_b+1} \dots \eqn{a} and \eqn{b} are \eqn{\in}
    \eqn{[1,R]}. \cr 
    -> Reasonable values need to be specified for each locus.
  }
}

\author{
  Celine Becquet - \email{celine.becquet@gmail.com}. 
}

\seealso{
  Text files in the format of \code{\link{info_loci}} are used as input for
  the functions \code{\link{simulate_data}} and
  \code{\link{estimate_IMc}} in cases where one or more regions are
  spanned by multiple loci. \cr
  The function \code{\link{check_reg}} checks the format of the
  files like \code{\link{info_loci}}.\cr
  Other data file examples used as input by the functions
  \code{\link{simulate_data}} and \code{\link{estimate_IMc}}: \cr 
  \code{\link{info_region}}, \code{\link{data_file}},
  \code{\link{param_sim}} and \code{\link{param_est}}. \cr
  Lists of definitions of the symbols and parameters mentioned in
  this file are found in \code{\link{Rmspack-package}}.
}

\examples{
### Create the data file.
data(simulation_files) # download the data
info_loci

## Write the file of information on the multiple loci in the local directory.
write.table(file="loci", x=info_loci, row.name=FALSE, col.names=FALSE, quote=FALSE)
# Creates the file "loci" containing: 
read.csv("loci",header =FALSE,sep="", comment="#")

# Clean up the directory.
unlink("loci")
}

\keyword{datasets}
