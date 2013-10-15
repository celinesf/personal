\name{info_region}
\alias{info_region}
\docType{data}

\title{ 
  Example of a file with information on regions. 
}
\description{
  The file \code{\link{info_region}} provides an example of the file
  \code{regfile}, which contains information on the independent genomic
  regions required to simulate data with the \R function
  \code{\link{simulate_data}} or estimate parameters of extensions of the
  Isolation-with-migration model with the \R function
  \code{\link{estimate_IMc}}. 
}

\usage{
data(info_region)
}

\format{
  The file \code{\link{info_region}} provides an example of the file
  \code{regfile}, which contains information on the independent genomic
  regions required to simulate data with the \R function
  \code{\link{simulate_data}} or estimate parameters of extensions of the
  Isolation-with-migration model with the \R function
  \code{\link{estimate_IMc}}. \cr 
  The first line of the file is ignored. The information on the
  independent genomic regions starts on the SECOND line. \cr 
  Each genomic region is described by ten values: 
   \tabular{ll}{
    \eqn{r} \tab 
    : The genomic region number, \eqn{r} \eqn{\in} \eqn{[1,R]}. \cr 
    \tab \eqn{R} is the number of independent genomic regions considered
    as specified with the keyword \emph{\code{nregions}} in the file
    \code{paramfile} (see \code{\link{param_sim}} and
    \code{\link{param_est}}). \cr 
    \tab ATTENTION: The independent genomic region numbers, \eqn{r}, need
    to be in order, i.e., region 1 starts on line 2, region 2 on line 3
    \dots region \eqn{R} on line \eqn{R+1}. \cr 
    \emph{\code{Region name}} \tab 
    : The name of the genomic region \eqn{r}, which should contain at least
    ONE non-numerical character. \cr 
    \eqn{x} \tab 
    : The inheritance scalar for the genomic region \eqn{r} (i.e.,
    \code{"1"} for autosomal region, \code{"0.75"} for X- and \code{"0.5"}
    for Y- and mtDNA-linked region). \cr 
    \eqn{v} \tab 
    : The mutation rate scalar for the genomic region \eqn{r} (which can be
    estimated e.g., from divergence data). \cr
    \eqn{w} \tab 
    : The recombination scalar for the genomic region \eqn{r}.\cr
    \tab -- Usually \eqn{w=\beta}, the ratio of the locus-specific
    population recombination rate per bp over \eqn{\rho=4N_1*c}. \cr 
    \tab -- If an estimate of the region-specific population recombination
    rate per bp is available for each region from linkage
    disequilibrium analysis, \eqn{\rho_o=4N_1*c_o}, set
    \eqn{w=\beta*\rho_o} to incorporate this knowledge in the
    simulation or estimation (with \code{"}\emph{\code{rho}} \code{1"} in the
    file \code{paramfile}, see \code{\link{param_sim}} and
    \code{\link{param_est}}). \cr
    \tab In this case, \eqn{w} is the scaled
    sex-averaged region-specific population recombination rate per
    bp, i.e., for an X-linked locus \eqn{c_o} is the female
    recombination rate and \eqn{\beta=0.5} so that
    \eqn{\beta*\rho_o=2N_1*c_o}. \cr 
    \tab -- If an estimate of the region-specific recombination rate per bp
    is available for each region from pedigree analysis, \eqn{c_o},
    set \eqn{w=\beta*c_o} to incorporate this knowledge in the simulation or estimation
    (with \code{"}\emph{\code{rho}} \code{2 }\eqn{\mu}\code{"} in the file
    \code{paramfile}, see \code{\link{param_sim}} and
    \code{\link{param_est}}).\cr
    \tab In this case, \eqn{w} is the scaled sex-averaged
    region-specific recombination rate per bp, i.e., for an X-linked locus,
    \eqn{c_o} is the estimated female recombination rate so the
    scaled sex-averaged recombination rate is
    \eqn{\beta*c_o=0.5c}. \cr 
    \eqn{n_1} \tab 
    : The sample size from population 1 for the genomic region \eqn{r}. \cr 
    \eqn{n_2} \tab 
    : The sample size from population 2 for the genomic region \eqn{r}. \cr 
    \eqn{z_s} \tab 
    : The start position of the genomic region \eqn{r} in bp. \cr 
    \eqn{z_e} \tab 
    : The end position of the genomic region \eqn{r} in bp. \cr 
    \eqn{Y} \tab 
    : The number of loci spanning the genomic region \eqn{r}. \cr
    \tab \eqn{Y=1} means that the whole genomic region \eqn{r} is
    considered as a single locus. \cr 
    \tab \eqn{Y>1} means that multiple loci are considered within the
    genomic region \eqn{r} (the information on the loci should be given in
    the file \code{locifile}, see \code{\link{info_loci}}). 
  }
  Where: 
  \tabular{ll}{ 
    \eqn{\rho=4N_1*c}\tab
    : The genomic average population intra-region recombination rate per
    bp per generation. \cr 
    \eqn{\rho_o=4N_1*c_o} \tab 
    : The estimate of the region-specific population recombination
    rate per bp per generation for the genomic region considered from
    linkage disequilibrium analysis. \cr
    \eqn{N_1} \tab
    : The effective population size in population 1 (the reference
    population). \cr
    \eqn{c} \tab 
    : The genomic generational cross-over rate per bp. \cr
    \eqn{c_o} \tab 
    : The estimate of the region-specific cross-over rate per bp
    per generation for the genomic region considered from pedigree
    analysis. \cr
    \eqn{\beta} \tab
    : The ratio of the region-specific population recombination
    rate per bp over \eqn{\rho=4N_1*c} for the genomic region
    considered.\cr
    \tab \eqn{\beta=}\code{"1"} (\code{"0.5"} in
    \emph{Drosophila}) for autosomal region, \code{"0.5"} for X- and
    \code{"0"} for Y- and mtDNA-linked region. 
  }
}

\note{
  \item{ ATTENTION: It is the user's responsibility to mind the following
    restrictions:}{
    -> The information on the first genomic region starts on the SECOND
    line of the file. \cr 
    -> The independent genomic region numbers, \eqn{r}, need to be in
    order, i.e., region 1 starts on line 2, region 2 on line 3 \dots region
    \eqn{R} on line \eqn{R+1}. \cr 
    -> Reasonable values need to be specified for each genomic region.
  }
}

\author{
  Celine Becquet - \email{celine.becquet@gmail.com}. 
}

\seealso{
  Text files in the format of \code{\link{info_region}} are used as input for
  the functions \code{\link{simulate_data}} and
  \code{\link{estimate_IMc}}. \cr
  The functions \code{\link{check_reg}} and \code{\link{get_rho}} check
  the format of the files like \code{\link{info_region}}.\cr 
  Other data file examples used as input by the functions
  \code{\link{simulate_data}} and \code{\link{estimate_IMc}}: \cr 
  \code{\link{info_loci}}, \code{\link{data_file}},
  \code{\link{param_sim}} and \code{\link{param_est}}. \cr
  Lists of definitions of the symbols and parameters mentioned in
  this file are found in \code{\link{Rmspack-package}}.
} 

\examples{
### Create the data file.
data(simulation_files) # download the data
info_region

## Write the file of information on the independent-genomic regions in the local directory.
write.table(file="regions", x=info_region, row.name=FALSE, col.names=FALSE, quote=FALSE)
# Creates the file "regions" containing: 
read.csv("regions",header =FALSE,sep="", comment="#")

# Clean up the directory.
unlink("regions")
}

\keyword{datasets}
