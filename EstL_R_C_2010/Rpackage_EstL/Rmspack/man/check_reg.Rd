\name{check_reg}
\alias{check_reg} 
\title{
  Function to check the format of the files with the information on the
  genomic regions
  and multiple loci. 
}
\description{
  The function \code{\link{check_reg}} checks that the information on the
  independent genomic regions described in the file \code{regfile} and
  the information on the loci for the multi-locus genomic regions in the
  file \code{locifile} are in the correct format to simulate data with
  the function \code{\link{simulate_data}} or estimate models of
  Isolation-with-migration and possible extensions with the function
  \code{\link{estimate_IMc}}. \cr 
  This function is called by the \R functions \code{\link{simulate_data}}
  and \code{\link{estimate_IMc}}. 
}

\usage{
check_reg(nregions, info_region, info_loci, locifile)
}

\arguments{
  \item{nregions}{
    The number of independent genomic regions considered, \eqn{R}. \cr 
    \eqn{R} is specified with the keyword \emph{\code{nregions}} in the
    file \code{paramfile} (see \code{\link{param_sim}} and
    \code{\link{param_est}}). \cr 
  }
  \item{info_region}{
    The matrix \eqn{R} \eqn{\times} \eqn{10} of information on the \eqn{R}
    independent genomic regions. \cr 
    Each genomic region is described by ten values:
    \tabular{ll}{ 
      \eqn{r} \tab 
      : The genomic region number, \eqn{r} \eqn{\in} \eqn{[1,R]}. \cr
      \tab \eqn{R} is the number of independent genomic regions
      specified with the argument \emph{\code{nregions}}. \cr
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
      \code{\link{param_est}}).\cr
      \tab In this case, \eqn{w} is the scaled 
      sex-averaged region-specific population recombination rate per
      bp, i.e., for an X-linked locus \eqn{c_o} is the female
      recombination rate and \eqn{\beta=0.5} so that
      \eqn{\beta*\rho_o=2N_1*c_o}. \cr 
      \tab -- If an estimate of the region-specific recombination rate per bp
      is available for each region from pedigree analysis, \eqn{c_o},
      set \eqn{w=\beta*c_o} to incorporate this knowledge in the simulation or estimation
      (with \code{"}\emph{\code{rho}} \code{2"} in the file
      \code{paramfile}, see \code{\link{param_sim}} and
      \code{\link{param_est}}). \cr
      \tab In this case, \eqn{w} is the scaled
      sex-averaged region-specific recombination rate per bp, i.e., for an
      X-linked locus, \eqn{c_o} is the estimated female recombination rate
      so the scaled sex-averaged recombination rate is \eqn{\beta*c_o=0.5c}. \cr 
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
    }Where: 
    \tabular{ll}{ 
      \eqn{\rho=4N_1*c}\tab
      : The genomic average population intra-region recombination rate per bp per
      generation. \cr  
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
      analysis.\cr
      \eqn{\beta} \tab
      : The ratio of the region-specific population recombination
      rate per bp over \eqn{\rho=4N_1*c} for the genomic region
      considered.\cr
      \tab \eqn{\beta=}\code{"1"} (\code{"0.5"} in
      \emph{Drosophila}) for autosomal region, \code{"0.5"} for X- and
      \code{"0"} for Y- and mtDNA-linked region.
    } See \code{\link{info_region}} for further details.
  }
  \item{info_loci}{
    The matrix \eqn{\Sigma{Y_r}} \eqn{\times} \eqn{6} of information on the
    loci for the multi-locus genomic regions (as specified with \eqn{Y_r>1}
    for the multi-locus region \eqn{r} in \code{info_region[}\eqn{r}\code{]$V10}). \cr
    This matrix is empty unless \eqn{\Sigma{Y_r}>0}. \cr 
    Each locus is described by six values:
    \tabular{ll}{
      \eqn{r} \tab 
      : The multi-locus genomic region number, \eqn{r} \eqn{\in} \eqn{[1,R]}, 
      that this locus is part of. \cr 
      \tab ATTENTION: The multi-locus genomic region numbers, \eqn{r}, need
      to be in order, i.e., the information for the loci for the first
      multi-locus genomic region (\eqn{a}) are from line 2 to \eqn{Y_a+1}, 
      for the second multi-locus genomic region (\eqn{b}), the loci
      information starts on line \eqn{Y_a+2} to \eqn{Y_a+Y_b+1} \dots \eqn{a}
      and \eqn{b} \eqn{\in} \eqn{[1,R]}. \cr 
      \eqn{y} \tab 
      : The locus number, \eqn{y} \eqn{\in} \eqn{[1,Y]}. \cr 
      \tab \eqn{Y} is the total number of loci spanning the multi-locus
      genomic region \eqn{r} as specified in the matrix \code{info_region}. \cr 
      \tab ATTENTION: The loci numbers, \eqn{y}, need to be in order, i.e., 
      information for locus 1 of the first multi-locus genomic region \eqn{r}
      is on line 2, locus 2 on line 3 \dots locus \eqn{Y} on line
      \eqn{Y+1}. \cr 
      \eqn{n_1y} \tab 
      : The sample size from population 1 for the locus \eqn{y}.\cr 
      \tab  \eqn{n_1y}\eqn{\le}\eqn{n_1}, where \eqn{n_1} for the multi-locus
      genomic region \eqn{r} is specified in the matrix \code{info_region}. \cr 
      \eqn{n_2y} \tab 
      : The sample size from population 2 for the locus \eqn{y}.\cr 
      \tab  \eqn{n_2y}\eqn{\le}\eqn{n_2}, where \eqn{n_2} for the multi-locus
      genomic region \eqn{r} is specified in the matrix \code{info_region}. \cr 
      \eqn{z_sy} \tab 
      : The start position of the locus \eqn{y} in bp.\cr 
      \tab  \eqn{z_s1}\eqn{\ge}\eqn{z_s}, where \eqn{z_s} is the start position of the
      multi-locus genomic region \eqn{r} specified in the matrix
      \code{info_region}. \cr 
      \eqn{z_ey} \tab 
      : The end position of the locus \eqn{y} in bp.\cr 
      \tab  \eqn{z_eY}\eqn{\le}\eqn{z_e} where \eqn{z_e} is the end position
      of the multi-locus genomic region \eqn{r} specified in the matrix
      \code{info_region}. 
    } See \code{\link{info_loci}} for further details. 
  }
  \item{locifile}{
    The name of the file with the information on the loci for the
    multi-locus genomic regions (as specified with \eqn{Y>1} in the file
    \code{regfile}). \cr
    By default the file name is \code{"info.loc"}. 
  }
}

\value{
  The function \code{\link{check_reg}} outputs \code{"1"} if the
  format of the matrices \code{info_region} (and \code{info_loci}) (from
  the input file \code{regfile} (and \code{locifile}))
  for the calling function (either \code{\link{simulate_data}} or
  \code{\link{estimate_IMc}}) have the format/information required. \cr
  Otherwise, \code{\link{check_reg}} informs the calling function that
  it should stop and outputs an error message explaining what part of
  the format/information is incorrect: \cr 
  \tabular{ll}{
    e.g.:\tab \code{[1] "Error message on independent (potentially multi-locus) genomic regions."}\cr
    \tab \code{[1] 0}
  }
}

\note{
  \item{ATTENTION: It is the user's responsibility to mind the following
    restrictions:}{ 
    \item{In the matrix \code{info_region}: }{
      -> The independent genomic region numbers, \eqn{r}, need to be in
      order, i.e., region 1 starts on line 2, region 2 on line 3 \dots region
      \eqn{R} on line \eqn{R+1}. \cr 
      -> Reasonable values need to be specified for each genomic region. 
    }
    \item{In the matrix \code{info_loci}:}{
      -> The loci numbers, \eqn{y}, need to be in order, i.e., information
      for locus 1 of the first multi-locus genomic region is on line 2, locus
      2 on line 3 \dots locus \eqn{Y} on line \eqn{Y+1}. \cr 
      -> The multi-locus genomic region numbers, \eqn{r}, need to be in
      order, i.e., the information for the loci for the first multi-locus
      genomic region (\eqn{a}) are from line 2 to \eqn{Y_a+1}, for the second
      multi-locus genomic region (\eqn{b}), the loci information starts on
      line \eqn{Y_a+2} to \eqn{Y_a+Y_b+1} \dots \eqn{a} and \eqn{b} are
      \eqn{\in} \eqn{[1,R]}. \cr 
      -> Reasonable values need to be specified for each locus. 
    }
  }
}

\author{
  Celine Becquet - \email{celine.becquet@gmail.com}. 
}

\seealso{
  The functions \code{\link{simulate_data}} and
  \code{\link{estimate_IMc}} call the function \code{\link{check_reg}}, 
  which in turn calls the function \code{\link{error_message}}. \cr 
  Other functions to check the format of input files: \cr
   \code{\link{get_rho}}  and \code{\link{check_param}}. \cr
  Lists of definitions of the symbols and parameters mentioned in
  this file are found in \code{\link{Rmspack-package}}.
}

\examples{
### Write the files of information in the local directory.
data(simulation_files) # download the data

write.table(file="regions", x=info_region, row.name=FALSE, col.names=FALSE, quote=FALSE)
# Creates the file "regions" containing: 
read.csv("regions",header =FALSE,sep="", comment="#")

write.table(file="loci", x=info_loci, row.name=FALSE, col.names=FALSE, quote=FALSE)
# Creates the file "loci" containing: 
read.csv("loci",header =FALSE,sep="", comment="#")

### Create the inputs for the function.
reg=read.table("regions", skip=1, fill=TRUE)
loci=read.table("loci", skip=1, fill=TRUE)

## Case with no errors.
check_reg(nregions=4, info_region=reg, info_loci=loci, locifile="loci")

## Case with error in the # of regions
reg[1, 1]=2
check_reg(nregions=4, info_region=reg, info_loci=loci, locifile="loci")

# Clean up the directory.
unlink(c("regions", "loci"))
}

\keyword{error}
\keyword{print}
