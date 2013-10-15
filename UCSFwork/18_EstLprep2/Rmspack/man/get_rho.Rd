\name{get_rho}
\alias{get_rho}

\title{
  Function to calculate the region-specific population recombination
  rate. 
}

\description{
  The function \code{\link{get_rho}} checks that the recombination
  scalar \eqn{w} specified for the recombining genomic region \eqn{r}
  in the file \code{regfile} (see \code{\link{info_region}})
  is in the correct format to simulate data with the function
  \code{\link{simulate_data}}.\cr 
  If so \code{\link{get_rho}} outputs the region-specific
  population recombination rate, \eqn{\rho_r=\beta*4N_1c*(Z-1)} for
  the recombining region \eqn{r}. \cr
  This function is called by the \R function \code{\link{simulate_data}}.
}

\usage{
get_rho(param, nregion, Z, W, theta=1)
}

\arguments{
  \item{param}{
    This vector of three value specifying the parameters on the intra-region recombination
    rate.\cr
    \code{param} can have the five following forms:
    \tabular{lllll}{
      \code{param[1]=}\tab \eqn{\rho} \tab \tab \tab
      : This specifies that the genomic average population
      intra-region recombination rate per bp per generation is fixed to the
      value \eqn{\rho=4N_1*c}. \cr
      \code{param[1]=} \tab \code{1} \tab \tab \tab
      : This specifies that an estimate of the region-specific population
      recombination rate per bp, \eqn{\rho_o=4N_1*c_o}, is KNOWN from
      linkage disequilibrium analysis and specified with the parameter
      \eqn{w=\beta*\rho_o} for each independent recombining genomic  region (with the
      argument \code{W}).\cr
      \code{param=c(}\tab \code{2}\code{,}\tab \eqn{\mu}\code{,}\tab \code{NA)}\tab
      : This specifies that an estimate of the region-specific
      recombination rate per bp, \eqn{c_o}, is KNOWN from pedigree
      analysis and specified with the parameter \eqn{w=\beta*c_o}
      for each recombining region (with the
      argument \code{W}). \cr
      \tab \tab \tab \tab This also specifies \eqn{\mu}, an
      independent estimate of the genomic generational mutation rate
      per bp.\cr 
      \code{param=c(}\tab \code{-}\code{1}\code{,}\tab \eqn{1/}\eqn{\lambda}\code{,}\tab \code{NA)} \tab
      : This specifies that the intra-region recombination rate is
      UNKNOWN and the ratio of recombination over mutation
      rate for the genomic region \eqn{r}, \eqn{\alpha=c_r/\mu}, is drawn
      from an exponential distribution with mean \eqn{1/}\eqn{\lambda}. \cr 
      \code{param=c(}\tab \code{-}\code{2}\code{,}\tab \eqn{\nu}\code{,}\tab \eqn{\sigma}\code{)} \tab
      : This specifies that the intra-region recombination rate is
      UNKNOWN and the ratio of recombination over mutation
      rate for the genomic region \eqn{r}, \eqn{\alpha=c_r/\mu}, is drawn
      from an normal distribution with mean
      \eqn{\nu} and standard deviation \eqn{\sigma}. 
    } Where:
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
      analysis. \cr 
      \eqn{\alpha=c_r/\mu} \tab : Drawn from a prior distribution. \cr
      \eqn{c_r} \tab
      : The generational region-specific cross-over rate per bp for the
      genomic region \eqn{r}. \cr
      \eqn{\mu} \tab
      : The genomic generational mutation rate per bp.\cr
      \eqn{w} \tab 
      : The recombination scalar for the genomic region considered
      (as specified with the argument \code{W}).\cr 
      \eqn{\beta} \tab
      : The ratio of the region-specific population recombination
      rate per bp over \eqn{\rho=4N_1*c} for the genomic region
      considered.\cr 
      \tab \eqn{\beta=}\code{"1"} (\code{"0.5"} in
      \emph{Drosophila}) for autosomal region, \code{"0.5"} for X- and
      \code{"0"} for Y- and mtDNA-linked region.
    }
    The information on the parameter \eqn{\rho} can be
    specified with the keyword \emph{\code{rho}} in the file
    \code{paramfile} (see \code{\link{param_sim}}).
  }
  \item{nregion}{
    The multi-locus genomic region number
    considered, \eqn{r} \eqn{\in} \eqn{[1,R]}. \cr
    \eqn{R} is the number of independent genomic regions to
    simulate as specified with the keyword \emph{\code{nregions}} in the
    file \code{paramfile} (see \code{\link{param_sim}}).
  }
  \item{Z}{
    The size \eqn{Z} in bp of the genomic region \eqn{r} (as
    specified with the argument \code{nregion}). \cr
    \eqn{z_s} and \eqn{z_e} (such as \code{Z=}\eqn{Z=z_e-z_s}) are
    specified in the file \code{regfile} (see \code{\link{info_region}}).
  }
  \item{W}{
    The recombination scalar for the genomic region \eqn{r} (as
    specified with the argument \code{nregion}) specified in the
    file \code{regfile}, \eqn{w} (see \code{\link{info_region}}).
    \tabular{l}{
      -- Usually \eqn{w=\beta}, the ratio of the locus-specific
      population recombination rate per bp over \eqn{\rho=4N_1*c}. \cr 
      -- If an estimate of the region-specific population recombination
      rate per bp is available for each region from linkage
      disequilibrium analysis, \eqn{\rho_o=4N_1*c_o}, set
      \eqn{w=\beta*\rho_o} to incorporate this knowledge in the
      simulation (with \code{"}\emph{\code{rho}} \code{1"} in the
      file \code{paramfile}, see \code{\link{param_sim}}). \cr
      In this case, \eqn{w} is the scaled
      sex-averaged region-specific population recombination rate per
      bp, i.e., for an X-linked locus \eqn{c_o} is the female
      recombination rate and \eqn{\beta=0.5} so that
      \eqn{\beta*\rho_o=2N_1*c_o}. \cr 
      -- If an estimate of the region-specific recombination rate per bp
      is available for each region from pedigree analysis, \eqn{c_o},
      set \eqn{w=\beta*c_o} to incorporate this knowledge in the simulation
      (with \code{"}\emph{\code{rho}} \code{2"} in the file
      \code{paramfile}, see \code{\link{param_sim}}).\cr
      In this case, \eqn{w} is the scaled sex-averaged
      region-specific recombination rate per bp, i.e., for an X-linked locus,
      \eqn{c_o} is the estimated female recombination rate so the
      scaled sex-averaged recombination rate is \eqn{\beta*c_o=0.5c}. 
    }
    ATTENTION: It is the user's responsibility to provide a reasonable value
    for the recombination scalar, \eqn{w} specified with the argument
    \code{W}.\cr
  }
  \item{theta}{
    The population mutation rate per bp per generation for
    population 1, \eqn{\theta_1=4N_1*\mu} (specified with the keyword
    \code{"}\emph{\code{theta}} \eqn{\theta_1}\code{"} in the file
    \code{paramfile}, see \code{\link{param_sim}}).\cr
    By default, \code{theta=1}: if this argument is not provided
    \code{\link{get_rho}} provides the unvariable portion of the
    intra-region recombination rate when \code{param[1]} \eqn{\in} \eqn{{-2,-1,2}}
    Where: 
    \tabular{ll}{
      \eqn{N_1} \tab
      : The effective population size in population 1 (the reference
      population so by default and unless specified, \eqn{N_1=N_2=N_A}). \cr 
      \eqn{\mu} \tab
      : The genomic generational mutation rate per bp. 
    } 
  }
}

\value{
  The function \code{\link{get_rho}} outputs error messages regarding
  the format of the recombination scalar \code{W=}\eqn{w} for the
  recombining region \code{nregion=}\eqn{r} as well as the following data
  frame: 
  \item{\code{$ok}}{
    \code{$ok} takes the value \code{"1"} if the
    recombination scalar \code{W=}\eqn{w} for the recombining region
    \code{nregion=}\eqn{r} has the format/information required. \cr 
    Otherwise, \code{$ok} takes the value \code{"0"} to inform the function
    \code{\link{simulate_data}} that it should stop. \cr
    In this later case, \code{\link{get_rho}} outputs an error
    message explaining what part of the format/information is
    incorrect.
  }
  \item{\code{$rho}}{
    The region-specific population recombination rate for the
    genomic region \eqn{r} (as specified with the
    argument \code{nregion}), \code{$rho=}\eqn{\rho_r=\beta*4N_1c*(Z-1)}.
    \tabular{llllll}{
      If\tab \code{param[1]=}\tab \eqn{\rho} \tab \tab \tab
      : The genomic average population
      intra-region recombination rate per bp per generation is fixed to the
      value \eqn{\rho=4N_1*c} . \cr
      \tab \tab \tab \tab \tab \code{$rho=}\eqn{\rho*w*(Z-1)}
      for the genomic region \eqn{r} (here, \code{W=}\eqn{w=\beta}). \cr 
      If\tab \code{param[1]=} \tab \code{1} \tab \tab \tab
      : An estimate of the region-specific population
      recombination rate per bp, \eqn{\rho_o=4N_1*c_o}, is KNOWN from
      linkage disequilibrium analysis and specified with
      \code{W=}\eqn{w=\beta*\rho_o}.\cr
      \tab \tab \tab \tab \tab \code{$rho=}\eqn{w*(Z-1)} for
      the genomic region for the genomic region \eqn{r}. \cr
      If\tab \code{param=c(}\tab \code{2}\code{,}\tab \eqn{\mu}\code{,}\tab \code{NA)}\tab
      : An estimate of the region-specific
      recombination rate per bp, \eqn{c_o}, is KNOWN from pedigree
      analysis and specified with \code{W=}\eqn{w=\beta*c_o}. \cr 
      \tab \tab \tab \tab \tab In this case,
      \code{$rho=}\eqn{\theta_1*w*(Z-1)/\mu} for the genomic region
      \eqn{r}. \cr 
      If\tab \code{param=c(}\tab \code{-}\code{1}\code{,}\tab \eqn{1/}\eqn{\lambda}\code{,}\tab \code{NA)} \tab
      : The intra-region recombination rate is
      UNKNOWN and the ratio of recombination over mutation
      rate for the genomic region \eqn{r}, \eqn{\alpha=c_r/\mu}, is drawn
      from an exponential distribution with mean \eqn{1/}\eqn{\lambda}. \cr 
      \tab \tab \tab \tab \tab
      \code{$rho=}\eqn{w*\alpha*(Z-1)*\theta_1} for the
      genomic region \eqn{r}. \cr 
      If\tab \code{param=c(}\tab \code{-}\code{2}\code{,}\tab \eqn{\nu}\code{,}\tab \eqn{\sigma}\code{)} \tab
      : The intra-region recombination rate is
      UNKNOWN and the ratio of recombination over mutation
      rate for the genomic region \eqn{r}, \eqn{\alpha=c_r/\mu}, is drawn
      from an normal distribution with mean
      \eqn{\nu} and standard deviation \eqn{\sigma}. \cr
      \tab \tab \tab \tab \tab
      \code{$rho=}\eqn{w*\alpha*(Z-1)*\theta_1} for
      the genomic region \eqn{r}. 
    } Where:
    \tabular{ll}{
      \eqn{\theta_1=4N_1*\mu}\tab
      : The population mutation rate per bp per generation for
      population 1 (as specified with the argument \code{theta}).\cr
      \eqn{\rho_r=\beta*(Z}\eqn{-}\eqn{1)*4N_1c}\tab 
      : The region-specific recombination rate per generation for the
      genomic region considered. \cr  
      \eqn{Z} \tab
      : The size in bp of the genomic region considered (as specified
      with the argument \code{Z}).
    }
    ATTENTION: \code{\link{get_rho}} requires that
    \eqn{rho_r/(\theta_1*(Z-1))\le100}. \cr
  }
  e.g.:\cr
  \code{[1] "PROBLEM: Error message on the recombination rate."}\cr
  \code{$ok} \cr
  \code{[1] 0}\cr
  \code{$rho}\cr
  \code{[1] 9999}\cr
}

\note{
  \item{ATTENTION: It is the user's responsibility to provide a reasonable value
    for the recombination scalar, \eqn{w} for a recombining
    genomic region (specified with the argument \code{W}).}{ 
    -> \code{\link{get_rho}} requires that
    \eqn{rho_r/(\theta_1*(Z-1))\le100}. \cr 
  }
}

\author{
  Celine Becquet - \email{celine.becquet@gmail.com}. 
}

\seealso{
  The function \code{\link{simulate_data}} calls the function
  \code{\link{get_rho}}, which in turn calls the function
  \code{\link{error_message}}.\cr
  The function \code{\link{get_rho}} checks the format of the
  files like \code{\link{param_sim}} and \code{\link{info_region}}. \cr 
  Other functions to check the format of input files:\cr
  \code{\link{check_param}} and \code{\link{check_reg}}. \cr
  Lists of definitions of the symbols and parameters mentioned in
  this file are found in \code{\link{Rmspack-package}}.
}

\examples{
### Write the files of information in the local directory.
data(simulation_files) # download the data

write.table(file="sim", x=param_sim, row.name=FALSE, col.names=FALSE, quote=FALSE)
# Creates the file "sim" containing (with comments):
read.csv("sim",header =FALSE,sep="", comment="#")

write.table(file="regions", x=info_region, row.name=FALSE, col.names=FALSE, quote=FALSE)
# Creates the file "regions" containing: 
read.csv("regions",header =FALSE,sep="", comment="#")

### Create the inputs for the function.
listparam=c("theta_1","M_present","theta_2","theta_A","T_split","T_change","M_change","rho", "nregions")
param=scan("sim", comment.char="#", what=c("c", "n","n","n"), allowEscapes=FALSE, fill=TRUE, sep="\n", strip.white=TRUE, quiet=TRUE)
vparam=order_param(param, listparam)
reg=read.table("regions", skip=1, fill=TRUE)

nreg=1
## Case with no errors.
get_rho(param=vparam[8,], nregion=nreg, Z=reg[nreg,9]-reg[nreg,8], W=reg[nreg,5],theta=vparam[1,1])

## Case with error: assumes w=4N_1*c_o
vparam[8,1]=1
get_rho(param=vparam[8,], nregion=nreg, Z=reg[nreg,9]-reg[nreg,8], W=reg[nreg,5],theta=vparam[1,1])

# Clean up the directory. 
unlink(c("sim", "regions"))
}

\keyword{error}
\keyword{print}

