\name{regfile}
\alias{regfile}
\docType{data}
\title{
  Example file of information on the genomic regions
}
\description{
  The file \code{\link{regfile}} provides an example of the file
  of information on the independent genomic
  regions considered in the functions
  \code{\link{simulate_data}} and 
  \code{\link{input_EstL}}.
}
\usage{
data(regfile)
}
\format{
  The first line of the file is ignored (in the example it contains the
  header). The information on the
  independent genomic regions starts on the SECOND line and
  separated by space or tabulations. \cr 
  Each genomic region is described by ten values:
  \describe{
    \item{\eqn{r}}{
      The genomic region number, \eqn{r} \eqn{\in} \eqn{[1,R]}. \cr 
      \eqn{R} is the number of independent genomic regions considered
      as specified with the keyword \emph{\code{nregions}} in the file
      \code{\link{paramfile}}.
    }
    \item{\code{Region name}}{
      The name of the genomic region, which should contain at least
      ONE non-numerical character. 
    }
    \item{\eqn{x}}{
      The inheritance scalar (i.e.,
      \code{"1"} for autosomal region, \code{"0.75"} for X- and \code{"0.5"}
      for Y- and mtDNA-linked region)
    }
    \item{\eqn{v}}{
      The mutation rate scalar (which can be
      estimated e.g., from divergence data).
    }
    \item{\eqn{w}}{
      The recombination scalar.
      \describe{
        \item{\eqn{w=\beta}}{
          The ratio of the region-specific population recombination
          rate per bp over the genomic average \eqn{\rho}
          per bp for the genomic region considered.\cr
          In this case, the region-specific population
          recombination rates are obtained with \eqn{\rho*w(Z-1)}, where
          \eqn{\rho} is specified in the file \code{\link{paramfile}} with
          "\emph{\code{rho}} \eqn{\rho}".
        }
        \item{\eqn{w=\beta\rho_o}}{
          The scaled sex-averaged region-specific population
          recombination rate per bp, assuming that an estimate of the
          region-specific population recombination rate per bp,\eqn{\rho_o}, is
          KNOWN from linkage disequilibrium analysis.\cr
          In this case, set  "\emph{\code{rho}} \code{1}"  in the
          file \code{\link{paramfile}} and the region-specific population
          recombination rates are obtained with and fixed to 
          \eqn{w*(Z-1)}.          
        }
        \item{\eqn{w=\beta*c_o}}{
          The scaled sex-averaged region-specific recombination rate
          per bp, assuming that an estimate of the region-specific
          recombination rate per bp, \eqn{c_o}, is KNOWN from pedigree
          analysis.  \cr
          In this case, set "\emph{\code{rho}} \code{2} \eqn{\mu}" in the
          file \code{\link{paramfile}} and recombination rates are
          obtained with \eqn{w(Z-1)\theta_1/\mu}.
        }
      }
    }% end w
    \item{\eqn{n_1}}{ 
      The sample size from population 1 for the genomic region.} 
    \item{\eqn{n_2}}{ 
      The sample size from population 2 for the genomic region.} 
    \item{\eqn{z_s}}{ 
      The start position of the genomic region in bp. } 
    \item{\eqn{z_e}}{  
      The end position of the genomic region  in bp. }
    \item{\eqn{Y=1}}{ 
      The number of loci spanning the genomic region.
    }
  }% end off describe
}

\note{
  ATTENTION: It is the user's responsibility to mind the following
  restrictions:
  \itemize{
    \item The information on the first genomic region starts on the SECOND
    line of the file. \cr 
    \item The independent genomic region numbers, \eqn{r}, need to be in
    order, i.e., region 1 starts on line 2, region 2 on line 3 \dots region
    \eqn{R} on line \eqn{R+1}. \cr 
    \item Reasonable values need to be specified for each information and genomic region.
  }
}
\seealso{
  Files in the format of \code{\link{regfile}} are used as input for
  the functions \code{\link{simulate_data}} and
  \code{\link{input_EstL}}. \cr
  The functions \code{\link{check_reg}} and \code{\link{get_rho}} check
  the format of the files like \code{\link{regfile}}.\cr 
 \code{\link{liststats}}, \code{\link{listparam}},  \code{\link{msout}}, \code{\link{paramfile}}
  are other data and file examples used by the functions
  \code{\link{simulate_data}} and/or
  \code{\link{input_EstL}}.\cr
  Lists of definitions of the symbols and parameters mentioned in
  this file are found in  \code{\link{EstLprep-package}}.
}
\references{
  See \code{References} section in \code{\link{EstLprep-package}}.
}
\examples{
### Create the data file.
data(regfile) # download the data

## Write the file of information on the independent-genomic regions in the local directory.
write.table(file="regions", x=regfile, row.name=FALSE, col.names=FALSE, quote=FALSE)

# Clean up the directory.
unlink("regions")
}
\keyword{datasets}
