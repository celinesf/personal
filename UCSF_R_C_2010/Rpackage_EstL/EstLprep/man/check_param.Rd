\name{check_param}
\alias{check_param}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Check of the parameter values
}
\description{
  The function \code{\link{check_param}} checks that the information on
  the parameter values and prior distributions are specified correctly
  in the file \code{\link{paramfile}} are in the correct format to simulate data
  under or to estimate parameters of IM or I2M models.
  \code{\link{simulate_data}} (or  \code{\link{input_EstL}}). \cr  
  This function is called within the \R functions
  \code{\link{simulate_data}} and \code{\link{input_EstL}}. 
}
\usage{
check_param(param, paramfile, listparam, type=NA, summaries=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{
    The matrix \eqn{9} (or \eqn{11}) \eqn{\times} \eqn{3} of parameter
    values for the IM or I2M model and other information
    recorded from the file \code{\link{paramfile}} and ordered with
    \code{\link{order_param}}. 
    \describe{
      \item{\strong{1. IM / I2M MODEL PARAMETERS}}{
        \code{param[1:7,]} contain the information on seven  I2M
        parameters listed below. \cr
        Either one value is recorded (fixed parameter) or two floating
        numbers and an integer (describing a grid of parameter values): \cr 
        "\eqn{\Theta_l}\code{ }\eqn{\Theta_u}\code{ }\eqn{\Theta_n}". 
        \describe{
          \item{\strong{a. REQUIRED PARAMETERS}}{
            \describe{
              \item{\code{param[1,]=}\eqn{\theta_1}}{ 
                records the generational population mutation
                rate per bp for the default population (required) as
                specified with the keyword \emph{\code{theta_1}} in the
                file \code{\link{paramfile}}.
              }% end theta 
              \item{\code{param[2,]=}\eqn{M_p}}{
                records the number of migrants exchanged
                each generation by the TWO populations at present  as
                specified with the keyword \emph{\code{M_present}} in the
                file \code{\link{paramfile}}.
              }% end M_p
            }% end descibed req arg
          }
          \item{\strong{b. OPTIONAL DEMOGRAPHIC PARAMETERS}}{
            \describe{
              \item{\code{param[3,]=}\eqn{\theta_2}}{
                records the generational population mutation rate 
                per bp for population 2 as
                specified with the keyword \emph{\code{theta_2}} in the
                file \code{\link{paramfile}}.
              }% end theta2 
              \item{\code{param[4,]=}\eqn{\theta_A}}{
                records the ancestral population mutation rate per bp per
                generation as
                specified with the keyword \emph{\code{theta_A}} in the
                file \code{\link{paramfile}}.
              }% end thetaA
              \item{\code{param[5,]=}\eqn{T_s}}{
                records the split time in unit of \eqn{4N_1}
                generations between the TWO populations as
                specified with the keyword \emph{\code{T_split}} in the
                file \code{\link{paramfile}}. 
              }% end T_split
              \item{\code{param[6,]=}\eqn{\epsilon}}{
                records the ratio \eqn{\epsilon=T_c/T_s} as
                specified with the keyword \emph{\code{T_change}} in the
                file \code{\link{paramfile}}. 
              }% end T_change
              \item{\code{param[7,]=}\eqn{M_c}}{
                records the number of migrants exchanged
                each generation by the TWO populations since the split until the time
                of gene flow rate change (i.e., \eqn{T_c<t<T_s}) as
                specified with the keyword \emph{\code{M_change}} in the
                file \code{\link{paramfile}}.
              }% end M_c
            }% en descibe optional demographic parameters
          }% end item optio demo
        }%end decribe demo
      }% end item demo
      \item{\strong{2. OPTIONAL NUISANCE PARAMETER }}{
        \describe{
          \item{\code{param[8,]=}\eqn{\rho}}{
            records the parameters on the intra-region
            recombination rate as
            specified with the keyword \emph{\code{rho}} in the
            file \code{\link{paramfile}}.
          }% end item rho
        } % end describe rho
      }% end item nuisance
      \item{\strong{3. OPTIONAL OTHER PARAMETERS}}{
        \describe{
          \item{\code{param[9,1]=}\eqn{R}}{
            records
            the number of independent genomic regions considered  as
            specified with the keyword \emph{\code{nregions}} in the
            file \code{\link{paramfile}}.
          }
        }% end extra
      }% end item other
      \item{\strong{4. OPTIONAL PARAMETERS SPECIFIC TO FUNCTION
          INPUT_ESTL}}{
        \describe{
          \item{\code{param[10,1]=}\eqn{H}}{
            recordsthe number of data sets to simulate per set of parameters
            (i.e., per grid point) to estimate the likelihood of the data given
            the set of model parameters  as
            specified with the keyword \emph{\code{howmany}} in the
            file \code{\link{paramfile}}.
          }% end howmany
          \item{\code{param[11,1]=}\eqn{J}}{  specifies
            records the number of jobs to run in parallel to perform the
            estimation of the posterior distribution of the parameters of the
            extension of the isolation-with-migration model  as
            specified with the keyword \emph{\code{parallel}} in the
            file \code{\link{paramfile}}.
          }% end parallel
        }% end dec
      }% end item stepc
    }% end descibe
  }
  \item{paramfile}{
    The name of the file containing the information in the format
    described in \code{\link{paramfile}}. 
  }
  \item{listparam}{
    The list of possible keywords that can be used in the file
    \code{\link{paramfile}}. See \code{\link{listparam}} for the full list.
  }
  \item{type}{
    The matrix \eqn{1} \eqn{\times} \eqn{4} recording the values
    specified with the keyword \emph{\code{type}} in
    \code{\link{paramfile}} and pre-processed with \code{\link{order_param}}. 
    \describe{
      \item{\code{type[1,1]}=0 or 1}{
        for haplotype and genotype data, respectively. }
      \item{\code{type[1,2]}=1 or 2}{
        for phased and unphased data, respectively. }
      \item{\code{type[1,3:4]}}{
        the mininum and maximum distance between
        pairs of SNPs used to calculate measures of LD and
        recombination, respectively.
      }
    }
  }% end type
  \item{summaries}{
    The string of character listing the summary statistics considered
    (used when called by \code{\link{input_EstL}} only) and pre-processed with
    \code{\link{order_param}}. \cr
    See the full list of summary statistics keywords in
    \code{\link{liststats}}. 
  }
}% end argument


\value{
  \code{\link{check_param}} may output "WARNING"s about errors in the
  file \code{\link{paramfile}} that were replaced by default values or
  other in the output structure. \cr
  If it outputs "PROBLEM", then \code{\link{check_param}} could not
  resolve the issue and output \code{$ok=0}.    
  \item{\code{$ok}}{
    \code{$ok} takes the value \code{"1"} if the format of the input
    file \code{paramfile}  has the format/information required. \cr 
    Otherwise, \code{$ok} takes the value \code{"0"} to inform the
    function \code{\link{simulate_data}} or \code{\link{input_EstL}}
    that it should stop. \cr 
    In this later case, \code{\link{check_param}} outputs an error
    message explaining what part of the format/information is
    incorrect.
  }
  \item{\code{$param}}{
    The same matrix as provided with argument \code{param} but with:
    \itemize{
      \item The values recorded from the file \code{\link{paramfile}} if
      correct.
      \item Different values from those  recorded from the file
      \code{\link{paramfile}} is there were some problems.
      \item Default values if not specified in the file
      \code{\link{paramfile}}.
    }
  }
  \item{\code{$type}}{
    The same matrix as provided with argument \code{type} if
    correct or the default values "0 1 5 10" otherwise.
  }% end type
  \item{\code{$summaries}}{
    The  same string of character listed with argument \code{summaries}
    if correct or the default list
    "S_1 S_2 S_sl S_sh F_st D_1 D_2 D_star1 D_star2" otherwise. 
  }
}% end Values
\references{
  See \code{References} section in \code{\link{EstLprep-package}}.
}
\author{
  Celine Becquet - \email{celine.becquet@gmail.com}.
}
\note{
  ATTENTION: It is the user's responsibility to mind the following
  restrictions. \code{\link{check_param}} requires that:
  \itemize{
    \item   Either the parameter value
    (i.e. \code{param[1,1]=}\eqn{\theta_1>0}) or the prior distribution
    for \eqn{theta_1} (see above for the restrictions) be specified in
    \code{param[1,]}. \cr
    \item  \eqn{M_p>0} when fixed or
    \eqn{0\le}\eqn{M_pl}\eqn{\le}\eqn{M_pu} and \eqn{M_pn>0} when
    variable (or estimated) be specified in \code{param[2,]} in case of an island
    model. \cr 
    \item  \code{param[j,i]=}\eqn{theta_i}\eqn{\ge}\eqn{0} for any \eqn{i} \eqn{\in}
    \eqn{{2,A}} and \eqn{j} \eqn{\in} \eqn{{2,3}}. \cr 
    \item   \eqn{0}\eqn{\le}\eqn{T_c}\eqn{\le}\eqn{T_s} when specified and fixed
    in \code{param[5,1]} and \code{param[6,1]}.\cr
    \item  \code{param[10,1]=}\eqn{H>0} when specified. \cr
    \item   \code{param[11,1]=}\eqn{J>0} when specified.\cr
  }
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  The functions  \code{\link{simulate_data}} and
  \code{\link{input_EstL}} calls the function
  \code{\link{check_param}} after calling the function \code{\link{order_param}}.
  \code{\link{check_param}} calls the function
  \code{\link{error_message}}.\cr
  See the file \code{\link{paramfile}} for an example of the format
  expected by \code{\link{check_param}}.\cr
  Other functions to check the format of input files: \cr
  \code{\link{get_rho}} and \code{\link{check_reg}}. \cr
  Lists of definitions of the symbols and parameters mentioned in
  this file are found in \code{\link{EstLprep-package}}.
}
\examples{
data(Example_files) # download the data

## The list of possible keywords in paramfile
listparam

## Prep for the function check_param
write.table(file="est", x=paramfile, row.name=FALSE, col.names=FALSE, quote=FALSE)
param=scan("est", comment.char="#", what=c("c", "n", "n", "n"), allowEscapes=FALSE, fill=TRUE, sep="\n", strip.white=TRUE, quiet=TRUE) 
vparam=order_param(param, listparam) 
vparam

## Case with no errors.
check_param(param=vparam$vparam, paramfile="est", listparam=listparam,type=vparam$type,summaries=vparam$summaries)

## Case with warning
vparam$vparam[3, 1]=0.001 # theta_2l==theta_2u
check_param(param=vparam$vparam, paramfile="est", listparam=listparam,type=vparam$type,summaries=vparam$summaries)

## Case with error in defining the prior of theta_2.
vparam$vparam[3, 1]=0.002 # theta_2l>theta_2u
check_param(param=vparam$vparam, paramfile="est", listparam=listparam,type=vparam$type,summaries=vparam$summaries)

# Clean up the directory.
unlink("est")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{error}
\keyword{print}
