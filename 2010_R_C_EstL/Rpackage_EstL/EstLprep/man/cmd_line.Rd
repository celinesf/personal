\name{cmd_line}
\alias{cmd_line}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a ms-like command line
}
\description{
  The function \code{\link{cmd_line}} uses information about the sample
  and a demographic model and generates the \code{ms}-like command line
  with the tags and values required as input by the \R function
  \code{\link{ms}}. \cr
}
\usage{
cmd_line(nsam, theta, howmany = 1, rho = 0, structure = 1, time = 0,
migration = 0, extra = "NA", seeds = "NA") 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{}{\strong{REQUIRED ARGUMENTS:}}
  \item{nsam}{
    \code{nsam=}\eqn{n} specifies the total number of chromosomes sampled
    for the genomic region considered. \cr 
    \code{\link{cmd_line}} initiates the \code{ms} command line with the
    string of characters: \cr
    "\code{./msR }\eqn{n}\code{ }\eqn{H}". \cr 
    By default \eqn{H=1} is the number of independent data sets to simulate with
    different gene genealogy samples for the genomic region
    considered. Can be specified with the argument \code{howmany=}\eqn{H>0}. \cr 
  } 
  \item{theta}{
    \code{theta=}\eqn{\theta_1} specifies the generational population
    mutation rate for the genomic region considered.\cr 
    Can be a vector of up to three floating numbers: 
    \code{theta=c(}\eqn{\theta_1}\code{, }\eqn{\theta_2}\code{,
    }\eqn{\theta_A}\code{)} to specify the region-specific generational
    mutation rate for population 1, 2 and their ancestral population,
    respectively.\cr 
    \code{\link{cmd_line}} adds the following strings of characters to the
    \code{ms} command line given each value:\cr
    \itemize{
      \item "\code{-t }\eqn{\theta_1}". \cr
      \item "\code{-n 2 } \eqn{\theta_2/\theta_1}". \cr
		\code{theta[2]} is ignored unless there is population
      structure (as specified with \code{structure[1]=}\eqn{2} and \eqn{n_1>0}
      and \eqn{n_2>0} with the argument \code{structure}). \cr 
      \item "\code{-eN }\eqn{T_s}\code{ }\eqn{N_A/N_1}". \cr
      \code{theta[3]} is ignored unless there is a population split (as
      specified with \eqn{T_s>0} with the argument \code{time}).
    }
  }% end theta
  \item{}{\strong{OPTIONAL ARGUMENTS:}}
  \item{howmany}{
    Use \code{howmany=}\eqn{H} to specify the number of independent gene
    genealogy samples to simulate for the genomic region
    considered.\cr
    \eqn{H=1} by default. \cr 
    \code{\link{cmd_line}} adds "\eqn{H}" to the string of
    characters "\code{./msR }\eqn{n}\code{ }\eqn{H}". \cr  
  }
  \item{rho}{
    Use \code{rho=c(}\eqn{\rho}\code{, }\eqn{Z}\code{)} to specify
    the generational population recombination rate for the genomic
    region considered (of length Z+1 bp). \cr
    By default there is no recombination. \cr 
    \code{\link{cmd_line}} adds the string of characters
    "\code{-r }\eqn{\rho}\code{ }\eqn{Z}" to the
    \code{ms} command line.\cr
  }
  \item{structure}{
    Use \code{structure=c(2, }\eqn{n_1}\code{, }\eqn{n_2}\code{)} to
    specify the sample sizes from population 1 and 2 in a model with
    population structure with at most TWO populations remaining at
    present. \cr 
    By default, \code{structure=1}, i.e., no structure. \cr 
    \code{\link{cmd_line}} adds the string of characters
    "\code{-I 2 }\eqn{n_1}\code{ }\eqn{n_2}\code{ }\eqn{M_p}"
    to the \code{ms} command line. \cr
    \eqn{M_p} is specified with the argument \code{migration}.
  }
  \item{time}{
    Use \code{time} to specify time of events. \cr
    Can take either one value (\code{time=}\eqn{T_s}) and up
    to TWO floating numbers:   \code{time=c(}\eqn{T_s}\code{, }\eqn{T_c}\code{)}.
    \describe{
       \item{\eqn{T_s}}{The split time in unit of \eqn{4N_1} generations
         between the TWO populations. \cr
         By default \eqn{T_s=0} to specify a model without
         population split. \cr 
         \code{\link{cmd_line}} adds the string of characters
         "\code{-ej }\eqn{T_s}\code{ 2 1}" to
         the \code{ms} command line, instructing the \R function
         \code{\link{ms}} to model a population split. \cr 
         \code{time[1]} is ignored unless there is population
           structure (as specified with the argument
           \code{structure=c(2,}\eqn{n_1},\eqn{n_2}\code{)}).
         }
       \item{\eqn{T_c}}{ 
         The time at which the rate of gene flow changed between the two
         populations in unit of \eqn{4N_1} generations. \cr
         By default there is no change of gene flow rate since the
         split. \cr 
         \code{\link{cmd_line}} adds "\eqn{T_c}" in the string of characters
         "\code{-eM }\eqn{T_c}\code{ }\eqn{M_c}". \cr
         \eqn{T_c} is ignored unless
         \itemize{
           \item  There is population structure
           (\code{structure=c(2,}\eqn{n_1},\eqn{n_2}\code{)}     
           specified)
           \item \eqn{0<T_c<T_s}.
           \item There is a population split (i.e., \eqn{T_s>0})
           \item \eqn{M_c \not=M_p} (as specified with the argument
           \code{migration}). 
         }
       }       
    }% end desc
  }% end T
  \item{migration}{
    Use \code{migration} to specify symmetrical gene flow rates between TWO
    populations in the model. \cr
    By default there is no gene flow. \cr 
    Can take either one value (\code{migration=}\eqn{M_p})
    and up to TWO floating numbers: \cr 
    \code{migration=c(}\eqn{M_p}\code{, }\eqn{M_c}\code{)}.\cr
    \describe{
      \item{\eqn{M_p}}{ The number of migrants exchanged each generation by the TWO
        populations at present. \cr
        \eqn{M_p} is ignored unless there is population
        structure (as specified with the argument
        \code{structure=c(2,}\eqn{n_1},\eqn{n_2}\code{)}). \cr        
        \code{\link{cmd_line}} adds "\eqn{M_p}" to the string of characters
        "\code{-I 2 }\eqn{n_1}\code{ }\eqn{n_2}\code{ }\eqn{M_p}". 
        \itemize{
          \item \eqn{M_p} is the symmetrical rate of gene flow between
          two population in an island model (i.e., \eqn{0<t<\inf}) when:\cr
          \itemize{
            \item A population structure is specified.
            \item \eqn{T_s=0} (i.e., specified with the argument
            \code{time} or unspecified).
          }
          \item \eqn{M_p} is the constant symmetrical rate of gene flow
          since the split until present (i.e., \eqn{0<t<T_s}) if there is a
          population split (as specified with \eqn{T_s>0} with the argument
          \code{time}).\cr 
          \item \eqn{M_p} is the constant symmetrical rate of gene flow
          since the time of gene flow rate change until present (i.e.,
          \eqn{0<t<T_c}) if a time at which the gene flow rate changed is
          specified with:\cr
          \itemize{
            \item \eqn{0<T_c<T_s} (as specified with the argument \code{time}). \cr 
            \item \eqn{M_c \not=M_p}. \cr 
          }
        }
      }%end desc Mp
      \item{\eqn{M_c}}{ The number of migrants exchanged each generation by the TWO
        populations since the split until the time of gene flow rate
        change (i.e., \eqn{T_c<t<T_s}). \cr 
        By default there is no change of gene flow rate since the split. \cr 
        \code{\link{cmd_line}} adds "\eqn{M_c}" to the
        string of characters
        "\code{-eM }\eqn{T_c}\code{ }\eqn{M_c}". \cr  
        \eqn{M_c} is ignored unless:
        \itemize{
          \item  There is a population split and \eqn{0<T_c<T_s} (as
          specified with the argument \code{time}). \cr 
          \item \eqn{M_c \not=M_p}.
        }
      }% end Mc
    }% end describ
    By default the migration is assumed symmetrical but asymmetrical
    migration may be specified with the argument \code{extra} instead. 
  }
  \item{extra}{
    Use \code{extra} to add more complications to the model. This argument
    takes a string of characters with the extra \code{ms} tags to specify
    any model not possible with the other arguments. To use this argument, 
    the user should be familiar with the documentation of Hudson's
    \code{ms} (\url{http://pps-spud.uchicago.edu/~cbecquet/msdoc.pdf}). 
  }
  \item{seeds}{
    Use \code{seeds=c(seed1, seed2, seed3)} to specify a random seed to run
    the \R function \code{\link{ms}}. \code{seeds} is a vector of THREE
    integers. \cr
    code{\link{cmd_line}} adds the
      string of characters
       "\code{-s seed1 seed2 seed3}" to the \code{ms} command line.
  }
}
\details{
  The function \code{\link{cmd_line}} allows for up to TWO populations to
  generate an IM model and and I2M model with a change in
  migration rate once since the split. \cr 
  The gene flow rates are assumed symmetrical and constant over time. \cr 
  With the argument \code{extra}, command lines for any possible
  \code{ms} model with any number of populations can be generated (see
  \url{http: //pps-spud.uchicago.edu/~cbecquet/msdoc.pdf}). 
}
\value{
  The function \code{\link{cmd_line}} outputs a string of characters with
  the tags and values required to run the \R function
  \code{\link{ms}}. \cr
  The argument \code{cmdline} of the \R function
  \code{\link{ms}} needs to have a
  minimum of FIVE "words": \cr 
  "\code{./msR }\eqn{n}\code{ }\eqn{H}\code{ -t }\eqn{\theta_1}".\cr
  (The first word \code{"./msR"} can be any string of characters  without white space).
}
\references{
  See \code{References} section in \code{\link{EstLprep-package}}.
}
\author{
 Celine Becquet - \email{celine.becquet@gmail.com}.
}
\note{
  ATTENTION: It is the user's responsibility to mind the following
  restrictions: 
  \itemize{
    \item The user needs to be comfortable with the use of the C program
    \code{ms} and needs to be familiar with the documentation at
    \url{http://pps-spud.uchicago.edu/~cbecquet/msdoc.pdf}. \cr
    \item The mutations and recombiantion rates are given for the
    genomic region of interest (i.e., not per bp).
    \item  \code{\link{cmd_line}} requires that \eqn{n>0} be
    specified with the argument \code{nsam} . \cr
    \item  \code{\link{cmd_line}} requires that \eqn{\theta_1*x*v*Z>0} be
    specified with the argument \code{theta}. \cr
    \item  \code{\link{cmd_line}} requires that \code{theta[i]}\eqn{>0} for any
    \eqn{i} \eqn{\in} \eqn{[1,3]} when specified. \cr 
    \item  If \eqn{H>1*10^5}, \code{\link{cmd_line}} requires that
    \code{howmany} be specified as a string of characters (e.g., 
    \code{howmany="100000"}). \cr 
    \item  \code{\link{cmd_line}} requires that \code{rho} be a vector of TWO
    values. \cr
    \item  If there are TWO populations (as specified with
    \code{structure[1]=}\eqn{2} and \eqn{n_1>0} and \eqn{n_2>0} with the
    argument \code{structure}), \code{\link{cmd_line}} requires that
    \eqn{n=n_1+n_2} be specified with the argument \code{nsam}. \cr
    \item  \code{\link{cmd_line}} requires that \eqn{M_p>0} be specified
    with the argument \code{migration} in case of an island model. \cr 
    \item  \code{\link{cmd_line}} requires
    \eqn{0}\eqn{\le}\eqn{T_c}\eqn{\le}\eqn{T_s} when specified with the
    argument \code{time}.
  }
}

\seealso{
  The functions \code{\link{simulate_data}} and
  \code{\link{input_EstL}} call the function \code{\link{cmd_line}}. \cr 
  The \R function \code{\link{ms}} uses the output of the function
  \code{\link{cmd_line}} as input. \cr
  Lists of definitions of the symbols and parameters mentioned in
  this file are found in \code{\link{EstLprep-package}}.
}
\examples{
### Create the inputs for the function.
np=c(2, 5, 6)      # npop n1 n2.
ns=np[2]+np[3]     # total number of samples.
h=2                # independent gene genealogies.
s=c("1", "2", "3") # Seed numbers.

th=c(5, 10)   # 4Ne*mu*Z for pop1 and pop2, here theta2/theta1=2.
rh=c(1, 1000) # 4Ne*c*Z and Z, recombination rate and size of locus in bp.

T=c(1, 0.5)   # Time of divergence and time of migration rate change in unit of 4Ne generations. 
M=c(5, 0)     # Migration rate until present and eM, the migration rate between T_div and T_change, forward in time. 

## Generate a command for an IM model with a change in gene flow rate 0.5/4N_1 generations ago.
cmd_line(nsam=ns, howmany=h, theta=th, rho=rh, structure=np, seeds=s,time=T, migration=M)  

## Generate a command with a change in effective population size .5/4Ne generations ago. 
cmd_line(nsam=ns, howmany=h, theta=th, extra=c("-eN 0.5 2"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
\keyword{IO}
