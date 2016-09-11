\name{cmd_line}
\alias{cmd_line}

\title{
  Function to generate ms-like command lines. 
}

\description{
  The function \code{\link{cmd_line}} takes values for the number of
  samples, the region-specific population mutation and recombination
  rates, the parameters of the demographic model to simulate... and
  generates the \code{ms}-like command line with the tags and values
  required as input by the \R function \code{\link{ms}}. \cr 
  This function is called by the \R functions \code{\link{simulate_data}}
  and \code{\link{estimate_IMc}}. 
}

\usage{
cmd_line(nsam, theta, howmany=1, rho=0, structure=1, time=0, migration=0, extra="NA", seeds="NA") 
}

\arguments{
  \item{}{\strong{REQUIRED ARGUMENTS:}}
  \item{nsam}{
    \code{nsam=}\eqn{n} specifies the total number of chromosomes sampled
    for the genomic region considered. \cr 
    ATTENTION: \code{\link{cmd_line}} requires that \eqn{n>0} be
    specified with the argument \code{nsam}. \cr 
    \code{\link{cmd_line}} initiates the \code{ms} command line with the
    string of characters: \cr
    \code{"./msR }\eqn{n}\code{ }\eqn{H}\code{"}. \cr 
    \eqn{H} is the number of independent data sets to simulate with
    different gene genealogy samples for the genomic region
    considered.\cr 
    By default \eqn{H=1} unless it is specified with the argument
    \code{howmany=}\eqn{H>0}. \cr 
    ATTENTION, if there are TWO populations (as specified with
    \code{structure[1]=2} and \eqn{n_1>0} and \eqn{n_2>0} with the argument
    \code{structure}), \code{\link{cmd_line}} requires
    \code{nsam=}\eqn{n=n_1+n_2}. 
  } 
  \item{theta}{
    ATTENTION: \code{\link{cmd_line}} requires that
    \eqn{\theta_1*x*v*Z>0} be specified with the argument
    \code{theta} (i.e., \code{theta}\eqn{\theta_1*x*v*Z}). \cr 
    \code{theta} can be a vector of up to three floating numbers: \cr 
    \code{theta=c(}\eqn{\theta_1*x*v*Z}\code{, }\eqn{\theta_2*x*v*Z}\code{, }\eqn{\theta_A*x*v*Z}\code{)}. 
    \tabular{ll}{
      \eqn{\theta_1*x*v*Z} \tab
      : The region-specific population mutation rate per generation for
      population 1. \cr
      \tab \code{\link{cmd_line}} adds the string of characters:
      \code{"}\code{-t }\eqn{\theta_1*x*v*Z}\code{"} to the
      \code{ms} command line, specifying the region-specific
      generational population mutation rate (required to run the \R
      function \code{\link{ms}}). \cr 
      \eqn{\theta_2*x*v*Z} \tab
      : The region-specific population mutation rate per generation for
      population 2. \cr
      \tab \code{\link{cmd_line}} adds the string of characters
      \code{"}\code{-n 2 }\eqn{N_2/N_1}\code{"} to the \code{ms}
      command line. \cr 
      \tab \code{theta[2]} is ignored unless there is population
      structure (as specified with \code{structure[1]=}\eqn{2} and \eqn{n_1>0}
      and \eqn{n_2>0} with the argument \code{structure}). \cr 
      \eqn{\theta_A*x*v*Z} \tab
      : The region-specific ancestral population mutation rate per
      generation. \cr 
      \tab \code{\link{cmd_line}} adds the string of characters
      \code{"}\code{-eN }\eqn{T_s}\code{ }\eqn{N_A/N_1}\code{"}
      to the \code{ms} command line. \cr 
      \tab \code{theta[3]} is ignored unless there is a population split (as
      specified with \eqn{T_s>0} with the argument \code{time}). 
    }
    Where: 
    \tabular{ll}{
      \eqn{\theta_1=4N_1*\mu} \tab
      : The population mutation rate per bp per generation for
      population 1 (required). \cr 
      \eqn{\theta_2=4N_2*\mu} \tab
      : The population mutation rate per bp per generation for
      population 2. \cr 
      \eqn{\theta_A=4N_A*\mu} \tab
      : The ancestral population mutation rate per bp per
      generation. \cr 
      \eqn{N_1} \tab
      : The effective population size in population 1 (the reference
      population so by default and unless specified, \eqn{N_1=N_2=N_A}). \cr 
      \eqn{N_2} \tab
      : The effective population size in population 2. \cr
      \eqn{N_A} \tab
      : The ancestral effective population size. \cr
      \eqn{\mu} \tab
      : The genomic generational mutation rate per bp. \cr
      \eqn{x} \tab
      : The inheritance scalar for the genomic region considered (i.e., 
      \code{"1"} for autosomal region, \code{"0.75"} for X- and \code{"0.5"}
      for Y- and mtDNA-linked region). \cr 
      \eqn{v} \tab
      : The mutation rate scalar for the genomic region considered
      (which can be estimated e.g., from divergence data). \cr 
      \eqn{Z} \tab
      : The size in bp of the genomic region considered. 
    } ATTENTION: \code{\link{cmd_line}} requires that \code{theta[i]}\eqn{>0}
    for any \eqn{i} \eqn{\in} \eqn{[1,3]} when specified. 
  }% end theta
  \item{}{\strong{OPTIONAL ARGUMENTS:}}
  \item{howmany}{
    Use \code{howmany=}\eqn{H} to specify the number of independent gene
    genealogy samples to simulate for the genomic region
    considered.\cr
    \code{howmany=}\eqn{H=1} by default. \cr 
    \code{\link{cmd_line}} adds \code{"}\eqn{H}\code{"} to the string of
    characters \code{"}\code{./msR }\eqn{n}\code{ }\eqn{H}\code{"}. \cr 
    ATTENTION: If \eqn{H>1*10^5}, \code{\link{cmd_line}} requires that
    \code{howmany} be specified as a string of characters (e.g., 
    \code{howmany="100000"}). 
  }
  \item{rho}{
    Use \code{rho=c(}\eqn{\rho*w*(Z-1)}\code{, }\eqn{Z}\code{)} to specify
    the region-specific population recombination rate per generation. \cr
    By default there is no recombination. \cr 
    \code{\link{cmd_line}} adds the string of characters
    \code{"}\code{-r }\eqn{\rho*w*(Z-1)}\code{ }\eqn{Z}\code{"} to the
    \code{ms} command line, instructing the \R function \code{\link{ms}}
    to generate ancestral recombination graphs instead of gene
    genealogies.\cr
    Where,
    \tabular{ll}{
      \eqn{\rho=4N_1*c} \tab
      : The genomic average population intra-region recombination rate
      per bp per generation. \cr 
      \eqn{c} \tab
      : The generational cross-over rate per bp. \cr
      \eqn{w} \tab
      : The recombination scalar for the genomic region considered.\cr
      \tab \eqn{w=\beta} is the ratio of the region-specific
      population recombination rate per bp over \eqn{\rho=4N_1*c}.\cr
      \tab \eqn{\beta=}\code{"1"} (\code{"0.5"} in
      \emph{Drosophila}) for autosomal region, \code{"0.5"} for X- and
      \code{"0"} for Y- and mtDNA-linked region.
    } ATTENTION: \code{\link{cmd_line}} requires that \code{rho} be a vector
    of TWO values.
  }
  \item{structure}{
    Use \code{structure=c(2, }\eqn{n_1}\code{, }\eqn{n_2}\code{)} to
    specify population structure in a model with at most TWO populations
    remaining at present. \cr
    By default, \code{structure=1}\eqn{1} to specify a model without
    structure. \cr 
    \code{\link{cmd_line}} adds the string of characters
    \code{"}\code{-I 2 }\eqn{n_1}\code{ }\eqn{n_2}\code{ }\eqn{M_p}\code{"}
    to the \code{ms} command line, instructing the \R function
    \code{\link{ms}} to model population structure. \cr 
    \tabular{ll}{
      \eqn{n_1} \tab
      : The sample size from population 1 for the genomic region
      considered. \cr 
      \tab \code{structure[2]} is ignored unless
      \code{structure[1]=}\eqn{2}. \cr 
      \eqn{n_2} \tab
      : The sample size from population 2 for the genomic region
      considered.\cr 
      \tab \code{structure[3]} is ignored unless
      \code{structure[1]=}\eqn{2}. 
    } ATTENTION: \cr
    -> When there is population structure (as specified with
    \code{structure[1]=}\eqn{2} and \eqn{n_1>0} and \eqn{n_2>0} with the
    argument \code{structure}), \code{\link{cmd_line}} requires
    \code{nsam=}\eqn{n=n_1+n_2}. \cr 
    -> \code{\link{cmd_line}} requires that \eqn{M_p>0} be specified
    with the argument \code{migration} in case of an island model (i.e.,
    \eqn{T_s=0} when specified with the argument \code{time} or unspecified). \cr 
  }
  \item{time}{
    Use \code{time} to specify time of events. \cr
    \code{time} can take either one value (\code{time=}\eqn{T_s}) and up
    to TWO floating numbers: \cr 
    \code{time=c(}\eqn{T_s}\code{, }\eqn{T_c}\code{)}. 
    \tabular{ll}{
      \eqn{T_s} \tab
      : The split time in unit of \eqn{4N_1} generations between the TWO
      populations. \cr
      \tab By default \eqn{T_s=0} to specify a model without
      population split. \cr 
      \tab \code{\link{cmd_line}} adds the string of characters
      \code{"}\code{-ej }\eqn{T_s}\code{ 2 1}\code{"} to
      the \code{ms} command line, instructing the \R function
      \code{\link{ms}} to model a population split. \cr 
      \tab \code{time[1]} is ignored unless there is population
      structure (as specified with \code{structure[1]=}\eqn{2} and \eqn{n_1>0}
      and \eqn{n_2>0} with the argument \code{structure}). \cr 
      \eqn{T_c} \tab
      : The time at which the rate of gene flow changed between the two
      populations in unit of \eqn{4N_1} generations. \cr
      \tab By default there is no change of gene flow rate since the
      split. \cr 
      \tab \code{\link{cmd_line}} adds the string of characters
      \code{"}\code{-eM }\eqn{T_c}\code{ }\eqn{M_c}\code{"} to
      the \code{ms} command line, instructing the \R function
      \code{\link{ms}} to change gene flow rate from \eqn{M_p} to \eqn{M_c}
      at time \eqn{T_c}, backward in time. \cr
      \tab * \eqn{T_c} is ignored unless: \cr
      \tab --- There is a population split and \eqn{0<T_c<T_s} (as
      specified with the argument \code{time}). \cr 
      \tab --- \eqn{M_c \not=M_p} (as specified with the argument 
      \code{migration}). 
    } ATTENTION: \code{\link{cmd_line}} requires
    \eqn{0}\eqn{\le}\eqn{T_c}\eqn{\le}\eqn{T_s} when specified. 
  }
  \item{migration}{
    Use \code{migration} to specify symmetrical gene flow rates between TWO
    populations in the model. \cr
    By default there is no gene flow. \cr 
    \code{migration} can take either one value (\code{migration=}\eqn{M_p})
    and up to TWO floating numbers: \cr 
    \code{migration=c(}\eqn{M_p}\code{, }\eqn{M_c}\code{)}. 
    \tabular{ll}{
      \eqn{M_p=4N_1*m_p} \tab 
      : The number of migrants exchanged each generation by the TWO
      populations at present. \cr 
      \tab \eqn{M_p} is ignored unless there is population
      structure (as specified with \code{structure[1]=}\eqn{2} and
      \eqn{n_1>0} and \eqn{n_2>0} with the argument
      \code{structure}). \cr
      \tab \code{\link{cmd_line}} adds \code{"}\eqn{M_p}\code{"} to the
      string of characters
      \code{"}\code{-I 2 }\eqn{n_1}\code{ }\eqn{n_2}\code{ }\eqn{M_p}\code{"}. \cr  
      \tab \eqn{M_p} has different definitions depending of the model
      considered:\cr 
      \tab * \eqn{M_p} is the symmetrical rate of gene flow between
      two population in an island model (i.e., \eqn{0<t<\inf}) as
      specified with:\cr 
      \tab -- A population structure (as specified with
      \code{structure[1]=}\eqn{2} and \eqn{n_1>0} and \eqn{n_2>0} with the
      argument \code{structure}). \cr
      \tab -- Either \eqn{T_s=0} when specified with the argument
      \code{time} or unspecified.\cr
      \tab ATTENTION: \code{\link{cmd_line}} requires that
      \eqn{M_p>0} be specified in case of an island model. \cr
      \tab * \eqn{M_p} is the constant symmetrical rate of gene flow
      since the split until present (i.e., \eqn{0<t<T_s}) if there is a
      population split (as specified with \eqn{T_s>0} with the argument
      \code{time}).\cr 
      \tab * \eqn{M_p} is the constant symmetrical rate of gene flow
      since the time of gene flow rate change until present (i.e.,
      \eqn{0<t<T_c}) if a time at which the gene flow rate changed is
      specified with:\cr 
      \tab -- \eqn{0<T_c<T_s} (as specified with \eqn{T_s>0} and
      \eqn{0<T_c<T_s} with the argument \code{time}). \cr 
      \tab -- \eqn{M_c \not=M_p} (as specified with the argument
      \code{migration}). \cr
      \eqn{M_c=4N_1*m_c} \tab 
      : The number of migrants exchanged each generation by the TWO
      populations since the split until the time of gene flow rate
      change (i.e., \eqn{T_c<t<T_s}). \cr 
      \tab By default there is no change of gene flow rate since the
      split. \cr 
      \tab * \eqn{M_c} is ignored unless: \cr 
      \tab -- There is a population split and \eqn{0<T_c<T_s} (as
      specified with the argument \code{time}). \cr 
      \tab -- \eqn{M_c \not=M_p} (as specified with the argument
      \code{migration}). 
    }
    Where: 
    \tabular{ll}{
      \eqn{m_p} \tab 
      : The generational fraction of migrant individuals at present. \cr 
      \eqn{m_c} \tab
      : The generational fraction of migrant individuals between
      \eqn{T_c<t<T_s}.
    }
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
    integers. 
  }
}% end arguments

\details{
  The function \code{\link{cmd_line}} allows up to TWO populations to
  generate an IM model and a more complex model with a change in
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
  The \code{ms}-like command line \code{cmdline} needs to have a
  minimum of FIVE "words" (The first word can be any string of characters
  without white space, but \code{\link{cmd_line}} outputs \code{"./msR"}): \cr 
  \code{"./msR }\eqn{n}\code{ }\eqn{H}\code{ -t }\eqn{\theta*x*v*Z}\code{"}.
}

\note{
  \item{ATTENTION: It is the user's responsibility to mind the following
    restrictions:}{
    -> The user needs to be comfortable with the use of the C program
    \code{ms} and needs to be familiar with the documentation at
    \url{http://pps-spud.uchicago.edu/~cbecquet/msdoc.pdf}. \cr 
    -> \code{\link{cmd_line}} requires that \eqn{n>0} be
    specified with the argument \code{nsam} . \cr
    -> \code{\link{cmd_line}} requires that \eqn{\theta_1*x*v*Z>0} be
    specified with the argument \code{theta}. \cr
    -> \code{\link{cmd_line}} requires that \code{theta[i]}\eqn{>0} for any
    \eqn{i} \eqn{\in} \eqn{[1,3]} when specified. \cr 
    -> If \eqn{H>1*10^5}, \code{\link{cmd_line}} requires that
    \code{howmany} be specified as a string of characters (e.g., 
    \code{howmany="100000"}). \cr 
    -> \code{\link{cmd_line}} requires that \code{rho} be a vector of TWO
    values. \cr
    -> If there are TWO populations (as specified with
    \code{structure[1]=}\eqn{2} and \eqn{n_1>0} and \eqn{n_2>0} with the
    argument \code{structure}), \code{\link{cmd_line}} requires that
    \eqn{n=n_1+n_2} be specified with the argument \code{nsam}. \cr
    -> \code{\link{cmd_line}} requires that \eqn{M_p>0} be specified
    with the argument \code{migration} in case of an island model. \cr 
    -> \code{\link{cmd_line}} requires
    \eqn{0}\eqn{\le}\eqn{T_c}\eqn{\le}\eqn{T_s} when specified with the
    argument \code{time}.
  }
}

\author{
  Celine Becquet - \email{celine.becquet@gmail.com}. 
}

\references{
  Hudson, R. R., 1983. Properties of a neutral allele model with
  intragenic recombination. Theor. Popul. Biol. 23:183-201. \cr 
  Hudson, R. R., 1990. Gene genealogies and the coalescent process, in
  D. Futuyma and J. Antonovics (eds), Oxford Surveys in Evolutionary
  Biology, Vol. 7:1-44. \cr 
  Hudson, R. R., 2002. Generating samples under a Wright-Fisher neutral
  model of genetic variation. Bioinformatics 18:337-338. \cr 
  \url{http://pps-spud.uchicago.edu/~cbecquet/msdoc.pdf}. 
}

\seealso{
  The functions \code{\link{simulate_data}} and
  \code{\link{estimate_IMc}} call the function \code{\link{cmd_line}}. \cr 
  The \R function \code{\link{ms}} uses the output of the function
  \code{\link{cmd_line}} as input. \cr
  Lists of definitions of the symbols and parameters mentioned in
  this file are found in \code{\link{Rmspack-package}}.
}

\examples{
### Create the inputs for the function.
np=c(2, 5, 6) # npop n1 n2.
ns=np[2]+np[3] # total number of samples.
h=2 # independent gene genealogies.
s=c("1", "2", "3") # Seed numbers.

th=c(5, 10) # 4Ne*mu*Z for pop1 and pop2, here theta2/Theta1=2.
rh=c(1, 1000) # 4Ne*c*Z and Z, recombination rate and size of locus in bp.

T=c(1, 0.5) # Time of divergence and time of migration rate change in unit of 4Ne generations. 
M=c(5, 0) # Migration rate until present and eM, the migration rate between T_div and T_change, forward in time.

## Generate a command for an IM model with a change in gene flow rate 0.5/4N_1 generations ago.
cmd_line(nsam=ns, howmany=h, theta=th, rho=rh, structure=np, seeds=s, time=T, migration=M) 

## Generate a command with a change in effective population size .5/4Ne generations ago.
cmd_line(nsam=ns, howmany=h, theta=th, extra=c("-eN 0.5 2"))
}

\keyword{misc}
\keyword{IO}
