\name{ms}
\alias{ms}

\title{
  Function to call the C program ms. 
}

\description{
  The function \code{\link{ms}} calls Hudson's C program \code{ms} to
  simulate independent data sets for a genomic region of interest.
  \code{\link{ms}} uses as input an \code{ms}-like command line defined
  by the user. \cr 
  This function is called by the \R functions \code{\link{simulate_data}}
  and \code{\link{estimate_IMc}}. 
}

\usage{
ms(cmdline)
}

\arguments{
  \item{cmdline}{
    The string of character with all the tags and values required by the C
    program \code{ms}
    (cf. \url{http://pps-spud.uchicago.edu/~cbecquet/msdoc.pdf}). \cr 
    This string is of the same format as output by the function
    \code{\link{cmd_line}}. 
  }
}

\details{
  ATTENTION: The \code{ms} command line \code{cmdline} needs to have a
  minimum of FIVE "words" (The first word can be any string of characters
  without white space): \cr 
  \code{"./ms }\eqn{n}\code{ }\eqn{H}\code{ -t }\eqn{\theta_1*x*v*Z}\code{"}. 
  \tabular{ll}{
    \eqn{n} \tab
    : The total number of chromosomes sampled for the genomic region
    considered. \cr
    \tab ATTENTION: \code{\link{ms}} requires that \eqn{n>0} be specified. \cr
    \eqn{H} \tab
    : The number of independent data sets to simulate with different gene
    genealogy samples for the genomic region considered. \cr 
    \eqn{\theta_1*x*v*Z} \tab
    : The region-specific population mutation rate per generation for
    population 1. \cr 
    
    \tab ATTENTION: \code{\link{ms}} requires that \eqn{\theta*x*v*Z>0} be
    specified. \cr
  }
  \item{}{Where:
    \tabular{ll}{
      \eqn{\theta_1=4N_1*\mu} \tab
      : The population mutation rate per bp
      per generation for population 1. \cr
      \eqn{N_1} \tab
      : The effective population size in population 1 (the reference
      population so by default and unless specified, \eqn{N_1=N_2=N_A}). \cr 
      \eqn{\mu} \tab
      : The genomic generational mutation rate per bp. \cr
      \eqn{x} \tab
      : The inheritance scalar for the genomic region considered (i.e.,
      \code{"1"} for autosomal region, \code{"0.75"} for X- and \code{"0.5"}
      for Y- and mtDNA-linked region). \cr 
      \eqn{v} \tab
      : The mutation rate scalar for the genomic region considered (which can
      be estimated e.g., from divergence data). \cr 
      \eqn{Z} \tab
      : The size in bp of the genomic region considered. 
    }
  }
}

\value{ 
  The function \code{\link{ms}} outputs a data frame with the
  information usually output by Hudson's C program \code{ms}.
  \item{$cmdline_seeds}{
    The vector with two entries: the \code{ms} command line and the three
    random seeds used to run the C program \code{ms}. 
  }
  \item{$segsites}{
    The vector of \eqn{H} integers \eqn{S_h}, the number of segregating
    sites found in the simulated data set \eqn{h}, \eqn{h} \eqn{\in}
    \eqn{[1,H]}. \cr
    \eqn{H} is specified in the \code{ms} command line \code{$cmdline}. 
  }
  \item{$positions}{
    The vector of \eqn{H} lists of positions of the segregating sites
    found in the \eqn{H} independent data sets. 
    \tabular{ll}{
      \code{$positions[[}\eqn{h}\code{]]} \tab
      : The list of positions of the \eqn{S_h} segregating sites found in the
      simulated data set \eqn{h} \eqn{\in} \eqn{[1,H]} along the genomic
      region considered. \cr
      \code{$positions[[}\eqn{h}\code{]][}\eqn{s}\code{]} \tab
      : The position of the segregating site \eqn{s} \eqn{\in}
      \eqn{[1,S_h]}, \eqn{p_s} \eqn{\in} \eqn{[0,1]} relative to \eqn{Z}, the
      size of the genomic region simulated. 
    }
  }
  \item{$haplotypes}{
    The matrix \eqn{H} \eqn{\times} \eqn{n} of the haplotypes simulated for
    the \eqn{H} independent data sets. \cr 
    \tabular{ll}{
      \code{$haplotypes[}\eqn{h,}\code{]} \tab
      : The \eqn{n} haplotypes simulated with the C program \code{ms} for the
      data set \eqn{h}, \eqn{h} \eqn{\in} \eqn{[1,H]}. \cr 
      \code{$haplotypes[}\eqn{h,i}\code{]} \tab
      : The string of \eqn{S_h} characters corresponding to the alleles found
      at the \eqn{S_h} segregating sites in the data set \eqn{h}, \eqn{i}
      \eqn{\in} \eqn{[1,n]}. \cr 
      \tab i.e., \code{"0"}, \code{"1"} or \code{"?"}, for ancestral,
      derived and missing allele, respectively. 
    }
  }
}

\note{
  \item{ Note that when \eqn{H>1}, \code{\link{ms}} assumes that each of the
    \eqn{H} data sets are simulated under the same model (i.e., same sample
    size \eqn{n}, population mutation rates and demographic history). \cr 
    ATTENTION: It is the user's responsibility to mind the following
    restrictions: \cr }{
    -> The user needs to be comfortable with the use of the C program
    \code{ms} and needs to be familiar with the documentation at
    \url{http://pps-spud.uchicago.edu/~cbecquet/msdoc.pdf}. \cr 
    -> The \code{ms}-like command line \code{cmdline} needs to have a
    minimum of FIVE "words" (The first word can be any string of characters
    without white space): \cr 
    \code{"./ms }\eqn{n}\code{ }\eqn{H}\code{ -t }\eqn{\theta*x*v*Z}\code{"}.
    \tabular{lll}{ 
      \tab --- \tab \code{\link{ms}} requires that \eqn{n>0} be specified.\cr
      \tab --- \tab \code{\link{ms}} requires that \eqn{\theta_1*x*v*Z>0} be
      specified. 
    }
    -> \code{\link{ms}} requires \eqn{N_i/N_1>0} for any
    \eqn{i} \eqn{\in} \eqn{{2,A}} when specified with
    \code{"-n 2 }\eqn{N_2/N_1}\code{"} and
    \code{"-eN }\eqn{T_s}\code{ }\eqn{N_A/N_1}\code{"} in the
    command line. \cr
    \tabular{llll}{ 
      \tab --- \tab Where: \tab \cr
      \tab \tab \eqn{N_2} \tab
      : The effective population size in population 2. \cr
      \tab \tab \eqn{N_A} \tab
      : The ancestral effective population size. \cr
      \tab \tab \eqn{T_s} \tab
      : The split time in unit of \eqn{4N_1} generations between
      the TWO populations. \cr
    }
    -> \code{\link{ms}} requires that the region-specific
    population recombination rate per generation be specified with
    \code{"-r }\eqn{\rho*w*(Z-1)}\code{ }\eqn{Z}\code{"} in the command
    line. 
    \tabular{llll}{
      \tab --- \tab Where: \tab \cr
      \tab \tab \eqn{\rho=4N_1*c} \tab
      : The genomic average population intra-region recombination rate
      per bp per generation. \cr
      \tab \tab \eqn{c} \tab
      : The generational cross-over rate per bp. \cr
      \tab \tab \eqn{w} \tab
      : The recombination scalar for the genomic region considered.\cr
      \tab \tab \tab \eqn{w=\beta} is the ratio of the locus-specific
      population recombination rate per bp over \eqn{\rho=4N_1*c}.\cr
      \tab \tab \tab \eqn{\beta=}\code{"1"} (\code{"0.5"} in \emph{Drosophila}) for autosomal
      region, \code{"0.5"} for X- and \code{"0"} for Y- and
      mtDNA-linked region. \cr 
      \tab \tab \eqn{Z} \tab
      : The size of the genomic region considered in bp. 
    }
    -> If there are TWO populations (as specified with
    \code{"-I 2 }\eqn{n_1}\code{ }\eqn{n_2}\code{ }\eqn{M_p}\code{"} 
    in the command line), \code{\link{ms}} requires \eqn{n=n_1+n_2} and
    \eqn{n_1>0} and \eqn{n_2>0}. \cr
    \tabular{llll}{
      \tab --- \tab Where: \tab \cr
      \tab \tab \eqn{n_1} \tab
      : The sample size from population 1 for the genomic region
      considered. \cr 
      \tab \tab \eqn{n_2} \tab
      : The sample size from population 2 for the genomic region
      considered.\cr 
      \tab \tab \eqn{M_p=4N_1*m_p} \tab
      : The number of migrants exchanged each generation by the
      TWO populations until present.\cr
      \tab \tab \tab \eqn{M_p} has different definitions depending of the model
      considered:\cr 
      \tab \tab \tab * \eqn{M_p} is the symmetrical rate of gene flow between
      two population in an island model (i.e., \eqn{0<t<\inf}) as
      specified with:\cr 
      \tab \tab \tab -- A population structure (as specified with
      \eqn{n_1>0} and \eqn{n_2>0} with
      \code{"-I 2 }\eqn{n_1}\code{ }\eqn{n_2}\code{ }\eqn{M_p}\code{"}
      in the command line . \cr
      \tab \tab \tab -- Either \eqn{T_s=0} when specified with
      \code{"-ej }\eqn{T_s}\code{ }\code{2 1}\code{"}) or unspecified.\cr
      \tab \tab \tab * \eqn{M_p} is the constant symmetrical rate of gene flow
      since the split until present (i.e., \eqn{0<t<T_s}) if there is a
      population split (as specified with \eqn{T_s>0} with
      \code{"-ej }\eqn{T_s}\code{ }\code{2 1}\code{"}).\cr 
      \tab \tab \tab * \eqn{M_p} is the constant symmetrical rate of gene flow
      since the time of gene flow rate change until present (i.e.,
      \eqn{0<t<T_c}) if a time at which the gene flow rate changed is
      specified with:\cr 
      \tab \tab \tab -- \eqn{0<T_c<T_s} (as specified with \eqn{T_s>0} and
      \eqn{0<T_c<T_s} with \code{"-ej }\eqn{T_s}\code{ 2 1}\code{"} and
      \code{"-eM }\eqn{T_c}\code{ }\eqn{M_c}\code{"} in the command
      line). \cr 
      \tab \tab \tab -- \eqn{M_c \not=M_p}. \cr
      \tab \tab \eqn{m_p} \tab 
      : The generational fraction of migrant individuals at present. \cr
    }
    -> In an island model (i.e., no population split specified with
    \eqn{T_s>0} in
    \code{"-ej }\eqn{T_s}\code{ }\code{2 1}\code{"}
    in the command line), \code{\link{ms}} requires that the rate of
    gene flow between the two structured populations be specified with
    \eqn{M_p>0} in 
    \code{"-I 2 }\eqn{n_1}\eqn{n_2}\code{ }\eqn{M_p}\code{"} in 
    the command line. \cr
    -> \code{\link{ms}} requires
    \eqn{0}\eqn{\le}\eqn{T_c}\eqn{\le}\eqn{T_s} when specified with
    \code{"-ej }\eqn{T_s}\code{ 2 1}\code{"} and
    \code{"-eM }\eqn{T_c}\code{ }\eqn{M_c}\code{"} in the command
    line. 
    \tabular{llll}{
      \tab --- \tab Where: \tab \cr
      \tab \tab \eqn{T_c} \tab
      : The time at which the rate of gene flow changed between
      the two populations in unit of \eqn{4N_1} generations. \cr
      \tab \tab \eqn{M_c=4N_1*m_c} \tab
      : The number of migrants exchanged each generation by the TWO
      populations since the split until the time of gene flow rate
      change (i.e., \eqn{T_c<t<T_s}). \cr
      \tab \tab \eqn{m_c} \tab
      : The generational fraction of migrant individuals between
      \eqn{T_c<t<T_s}.
    }
    -> In general, mind the restrictions listed in
    \code{\link{cmd_line}}. 
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
  The \R function \code{\link{ms}} calls Dick Hudson's' program \code{ms}
  (see documentation at
  \url{http://pps-spud.uchicago.edu/~cbecquet/msdoc.pdf}). \cr 
  The functions \code{\link{simulate_data}} and
  \code{\link{estimate_IMc}} call the function \code{\link{ms}}. \cr 
  The following functions input, output and/or manipulate the object
  output by the function \code{\link{ms}}: \cr 
  \code{\link{get_loci}}, \code{\link{get_ms_output}},
  \code{\link{out_ms_loci}}, \code{\link{output_ms}},
  \code{\link{stats_pop}} and
  \code{\link{write_ms_output}}. \cr 
  The function \code{\link{cmd_line}} generates inputs for the \R function
  \code{\link{ms}}. \cr
  Lists of definitions of the symbols and parameters mentioned in
  this file are found in \code{\link{Rmspack-package}}.
}

\examples{
### Create the input for the function.
msline="./ms 8 2 -t 1 -r 0.1 9999 -I 2 4 4 10 -ej 1 1 2 -eN 1 3"

ms(msline)
}

\keyword{programming}