\name{get_loci}
\alias{get_loci}

\title{
  Extract the data for a locus from the  data frame of the kind
  output by the \R function ms for a multi-locus region.
}

\description{
  The function  \code{\link{get_loci}} extracts the data for the locus
  \code{nlocus} from a data frame \code{msframe} of the kind output by
  the \R function 
  \code{\link{ms}} for a multi-locus region.\cr
  This function is called within the \R functions
  \code{\link{simulate_data}} and \code{\link{estimate_IMc}}. 
  
}
\usage{
get_loci(msframe, info_region, info_loci, nlocus)
}

\arguments{
  \item{msframe}{
    The data frame of the kind output by the \R function \code{\link{ms}}
    for the multi-locus genomic region considered (with some
    variations from Hudson's program). \cr 
    \item{\code{msframe$cmdline_seeds}}{
      The vector with two entries: the \code{ms} command line and the
      three random seeds used (or that could have been used) to run the
      C program \code{ms} for the multi-locus genomic region
      considered.\cr 
    }
    \item{\code{msframe$segsites}}{
      The integer \eqn{S}, the number of segregating sites found in the
      simulated data set for the multi-locus genomic region considered. 
    }
    \item{\code{msframe$positions}}{
      The list of positions in base pairs of the \eqn{S} segregating sites
      found along the multi-locus genomic region considered. \cr
      \code{msframe$positions[[1]][}\eqn{j}\code{]} \eqn{\in}
      \eqn{[z_s,z_e]} is the position of the segregating site \eqn{j}
      \eqn{\in} \eqn{1,S}.\cr 
      i.e., the value output by the \R function \code{\link{ms}}
      \eqn{\times} \eqn{Z} for the multi-locus genomic region considered.\cr 
      \eqn{Z=z_e-z_s} is the size in bp of the multi-locus genomic region
      considered as specified in the vector \code{info_region}.
    }
    \item{\code{msframe$haplotypes}}{
      The vector of \eqn{n=n_1+n_2} haplotypes simulated for the
      multi-locus genomic region considered. \cr
      \eqn{n_1} and \eqn{n_2} are the sample sizes in the two
      populations for the multi-locus genomic region considered
      specified in the vector \code{info_region}. \cr 
      \code{msframe$haplotypes[}\eqn{1,i}\code{]} is the string of \eqn{S}
      characters corresponding to the alleles found at the \eqn{S}
      segregating sites for the multi-locus genomic region considered,
      \eqn{i} \eqn{\in} \eqn{[1,n]}. \cr
      i.e., \code{"0"}, \code{"1"} or \code{"?"}, for ancestral,
      derived and missing allele, respectively. 
    }
    See \code{\link{ms}} for further details. 
  }
  \item{info_region}{
    The vector of \eqn{10} items with the information for the
    multi-locus genomic region \eqn{r}.\cr
    \tabular{ll}{ 
      \eqn{r} \tab 
      : The multi-locus genomic region number. \cr 
      \emph{\code{Region name}} \tab 
      : The name of the multi-locus genomic region \eqn{r}, which should
      contain at least ONE non-numerical character. \cr 
      \eqn{x} \tab
      : The inheritance scalar for the multi-locus genomic region
      \eqn{r} (i.e., \code{"1"} for autosomal region, \code{"0.75"} for X- and
      \code{"0.5"} for Y- and mtDNA-linked region). \cr 
      \eqn{v} \tab 
      : The mutation rate scalar for the multi-locus genomic region
      \eqn{r} (which can be estimated e.g., from divergence data). \cr 
      \eqn{w} \tab 
      : The recombination scalar for the multi-locus genomic region \eqn{r}.\cr
      \tab -- Usually \eqn{w=\beta}, the ratio of the locus-specific
      population recombination rate per bp over \eqn{\rho=4N_1*c}. \cr 
      \tab -- If an estimate of the region-specific population recombination
      rate per bp is available for each of  the independent recombining genomic
  regions from linkage 
      disequilibrium analysis, \eqn{\rho_o=4N_1*c_o}, set
      \eqn{w=\beta*\rho_o} to incorporate this knowledge in the
      simulation (with \code{"}\emph{\code{rho}} \code{1"} in the
      file \code{paramfile}, see \code{\link{param_sim}}).\cr
      \tab In this case,
      \eqn{w} is the scaled sex-averaged region-specific population
      recombination rate per bp, i.e., for an X-linked locus \eqn{c_o} is
      the female recombination rate and \eqn{\beta=0.5} so that
      \eqn{\beta*\rho_o=2N_1*c_o}. \cr 
      \tab -- If an estimate of the region-specific recombination rate per bp
      is available for each of the independent recombining genomic regions from
      pedigree analysis, \eqn{c_o}, 
      set \eqn{w=\beta*c_o} to incorporate this knowledge in the simulation
      (with \code{"}\emph{\code{rho}} \code{2"} in the file
      \code{paramfile}, see \code{\link{param_sim}}).\cr
      \tab In this case,
      \eqn{w} is the scaled sex-averaged region-specific recombination rate
      per bp, i.e., for an X-linked locus, \eqn{c_o} is the estimated female
      recombination rate so the scaled sex-averaged recombination rate is
      \eqn{\beta*c_o=0.5c}. \cr 
      \eqn{n_1} \tab 
      : The sample size from population 1 for the multi-locus genomic
      region \eqn{r}. \cr 
      \eqn{n_2} \tab 
      : The sample size from population 2 for the multi-locus genomic
      region \eqn{r}. \cr
      \eqn{z_s} \tab 
      : The start position of the multi-locus genomic region \eqn{r} in
      bp. \cr 
      \eqn{z_e} \tab 
      : The end position of the multi-locus genomic region \eqn{r} in
      bp. \cr 
      \eqn{Y} \tab 
      : The number of loci spanning the multi-locus genomic region
      \eqn{r}. \cr 
    }Where: 
    \tabular{ll}{
      \eqn{\rho=4N_1*c}\tab
      : The genomic average population intra-region recombination rate per bp per
      generation. \cr
      \eqn{\rho_o=4N_1*c_o} \tab 
      : The estimate of the region-specific population recombination
      rate per bp per generation for the multi-locus genomic region considered from
      linkage disequilibrium analysis. \cr
      \eqn{N_1} \tab
      : The effective population size in population 1 (the reference
      population). \cr
      \eqn{c} \tab 
      : The genomic generational cross-over rate per bp. \cr
      \eqn{c_o} \tab 
      : The estimate of the region-specific cross-over rate per bp
      per generation for the multi-locus genomic region considered from pedigree
      analysis. \cr
      \eqn{\beta} \tab
      : The ratio of the region-specific population recombination
      rate per bp over \eqn{\rho=4N_1*c} for the multi-locus genomic region
      considered.\cr
      \tab \eqn{\beta=}\code{"1"} (\code{"0.5"} in
      \emph{Drosophila}) for autosomal region, \code{"0.5"} for X- and
      \code{"0"} for Y- and mtDNA-linked region.
    } 
    See \code{\link{info_region}} for further details. 
  }
  \item{info_loci}{
    The vector of \eqn{6} items with the information on the
    loci for the multi-locus genomic region \eqn{r}.
    \eqn{Y>1} is specified in \code{info_region[10]}. \cr 
    This matrix is empty unless \eqn{Y>0}. \cr 
    The six values are:
    \tabular{ll}{
      \eqn{r} \tab 
      : The multi-locus genomic region number, \eqn{r}, that this locus
      is part of. \cr 
      \eqn{y} \tab 
      : The locus number, \eqn{y} \eqn{\in} \eqn{[1,Y]}. \cr
      \eqn{n_1y} \tab 
      : The sample size from population 1 for the locus \eqn{y}.\cr
      \tab \eqn{n_1y}\eqn{\le}\eqn{n_1}, where \eqn{n_1} for the multi-locus
      genomic region \eqn{r} is specified in the vector \code{info_region}. \cr 
      \eqn{n_2y} \tab 
      : The sample size from population 2 for the locus \eqn{y}.\cr
      \tab \eqn{n_2y}\eqn{\le}\eqn{n_2}, where \eqn{n_2} for the multi-locus
      genomic region \eqn{r} is specified in the vector
      \code{info_region}. \cr 
      \eqn{z_sy} \tab 
      : The start position of the locus \eqn{y} in bp.\cr
      \tab \eqn{z_s1}\eqn{\ge}\eqn{z_s}, where \eqn{z_s} is the start position of the
      multi-locus genomic region \eqn{r} specified in the vector
      \code{info_region}. \cr 
      \eqn{z_ey} \tab 
      : The end position of the locus \eqn{y} in bp.\cr
      \tab \eqn{z_eY}\eqn{\le}\eqn{z_e} where \eqn{z_e} is the end position of
      the multi-locus genomic region \eqn{r} specified in the vector
      \code{info_region}. 
    }
    See \code{\link{info_loci}} for further details. 
  }  
  \item{nlocus}{
    \code{nlocus=}\eqn{y}, the locus number, \eqn{y} \eqn{\in} \eqn{[1,Y]}. 
  }
}

\details{
  The function \code{\link{get_loci}} takes \code{msframe}, the data
  frame output by the \R function \code{\link{ms}} for a multi-locus
  genomic region, and updates to include only the segregating sites located
  within  the locus \eqn{y} (specified with the argument \code{nlocus})
  in the multi-locus genomic region
  considered. \cr  
}

\value{
  The function \code{\link{get_loci}} outputs \code{msframe}, the
  data frame of the kind output by the \R function \code{\link{ms}}, but
  updated for the locus \eqn{y} in the multi-locus genomic region considered:
  \item{}{\strong{Matrices from \code{msframe}:}
    \item{\code{$cmdline_seeds}}{
      The vector with two entries: the \code{ms} command line and the three
      random seeds if the data where simulated withthe C program
  \code{ms} for the multi-locus genomic region considered.\cr
      i.e., \code{msframe$cmdline_seeds}.
    }
    \item{\code{$segsites}}{
      The integer \eqn{S_y}, the number of segregating sites found in the
      data set for locus \eqn{y}.
      i.e., the number of segregating sites found within the
      positions \eqn{z_sy} and \eqn{z_ey}. 
    }
    \item{\code{$positions}}{
      The list of positions in base pairs of the \eqn{S_y} segregating sites found
      along the locus \eqn{y}. \cr
      i.e.,
      \code{$positions[[1]][}\eqn{j}\code{]=msframe$positions[[1]][$ns[}\eqn{j}\code{]]}
      \eqn{in} \eqn{[z_sy,z_ey]}, \eqn{j} \eqn{\in} \eqn{[1:S_y]}.
    }
    \item{\code{$haplotypes}}{
      The vector of \eqn{n=n_1+n_2} haplotypes simulated for the multi-locus
      genomic region considered. \cr 
      \eqn{n_1} and \eqn{n_2} are the sample sizes for the multi-locus
      genomic region considered specified in the vector \code{info_region}.\cr 
      \code{$haplotypes[}\eqn{1,i}\code{]} is the string of \eqn{S_y}
      characters corresponding to the alleles found at the \eqn{S_y}
      segregating sites within locus \eqn{y}, \eqn{i}
      \eqn{\in} \eqn{[1,n]}. \cr
       i.e., \code{$haplotypes[}\eqn{1,i}\code{]} corresponds to the
      \eqn{S_y} characters at positions \eqn{s} reported in \code{$ns}
      withing the string of character
      \code{msframe$haplotypes[}\eqn{1,i}\code{]}. \cr
      If \eqn{n_py<n_p} (where \eqn{n_1y} and \eqn{n_2y} are the sample
      sizes for the locus \eqn{y} specified in the matrix \code{info_loci}),
      the segregating sites spanning the positions of locus
      \eqn{y} at \eqn{n_p-n_py} haplotypes are noted as
      missing data with \code{"?"}.
    }
  }
  \item{}{\strong{Extra vector:}
    \item{\code{$ns}}{
      The list of \eqn{S_y} positions \eqn{s} within the string of characters
      \code{msfram$haplotypes} and the vector \code{msfram$positions} of the
      \eqn{S_y} segregating sites found along the locus \eqn{y}. \cr 
      \eqn{s} \eqn{\in} \eqn{[1,S]}, where \eqn{S} is the original number of
      segregating sites found for the multi-locus genomic region considered
      in \code{msframe}. 
    }
  }
}

\author{
  Celine Becquet - \email{celine.becquet@gmail.com}. 
}

\seealso{
  The functions \code{\link{simulate_data}} and
  \code{\link{estimate_IMc}} call the function \code{\link{get_loci}}. \cr 
  \code{\link{out_ms_loci}} is another functions to manipulate data
  frame of the kind output by the \R function \code{\link{ms}}. \cr
  Lists of definitions of the symbols and parameters mentioned in
  this file are found in \code{\link{Rmspack-package}}.
}

\examples{
### Create the inputs for the function.
data(simulation_files) # download the data

## Information on region
write.table(file="regions", x=info_region, row.name=FALSE, col.names=FALSE, quote=FALSE)
# Creates the file "regions" containing: 
info_region=read.table("regions",header=FALSE,skip=1)[3,]

## Information on loci
write.table(file="loci", x=info_loci, row.name=FALSE, col.names=FALSE, quote=FALSE)
# Creates the file "loci" containing:
info_loci=read.table("loci",header=FALSE,skip=1)
info_loci=info_loci[info_loci$V1==info_region$V1,]

msline="./ms 8 1 -t 1 -r 0.1 9999 -I 2 4 4 10 -ej 1 1 2 -eN 1 3" 
msframe=ms(msline)
msframe$positions[[1]]=msframe$positions[[1]]*(info_region$V9-info_region$V8)
msframe

get_loci(msframe=msframe, info_region=info_region, info_loci=info_loci, nlocus=1)

get_loci(msframe=msframe, info_region=info_region, info_loci=info_loci, nlocus=2)

# Clean up the directory.
unlink(c("regions","loci"))
}

\keyword{ misc }