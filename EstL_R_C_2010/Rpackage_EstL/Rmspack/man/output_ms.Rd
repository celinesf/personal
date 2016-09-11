\name{output_ms}
\alias{output_ms}

\title{
  Function to write the output of the \R function ms for a (potentially
  multi-locus) genomic region into a text file with ms-like format. 
}
\description{
  The function \code{\link{output_ms}} outputs and writes the data frame
  \code{msframe} output by the \R function \code{\link{ms}} updated for
  multiple-loci if applicable into an
  \code{ms}-like text file named \code{"}\code{datafile}\eqn{r}\code{"},
  where \eqn{r} is the number of the genomic region considered specified
  by the argument \code{nregion}. \cr 
  This function is called by the \R function \code{\link{simulate_data}}. 
}

\usage{
output_ms(datafile, nregion, msframe, info_region, info_loci)
}

\arguments{
  \item{datafile}{
    The argument \code{datafile} specifies the location and the string
    of characters that initiates the name of the text file output by
    the function \code{\link{output_ms}}. 
  }
  \item{nregion}{
    The number of the genomic region considered, \eqn{r} \eqn{\in}
    \eqn{[1,R]}.\cr 
    \eqn{R} is the number of independent genomic regions to simulate as
    specified with the keyword \emph{\code{nregions}} in the file
    \code{paramfile}, see \code{\link{param_sim}}). 
  }
  \item{msframe}{
    The data frame of the kind output by the \R function \code{\link{ms}}
    for the genomic region \eqn{r}. \cr
    \item{\code{msframe$cmdline_seeds}}{
      The vector with two entries: the \code{ms} command line and the three
      random seeds used to run the C program \code{ms} for the genomic
      region \eqn{r}.\cr 
    }
    \item{\code{msframe$segsites}}{
      The integer \eqn{S}, the number of segregating sites found in the
      simulated data set for the genomic region \eqn{r}. 
    }
    \item{\code{msframe$positions}}{
      The list of positions of the \eqn{S} segregating sites
      found along the genomic region \eqn{r}. \cr
      \code{msframe$positions[[1]][}\eqn{j}\code{]} \eqn{\in}
      \eqn{[0,1]} relative to \eqn{Z} is the position of the segregating
      site \eqn{j} \eqn{\in} \eqn{[1,S]}.\cr 
      \eqn{Z=z_e-z_s} is the size in bp of the genomic region \eqn{r} as
      specified in the vector \code{info_region}. 
    }
    \item{\code{msframe$haplotypes}}{
      The vector of \eqn{n=n_1+n_2} haplotypes simulated for the genomic
      region \eqn{r}. \cr
      \eqn{n_1} and \eqn{n_2} are the sample sizes in the two
      populations for the genomic region \eqn{r} specified in the vector
      \code{info_region}. 
      \code{msframe$haplotypes[}\eqn{1,i}\code{]} is the string of \eqn{S}
      characters corresponding to the alleles found at the \eqn{S}
      segregating sites for the genomic region \eqn{r}, \eqn{i} \eqn{\in}
      \eqn{[1,n]}.\cr 
      i.e., \code{"0"}, \code{"1"} or \code{"?"}, for ancestral),
      derived and missing allele, respectively. 
    } See \code{\link{ms}} for further details. 
  }
  \item{info_region}{
    The vector of \eqn{10} items with the information for the genomic
    region \eqn{r}.\cr
    \tabular{ll}{ 
      \eqn{r} \tab 
      : The genomic region number. \cr 
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
      rate per bp is available for each of the independent recombining genomic
      regions from linkage disequilibrium analysis, \eqn{\rho_o=4N_1*c_o}, set
      \eqn{w=\beta*\rho_o} to incorporate this knowledge in the
      simulation (with \code{"}\emph{\code{rho}} \code{1"} in the
      file \code{paramfile}, see \code{\link{param_sim}}).\cr
      \tab In this case,
      \eqn{w} is the scaled sex-averaged region-specific population
      recombination rate per bp, i.e., for an X-linked locus \eqn{c_o} is
      the female recombination rate and \eqn{\beta=0.5} so that
      \eqn{\beta*\rho_o=2N_1*c_o}. \cr 
      \tab -- If an estimate of the region-specific recombination rate per bp
      is available for each of the independent recombining genomic
      regions from pedigree analysis, \eqn{c_o}, 
      set \eqn{w=\beta*c_o} to incorporate this knowledge in the simulation
      (with \code{"}\emph{\code{rho}} \code{2"} in the file
      \code{paramfile}, see \code{\link{param_sim}}).\cr
      \tab In this case,
      \eqn{w} is the scaled sex-averaged region-specific recombination rate
      per bp, i.e., for an X-linked locus, \eqn{c_o} is the estimated female
      recombination rate so the scaled sex-averaged recombination rate is
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
    }Where: 
    \tabular{ll}{  
      \eqn{\rho=4N_1*c} \tab
      : The genomic average population intra-region recombination rate
      per bp per generation. \cr
      \eqn{\rho_o=4N_1*c_o} \tab 
      : The estimate of the region-specific population recombination
      rate per bp per generation for the genomic region \eqn{r} from
      linkage disequilibrium analysis. \cr  
      \eqn{N_1} \tab
      : The effective population size in population 1 (the reference
      population). \cr
      \eqn{c} \tab 
      : The genomic generational cross-over rate per bp. \cr
      \eqn{c_o} \tab 
      : The estimate of the region-specific cross-over rate per bp
      per generation for the genomic region \eqn{r} from pedigree
      analysis. \cr
      \eqn{\beta} \tab
      : The ratio of the region-specific population recombination
      rate per bp over \eqn{\rho=4N_1*c} for the genomic region
      \eqn{r}.\cr
      \tab \eqn{\beta=}\code{"1"} (\code{"0.5"} in
      \emph{Drosophila}) for autosomal region, \code{"0.5"} for X- and
      \code{"0"} for Y- and mtDNA-linked region. 
    } See \code{\link{info_region}} for further details. 
  }
  \item{info_loci}{
    The matrix \eqn{Y} \eqn{\times} \eqn{6} of information on the
    loci for the multi-locus genomic region \eqn{r} (as specified with
    \eqn{Y>1} in the vector \code{info_region}). \cr 
    This matrix is empty unless \eqn{Y>0}. \cr 
    Each locus is described by six values:
    \tabular{ll}{
      \eqn{r} \tab 
      : The multi-locus genomic region number, \eqn{r}, that this locus
      is part of. \cr 
      \eqn{y} \tab 
      : The locus number, \eqn{y} \eqn{\in} \eqn{[1,Y]}. \cr
      \tab \eqn{Y} is the total number of loci spanning the multi-locus
      genomic region \eqn{r} as specified in the vector \code{info_region}. \cr 
      \tab ATTENTION: The loci numbers, \eqn{y}, need to be in order,
      i.e., information for locus 1 of the first multi-locus genomic region
      \eqn{r} is on line 2, locus 2 on line 3 \dots locus \eqn{Y} on line
      \eqn{Y+1}. \cr 
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
    } See \code{\link{info_loci}} for further details. 
  }
}

\details{
  The function \code{\link{output_ms}} creates a text file named
  \code{"}\code{datafile}\eqn{r}\code{"} in the local directory. \cr
  \eqn{r} is the number of the genomic region \eqn{r} specified with the
  argument \code{nregion}. \cr 
  In the text file \code{"}\code{datafile}\eqn{r}\code{"},
  \code{\link{output_ms}} writes the data frame \code{msframe} output by
  the \R function \code{\link{ms}} for the genomic region \eqn{r} in an
  \code{ms}-like format. \cr 
  If there are multiple loci in the genomic region \eqn{r},
  \code{\link{output_ms}} calls the function \code{\link{out_ms_loci}} to
  update \code{msframe} accounting for the multiple loci information
  before writing in the text file. 
}

\value{  
  \item{\code{"}\code{datafile}\eqn{r}\code{"}}{
    The function \code{\link{output_ms}} creates the text file named
  \code{"}\code{datafile}\eqn{r}\code{"}, where \eqn{r} is the number of the
  genomic region considered specified with the argument
  \code{nregions}.\cr
    In the file \code{"}\code{datafile}\eqn{r}\code{"},
    \code{\link{output_ms}} writes the updated data frame \code{msframe}
  (see below) 
    for the genomic region \eqn{r} in the usual \code{ms} format (with
    some variations from Hudson's program):
    \item{Line 1}{
      The \code{ms} command line used to run the C program \code{ms}. \cr
      i.e., \code{msframe$cmdline_seeds[1]}.
    } 
    \item{Line 2}{   
      The three random seeds used to run the C program \code{ms}. \cr
      i.e., \code{msframe$cmdline_seeds[2]}.
    } 
    \item{Line 3}{   
      An empty line.
    }
    \item{Line 4}{
      \code{"\\\\"} to specify that the data set simulated for the genomic region
      \eqn{r} starts below.
    } 
    \item{Line 5}{  
      \code{"}\code{segsites: }\eqn{S}\code{"}, the number of segregating
      sites found in the data set simulated for the genomic region
      \eqn{r}.\cr
      i.e., \code{msframe$segsites[1]}. \cr 
      -- In case of a multi-locus genomic region (as specified with
      \eqn{Y>1} in the vector \code{info_region}), \eqn{S=S_n},
      the number of segregating sites left once gaps and missing data for
      the multiple loci are accounted for. \cr  
      \eqn{S_n=\Sigma(S_y)}, where \eqn{S_y} is the
      number of segregating sites found in locus \eqn{y}, \eqn{y} \eqn{\in}
      \eqn{[1,Y]} for the multi-locus genomic region \eqn{r}.
    }  
    \item{Line 6}{  
      \code{"}\code{position:}\code{"} followed by the \eqn{S} positions
      in base pairs of the \eqn{S} segregating sites found in the data set
      along the genomic region \eqn{r}. \cr
      The position of the segregating site \eqn{j} \eqn{\in}
      \eqn{[1,S]} is \eqn{p_j} \eqn{\in} \eqn{[z_s,z_e]}. \cr
      i.e., \eqn{p_j=}\code{msframe$positions[[1]][}\eqn{j}\code{]}
      \eqn{\times} \eqn{Z}. \cr
      \eqn{Z=z_e-z_s} is the size of the genomic region \eqn{r} as
      specified in the vector \code{info_region}. \cr
      -- In case of a multi-locus genomic region (as specified with
      \eqn{Y>1} in the vector \code{info_region}), the \eqn{S_n} 
      \eqn{p_j} values are all \eqn{\in} \eqn{[z_sy, z_ey]}, \eqn{y}
      \eqn{\in} \eqn{[1,Y]}.
    }
    \item{Lines 7  to 7\emph{+n}}{  
      \eqn{n=n_1+n_2} lines follow corresponding to the \eqn{n} haplotypes
      simulated with the C program \code{ms} for the genomic region
      \eqn{r}. \cr 
      \eqn{n_1} and \eqn{n_2} are the sample sizes in the two
      populations for the genomic region \eqn{r} specified in the vector
      \code{info_region}. \cr 
      Each line \eqn{j=i+7} \eqn{\in} \eqn{[7,7+n]} of
      haplotype is a string of \eqn{S} characters corresponding to the
      alleles found at the \eqn{S} segregating sites in the data set
      simulated for the genomic region \eqn{r}. \cr 
      i.e., \code{msframe$haplotypes[}\eqn{1,i}\code{]}, \eqn{i}
      \eqn{\in} \eqn{[1,n]}, with the values \code{"0"}, \code{"1"} or
      \code{"?"}, for ancestral, derived and missing allele,
      respectively. \cr
      -- In case of a multi-locus genomic region (as specified with
      \eqn{Y>1} in the vector \code{info_region}),
      \code{msframe$haplotypes[}\eqn{1,i}\code{]} contains only \eqn{S_n}
      characters. \cr 
      -- If for any locus \eqn{y} \eqn{n_py<n_p} (where \eqn{n_1y} and
      \eqn{n_2y} are the sample sizes for the locus \eqn{y} specified in the
      matrix \code{info_loci}), the segregating sites spanning the positions
      of locus \eqn{y} at \eqn{n_p-n_py} randomly chosen haplotypes are noted as
      missing data with \code{"?"}.
    }
  }
   \item{Updated msframe}{
    -- In case of a uni-locus genomic region 
    \code{\link{output_ms}} also outputs the data frame
  \code{msframe}, but with \code{$positions} in bp.\cr
    -- In case of a multi-locus genomic region (as specified with
    \eqn{Y>1} in the vector \code{info_region}),
    \code{\link{output_ms}} outputs the data frame \code{msframe} with
  \code{$positions} in bp and updated
    to account for the multiple loci as output by \code{\link{out_ms_loci}}:
    \item{\code{$cmdline_seeds}}{
      The vector with two entries: the \code{ms} command line and the three
      random seeds used to run the C program \code{ms} for the multi-locus genomic
      region considered.\cr
      i.e., \code{msframe$cmdline_seeds}. \cr 
    }
    \item{\code{$segsites}}{
      The integer \eqn{S}, the number of segregating sites found in the
      simulated data set for  genomic region \eqn{r}.\cr
      -- In case of a multi-locus genomic region (as specified with
      \eqn{Y>1} in the vector \code{info_region}),
      \code{$segsites}\eqn{=S_n=\Sigma(S_y)\le S}, where \eqn{S_y} is the
      number of segregating sites found in locus \eqn{y}, \eqn{y} \eqn{\in}
      \eqn{[1,Y]} for the multi-locus genomic region \eqn{r}. 
    }
    \item{\code{$positions}}{
      The list of positions in base pairs of the \eqn{S} segregating
      sites found along the genomic region \eqn{r}.\cr
       The position of the segregating site \eqn{j} \eqn{\in}
      \eqn{[1,S]} is \eqn{p_j} \eqn{\in} \eqn{[z_s,z_e]}. \cr
      i.e., \eqn{p_j=}\code{msframe$positions[[1]][}\eqn{j}\code{]}
      \eqn{\times} \eqn{Z}. \cr
      \eqn{Z=z_e-z_s} is the size of the genomic region \eqn{r} as
      specified in the vector \code{info_region}. \cr
      -- In case of a multi-locus genomic region (as specified with
      \eqn{Y>1} in the vector \code{info_region}), 
      the \eqn{S_n} values in \code{$positions[[1]]} are all \eqn{\in}
      \eqn{[z_sy, z_ey]}, \eqn{y} \eqn{\in} \eqn{[1,Y]}.
    }
    \item{\code{$haplotypes}}{
      The vector of \eqn{n=n_1+n_2} haplotypes simulated for the multi-locus
      genomic region \eqn{r}. \cr 
      \eqn{n_1} and \eqn{n_2} are the sample sizes for the multi-locus
      genomic region \eqn{r} specified in the vector \code{info_region}.\cr 
      \code{$haplotypes[}\eqn{1,i}\code{]} is the string of \eqn{S}
      characters corresponding to the alleles found at the \eqn{S} segregating sites
      for the genomic region \eqn{r}, \eqn{i}
      \eqn{\in} \eqn{[1,n]}. \cr
       -- In case of a multi-locus genomic region (as specified with
      \eqn{Y>1} in the vector \code{info_region}),
      \code{msframe$haplotypes[}\eqn{1,i}\code{]} contains only \eqn{S_n}
      characters. \cr 
      -- If for any locus \eqn{y} \eqn{n_py<n_p} (where \eqn{n_1y} and
      \eqn{n_2y} are the sample sizes for the locus \eqn{y} specified in the
      matrix \code{info_loci}), the segregating sites spanning the positions
      of locus \eqn{y} at \eqn{n_p-n_py} randomly chosen haplotypes are noted as
      missing data with \code{"?"}.
    }
    \item{\code{$ns}}{
      The vector of \eqn{S} (or \eqn{S_n})  positions \eqn{s} \eqn{\in}
  \eqn{[1,S]} within the string of characters \code{msframe$haplotypes}
  and the vector \code{msframe$positions} of the \eqn{S} (or \eqn{S_n})
  segregating sites found along the uni- (or multi-)locus genomic region
  \eqn{r}. 
    }
  }
}

\note{
  \item{ ATTENTION: It is the user's responsibility to mind the following
    restrictions:}{
    \item{In the vector \code{info_region}: }{
      -> Reasonable values need to be specified for each genomic region.
    }
    \item{In the matrix \code{info_loci}:}{
      -> The loci numbers, \eqn{y}, need to be in order, i.e., information
      for locus 1 of the first multi-locus genomic region is on line 2, locus
      2 on line 3 \dots locus \eqn{Y} on line \eqn{Y+1}. \cr 
      -> Reasonable values need to be specified for each locus. 
    }
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
  The function \code{\link{simulate_data}} calls the function
  \code{\link{output_ms}}, which in turn calls the function
  \code{\link{out_ms_loci}}. \cr
  \code{\link{data_file}} is the data file example output by the
  function \code{\link{output_ms}}. \cr 
  \code{\link{write_ms_output}} is another function that outputs the
  results from the \R function \code{\link{ms}} into an \code{ms}-like
  text file. \cr
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
msframe

## Create the text file "reg5" is the current directory.
output_ms(datafile="reg",nregion=5, msframe=msframe,
info_region=info_region, info_loci=info_loci)
read.csv("reg5",header =FALSE,sep="\n", comment="#")

# Clean up the directory.
unlink(c("reg5","regions","loci"))
}

\keyword{IO}
\keyword{print}
