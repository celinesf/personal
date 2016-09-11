\name{data_file}
\alias{data_file}
\docType{data}

\title{
  Example of a test file with the data for a region in an ms-like format.
}
\description{
  The file \code{\link{data_file}} provides an example of the file
  \code{datafile}, which contains the data for a genomic region in the
  format as output by Hudson's C program \code{ms} required to estimate
  parameters of extensions of the Isolation-with-migration model with
  the \R function \code{\link{estimate_IMc}} or to transform into the data
  frame output by the \R function \code{\link{ms}} with the \R function
  \code{\link{get_ms_output}}.\cr 
  Text files in the format of \code{\link{data_file}} are output by the \R
  functions:\cr
  \code{\link{output_ms}}, \code{\link{simulate_data}} and
  \code{\link{write_ms_output}}. 
}

\usage{
data(data_file)
}

\format{
  The file \code{\link{data_file}} contains an example of the parameter file
  \code{datafile} (or \code{msoutput}), which contains the data for
  a genomic region in the usual \code{ms} format (with some variations
  from Hudson's program) required to estimate parameters of extensions
  of the Isolation-with-migration model with the \R function
  \code{\link{estimate_IMc}} (or to transform into the data frame output
  by the \R function \code{\link{ms}} with the \R function
  \code{\link{get_ms_output}}):
  \item{Line 1}{
    The \code{ms} command line used (or that could have been used) to
    run the C program \code{ms}.\cr 
    ATTENTION: This string of character needs to
    contain at minimum these SEVEN "words" to use the function
    \code{\link{estimate_IMc}} (only the THREE first "words" are
    required for the function \code{\link{get_ms_output}}, the first
    word can be any string of characters without white space): \cr
    \code{"./ms }\eqn{n}\code{ 1 -I 2 }\eqn{n_1}\code{ }\eqn{n_2}\code{"}.
    \tabular{ll}{
      \eqn{n} \tab 
      : The total number of chromosomes sampled for the genomic region
      considered. \eqn{n=n_1+n_2}. \cr 
      \eqn{n_1} \tab
      : The sample size from population 1 for the genomic region
      considered. \cr 
      \eqn{n_2} \tab
      : The sample size from population 2 for the genomic region
      considered.
    }
  }
  \item{Line 2}{
    The three random seeds used to run the C program \code{ms} if the
    data were simulated. \cr 
    This line is not used and can be empty. 
  }
  \item{Line 3}{
    An empty line.
  }
  \item{Line 4}{
    \code{"\\\\"} to specify that the data set for the genomic region
    considered starts below.
  }
  \item{Line 5}{
    \code{"}\code{segsites: } \eqn{S}\code{"}, the number of
    segregating sites found in the data set for the genomic region
    considered.
  }
  \item{Line 6}{
    \code{"}\code{position: }\code{"} followed by the \eqn{S}
    positions in base pairs of the \eqn{S} segregating sites found in the
    data set along the genomic region considered. \cr 
    The position of the segregating site \eqn{s}
    \eqn{\in} \eqn{[1,S]} is \eqn{p_s} \eqn{\in} \eqn{[0,Z]}. \cr 
    \eqn{Z} is the size of the genomic region
    considered.
  }
  \item{Line 7 to 7\emph{+n}}{
    \eqn{n=n_1+n_2} lines follow corresponding to the \eqn{n}
    haplotypes sampled for the genomic 
    region considered. \cr 
     \eqn{n_1} and \eqn{n_2} are the sample sizes in the
    two populations for the genomic region considered specified in the file
    \code{regfile} (see \code{\link{info_region}}). \cr 
     Each line \eqn{j=i+7} \eqn{\in} \eqn{[7,7+n]} of
    haplotype is a string of \eqn{S} characters corresponding to the
    alleles found at the \eqn{S} segregating sites in the data set
    for the genomic region considered. \cr 
     i.e., \code{"0"}, \code{"1"} or \code{"?"}, for
    ancestral (when known), derived and missing allele, respectively. 
  }
}

\note{
  \item{ ATTENTION: It is the user's responsibility to mind the following
    restrictions:}{ 
    -> The user needs to be comfortable with the use of the C program
    \code{ms} and needs to be familiar with the documentation at
    \url{http://pps-spud.uchicago.edu/~cbecquet/msdoc.pdf}. \cr 
    -> The \code{ms}-like command line (i.e., line 1 of the file
    \code{datafile}) needs to contain at minimum these SEVEN "words" to
    use the function \code{\link{estimate_IMc}} (only the THREE first
    "words" are required for the function \code{\link{get_ms_output}},
    the first word can be any string of characters without white space):
    \cr
    \code{"./ms }\eqn{n}\code{ 1 -I 2 }\eqn{n_1}\code{ }\eqn{n_2}\code{"}.
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
  Text files in the format of \code{\link{data_file}} are used as input for
  the functions \code{\link{estimate_IMc}} and \code{\link{get_ms_output}}.\cr 
  Text files in the format of \code{\link{data_file}} are output by the \R
  functions:\cr
  \code{\link{output_ms}}, \code{\link{simulate_data}} and
  \code{\link{write_ms_output}}.\cr 
  Other data file examples used as input by the function
  \code{\link{estimate_IMc}}: \cr 
  \code{\link{info_region}}, \code{\link{info_loci}},
  and \code{\link{param_est}}. \cr
  Lists of definitions of the symbols and parameters mentioned in
  this file are found in \code{\link{Rmspack-package}}.
}

\examples{
### Create the data file.
data(simulation_files) # download the data
data_file

## Write the file of information on the parameters in the local directory.
write.table(file="msout", x=data_file, row.name=FALSE, col.names=FALSE, quote=FALSE)
# Creates the file "msout" containing:
scan("msout",sep="\n", blank.lines.skip=FALSE,what=("c"))

# Clean up the directory.
unlink("msout")
}

\keyword{datasets}
