\name{get_ms_output}
\alias{get_ms_output}

\title{
  Function to transform a text file with data in ms-like format into a
  data frame in the format as output by the \R function ms. 
}

\description{
  The function \code{\link{get_ms_output}} transforms text file in the
  format as output by Hudson's C program \code{ms} into the data frame
  output by the \R function \code{\link{ms}}. \cr 
  This function is called by the \R function \code{\link{estimate_IMc}}. 
}

\usage{
get_ms_output(msoutput="msout")
}

\arguments{
  \item{msoutput}{
    The name of the text file containing the data (either simulated or from
    resequencing) in \code{ms}-like format. \cr 
    By default the file name is \code{"msout"}. \cr
    The file \code{msoutput} is in the same format as output by the \R functions
    \code{\link{simulate_data}}, \code{\link{write_ms_output}},
    \code{\link{output_ms}} and by Hudson's C program
    \code{ms}.\cr
    \strong{FILE HEADER:}
    \item{Line 1}{ The \code{ms} command line used (or that could have been used) to
      run the C program \code{ms}. \cr
      ATTENTION: This string of character needs to have a
      minimum of THREE "words" (The first word can be any string of
      characters without white space): \cr 
      \code{"./ms }\eqn{n}\code{ }\eqn{H}\code{"}. \cr
      \tab \tab \eqn{n} \tab
      : The total number of chromosomes sampled for the genomic region
      considered. \cr 
      \tab \tab \eqn{H} \tab
      : The number of independent data sets described in the file for the
      genomic region considered.
    }
    \item{ Line 2}{
      The three random seeds used to run the C program \code{ms} if the
      data were simulated. \cr 
      This line is not used and can be empty.
    }
    \item{Line 3}{
      An empty line separates the information of the \eqn{H} independent
      data sets.
    }
    \strong{INFORMATION FOR EACH INDEPENDENT DATA SET: }\cr 
    The line count below corresponds to the lines within a block of
    information for the data set \eqn{h}:
    \item{Line 1}{\code{"\\\\"} to specify that the data set \eqn{h} \eqn{\in}
      \eqn{[1,H]} starts below.
    }
    \item{Line 2}{
      \code{"}\code{segsites: }\eqn{S_h}\code{"}, the number of segregating
      sites found in the data set \eqn{h}.
    }
    \item{Line 3}{
      \code{"}\code{position: }\code{"} followed by the list of positions
      of the \eqn{S_h} segregating sites found in the data set
      \eqn{h} along the genomic region considered. \cr
      \tab \tab \tab The position of the segregating site \eqn{s} \eqn{\in}
      \eqn{[1,S_h]} is \eqn{p_s}:
      \tabular{ll}{
        \eqn{0}\eqn{\le}\eqn{p_s}\eqn{\le}\eqn{1}   \tab 
        : If the file \code{msoutput} was output by the C program \code{ms} or by
        the \R function \code{\link{write_ms_output}} (i.e., from a data frame
        output by the \R function \code{\link{ms}}), \eqn{p_s} \eqn{\in}
        \eqn{[0,1]} relative to \eqn{Z}, the size of the genomic region
        considered. \cr 
        \eqn{0}\eqn{\le}\eqn{p_s}\eqn{\le}\eqn{Z} \tab 
        If the file \code{msoutput} was output by the \R functions
        \code{\link{simulate_data}} or \code{\link{output_ms}} or generated
        from experimental data, \eqn{p_s} \eqn{\in} \eqn{[0,Z]} in base
        pair.
      }
    }
    \item{Lines 4  to 4\emph{+n}}{
      \eqn{n} lines follow corresponding to the \eqn{n} haplotypes
      sampled in the data set \eqn{h}. \cr 
      Each line \eqn{j=i+4} \eqn{\in} \eqn{[4,4+n]} of
      haplotype is a string of \eqn{S_h} characters corresponding to the
      alleles found at the \eqn{S_h} segregating sites in the data set
      \eqn{h}.\cr 
      i.e., \code{"0"}, \code{"1"} or \code{"?"}, for
      ancestral (when known), derived and missing allele, respectively.
    }
  }
}

\value{
  The function \code{\link{get_ms_output}} outputs a data frame of the
  same format as output by the \R function \code{\link{ms}}: 
  \item{$cmdline_seeds}{
    The vector with two entries: the \code{ms} command line and the
    three random seeds used (or that could have been used) to run the C
    program \code{ms}. \cr
    i.e., lines 1 and 2 of the file \code{msoutput}. 
  }
  \item{$segsites}{
    The vector of \eqn{H} integers \eqn{S_h}, the number of segregating
    sites found in the data set \eqn{h}, \eqn{h} \eqn{\in} \eqn{[1,H]}. \cr 
    i.e., line 2 of the block of information for the data set
    \eqn{h}. \cr
    \eqn{H} is specified in the \code{ms} command line.
  }
  \item{$positions}{
    The vector of \eqn{H} lists of positions of the segregating sites
    found in the \eqn{H} independent data sets.
    \tabular{ll}{
      \code{$positions[[}\eqn{h}\code{]]} \tab
      : The list of positions of the \eqn{S_h} segregating sites found in the
      data set \eqn{h} \eqn{\in} \eqn{[1,H]} along the genomic
      region considered. \cr
      \code{$positions[[}\eqn{h}\code{]][}\eqn{s}\code{]} \tab      
      : The position of the segregating site \eqn{s} \eqn{\in} \eqn{[1,S_h]}. \cr 
      \tab i.e., the value \eqn{p_s} in the line 3 of the block of
      information for the data set \eqn{h} of the file \code{msoutput}.
    }
  }
  \item{$haplotypes}{
    The matrix \eqn{H} \eqn{\times} \eqn{n} of the haplotypes for the
    \eqn{H} independent data sets.
     \tabular{ll}{
     \code{$haplotypes[}\eqn{h,}\code{]}\tab
      : The \eqn{n} haplotypes for the data set \eqn{h}, \eqn{h} \eqn{\in}
      \eqn{[1,H]}.\cr 
     \code{$haplotypes[}\eqn{h,i}\code{]}\tab
      : The string of \eqn{S_h} characters corresponding to the alleles found
      at the \eqn{S_h} segregating sites in the data set \eqn{h},
      \eqn{i} \eqn{\in} \eqn{[1,n]}. \cr 
      \tab i.e., line 4+\eqn{i} of the block of information for the data set
      \eqn{h} of the file \code{msoutput}, with the values \code{"0"},
      \code{"1"} or \code{"?"}, for ancestral (when known), derived and
      missing allele, respectively. 
    }
  }
}

\note{
  \item{ ATTENTION: It is the user's responsibility to mind the following
    restrictions:}{
    -> The user needs to be comfortable with the use of the C program
    \code{ms} and needs to be familiar with the documentation at
    \url{http://pps-spud.uchicago.edu/~cbecquet/msdoc.pdf}.\cr
    -> The \code{ms}-like command line (i.e., line 1 of the file
    \code{msoutput}) needs to have a minimum of THREE "words" (The first
    word can be any string of characters without white space): \cr
    \code{"./ms }\eqn{n}\code{ }\eqn{H}\code{"}. 
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
  The function \code{\link{estimate_IMc}} calls the function
  \code{\link{get_ms_output}}, which in turn calls the function
  \code{\link{error_message}}. \cr
  \code{\link{data_file}} is the data file example used as input by the
  function \code{\link{get_ms_output}}. \cr 
  The following functions write data frame output by the \R function
  \code{\link{ms}} into text file in \code{ms}-like format: \cr 
  \code{\link{write_ms_output}} and \code{\link{output_ms}}.\cr
  Lists of definitions of the symbols and parameters mentioned in
  this file are found in \code{\link{Rmspack-package}}. 
}

\examples{
### Create the input for the function.
data(simulation_files) # download the data

## Write the file of information on the parameters in the local directory.
write.table(file="msout", x=data_file, row.name=FALSE, col.names=FALSE, quote=FALSE) 
# Creates the file "msout" containing
scan("msout",sep="\n", blank.lines.skip=FALSE,what=("c"))

get_ms_output(msoutput="msout")

# Clean up the directory.
unlink("msout")
}

\keyword{IO}
\keyword{misc}
