\name{write_ms_output}
\alias{write_ms_output}

\title{
  Function to write the data frame output by the \R function ms into a
  text file with ms-like format. 
}

\description{
  The function \code{\link{write_ms_output}} writes the data frame
  \code{msframe} output by the \R function \code{\link{ms}} into an
  \code{ms}-like text file named \code{msoutput}. 
}

\usage{
write_ms_output(msframe, msoutput="msout")
}

\arguments{
  \item{msframe}{
    The data frame of the kind output by the \R function
    \code{\link{ms}}. \cr 
    \item{\code{msframe$cmdline_seeds}}{
      The vector with two entries: the \code{ms} command line and the three
      random seeds used (or that could have been used) to run the
      \R function \code{\link{ms}}. \cr 
      ATTENTION: The \code{ms} command line, \code{msframe$cmdline_seeds[1]}, needs
      to have a minimum of THREE "words" (The first word can be any
      string of characters without white space): \cr 
      \code{"./ms }\eqn{n}\code{ }\eqn{H}\code{"}. 
      \tabular{ll}{
        \eqn{n} \tab
        : The total number of chromosomes sampled for the genomic region
        considered. \cr 
        \eqn{H} \tab
        : The number of independent data sets described for the genomic region
        considered. \cr 
      }
      \code{msframe$cmdline_seeds[2]} can be either empty or a string of characters
      with the three random seeds. 
    }
    \item{\code{msframe$segsites}}{
      The vector of \eqn{H} integers \eqn{S_h}, the number of segregating
      sites found in the data set \eqn{h}, \eqn{h} \eqn{\in}
      \eqn{[1,H]}. 
    }
    \item{\code{msframe$positions}}{
      The vector of \eqn{H} lists of positions of the segregating sites
      found in the \eqn{H} independent data sets. \cr 
      \tabular{ll}{
        \code{msframe$positions[[}\eqn{h}\code{]]} \tab
        : The list of positions of the \eqn{S_h} segregating sites found
        in the data set \eqn{h} \eqn{\in} \eqn{[1,H]} along the genomic region
        considered. \cr 
        \code{msframe$positions[[}\eqn{h}\code{]][}\eqn{s}\code{]} \tab 
        : The position of the segregating site \eqn{s} \eqn{\in}
        \eqn{[1,S_h]}, \eqn{p_s}.
      }
      \tabular{ll}{
        \eqn{0}\eqn{\le}\eqn{p_s}\eqn{\le}\eqn{1} \tab
        : If the data frame \code{msframe} was output by the \R function
        \code{\link{ms}}, \eqn{p_s} \eqn{\in} \eqn{[0,1]} relative to \eqn{Z},
        the size of the genomic region considered. \cr 
        \eqn{0}\eqn{\le}\eqn{p_s}\eqn{\le}\eqn{Z} \tab 
        : If the data frame \code{msframe} was output by the \R functions
        \code{\link{simulate_data}} or \code{\link{output_ms}} or generated
        from experimental data, \eqn{p_s} \eqn{\in} \eqn{[0,Z]} in base pair. \cr 
      }
    }
    \item{\code{msframe$haplotypes}}{
      The matrix \eqn{H} \eqn{\times} \eqn{n} of the haplotypes for
      the \eqn{H} independent data sets. \cr 
      \tabular{ll}{
        \code{msframe$haplotypes[}\eqn{h,}\code{]} \tab
        : The \eqn{n} haplotypes for the data set \eqn{h}, \eqn{h}
        \eqn{\in} \eqn{[1,H]}. \cr 
        \code{msframe$haplotypes[}\eqn{h,i}\code{]} \tab
        : The string of \eqn{S_h} characters corresponding to the alleles found
        at the \eqn{S_h} segregating sites in the data set \eqn{h},
        \eqn{i} \eqn{\in} \eqn{[1,n]}. \cr 
        \tab i.e., \code{"0"}, \code{"1"} or \code{"?"}, for ancestral (when known),
        derived and missing allele, respectively. 
      }
    }
  }
  \item{msoutput}{
    The name of the text file in which to write the data frame
    \code{msframe} into an \code{ms}-like format. \cr 
    The default name is \code{"msout"}. 
  }
}

\value{
  The function \code{\link{write_ms_output}} creates a text file named as
  specified with the argument \code{msoutput} in the local directory.
  In this text file, \code{\link{write_ms_output}} writes the information
  from the data frame \code{msframe} in the same format as output by the \R
  functions \code{\link{simulate_data}}, \code{\link{write_ms_output}},
  \code{\link{output_ms}} and by Hudson's C program \code{ms}. \cr 
  \item{}{\strong{FILE HEADER:}
    \item{Line 1}{  
      The \code{ms} command line used (or that could have been used) to
      run the C program \code{ms}. \cr 
      i.e., \code{msframe$cmdline_seeds[1]}.} 
    \item{Line 2}{  
      The three random seeds used to run the C program \code{ms} if the
      data were simulated. \cr
      i.e., \code{msframe$cmdline_seeds[2]}.} 
    \item{Line 3}{
      An empty line separates the information of the \eqn{H} independent
      data sets.
    }
  }
  \item{}{\strong{INFORMATION FOR EACH INDEPENDENT DATA SET: } \cr
    The line count below corresponds to the lines within a block of
    information for the data set \eqn{h}: 
    \item{Line 1}{  
      \code{"\\\\"} to specify that the data set \eqn{h} \eqn{\in}
      \eqn{[1,H]} starts below.
    } 
    \item{Line 2}{   
      \code{"}\code{segsites: }\eqn{S_h}\code{"}, the number of segregating
      sites found in the data set \eqn{h}. \cr
      i.e., \code{msframe$segsites[}\eqn{h}\code{]}.
    }
    \item{Line 3}{   
      \code{"}\code{position: }\code{"} followed by the list of positions
      of the \eqn{S_h} segregating sites found in the data set
      \eqn{h} along the genomic region considered. \cr
      The position of the segregating site
      \eqn{s} \eqn{\in} \eqn{[1,S_h]}, \eqn{p_s}, is the value
      \code{msframe$positions[[}\eqn{h}\code{]][}\eqn{s}\code{]}.
    }
    \item{Lines 4  to 4\emph{+n}}{  
      \eqn{n} lines follow corresponding to the \eqn{n} haplotypes in the
      data set \eqn{h}. \cr 
      Each line \eqn{j=i+4} \eqn{\in} \eqn{[4,4+n]} of
      haplotype is a string of \eqn{S_h} characters corresponding to the
      alleles found at the \eqn{S_h} segregating sites in the data set
      \eqn{h}.\cr 
      i.e., \code{msframe$haplotypes[}\eqn{1,i}\code{]}, \eqn{i}
      \eqn{\in} \eqn{[1,n]}, with the values \code{"0"}, \code{"1"} or
      \code{"?"}, for ancestral (when known), derived and missing allele,
      respectively. 
    }
  }
}

\note{
  \item{ ATTENTION: It is the user's responsibility to mind the following
    restrictions:}{
    -> The data frame \code{msframe} needs to be of the same format as
    output by the \R functions \code{\link{ms}} and
    \code{\link{get_ms_output}}. \cr
    -> The \code{ms}-like command line, \code{cmdline_seeds[1]}, needs to
    have a minimum of THREE "words" (The first word can be any string of
    characters without white space): \cr 
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
  \code{\link{output_ms}} is another function that writes the data frame
  from the \R function \code{\link{ms}} into a \code{ms}-like text
  file. \cr
  \code{\link{data_file}} is the data file example output by the
  function \code{\link{write_ms_output}}. \cr 
  \code{\link{get_ms_output}}, in turn, is a function that reads from an
  \code{ms}-like text file and created a data frame of the kind output by
  the \R function \code{\link{ms}}. \cr
  Lists of definitions of the symbols and parameters mentioned in
  this file are found in \code{\link{Rmspack-package}}.
}

\examples{
### Create the input for the function.
msline="./ms 8 2 -t 1 -r 0.1 9999 -I 2 4 4 10 -ej 1 1 2 -eN 1 3"
msframe=ms(msline)
msframe

write_ms_output(msframe=msframe, msoutput="msout")
read.csv("msout",header =FALSE,sep="\n", comment="#")

# Clean up the directory.
unlink("msout")
}

\keyword{IO}
\keyword{print}
