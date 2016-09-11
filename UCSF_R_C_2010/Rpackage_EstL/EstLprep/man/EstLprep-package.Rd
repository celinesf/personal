\name{EstLprep-package}
\alias{EstLprep-package}
\alias{EstLprep}
\docType{package}
\title{
  EstLprep: A package to prepare input files required by EstL
}
\description{
  \code{\link{EstLprep}} allows to generate the input files used by the C program of estimation \code{EstL}.
}

\details{
  \tabular{ll}{
 	 Package: \tab EstLprep\cr
	 Type: \tab Package\cr
	 Version: \tab 1.0\cr
	 Date: \tab 2011-03-29\cr
	 License: \tab \cr
	 LazyLoad: \tab no\cr
  }
  \code{\link{liststats}} lists the summary statistics that can be considered in the
  estimation program \code{EstL} and that can be used as input by
  \code{\link{input_EstL}}.  
}

\section{Symbols}{
  The following tables list the definitions of the symbols mentioned in
  this file.
  \describe{
    \item{\bold{\strong{List of recuring symbols below}}}{
      \describe{
        \item{\eqn{n} and \eqn{n_i}}{The total sample size and the
          sample size from population \eqn{i} for the genomic region considered,
          \eqn{i \in [1,2]}, respectively.
        }
        \item{\eqn{N_i}}{ The effective population size in population
          \eqn{i} (by default and unless specified, \eqn{N_1=N_2=N_A}).
        } 
        \item{\eqn{R}}{The number of independent genomic regions
          considered. 
        }
        \item{\eqn{Z}}{
          The size in bp of the genomic region considered.
        }
        \item{\eqn{Y}}{
          The number of loci spanning the genomic region considered.
        }
        \item{\eqn{p}}{The frequency of an allele (usually the derived
          allele) at a segregating sites in a specific sample (this
          number may include missing data).}
      }
    }% end of recuring symbols
    
    \item{\bold{\strong{Population mutation rate parameters, \eqn{\theta}}}}{
      \itemize{
        \item \describe{
          \item{\eqn{\theta_i = 4N_i \mu}}{ 
            The generational population
            mutation rate when specified per base pair.}
        }
        \item \describe{
          \item{\eqn{\theta_i = 4N_i \mu xvZ}}{
            The region-specific generational population mutation rate when
            specified for the genomic region considerd.
          }
          \item{Where:}{
            \describe{
              \item{\eqn{\mu}}{ The genomic generational mutation rate per bp.
              }
              \item{\eqn{x}}{   The inheritance scalar for the genomic region considered (i.e.,
                \code{"1"} for autosomal region, \code{"0.75"} for X- and \code{"0.5"}
                for Y- and mtDNA-linked region).
              }
              \item{ \eqn{v}}{  The mutation rate scalar for the genomic region considered (which
                can be estimated e.g., from divergence data).
              }
            }
          }% end where
        }% end desc
      }% end item
    }% end mut
    \item{\bold{\strong{Population recombination rate parameters,
          \eqn{\rho}}}}{
      \itemize{
        \item \describe{
          \item{\eqn{\rho = 4N_1 c}}{
            The generational genomic average  intra-region population
            recombination rate when specified per base pair.
        }}
        \item \describe{
          \item{\eqn{\rho = 4N_1 c w(Z-1)}}{
            The generational region-specific intra-region population
            recombination rate when specified for the genomic region
            considerd.
        }}
        \item \describe{
         \item{\eqn{\rho_o}}{
           The estimate from linkage disequilibrium analysis of the
           generational region-specific population recombination
           rate per bp for the genomic region considered .}}
         \item \describe{
         \item{\eqn{c_o}}{ 
           The estimate from pedigree analysis of the generational
           region-specific cross-over rate per bp 
           for the genomic region considered. }}
       \item \describe{
         \item{\eqn{\aplha=c_r/\mu}}{
           The ratio of recombination over mutation
           rate per bp for the genomic region \eqn{r}.
           \eqn{\alpha} is typically drawn from a distribution.
         }
         \item{Where:}{
          \describe{
            \item{\eqn{c}}{
              The generational cross-over rate per bp.
            }
            \item{\eqn{c_r}}{
              The generational region-specific cross-over rate per bp
              for the genomic region \eqn{r}. 
            }
            \item{\eqn{w}}{
              The recombination scalar for the genomic region
              considered. Usually \eqn{w=\beta}.
            }
            \item{\eqn{\beta}}{
              The ratio of the region-specific population recombination
              rate per bp over the genomic average \eqn{\rho}
              per bp for the genomic region considered.\cr
              I.e., \code{"1"} (\code{"0.5"} in \emph{Drosophila}) for
              autosomal region, \code{"0.5"} for X- and \code{"0"} for Y- and
              mtDNA-linked region.
            }
          }% end descb rho symbols
        }% end where
        }% end desc
      }% end itemize rho
    }% end RHO
    
    \item{ \bold{\strong{Parameters of time of events}}}{
      \describe{
        \item{\eqn{T_s}}{
          The split time in unit of \eqn{4N_1} generations between the
          TWO populations.}
        \item{\eqn{T_c}}{
          The time at which the rate of gene flow changed between the
          two populations in unit of \eqn{4N_1} generations.}
        \item{\eqn{\epsilon=T_c/T_s}}{} 
      }% end describe
    }% end time
    
    \item{ \bold{\strong{Parameters of gene flow rates}}}{
      \describe{
        \item{\eqn{M_p=4N_1 m_p}}{
          The number of migrants exchanged each generation by the TWO
          populations at present. \cr 
          \itemize{
            \item \eqn{M_p} is the symmetrical rate of gene flow between
            two population in an island model (i.e., \eqn{0<t<\inf}) when:\cr
            \itemize{
              \item A population structure is specified.
              \item \eqn{T_s=0}.
            }
            \item \eqn{M_p} is the constant symmetrical rate of gene flow
            since the split until present (i.e., \eqn{0<t<T_s}) if there is a
            population split (as specified with \eqn{T_s>0}).\cr 
            \item \eqn{M_p} is the constant symmetrical rate of gene flow
            since the time of gene flow rate change until present (i.e.,
            \eqn{0<t<T_c}) if a time at which the gene flow rate changed is
            specified with:\cr
            \itemize{
              \item \eqn{0<T_c<T_s}. \cr 
              \item \eqn{M_c \not=M_p}. \cr 
            }
          }% end Mp
        }%end desc Mp
        \item{\eqn{M_c=4N_1 m_c}}{
          The number of migrants exchanged each generation by the TWO
          populations since the split until the time of gene flow rate
          change (i.e., \eqn{T_c<t<T_s}). \cr 
        }
        \item{Where}{
          \describe{
            \item{\eqn{m_p}}{ 
              The generational fraction of migrant individuals at present.}
            \item{\eqn{m_c}}{
              The generational fraction of migrant individuals between
              \eqn{T_c<t<T_s}.}   
          }
        }
      }% end describ M<>
    }
  }% end itemize
}% end symbols


\author{
  Maintainer:  Celine Becquet - \email{celine.becquet@gmail.com}. 
}

\references{
  Becquet, C., and M. Przeworski. 2007. A new approach to estimate
  parameters of speciation models with application to apes. Genome
  Res. 17:1505-1519. \cr 
  Fay, J. C., and C. -I. Wu. 2000. Hitchhiking Under Positive Darwinian
  Selection. Genetics 155:1405-1413. \cr 
  Fu, Y. -X. and W. -H. Li. 1993. Statistical tests of neutrality of
  mutations. Genetics 133:693-709. \cr 
  Hudson, R. R., 1983. Properties of a neutral allele model with
  intragenic recombination. Theor. Popul. Biol. 23:183-201. \cr 
  Hudson, R. R., 1990. Gene genealogies and the coalescent process, in
  D. Futuyma and J. Antonovics (eds), Oxford Surveys in Evolutionary
  Biology, Vol. 7:1-44. \cr 
  Hudson, R. R., Slatkin, M., and Maddison, W. P., 1992. Estimation of
  levels of gene flow from DNA sequence data. Genetics, 132:583-589.\cr
  Hudson, R. R., 2000. A New Statistic for Detecting Genetic
  Differentiation. Genetics 155:2011-2014. \cr 
  Hudson, R. R., 2002. Generating samples under a Wright-Fisher neutral
  model of genetic variation. Bioinformatics 18:337-338. \cr 
  \url{http://pps-spud.uchicago.edu/~cbecquet/msdoc.pdf}. \cr 
  Nei, M., and W. H. Li. 1979. Mathematical model for studying genetic
  variation in terms of restriction
  endonucleases. Proc. Natl. Acad. Sci. USA 76:5269-5273. \cr 
  Tajima, F., 1989. Statistical method for testing the neutral mutation
  hypothesis by DNA polymorphism. Genetics 123:585-595. \cr 
  Wakeley, J., and J. Hey. 1997. Estimating ancestral population
  parameters. Genetics 145:847-855. \cr 
  Watterson, G. A., 1975. On the number of segregating sites in
  genetical models without
  recombination. Theor. Popul. Biol. 7:256-276. \cr 
  Wright, S., 1931. Evolution in Mendelian populations. Genetics
  16:97-159. 
}

\keyword{ package }
\seealso{
  \describe{
    \item{\bold{\strong{Functions in \code{\link{EstLprep}}}}}{
      \describe{
        \item{\code{\link{cmd_line}}}{
        }
      }% end describe functions
    }% end item functions

    \item{\bold{\strong{Example of data files for \code{\link{EstLprep}}}}}{
      \describe{
        \item{\code{\link{regfile}}}{
        }
      }% end describe functions
    }% end item functions
  }% end describe Calso
}% end c also
\examples{
print("hello")
}
