\name{liststats}
\alias{liststats}
\docType{data}
\title{
  List of possible summary statistics keywords
}
\description{
  Lists the 59 different keywords for the summary statistics that can be used as
  input by the function \code{\link{input_EstL}} and that are output by the
  function \code{\link{stats_pop}}.
}
\usage{data(liststats)}
\format{
  "pi"       
  [31] "pi_1"      "pi_2"      "theta_W"   "theta_W1"  "theta_W2" 
  [36] "theta_H"   "theta_H1"  "theta_H2"  "D"         "D_1"      
  [41] "D_2"       "H"         "H_1"       "H_2"       "D_star"   
  [46] "D_star1"   "D_star2"   "D_prime"   "D_prime1"  "D_prime2" 
  [51] "r_square"  "r_square1" "r_square2" "nH"        "nH_1"     
  [56] "nH_2"      "Rm"        "Rm_1"      "Rm_2"  
  \describe{
    The statistics are defined assuming that the ancestral allele in
    known. However these statistics are calculated similarly when the
    segregating sites can not be polarized, replacing
    "derived" by "minor allele".
    \item{\strong{Summaries of the polymorphisms}}{
      \describe{
        \item{\emph{\code{S}}}{
          \eqn{S}, the total number of segregating sites found in the data
          set.\cr
          \eqn{S=S_1+S_2+S_s+S_f}.
        }
        \item{\emph{\code{S_1}} and \emph{\code{S_2}}}{
          \eqn{S_1} (\eqn{S_2}), the number of segregating sites with the derived alleles
          polymorphic uniquely in the sample from population 1 (2).
        }
        \item{\emph{\code{S_s}}}{
          \eqn{S_s}, the number of segregating sites with the derived allele shared
          by both population samples.\cr
          \eqn{S_s=S_ss+S_sf1+S_sf2} and \eqn{S_s=S_sl+S_sh}.
        }
        \item{\emph{\code{S_ss}}}{
          \eqn{S_ss}, the number of segregating sites with the derived allele shared
          AND polymorphic in both population samples.
        }
        \item{\emph{\code{S_sf1}} and \emph{\code{S_sf2}}}{
          \eqn{S_sf1} (\eqn{S_sf2}), the number of segregating sites with the derived allele shared
          AND polymorphic ONLY in the sample from population 2 (or 1) (i.e., fixed in
          the sample from population 1 or 2). 
        }
        \item{\emph{\code{S_sl}} and \emph{\code{S_sh}}}{
          \eqn{S_sl} (\eqn{S_sh}), the number of segregating sites with the derived allele shared
          AND at low (high) frequency (i.e., \eqn{\le 10\%} (\eqn{> 10\%})).
        }
        \item{\emph{\code{S_f}}}{
          \eqn{S_f}, the number of segregating sites with the alleles fixed between
          the population sampled. \cr 
          \eqn{S_f=S_f1+S_f2}. 
        }
        \item{\emph{\code{S_f1}} and \emph{\code{S_f2}}}{
          \eqn{S_f1} (\eqn{S_f2}), the number of segregating sites with the derived allele fixed in
          the sample from population 1 (2).
        }
        \item{\emph{\code{S_o}}}{
          \eqn{S_o}, the number of singletons found in the data set. \cr 
          I.e.: \eqn{S_o=S_o1+S_o2} \eqn{\in} \eqn{[0,S]} where
          the frequency of the derived allele is \eqn{=1/(n)}.
        }
        \item{\emph{\code{S_o1}} and \emph{\code{S_o2}}}{
          \eqn{S_o1} (\eqn{S_o2}), the number of singletons found in the sample from population
          1 (2). \cr
          I.e., \eqn{S_o1 \in [0,S_1]} where the
          frequency of the derived allele is \eqn{1/n_1} (\eqn{S_o2 \in [0,S_2]} where the
          frequency of the derived allele is \eqn{1/n_2}). 
        }
      }% end describe S stats
    }% end item sum pol
    \item{\strong{Frequencies of the polymorphisms}}{
      The mean derived allele frequencies for seven types of
      segregating sites described above.
      \describe{
        \item{\emph{\code{F(S)}}}{
          The mean derived allele frequency over \eqn{n} haplotypes for
          the \eqn{S} segregating sites in the data set.
        }
        \item{\emph{\code{F(S_1)}} and \emph{\code{F(S_2)}}}{
          The mean derived allele frequency over \eqn{n_1} (\eqn{n_2}) haplotypes for
          the segregating sites in \eqn{S_1} (\eqn{S_2}).
        }
        \item{\emph{\code{F(S_s)}},
          \emph{\code{F(S_ss)}},\emph{\code{F(S_sf1)}},\emph{\code{F(S_sf2)}},
          \emph{\code{F(S_sl)}} and \emph{\code{F(S_sh)}}}{
          The mean derived allele frequency over \eqn{n} haplotypes for
          the segregating sites in \eqn{S_s}, \eqn{S_ss}, \eqn{S_sf1},
          \eqn{S_sf2}, \eqn{S_sl} and \eqn{S_sh}, respectively. 
        }
      }
    }% end item frequencies
    \item{\strong{Measures of differentiation}}{
      \describe{
        \item{\emph{\code{F_st}}}{
          \eqn{F_ST}, the level of differentiation between the two population
          samples. \cr 
          \eqn{F_ST=1-(H_w1+H_w2)/H_b} (Wright 1931;
          Hudson et al. 1992).
        }
        \item{\emph{\code{H_b}}}{
          \eqn{H_b}, the mean pairwise differences between both population
          samples.
        }
        \item{\emph{\code{H_w1}} and \emph{\code{H_w2}}}{
          \eqn{H_w1} (\eqn{H_w2}),  the mean pairwise differences within the sample from population
          1 (2).
        }
        \item{\emph{\code{Snn}}}{
          \eqn{S_nn}, the nearest-neighbor statistic, which measures of how often the
          nearest neighbors of haplotypes are found in the same population
          sample. \cr 
          \eqn{S_nn=\Sigma[X_i/n]=\Sigma[W_i/(T_i*n)]}, \eqn{i}
          \eqn{\in} \eqn{[1,n]} (Hudson 2000). \cr
          \describe{
            \item{\eqn{W_i}}{ The number of nearest neighbors to haplotype
              \eqn{i} that are from the same population sample as haplotype
              \eqn{i}.} 
            \item{\eqn{T_i}}{The number of nearest neighbors of haplotype
              \eqn{i}.}
            \item{\eqn{X_i}}{The fraction of nearest neighbors of haplotype
              \eqn{i} that are from the same population sample as haplotype
              \eqn{i}. 
            }
          }% end describe snn
        }% end snn
      } % enbd desc dif
    }% end differentiation
    \item{\strong{Estimators of the population mutation rate}}{
      \describe{
        \item{\emph{\code{pi}}, \emph{\code{pi_1}} and
          \emph{\code{pi_2}}}{
          The unbiased estimator of
          the population mutation rate calculated from the average nucleotide
          diversity (or heterozygosity) for the full sample and
          samples from population 1 and 2, \eqn{\pi} (Nei and Li 1979, but
          equations 10 and 12 of Tajima 1989). 
          \itemize{
            \item  \eqn{\pi=\Sigma[2*p*(1-p)*n/(n-1)]} for the segregating
            sites in \eqn{S}. \cr 
            \item \eqn{\pi_1=\Sigma[2*p*(1-p)*n_1/(n_1-1)]} for the segregating
            sites in \eqn{S_1}. \cr 
            \item \eqn{\pi_2=\Sigma[2*p*(1-p)*n_2/(n_2-1)]} for the segregating
            sites in \eqn{S_2}.
          }% end itemize
        }
        \item{\emph{\code{theta_w}}, \emph{\code{theta_W1}} and
          \emph{\code{theta_W2}}}{
          The unbiased estimator of
          the population mutation rate calculated from the number of segregating
          sites, \eqn{\theta_W} (Watterson 1975).
          \itemize{
            \item  \eqn{\theta_W=\Sigma[1/(1/n)]} for the segregating
            sites in  \eqn{S}.\cr 
            \item \eqn{\theta_W1=\Sigma[1/(1/n_1)]} for the segregating
            sites in  \eqn{S_1}. \cr 
            \item \eqn{\theta_W2=\Sigma[1/(1/n_2)]} for the segregating
            sites in  \eqn{S_2}.
          }% end itemize
        }
        \item{\emph{\code{theta_H}}, \emph{\code{theta_H1}} and
          \emph{\code{theta_H2}}}{
          The unbiased estimator of
          the population mutation rate weighted by the homozygosity of the
          derived allele, \eqn{\theta_H} (Fay and Wu 2000).
          \itemize{
            \item  \eqn{\theta_H=\Sigma[2*(p^2)/[n(n-1)]]} for the segregating
            sites in  \eqn{S}. \cr
            sites in  \eqn{S}.\cr 
            \item \eqn{\theta_H1=\Sigma[2*(p^2)/[n_1(n_1-1)]]} for the segregating
            sites in  \eqn{S_1}. \cr 
            \item \eqn{\theta_H2=\Sigma[2*(p^2)/[n_2(n_2-1)]]} for the segregating
            sites in  \eqn{S_2}.
          }% end itemize
        }
      }% end desc estim
    }% end item esti
    \item{\strong{Tests of neutrality}}{
      \describe{
        \item{\emph{\code{D}}, \emph{\code{D_1}} and
          \emph{\code{D_2}}}{
          The test of neutrality
          from the difference between \eqn{\pi} and \eqn{\theta_W}, Tajima's
          \eqn{D} (Tajima 1989).. 
          \itemize{
            \item \eqn{D=(\pi-\theta_W)/\sqrt(var(\pi-\theta_W))}. \cr 
            \item \eqn{D_1=(\pi_1-\theta_W1)/\sqrt(var(\pi_1-\theta_W1))}. \cr 
            \item \eqn{D_2=(\pi_2-\theta_W2)/\sqrt(var(\pi_2-\theta_W2))}.
          }% end itemize
        }
        \item{\emph{\code{H}}, \emph{\code{H_1}} and
          \emph{\code{H_2}}}{
          The test of neutrality
          from the difference between \eqn{\pi} and \eqn{\theta_H}, Fay &
          Wu's \eqn{H} (Fay and Wu 2000). 
          \itemize{
            \item   \eqn{H=\pi-\theta_H=-2*p(2*p-1)n/(n-1)} for the segregating
            sites in  \eqn{S}. \cr
            \item \eqn{H_1=\pi_1-\theta_H1=-2*p(2*p-1)n_1/(n_1-1)} for the segregating
            sites in  \eqn{S_1}. \cr 
            \item \eqn{H_2=\pi_2-\theta_H2=-2*p(2*p-1)n_2/(n_2-1)}
            for the segregating sites in  \eqn{S_2}.
          }% end itemize
        }
        \item{\emph{\code{D_star}}, \emph{\code{D_star1}} and
          \emph{\code{D_star2}}}{
          The  test of neutrality
          from the difference between the number of singletons and the number
          of segregating sites, Fu & Li's \eqn{D*} (Fu & Li 1993). 
          \itemize{
            \item  \eqn{D*=(S-S_o)/\sqrt(var(S-S_o))}.\cr
            \item \eqn{D_1*=(S_1-S_o1)/\sqrt(var(S_1-S_o1))}.\cr 
            \item \eqn{D_2*=(S_2-S_o2)/\sqrt(var(S_2-S_o2))}.
          }% end itemize
        }
      }% end decribe test
    }% end item test
    \item{\strong{Tests of neutrality}}{
      \describe{ 
   \item{\emph{\code{D_prime}}, \emph{\code{D_prime1}} and
          \emph{\code{D_prime2}}}{
        
          \itemize{
            \item  \cr
            \item \cr 
            \item 
          }% end itemize
        }
   \item{\emph{\code{r_square}}, \emph{\code{r_square1}} and
          \emph{\code{r_square2}}}{
        
          \itemize{
            \item  \cr
            \item \cr 
            \item 
          }% end itemize
        }
   \item{\emph{\code{nH}}, \emph{\code{nH_1}} and
          \emph{\code{nH_2}}}{
        
          \itemize{
            \item  \cr
            \item \cr 
            \item 
          }% end itemize
        }
   \item{\emph{\code{Rm}}, \emph{\code{Rm_1}} and
          \emph{\code{Rm_2}}}{
        
          \itemize{
            \item  \cr
            \item \cr 
            \item 
          }% end itemize
        }
        
      }% end decsrib LD
    }% end item ld measure
  }% end describe all stats

}
\details{
  %%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
  %%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
  See \code{References} section in \code{\link{EstLprep-package}}.
}
\examples{
    data(liststats)
    ## The list of possible summary statistics one can use with the keyword stats
    liststats
    }
\keyword{datasets}
