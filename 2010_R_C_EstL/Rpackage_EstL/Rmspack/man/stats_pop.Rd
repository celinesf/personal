\name{stats_pop}
\alias{stats_pop}

\title{
  Function to calculate summary statistics from a data frame as output
  by the \R function ms for a model with TWO populations. 
}

\description{
  The function \code{\link{stats_pop}} calculates various summary
  statistics for the data frame \code{msframe} of the same format as
  output by the \R function \code{\link{ms}} for a model with TWO
  populations.\cr
  This function is called by the \R functions
  \code{\link{simulate_data}} and \code{\link{estimate_IMc}}. 
}

\usage{
stats_pop(msframe, type="anc")
}

\arguments{
  \item{msframe}{
    The data frame of the kind output by the \R function \code{\link{ms}}
    for a model with TWO populations. \cr 
    The following information are required in the data frame
    \code{msframe}: 
    \item{\code{msframe$cmdline_seeds}}{
      The vector with two entries: the \code{ms} command line and the
      three random seeds used (or that could have been used) to run the
      \R function \code{\link{ms}}. \cr 
      ATTENTION: The \code{ms} command line, \code{msframe$cmdline_seeds[1]},
      needs to contain at minimum these SEVEN "words" (The first word
      can be any string of characters without white space): \cr 
      \code{"./ms }\eqn{n}\code{ }\eqn{H}\code{ -I 2 }\eqn{n_1}\code{ }\eqn{n_2}\code{ "}. 
      \tabular{ll}{
        \eqn{n} \tab
        : The total number of chromosomes sampled for the genomic region
        considered. \eqn{n=n_1+n_2}. \cr
        \eqn{n_1} \tab
        : The sample size from population 1 for the genomic region
        considered. \cr 
        \eqn{n_2} \tab
        : The sample size from population 2 for the genomic region
        considered.\cr 
        \eqn{H} \tab
        : The number of independent data sets described for the genomic
        region considered. 
      }
      \code{msframe$cmdline_seeds[2]} can be either empty or a string of characters
      with the three random seeds. 
    }
    \item{\code{msframe$segsites}}{
      The vector of \eqn{H} integers \eqn{S_h}, the number of
      segregating sites found in the data set \eqn{h}, \eqn{h} \eqn{\in}
      \eqn{[1,H]}. 
    }
    \item{\code{msframe$haplotypes}}{
      The matrix \eqn{H} \eqn{\times} \eqn{n} of the haplotypes for the
      \eqn{H} independent data sets. \cr 
      \tabular{ll}{
        \code{msframe$haplotypes[}\eqn{h,}\code{]} \tab
        : The \eqn{n} haplotypes for the data set \eqn{h}, \eqn{h}
        \eqn{\in} \eqn{[1,H]}. \cr 
        \code{msframe$haplotypes[}\eqn{h,i}\code{]} \tab
        : The string of \eqn{S_h} characters corresponding to the
        alleles found at the \eqn{S_h} segregating sites in the data set
        \eqn{h}, \eqn{i} \eqn{\in} \eqn{[1,n]}. \cr
        \tab i.e., \code{"0"}, \code{"1"} or \code{"?"}, for ancestral
        (when known), derived and missing allele, respectively. 
      }
    }
  }% end msframe 
  \item{type}{ 
    The string of three characters specifying whether the ancestral
    alleles of the genomic region considered is known when the function
    \code{\link{stats_pop}} calculates the summary statistics for
    each data set. \cr 
    \code{type} can take the two following values. 
    \tabular{ll}{
      \code{"anc"} \tab
      : By default \code{type="anc"}. \cr
      \tab \code{\link{stats_pop}} calculates the summary statistics for the
      data set \eqn{h} considering the character \code{"0"} at segregating
      site \eqn{s} as the ancestral allele (as specified in
      \code{msframe$haplotypes[}\eqn{h,}\code{]}), \eqn{s} \eqn{\in} \eqn{[1,S_h]},
      \eqn{h} \eqn{\in} \eqn{[1,H]}. \cr 
      \tab In turn, \code{\link{stats_pop}} considers the character
      \code{"1"} as the derived allele at segregating site \eqn{s}. \cr 
      \code{"maf"} \tab
      : \code{\link{stats_pop}} considers the allele with the minimum
      frequency at segregating site \eqn{s} as the derived allele
      (regardless of the character specified in
      \code{msframe$haplotypes[}\eqn{h,}\code{]}). 
    }
  }
}

\details{
  The function \code{\link{stats_pop}} takes into account and ignores
  missing data and makes sure that each site is indeed segregating in
  the full sample and the population samples. 
}

\value{
  The function \code{\link{stats_pop}} outputs vectors with \eqn{H}
  items for various summary statistics calculated for each of the
  \eqn{H} independent data sets. \cr 
  For the data set \eqn{h} \eqn{\in} \eqn{[1,H]}, the function
  \code{\link{stats_pop}} provides the following summary statistics: 
  \item{$type}{
    The sentence describing the type of summary statistics calculated:
    \cr 
    \tabular{ll}{
      \code{"Known ancestral allele"}\tab
      : Specifies that \code{\link{stats_pop}} calculated the summary statistics
      for the data sets considering the character \code{"0"} at each
      segregating site (as specified in \code{msframe$haplotypes}) as the
      ancestral allele. \cr 
      \tab  In turn, \code{\link{stats_pop}} considered the character
      \code{"1"} as the derived allele at each segregating site. \cr 
      \code{"Unknown ancestral allele"} \tab     
      : Specifies that \code{\link{stats_pop}} considered the allele with
      the minimum frequency at each segregating site as the derived allele
      (regardless of the character specified in \code{msframe$haplotypes}). 
    }
  } 
  \item{$segsites}{
    The vector of twelve summary statistics on the number of segregating sites in
    the full sample and the populations samples (Wakeley and Hey 1997;
    Becquet and Przeworski 2007): 
    \tabular{lll}{
      \code{$segsites[1]} \tab \eqn{S} \tab 
      : The total number of segregating sites found in the data set 
      \eqn{h}. \cr
      \tab \tab i.e., \code{msframe$segsites[h]}. \cr
      \tab \tab \eqn{S=S_1+S_2+S_s+S_f}. \cr
      \code{$segsites[2]=} \tab \eqn{S_1} \tab 
      : The number of segregating sites with the derived alleles
      polymorphic uniquely in the sample from population 1. \cr 
      \code{$segsites[3]=} \tab \eqn{S_2} \tab 
      : The number of segregating sites with the derived alleles
      polymorphic uniquely in the sample from population 2. \cr 
      \code{$segsites[4]=} \tab \eqn{S_s} \tab 
      : The number of segregating sites with the derived allele shared
      by both population samples. \cr 
      \tab \tab \eqn{S_s=S_ss+S_sfs+S_sf2} and \eqn{S_s=S_sl+S_sh}. \cr
      \code{$segsites[5]=} \tab \eqn{S_ss} \tab 
      : The number of segregating sites with the derived allele shared
      AND polymorphic in both population samples. \cr 
      \code{$segsites[6]=} \tab \eqn{S_sf1} \tab 
      : The number of segregating sites with the derived allele shared
      AND polymorphic ONLY in the sample from population 2 (i.e., fixed in
      the sample from population 1). \cr 
      \code{$segsites[7]=} \tab \eqn{S_sf2} \tab
      : The number of segregating sites with the derived allele shared
      AND polymorphic ONLY in the sample from population 1 (i.e., fixed in
      the sample from population 2). \cr
       \code{$segsites[8]=} \tab \eqn{S_sl} \tab 
      : The number of segregating sites with the derived allele shared
      AND  at low frequency (i.e., \eqn{p_s \le 10\%}). \cr 
      \code{$segsites[9]=} \tab \eqn{S_sh} \tab
      : The number of segregating sites with the derived allele shared
      AND  at high frequency (i.e., \eqn{p_s > 10\%}. \cr 
      \code{$segsites[10]=} \tab \eqn{S_f} \tab 
      : The number of segregating sites with the alleles fixed between
      the population sampled. \cr 
      \tab \tab \eqn{S_f=S_f1+S_f2}. \cr
      \code{$segsites[11]=} \tab \eqn{S_f1} \tab 
      : The number of segregating sites with the derived allele fixed in
      the sample from population 1. \cr 
      \code{$segsites[12]=} \tab \eqn{S_f2} \tab 
      : The number of segregating sites with the derived allele fixed in
      the sample from population 2. 
    }
  }
  \item{$singletons}{
    The vector of the three summary statistics on the number of singletons in
    the full sample and the populations samples: 
    \tabular{ll}{
      \code{$singletons[1]=}\eqn{S_o} \tab
      : The number of singletons found in the data set \eqn{h}. \cr 
      \tab i.e.: \eqn{S_o=S_o1+S_o2} \eqn{\in} \eqn{[0,S]}, where
      the frequency of the derived allele is \eqn{p_s=1/n_s}. \cr 
      \code{$singletons[2]=}\eqn{S_o1} \tab
      : The number of singletons found in the sample from population
      1. \cr 
      \tab i.e.: \eqn{S_o1} \eqn{\in} \eqn{[0,S_1]}, where the
      frequency of the derived allele is \eqn{q_s=1/n_1s}. \cr 
      \code{$singletons[3]=}\eqn{S_o2} \tab
      : The number of singletons found in the sample from population
      2. \cr 
      \tab i.e.: \eqn{S_o2} \eqn{\in} \eqn{[0,S_2]}, where the
      frequency of the derived allele is \eqn{r_s=1/n_2s}. 
    }
  } 
  \item{$frequencies}{
    The vector of nine summary statistics on the mean derived allele
    frequencies for the seven types of segregating sites (as specified in
    \code{$segsites[}\eqn{1: 9}\code{]}): 
    \tabular{ll}{
      \code{$frequencies[1]=}\eqn{F(S)} \tab 
      : The mean derived allele frequency over \eqn{n} haplotypes for
      the \eqn{S} segregating sites in the data set \eqn{h}. \cr 
      \code{$frequencies[2]=}\eqn{F(S_1)} \tab 
      : The mean derived allele frequency over \eqn{n_1} haplotypes for
      the \eqn{S_1} segregating sites with the derived alleles polymorphic
      uniquely in the sample from population 1. \cr 
      \code{$frequencies[3]=}\eqn{F(S_2)} \tab 
      : The mean derived allele frequency over \eqn{n_2} haplotypes for
      the \eqn{S_2} segregating sites with the derived alleles polymorphic
      uniquely in the sample from population 2. \cr 
      \code{$frequencies[4]=}\eqn{F(S_s)} \tab
      : The mean derived allele frequency over \eqn{n} haplotypes for
      the \eqn{S_s} segregating sites with the derived allele shared by both
      population samples. \cr 
      \code{$frequencies[5]=}\eqn{F(S_ss)} \tab
      : The mean derived allele frequency over \eqn{n} haplotypes for
      the \eqn{S_ss} segregating sites with the derived allele shared AND
      polymorphic in both population samples. \cr 
      \code{$frequencies[6]=}\eqn{F(S_sf1)} \tab 
      : The mean derived allele frequency over \eqn{n} haplotypes for
      the \eqn{S_sf1} segregating sites with the derived allele shared AND
      polymorphic ONLY in the sample from population 2 (i.e., fixed in the
      sample from population 1). \cr 
      \code{$frequencies[7]=}\eqn{F(S_sf1)} \tab 
      : The mean derived allele frequency over \eqn{n} haplotypes for
      the \eqn{S_sf2} segregating sites with the derived allele shared AND
      polymorphic ONLY in the sample from population 1 (i.e., fixed in the
      sample from population 2).\cr
       \code{$frequencies[8]=}\eqn{F(S_sl)} \tab 
      : The mean derived allele frequency over \eqn{n} haplotypes for
      the \eqn{S_sl} segregating sites with the derived allele shared
      AND at low frequency (i.e., \eqn{p_s \le 10\%}). \cr 
      \code{$frequencies[9]=}\eqn{F(S_sh)} \tab 
      : The mean derived allele frequency over \eqn{n} haplotypes for
      the \eqn{S_sh} segregating sites with the derived allele shared
      AND at high frequency (i.e., \eqn{p_s > 10\%}).
    }
  }
  \item{$fst}{
    The vector of five summary statistics on the levels of
    differentiation and nucleotide diversity between the populations: 
    \tabular{ll}{
      \code{$fst[1]=}\eqn{F_st} \tab 
      : The level of differentiation between the two population
      samples. \cr 
      \tab \eqn{F_st=1-H_w/H_b}, where \eqn{H_w=H_w1+H_w2} (Wright 1931;
      Hudson et al. 1992). \cr 
      \code{$fst[2]=}\eqn{H_b} \tab 
      : The mean pairwise differences between both population
      samples. \cr 
      \code{$fst[3]=}\eqn{H_w1} \tab 
      : The mean pairwise differences within the sample from population
      1. \cr 
      \code{$fst[4]=}\eqn{H_w2} \tab 
      : The mean pairwise differences within the sample from population
      2. \cr 
      \code{$fst[5]=}\eqn{S_nn} \tab 
      : The nearest-neighbor statistic, which measures of how often the
      nearest neighbors of haplotypes are found in the same population
      sample. \cr 
      \tab \eqn{S_nn=\Sigma[X_i/n]=\Sigma[W_i/(T_i*n)]}, \eqn{i}
      \eqn{\in} \eqn{[1,n]} (Hudson 2000). \cr 
      \tab \eqn{W_i} is the number of nearest neighbors to haplotype
      \eqn{i} that are from the same population sample as haplotype
      \eqn{i}. \cr 
      \tab \eqn{T_i} is the number of nearest neighbors of haplotype
      \eqn{i}. \cr 
      \tab \eqn{X_i} is the fraction of nearest neighbors of haplotype
      \eqn{i} that are from the same population sample as haplotype \eqn{i}. 
    }
  }
  \item{$pi}{
    The vector of three summary statistics on the unbiased estimator of
    the population mutation rate calculated from the average nucleotide
    diversity (or heterozygosity), \eqn{\pi} (Nei and Li 1979, but
    equations 10 and 12 of Tajima 1989). 
    \tabular{ll}{
      \code{$pi[1]=}\eqn{\pi} \tab
      : \eqn{\pi} for the data set \eqn{h}. \cr
      \tab \eqn{\pi=\Sigma[2*p_s*(1-p_s)*n_s/(n_s-1)]}, \eqn{s}
      \eqn{\in} \eqn{[1,S]}. \cr 
      \code{$pi[2]=}\eqn{\pi_1} \tab 
      : \eqn{\pi} for the sample from population 1. \cr
      \tab \eqn{\pi_1=\Sigma[2*q_s*(1-q_s)*n_1s/(n_1s-1)]}, \eqn{s}
      \eqn{\in} \eqn{[1,S_1]}. \cr 
      \code{$pi[3]=}\eqn{\pi_2} \tab 
      : \eqn{\pi} for the sample from population 2. \cr
      \tab \eqn{\pi_2=\Sigma[2*r_s*(1-r_s)*n_2s/(n_2s-1)]}, \eqn{s}
      \eqn{\in} \eqn{[1,S_2]}. 
    }
  } 
  \item{$thetaW}{
    The vector of three summary statistics on the unbiased estimator of
    the population mutation rate calculated from the number of segregating
    sites, \eqn{\theta_W} (Watterson 1975). \cr 
    \tabular{ll}{
      \code{$thetaW[1]=}\eqn{\theta_W} \tab 
      : \eqn{\theta_W} for the data set \eqn{h}. \cr
      \tab \eqn{\theta_W=\Sigma[1/(1/n_s)]}, \eqn{s} \eqn{\in}
      \eqn{[1,S]}. \cr 
      \code{$thetaW[2]=}\eqn{\theta_W1} \tab 
      : \eqn{\theta_W} for the sample from population 1. \cr
      \tab \eqn{\theta_W1=\Sigma[1/(1/n_1s)]}, \eqn{s} \eqn{\in}
      \eqn{[1,S_1]}. \cr 
      \code{$thetaW[3]=}\eqn{\theta_W2} \tab 
      : \eqn{\theta_W} for the sample from population 2. \cr
      \tab \eqn{\theta_W2=\Sigma[1/(1/n_2s)]}, \eqn{s} \eqn{\in}
      \eqn{[1,S_2]}. 
    }
  }
  \item{$thetaH}{
    The vector of three summary statistics on the unbiased estimator of
    the population mutation rate weighted by the homozygosity of the
    derived allele, \eqn{\theta_H} (Fay and Wu 2000). 
    \tabular{ll}{
      \code{$thetaH[1]=}\eqn{\theta_H} \tab 
      : \eqn{\theta_H} for the data set \eqn{h}. \cr
      \tab \eqn{\theta_H=\Sigma[2*(p_s^2)/[n_s(n_s-1)]]}, \eqn{s}
      \eqn{\in} \eqn{[1,S]}. \cr 
      \code{$thetaH[2]=}\eqn{\theta_H1} \tab 
      : \eqn{\theta_H} for the sample from population 1. \cr
      \tab \eqn{\theta_H1=\Sigma[2*(q_s^2)/[n_1s(n_1s-1)]]}, \eqn{s}
      \eqn{\in} \eqn{[1,S_1]}. \cr 
      \code{$thetaH[3]=}\eqn{\theta_H2} \tab 
      : \eqn{\theta_H} for the sample from population 2. \cr
      \tab \eqn{\theta_H2=\Sigma[2*(r_s^2)/[n_2s(n_2s-1)]]}, \eqn{s}
      \eqn{\in} \eqn{[1,S_2]}. 
    }
  }
  \item{$D}{
    The vector of three summary statistics on the test of neutrality
    from the difference between \eqn{\pi} and \eqn{\theta_W}, Tajima's
    \eqn{D} (Tajima 1989). 
    \tabular{ll}{
      \code{$D[1]=}\eqn{D} \tab 
      : Tajima's \eqn{D} for the data set \eqn{h}. \cr
      \tab \eqn{D=(\pi-\theta_W)/\sqrt(var(\pi-\theta_W))}. \cr 
      \code{$D[2]=}\eqn{D_1} \tab 
      : Tajima's \eqn{D} for the sample from population 1. \cr
      \tab \eqn{D_1=(\pi_1-\theta_W1)/\sqrt(var(\pi_1-\theta_W1))}. \cr 
      \code{$D[3]=}\eqn{D_2} \tab 
      : Tajima's \eqn{D} for the sample from population 2. \cr 
      \tab \eqn{D_2=(\pi_2-\theta_W2)/\sqrt(var(\pi_2-\theta_W2))}. 
    }
  }
  \item{$H}{
    The vector of three summary statistics on the test of neutrality
    from the difference between \eqn{\pi} and \eqn{\theta_H}, Fay &
    Wu's \eqn{H} (Fay and Wu 2000). 
    \tabular{ll}{
      \code{$H[1]=}\eqn{H} \tab 
      : Fay & Wu's \eqn{H} for the data set \eqn{h}. \cr
      \tab \eqn{H=\pi-\theta_H=-2*p_s(2*p_s-1)n_s/(n-1)}, \eqn{s}
      \eqn{\in} \eqn{[1,S]}. \cr 
      \code{$H[2]=}\eqn{H_1} \tab 
      : Fay & Wu's \eqn{H} for the sample from population 1. \cr
      \tab \eqn{H_1=\pi_1-\theta_H1=-2*q_s(2*q_s-1)n_1s/(n_1s-1)},
      \eqn{s} \eqn{\in} \eqn{[1,S_1]}. \cr 
      \code{$H[3]=}\eqn{H_2} \tab 
      : Fay & Wu's \eqn{H} for the sample from population 2. \cr 
      \tab \eqn{H_2=\pi_2-\theta_H2=-2*r_s(2*r_s-1)n_2s/(n_2s-1)},
      \eqn{s} \eqn{\in} \eqn{[1,S_2]}. 
    }
  }
  \item{$D_star}{
    The vector of three summary statistics on the test of neutrality
    from the difference between the number of singletons and the number
    of segregating sites, Fu & Li's \eqn{D*} (Fu & Li 1993). 
    \tabular{ll}{
      \code{$D_star[1]=}\eqn{D*} \tab 
      : Fu & Li's \eqn{D*} for the data set \eqn{h}. \cr
      \tab \eqn{D*=(S-S_o)/\sqrt(var(S-S_o))}. \cr
      \code{$D_star[2]=}\eqn{D*_1} \tab 
      : Fu & Li's \eqn{D*} for the sample from population 1. \cr
      \tab \eqn{D_1*=(S_1-S_o1)/\sqrt(var(S_1-S_o1))}. \cr
      \code{$D_star[3]=}\eqn{D*_2} \tab 
      : Fu & Li's \eqn{D*} for the sample from population 2. \cr
      \tab \eqn{D_2*=(S_2-S_o2)/\sqrt(var(S_2-S_o2))}.
    }
  } Where: \tabular{ll}{
    \eqn{n_s} \tab
    : The number of chromosomes without missing data in the data set
    \eqn{h} (as specified with \code{"?"} in
    \code{msframe$haplotypes[}\eqn{h,}\code{]}) at the segregating site
    \eqn{s}, \eqn{n_s} \eqn{\in} \eqn{[1,n]} and \eqn{s} \eqn{\in}
    \eqn{[1,S]}. \cr 
    \eqn{n_1s} \tab
    : The number of chromosomes without missing data in the sample from
    population 1 (as specified with \code{"?"} in
    \code{msframe$haplotypes[}\eqn{h,n_1s}\code{]}) at the segregating
    site \eqn{s}, \eqn{n_1s} \eqn{\in} \eqn{[1,n_1]} and \eqn{s}
    \eqn{\in} \eqn{[1,S_1]}. \cr 
    \eqn{n_2s} \tab
    : The number of chromosomes without missing data in the sample from
    population 2 (as specified with \code{"?"} in
    \code{msframe$haplotypes[}\eqn{h,n_2s}\code{]}) at the segregating
    site \eqn{s}, \eqn{n_2s} \eqn{\in} \eqn{[1,n_2]} and \eqn{s}
    \eqn{\in} \eqn{[1,S_2]}. \cr 
    \eqn{p_s} \tab
    : The frequency of the derived allele at the segregating site
    \eqn{s} in the data set \eqn{h}, \eqn{s} \eqn{\in} \eqn{[1,S]}. \cr 
    \eqn{q_s} \tab 
    : The frequency of the derived allele at the segregating site
    \eqn{s} in the sample from population 1, \eqn{s} \eqn{\in}
    \eqn{[1,S_1]}. \cr 
    \eqn{r_s} \tab
    : The frequency of the derived allele at the segregating site
    \eqn{s} in the sample from population 2, \eqn{s} \eqn{\in}
    \eqn{[1,S_2]}. 
  }
}

\note{
  \item{ATTENTION: It is the user's responsibility to mind the following
    restrictions:}{
    -> The data frame \code{msframe} needs to be of the same format as
    output by the \R functions \code{\link{ms}} and
    \code{\link{get_ms_output}} for a model with TWO populations (although
    the vector \code{msframe$positions} is ignored). \cr 
    -> The \code{ms} command line, \code{cmdline_seeds[1]}, needs to
    contain at minimum these SEVEN "words" (The first word
    can be any string of characters without white space): \cr 
    \code{"./ms }\eqn{n}\code{ }\eqn{H}\code{ -I 2 }\eqn{n_1}\code{ }\eqn{n_2} \code{"}. 
  }
}

\author{
  Celine Becquet - \email{celine.becquet@gmail.com}. 
}

\references{
  Becquet, C., and M. Przeworski. 2007. A new approach to estimate
  parameters of speciation models with application to apes. Genome
  Res. 17:1505-1519. \cr 
  Fay, J. C., and C.-I. Wu. 2000. Hitchhiking Under Positive Darwinian
  Selection. Genetics 155:1405-1413. \cr 
  Fu, Y.-X. and W.-H. Li. 1993. Statistical tests of neutrality of
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

\seealso{
  The function \code{\link{stats_pop}} uses the data frame output by
  \code{\link{ms}} or \code{\link{get_ms_output}} to calculate summary
  statistics. \cr
  \code{\link{stats_pop}} is called by the functions
  \code{\link{simulate_data}} and \code{\link{estimate_IMc}}.\cr
  The function \code{\link{get_stats}} uses the data frame output by
  \code{\link{stats_pop}} as input.
  Lists of definitions of the symbols and parameters mentioned in
  this file are found in \code{\link{Rmspack-package}}.
}

\examples{
### Create the input for the function. 
msline="./ms 8 2 -t 1 -I 2 4 4 0 -ej 0.5 1 2 -seeds 1 2 3"
msframe=ms(msline)

## Calculate summary statistics with known ancestral allele. 
stats_pop(msframe=msframe)

## Calculate summary statistics with unknown ancestral allele. 
stats_pop(msframe=msframe, type="maf")

## Add missing data. 
msframe$haplotypes[1,2]=c("??????")
msframe$haplotypes[2,]=c("0?001100", "0?001100", "?0110001", "?0110000", "10110000", "101100?0", "10110000", "10110000")
stats_pop(msframe=msframe)
}

\keyword{programming}
