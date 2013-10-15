\name{get_stats}
\alias{get_stats}

\title{
  Function to extrat the considered summary statistics at a genomic
  region or locus.
}

\description{
  The function  \code{\link{get_stats}} extracts the considered summary
  statistics (i.e., corresponding to the positions specified in the argument
  \code{summaries}) of the genomic region or locus \code{nlocus} from
  the data frame \code{statslocus} output by the function
  \code{\link{stats_pop}}. \cr
  It adds the considered summary statistics into the
   row number \code{nlocus} of the matrices with \code{nloci} rows of the
  data frame \code{sumstats}.
  This function is called within the \R functions
  \code{\link{simulate_data}} and \code{\link{estimate_IMc}}. 
}
\usage{
get_stats(nloci, nlocus, sumstats, statslocus, info, summaries)
}

\arguments{
  \item{nloci}{
    The  number of independent genomic region considered \eqn{R}
    (specified in  specified with the keyword \emph{\code{nloci}} in the file
    \code{paramfile},see \code{\link{param_sim}} and
    \code{\link{param_est}}).\cr
    Alternatively, the number of genomic loci considered, \eqn{\Sigma(Y_r)},
    \eqn{r} \eqn{\in} \eqn{[1,R]} (specified in the file \code{regfile}, see
    \code{\link{info_region}}).\cr
    \code{nloci} determines the number of rows of the matrices in the
    data frame \code{sumstats}.
  }
  \item{nlocus}{
    The  number of the genomic region
    considered, \code{nlocus=}\eqn{r} \eqn{\in} \eqn{[1,R]}.\cr
    Alternatively, the  number of the locus  considered,
    \code{nlocus=}\eqn{l} \eqn{\in} \eqn{[1,\Sigma(Y_r)]}. \cr
  }
  \item{sumstats}{
    The data frame of matrices with \code{nloci} rows output by
    \code{\link{get_stats}}. See section  \code{value} below.
  }
  \item{statslocus}{
    The data frame output by the function
    \code{\link{stats_pop}} for the genomic region  or 
    locus \code{nlocus}.
  }
  \item{info}{
    The vector of information on the  genomic region or 
    locus  \code{nlocus}.\cr
    -- If the genomic region \code{nlocus=}\eqn{r} is considered, the vector is of size
    \eqn{10} and \code{info$V1=}\eqn{r}.\cr
    The other nine information are irrelevant for the function
    \code{\link{get_stats}}), but see \code{\link{info_region}}.\cr
    -- If the genomic locus \code{nlocus=}\eqn{l} is considered, the vector is of size
    \eqn{<10} (usually \eqn{6}) and \code{info$V1=}\eqn{r} and
    \code{info$V2=}\eqn{y}.\cr
    The other four information are irrelevant for the function
    \code{\link{get_stats}}), but see \code{\link{info_loci}}.\cr     
    Where
    \tabular{ll}{
      \eqn{r}\tab 
      : The uni- or multi-locus genomic region number, \eqn{r} \eqn{\in}
      \eqn{[1,R]}. \cr 
      \eqn{R}\tab
      : The number of independent genomic regions considered
      as specified with the keyword \emph{\code{nloci}} in the file
      \code{paramfile} (see \code{\link{param_sim}} and
      \code{\link{param_est}}). \cr  
      \eqn{y}\tab 
      : The locus number, \eqn{y} \eqn{\in} \eqn{[1,Y_r]}. \cr 
      \eqn{Y_r}\tab
      : The total number of loci spanning the multi-locus
      genomic region \eqn{r} as specified in the file \code{regfile} (see
      \code{\link{info_region}}). 
    }
  }
  \item{summaries}{
    The data frame \code{$stats} of the positions of the considered summary statistics
    in the data frame \code{statslocus}  output by the function
    \code{\link{get_liststats}}. 
  }
}

\details{
  If \code{nlocus=1}, the data frame  \code{sumstats} is initialized.\cr
  i.e., the equivalent matrices for which
  \code{length(summaries$name)>0} are created with \code{nloci} rows and
  \code{length(summaries$name)} columns
}

\value{
  The function \code{\link{get_stats}} outputs a data frame similar to
  ones output by the functions \code{\link{stats_pop}} and
  \code{\link{get_liststats}}.\cr
  The data frame contains matrices with \code{nloci} rows and the
  same numbers of columns as in the equivalent vectors in the data frame
  \code{summaries}.\cr
  \strong{Extra matrix:}
  \item{$loci}{
    The matrix of \eqn{R} or \eqn{\Sigma{Y_r}}  \eqn{\times} \eqn{2}:
    \tabular{ll}{
      \code{$loci[,1]}\tab
      : The vector of \code{nloci} number of independent genomic regions
      or loci considered.\cr
      \tab -- If \code{nloci=}\eqn{R} independent genomic regions \eqn{r}
      are considered:\cr
      \tab \code{$loci[}\eqn{r}\code{,1]=}\eqn{r}, the genomic region number.\cr
      \tab -- If \code{nloci=}\eqn{\Sigma(Y_r)} loci are considered:\cr
      \tab \code{$loci[}\eqn{l}\code{,1]=}\eqn{r}, the genomic region number
      that locus \eqn{l} is part of.\cr
      \code{$loci[,2]}\tab
      :  The vector of \code{nloci} number of independent genomic regions
       considered.\cr
      \tab -- If \code{nloci=}\eqn{R} independent genomic regions \eqn{r}
      are considered:\cr
      \tab \code{$loci[}\eqn{r}\code{,2]=}\eqn{r}, the genomic region number.\cr
      \tab -- If \code{nloci=}\eqn{\Sigma(Y_r)} loci are considered:\cr
      \tab \code{$loci[}\eqn{l}\code{,2]=}\eqn{y}, the locus number
      \eqn{\in} \eqn{[1,Y_r]} corresponding to the locus \eqn{l}.
    }
  }
  \strong{Matrices from the data frames output by
    \code{\link{stats_pop}} and  \code{$stats} output by
    \code{\link{get_liststats}}:}
   \item{$type}{
     The sentence describing the type of summary statistics calculated:
     \cr 
     \tabular{ll}{
       \code{"Known ancestral allele"}\tab
       : Specifies that the summary statistic in \code{statslocus} were calculated 
       for the data sets assuming a known ancestral state at each
       segregating site. \cr 
       \code{"Unknown ancestral allele"} \tab     
       : Specifies that  the allele with
       the minimum frequency at each segregating site was considered  the derived allele. 
     }
   } 
  \item{$segsites}{
    The matrix \eqn{R} or \eqn{\Sigma{Y_r}}
     \eqn{\times} \code{length(summaries$stats$segsites)} 
    of summary statistics on the number of segregating sites in the full sample and
    the populations samples (Wakeley and Hey 1997;  Becquet and
    Przeworski 2007). \cr
    If \code{summaries$stats$segsites[}\eqn{i}\code{]=}\eqn{j},
    \eqn{i} and \eqn{j} \eqn{\in} \eqn{[1,12]},\cr
    \code{$segsites[}\eqn{l}\code{,}\eqn{i}\code{]=statslocus$segsites[}\eqn{j}\code{]},\cr
    where: 
    \tabular{lll}{
      \eqn{j=1}:\tab \code{"S"}\tab 
      : The total number of segregating sites found in a data set. \cr
      \tab \tab \eqn{S=S_1+S_2+S_s+S_f}. \cr
      \eqn{j=2}:\tab \code{"S_1"}\tab 
      : The number of segregating sites with the derived alleles
      polymorphic uniquely in the sample from population 1. \cr 
      \eqn{j=3}:\tab \code{"S_2"}\tab 
      : The number of segregating sites with the derived alleles
      polymorphic uniquely in the sample from population 2. \cr 
      \eqn{j=4}:\tab \code{"S_s"}\tab 
      : The number of segregating sites with the derived allele shared
      by both population samples. \cr 
      \tab \tab \eqn{S_s=S_ss+S_sfs+S_sf2} and \eqn{S_s=S_sl+S_sh}. \cr
      \eqn{j=5}:\tab \code{"S_ss"}\tab 
      : The number of segregating sites with the derived allele shared
      AND polymorphic in both population samples. \cr 
      \eqn{j=6}:\tab \code{"S_sf1"}\tab 
      : The number of segregating sites with the derived allele shared
      AND polymorphic ONLY in the sample from population 2 (i.e., fixed in
      the sample from population 1). \cr 
      \eqn{j=7}:\tab \code{"S_sf2"}\tab
      : The number of segregating sites with the derived allele shared
      AND polymorphic ONLY in the sample from population 1 (i.e., fixed in
      the sample from population 2). \cr
       \eqn{j=8} \tab \code{"S_sl"} \tab 
        : The number of segregating sites with the derived allele shared
        AND  at low frequency (i.e., \eqn{p_s \le 10\%}). \cr 
        \eqn{j=9} \tab \code{"S_sh"} \tab
        : The number of segregating sites with the derived allele shared
        AND  at high frequency (i.e., \eqn{p_s > 10\%}. \cr 
        \eqn{j=10}:\tab \code{"S_f"}\tab 
        : The number of segregating sites with the alleles fixed between
        the population sampled. \cr 
        \tab \tab \eqn{S_f=S_f1+S_f2}. \cr
        \eqn{j=11}:\tab \code{"S_f1"}\tab 
        : The number of segregating sites with the derived allele fixed in
        the sample from population 1. \cr 
        \eqn{j=12}:\tab \code{"S_f2"}\tab 
        : The number of segregating sites with the derived allele fixed in
        the sample from population 2. 
    }
  }
  \item{$singletons}{
    The matrix \eqn{R} or \eqn{\Sigma{Y_r}}
     \eqn{\times} \code{length(summaries$stats$singletons)} 
    of summary statistics on the number of singletons in the full sample
    and the populations samples. \cr
    If \code{summaries$stats$singletons[}\eqn{i}\code{]=}\eqn{j},
    \eqn{i} and \eqn{j} \eqn{\in} \eqn{[1,3]},\cr
    \code{$singletons[}\eqn{l}\code{,}\eqn{i}\code{]=statslocus$singletons[}\eqn{j}\code{]},\cr
    where: 
    \tabular{lll}{
      \eqn{j=1}:\tab \code{"S_o"}\tab
      : The number of singletons found in a data set. \cr 
      \tab \tab i.e.: \eqn{S_o=S_o1+S_o2} \eqn{\in} \eqn{[0,S]}, where
      the frequency of the derived allele is \eqn{p_s=1/n_s}. \cr 
      \eqn{j=2}:\tab \code{"S_o1"}\tab
      : The number of singletons found in the sample from population
      1. \cr 
      \tab \tab i.e.: \eqn{S_o1} \eqn{\in} \eqn{[0,S_1]}, where the
      frequency of the derived allele is \eqn{q_s=1/n_1s}. \cr 
      \eqn{j=3}:\tab \code{"S_o2"}\tab
      : The number of singletons found in the sample from population
      2. \cr 
      \tab \tab i.e.: \eqn{S_o2} \eqn{\in} \eqn{[0,S_2]}, where the
      frequency of the derived allele is \eqn{r_s=1/n_2s}. 
    }
  } 
  \item{$frequencies}{    
    The matrix \eqn{R} or \eqn{\Sigma{Y_r}}
     \eqn{\times} \code{length(summaries$stats$frequencies)} 
    of summary statistics on the mean derived allele
    frequencies for the seven types of segregating sites (as specified in
    \code{$segsites} with positions \eqn{j} \eqn{\in} \eqn{[1: 9]}). \cr
    If \code{summaries$stats$frequencies[}\eqn{i}\code{]=}\eqn{j},
    \eqn{i} and \eqn{j} \eqn{\in} \eqn{[1,9]},\cr
    \code{$frequencies[}\eqn{l}\code{,}\eqn{i}\code{]=statslocus$frequencies[}\eqn{j}\code{]},\cr
    where: 
    \tabular{lll}{
      \eqn{j=1}:\tab \code{"F(S)"}\tab 
      : The mean derived allele frequency over \eqn{n} haplotypes for
      the \eqn{S} segregating sites in a data set. \cr 
      \eqn{j=2}:\tab \code{"F(S_1)"}\tab 
      : The mean derived allele frequency over \eqn{n_1} haplotypes for
      the \eqn{S_1} segregating sites with the derived alleles polymorphic
      uniquely in the sample from population 1. \cr 
      \eqn{j=3}:\tab \code{"F(S_2)"}\tab 
      : The mean derived allele frequency over \eqn{n_2} haplotypes for
      the \eqn{S_2} segregating sites with the derived alleles polymorphic
      uniquely in the sample from population 2. \cr 
      \eqn{j=4}:\tab \code{"F(S_s)"}\tab
      : The mean derived allele frequency over \eqn{n} haplotypes for
      the \eqn{S_s} segregating sites with the derived allele shared by both
      population samples. \cr 
      \eqn{j=5}:\tab \code{"F(S_ss)"}\tab
      : The mean derived allele frequency over \eqn{n} haplotypes for
      the \eqn{S_ss} segregating sites with the derived allele shared AND
      polymorphic in both population samples. \cr 
      \eqn{j=6}:\tab \code{"F(S_sf1)"}\tab 
      : The mean derived allele frequency over \eqn{n} haplotypes for
      the \eqn{S_sf1} segregating sites with the derived allele shared AND
      polymorphic ONLY in the sample from population 2 (i.e., fixed in the
      sample from population 1). \cr 
      \eqn{j=7}:\tab \code{"F(S_sf1)"}\tab 
      : The mean derived allele frequency over \eqn{n} haplotypes for
      the \eqn{S_sf2} segregating sites with the derived allele shared AND
      polymorphic ONLY in the sample from population 1 (i.e., fixed in the
      sample from population 2). \cr
      \eqn{j=8}: \tab \code{"F(S_sl)"} \tab 
      : The mean derived allele frequency over \eqn{n} haplotypes for
      the \eqn{S_sl} segregating sites with the derived allele shared
      AND at low frequency (i.e., \eqn{p_s \le 10\%}). \cr 
      \eqn{j=9}: \tab \code{"F(S_sh)"} \tab 
      : The mean derived allele frequency over \eqn{n} haplotypes for
      the \eqn{S_sh} segregating sites with the derived allele shared
      AND at high frequency (i.e., \eqn{p_s > 10\%}).
    }
  }
  \item{$fst}{
    The matrix \eqn{R} or \eqn{\Sigma{Y_r}}
     \eqn{\times} \code{length(summaries$stats$fst)} 
    of summary statistics on the levels of
    differentiation and nucleotide diversity between the populations. \cr
    If \code{summaries$stats$fst[}\eqn{i}\code{]=}\eqn{j},
    \eqn{i} and \eqn{j} \eqn{\in} \eqn{[1,5]},\cr
    \code{$fst[}\eqn{l}\code{,}\eqn{i}\code{]=statslocus$fst[}\eqn{j}\code{]},\cr
    where: 
    \tabular{lll}{
      \eqn{j=1}:\tab \code{"F_st"}\tab 
      : The level of differentiation between the two population
      samples. \cr 
      \tab \tab \eqn{F_st=1-H_w/H_b}, where \eqn{H_w=H_w1+H_w2} (Wright 1931;
      Hudson et al. 1992). \cr 
      \eqn{j=2}:\tab \code{"H_b"}\tab 
      : The mean pairwise differences between both population
      samples. \cr 
      \eqn{j=3}:\tab \code{"H_w1"}\tab 
      : The mean pairwise differences within the sample from population
      1. \cr 
      \eqn{j=4}:\tab \code{"H_w2"}\tab 
      : The mean pairwise differences within the sample from population
      2. \cr 
      \eqn{j=5}:\tab \code{"S_nn"}\tab 
      : The nearest-neighbor statistic, which measures of how often the
      nearest neighbors of haplotypes are found in the same population
      sample. \cr 
      \tab \tab \eqn{S_nn=\Sigma[X_i/n]=\Sigma[W_i/(T_i*n)]}, \eqn{i}
      \eqn{\in} \eqn{[1,n]} (Hudson 2000). \cr 
      \tab \tab \eqn{W_i} is the number of nearest neighbors to haplotype
      \eqn{i} that are from the same population sample as haplotype
      \eqn{i}. \cr 
      \tab \tab \eqn{T_i} is the number of nearest neighbors of haplotype
      \eqn{i}. \cr 
      \tab \tab \eqn{X_i} is the fraction of nearest neighbors of haplotype
      \eqn{i} that are from the same population sample as haplotype \eqn{i}. 
    }
  }
  \item{$pi}{
    The matrix \eqn{R} or \eqn{\Sigma{Y_r}}
     \eqn{\times} \code{length(summaries$stats$pi)} 
    of summary statistics on the unbiased estimator of
    the population mutation rate calculated from the average nucleotide
    diversity (or heterozygosity), \eqn{\pi} (Nei and Li 1979, but
    equations 10 and 12 of Tajima 1989). \cr
    If \code{summaries$stats$pi[}\eqn{i}\code{]=}\eqn{j},
    \eqn{i} and \eqn{j} \eqn{\in} \eqn{[1,3]},\cr
    \code{$pi[}\eqn{l}\code{,}\eqn{i}\code{]=statslocus$pi[}\eqn{j}\code{]},\cr
    where: 
    \tabular{lll}{
      \eqn{j=1}:\tab \code{"pi"}\tab
      : \eqn{\pi} for a data set. \cr
      \tab \tab \eqn{\pi=\Sigma[2*p_s*(1-p_s)*n_s/(n_s-1)]}, \eqn{s}
      \eqn{\in} \eqn{[1,S]}. \cr 
      \eqn{j=2}:\tab \code{"pi_1"}\tab 
      : \eqn{\pi} for the sample from population 1. \cr
      \tab \tab \eqn{\pi_1=\Sigma[2*q_s*(1-q_s)*n_1s/(n_1s-1)]}, \eqn{s}
      \eqn{\in} \eqn{[1,S_1]}. \cr 
      \eqn{j=3}:\tab \code{"pi_2"}\tab 
      : \eqn{\pi} for the sample from population 2. \cr
      \tab \tab \eqn{\pi_2=\Sigma[2*r_s*(1-r_s)*n_2s/(n_2s-1)]}, \eqn{s}
      \eqn{\in} \eqn{[1,S_2]}. 
    }
  } 
  \item{$thetaW}{
    The matrix \eqn{R} or \eqn{\Sigma{Y_r}}
     \eqn{\times} \code{length(summaries$stats$theta_W)} 
    of summary statistics on the unbiased estimator of
    the population mutation rate calculated from the number of segregating
    sites, \eqn{\theta_W} (Watterson 1975). \cr
    If \code{summaries$stats$theta_W[}\eqn{i}\code{]=}\eqn{j},
    \eqn{i} and \eqn{j} \eqn{\in} \eqn{[1,3]},\cr
    \code{$theta_W[}\eqn{l}\code{,}\eqn{i}\code{]=statslocus$theta_W[}\eqn{j}\code{]},\cr
    where: 
    \tabular{lll}{
      \eqn{j=1}:\tab \code{"theta_W"}\tab 
      : \eqn{\theta_W} for a data set. \cr
      \tab \tab \eqn{\theta_W=\Sigma[1/(1/n_s)]}, \eqn{s} \eqn{\in}
      \eqn{[1,S]}. \cr 
      \eqn{j=2}:\tab \code{"theta_W1"}\tab 
      : \eqn{\theta_W} for the sample from population 1. \cr
      \tab \tab \eqn{\theta_W1=\Sigma[1/(1/n_1s)]}, \eqn{s} \eqn{\in}
      \eqn{[1,S_1]}. \cr 
      \eqn{j=3}:\tab \code{"theta_W2"}\tab 
      : \eqn{\theta_W} for the sample from population 2. \cr
      \tab \tab \eqn{\theta_W2=\Sigma[1/(1/n_2s)]}, \eqn{s} \eqn{\in}
      \eqn{[1,S_2]}. 
    }
  }
  \item{$thetaH}{
    The matrix \eqn{R} or \eqn{\Sigma{Y_r}}
     \eqn{\times} \code{length(summaries$stats$theta_H)} 
    of summary statistics on the unbiased estimator of
    the population mutation rate weighted by the homozygosity of the
    derived allele, \eqn{\theta_H} (Fay and Wu 2000). \cr
    If \code{summaries$stats$theta_H[}\eqn{i}\code{]=}\eqn{j},
    \eqn{i} and \eqn{j} \eqn{\in} \eqn{[1,3]},\cr
    \code{$theta_H[}\eqn{l}\code{,}\eqn{i}\code{]=statslocus$theta_H[}\eqn{j}\code{]},\cr
    where: 
    \tabular{lll}{
      \eqn{j=1}:\tab \code{"theta_H"}\tab 
      : \eqn{\theta_H} for a data set. \cr
      \tab \tab \eqn{\theta_H=\Sigma[2*(p_s^2)/[n_s(n_s-1)]]}, \eqn{s}
      \eqn{\in} \eqn{[1,S]}. \cr 
      \eqn{j=2}:\tab \code{"theta_H1"}\tab 
      : \eqn{\theta_H} for the sample from population 1. \cr
      \tab \tab \eqn{\theta_H1=\Sigma[2*(q_s^2)/[n_1s(n_1s-1)]]}, \eqn{s}
      \eqn{\in} \eqn{[1,S_1]}. \cr 
      \eqn{j=3}:\tab \code{"theta_H2"}\tab 
      : \eqn{\theta_H} for the sample from population 2. \cr
      \tab \tab \eqn{\theta_H2=\Sigma[2*(r_s^2)/[n_2s(n_2s-1)]]}, \eqn{s}
      \eqn{\in} \eqn{[1,S_2]}. 
    }
  }
  \item{$D}{
    The matrix \eqn{R} or \eqn{\Sigma{Y_r}}
     \eqn{\times} \code{length(summaries$stats$D)} 
    of summary statistics on the test of neutrality
    from the difference between \eqn{\pi} and \eqn{\theta_W}, Tajima's
    \eqn{D} (Tajima 1989. \cr
    If \code{summaries$stats$D[}\eqn{i}\code{]=}\eqn{j},
    \eqn{i} and \eqn{j} \eqn{\in} \eqn{[1,3]},\cr
    \code{$D[}\eqn{l}\code{,}\eqn{i}\code{]=statslocus$D[}\eqn{j}\code{]},\cr
    where:
    \tabular{lll}{
      \eqn{j=1}:\tab \code{"D"}\tab 
      : Tajima's \eqn{D} for a data set. \cr
      \tab \tab \eqn{D=(\pi-\theta_W)/\sqrt(var(\pi-\theta_W))}. \cr 
      \eqn{j=2}:\tab \code{"D_1"}\tab 
      : Tajima's \eqn{D} for the sample from population 1. \cr
      \tab \tab \eqn{D_1=(\pi_1-\theta_W1)/\sqrt(var(\pi_1-\theta_W1))}. \cr 
      \eqn{j=3}:\tab \code{"D_2"}\tab 
      : Tajima's \eqn{D} for the sample from population 2. \cr 
      \tab \tab \eqn{D_2=(\pi_2-\theta_W2)/\sqrt(var(\pi_2-\theta_W2))}. 
    }
  }
  \item{$H}{
    The matrix \eqn{R} or \eqn{\Sigma{Y_r}}
     \eqn{\times} \code{length(summaries$stats$H)} 
    of summary statistics on the test of neutrality
    from the difference between \eqn{\pi} and \eqn{\theta_H}, Fay &
    Wu's \eqn{H} (Fay and Wu 2000). \cr
    If \code{summaries$stats$H[}\eqn{i}\code{]=}\eqn{j},
    \eqn{i} and \eqn{j} \eqn{\in} \eqn{[1,3]},\cr
    \code{$H[}\eqn{l}\code{,}\eqn{i}\code{]=statslocus$H[}\eqn{j}\code{]},\cr
    where:
    \tabular{lll}{
      \eqn{j=1}:\tab \code{"H"}\tab 
      : Fay & Wu's \eqn{H} for a data set. \cr
      \tab \tab \eqn{H=\pi-\theta_H=-2*p_s(2*p_s-1)n_s/(n-1)}, \eqn{s}
      \eqn{\in} \eqn{[1,S]}. \cr 
      \eqn{j=2}:\tab \code{"H_1"}\tab 
      : Fay & Wu's \eqn{H} for the sample from population 1. \cr
      \tab \tab \eqn{H_1=\pi_1-\theta_H1=-2*q_s(2*q_s-1)n_1s/(n_1s-1)},
      \eqn{s} \eqn{\in} \eqn{[1,S_1]}. \cr 
      \eqn{j=3}:\tab \code{"H_2"}\tab 
      : Fay & Wu's \eqn{H} for the sample from population 2. \cr 
      \tab \tab \eqn{H_2=\pi_2-\theta_H2=-2*r_s(2*r_s-1)n_2s/(n_2s-1)},
      \eqn{s} \eqn{\in} \eqn{[1,S_2]}. 
    }
  }
  \item{$D_star}{
    The matrix \eqn{R} or \eqn{\Sigma{Y_r}}
     \eqn{\times} \code{length(summaries$stats$D_star)} 
    of summary statistics on the test of neutrality
    from the difference between the number of singletons and the number
    of segregating sites, Fu & Li's \eqn{D*} (Fu & Li 1993). \cr
    If \code{summaries$stats$D_star[}\eqn{i}\code{]=}\eqn{j},
    \eqn{i} and \eqn{j} \eqn{\in} \eqn{[1,3]},\cr
    \code{$D_star[}\eqn{l}\code{,}\eqn{i}\code{]=statslocus$D_star[}\eqn{j}\code{]},\cr
    where:
    \tabular{lll}{
      \eqn{j=1}:\tab \code{"D_star"}\tab 
      : Fu & Li's \eqn{D*} for a data set. \cr
      \tab \tab \eqn{D*=(S-S_o)/\sqrt(var(S-S_o))}. \cr
      \eqn{j=2}:\tab \code{"D_star1"}\tab 
      : Fu & Li's \eqn{D*} for the sample from population 1. \cr
      \tab \tab \eqn{D_1*=(S_1-S_o1)/\sqrt(var(S_1-S_o1))}. \cr
      \eqn{j=3}:\tab \code{"D_star2"}\tab 
      : Fu & Li's \eqn{D*} for the sample from population 2. \cr
      \tab \tab \eqn{D_2*=(S_2-S_o2)/\sqrt(var(S_2-S_o2))}.
    }
  }
    \item{}{
    Where: \tabular{ll}{
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

\author{
  Celine Becquet - \email{celine.becquet@gmail.com}. 
}

\seealso{
  The functions  \code{\link{simulate_data}} and
  \code{\link{estimate_IMc}} calls the function
  \code{\link{get_stats}}.\cr
  The functions \code{\link{stats_pop}} and  \code{\link{get_liststats}}
  output a data frame similar to 
  the one output by  \code{\link{get_stats}}. \cr
  The data frame output by \code{\link{get_stats}} is
  used as input for the function
  \code{\link{output_stats}}.
}

\examples{
### Create the inputs for the function.
data(simulation_files) # download the data

## Information on region
write.table(file="regions", x=info_region, row.name=FALSE, col.names=FALSE, quote=FALSE)
# Creates the file "regions" containing: 
info_region=read.table("regions",header =FALSE,skip=1)

## summary statistics considered
summaries=list()
summaries[[1]]=c("S_1","S_2","S_sl","S_sh","F_st","D_1","D_2","D_star1","D_star2")
liststats=list_stats
nstats=get_liststats(liststats=liststats,summaries=summaries)

## summary statistics values
msline="./ms 8 1 -t 1 -I 2 4 4 0 -ej 0.5 1 2 -seeds 1 2 3"
msframe=ms(msline) 
statslocus=stats_pop(msframe=msframe)

sumstats=list()
nloci=4
for(nlocus in 1:nloci)
{
## summary statistics values
msline="./ms 8 1 -t 1 -I 2 4 4 0 -ej 0.5 1 2"
msframe=ms(msline) 
statslocus=stats_pop(msframe=msframe)
sumstats=get_stats(nloci=nloci, nlocus=nlocus,sumstats=sumstats,statslocus=statslocus,info=info_region[nlocus,],summaries=nstats$stats)
}
sumstats

# Clean up the directory.
unlink("regions")
}

\keyword{misc}
