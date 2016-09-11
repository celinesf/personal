\name{list_stats}
\alias{list_stats}
\docType{data}

\title{
  The complete list of names of the summary statistics of the polymorphism.
}
\description{
  The file \code{\link{list_stats}} provides the list of \eqn{47} strings of
  characters describing the  summary statistics in the order output by
  \code{\link{stats_pop}}. 
}

\usage{
data(simulation_files) # download the data
list_stats
}

\format{
  The vector \code{\link{list_stats}} contains the list of \eqn{47} strings of
  characters describing the summary statistics in the order output by
  \code{\link{stats_pop}}:
  \tabular{llllllllll}{
    \code{1}\tab \code{"S"}\tab \code{"S_1"}\tab \code{"S_2"}\tab
    \code{"S_s"}\tab \code{"S_ss"} \tab \code{"S_sf1"}\tab \code{"S_sf2"}\tab
    \code{"S_sl"}\tab \code{"S_sh"} \cr 
    \code{10}\tab \code{"S_f"}\tab \code{"S_f1"}\tab \code{"S_f2"}\tab
    \code{"S_o"}\tab \code{"S_o1"}\tab \code{"S_o2"}\tab \code{"F(S)"}\tab
    \code{"F(S_1)"}\tab \code{"F(S_2)"}\cr
    \code{19} \tab \code{"F(S_s)"}\tab \code{"F(S_ss)"}\tab \code{"F(S_sf1)"}\tab
    \code{"F(S_sf2)"}\tab \code{"F(S_sl)"}\tab \code{"F(S_sh)"}\tab
    \code{"F_st"}\tab \code{"H_b"}\tab \code{"H_w1"} \cr 
    \code{28}\tab \code{"H_w2"}\tab \code{"Snn"} \tab \code{"pi"}\tab
    \code{"pi_1"}\tab \code{"pi_2"}\tab \code{"theta_W"}\tab
    \code{"theta_W1"}\tab \code{"theta_W2"}\tab \code{"theta_H"}\cr
    \code{37}\tab \code{"theta_H1"}\tab \code{"theta_H2"}\tab \code{"D"}\tab
    \code{"D_1"}\tab \code{"D_2"}\tab \code{"H"}\tab \code{"H_1"}\tab
    \code{"H_2"}\tab \code{"D_star"}\cr 
    \code{46} \tab \code{"D_star1"}\tab \code{"D_star2"} \tab \tab
    \tab \tab \tab \tab \tab 
  } 
}

\details{
  The strings of characters correspond to the following summary statistics:
  \tabular{ll}{
    \code{"S"} \tab 
    : The total number of segregating sites found in the data set for
    a genomic region or locus. \cr
    \tab \eqn{S=S_1+S_2+S_s+S_f}. \cr
    \code{"S_1"} \tab 
    : The number of segregating sites with the derived alleles
    polymorphic uniquely in the sample from population 1. \cr 
    \code{"S_2"} \tab 
    : The number of segregating sites with the derived alleles
    polymorphic uniquely in the sample from population 2. \cr 
    \code{"S_s"} \tab 
    : The number of segregating sites with the derived allele shared
    by both population samples. \cr 
    \ \tab \eqn{S_s=S_ss+S_sfs+S_sf2} and \eqn{S_s=S_sl+S_sh}. \cr
    \code{"S_ss"} \tab 
    : The number of segregating sites with the derived allele shared
    AND polymorphic in both population samples. \cr 
    \code{"S_sf1"} \tab 
    : The number of segregating sites with the derived allele shared
    AND polymorphic ONLY in the sample from population 2 (i.e., fixed in
    the sample from population 1). \cr 
    \code{"S_sf2"} \tab
    : The number of segregating sites with the derived allele shared
    AND polymorphic ONLY in the sample from population 1 (i.e., fixed in
    the sample from population 2). \cr
    \code{"S_sl"} \tab 
    : The number of segregating sites with the derived allele shared
    AND  at low frequency (i.e., \eqn{p_s \le 10\%}). \cr 
    \code{"S_sh"} \tab
    : The number of segregating sites with the derived allele shared
    AND  at high frequency (i.e., \eqn{p_s > 10\%}. \cr
    \code{"S_f"} \tab 
    : The number of segregating sites with the alleles fixed between
    the population sampled. \cr 
    \tab \eqn{S_f=S_f1+S_f2}. \cr
    \code{"S_f1"} \tab 
    : The number of segregating sites with the derived allele fixed in
    the sample from population 1. \cr 
    \code{"S_f2"} \tab 
    : The number of segregating sites with the derived allele fixed in
    the sample from population 2. \cr
    \code{"S_o"} \tab
    : The number of singletons found in the data set for
    a genomic region or locus. \cr 
    \tab i.e.: \eqn{S_o=S_o1+S_o2} \eqn{\in} \eqn{[0,S]}, where
    the frequency of the derived allele is \eqn{p_s=1/n_s}. \cr 
    \code{"S_o1"} \tab
    : The number of singletons found in the sample from population
    1. \cr 
    \tab i.e.: \eqn{S_o1} \eqn{\in} \eqn{[0,S_1]}, where the
    frequency of the derived allele is \eqn{q_s=1/n_1s}. \cr 
    \code{"S_o2"} \tab
    : The number of singletons found in the sample from population
    2. \cr 
    \tab i.e.: \eqn{S_o2} \eqn{\in} \eqn{[0,S_2]}, where the
    frequency of the derived allele is \eqn{r_s=1/n_2s}. \cr
    \code{"F(S)"} \tab 
    : The mean derived allele frequency over \eqn{n} haplotypes for
    the \eqn{S} segregating sites in the data set for
    a genomic region or locus. \cr 
    \code{"F(S_1)"} \tab 
    : The mean derived allele frequency over \eqn{n_1} haplotypes for
    the \eqn{S_1} segregating sites with the derived alleles polymorphic
    uniquely in the sample from population 1. \cr 
    \code{"F(S_2)"} \tab 
    : The mean derived allele frequency over \eqn{n_2} haplotypes for
    the \eqn{S_2} segregating sites with the derived alleles polymorphic
    uniquely in the sample from population 2. \cr 
    \code{"F(S_s)"} \tab
    : The mean derived allele frequency over \eqn{n} haplotypes for
    the \eqn{S_s} segregating sites with the derived allele shared by both
    population samples. \cr 
    \code{"F(S_ss)"} \tab
    : The mean derived allele frequency over \eqn{n} haplotypes for
    the \eqn{S_ss} segregating sites with the derived allele shared AND
    polymorphic in both population samples. \cr 
    \code{"F(S_sf1)"} \tab 
    : The mean derived allele frequency over \eqn{n} haplotypes for
    the \eqn{S_sf1} segregating sites with the derived allele shared AND
    polymorphic ONLY in the sample from population 2 (i.e., fixed in the
    sample from population 1). \cr 
    \code{"F(S_sf2)"} \tab 
    : The mean derived allele frequency over \eqn{n} haplotypes for
    the \eqn{S_sf2} segregating sites with the derived allele shared AND
    polymorphic ONLY in the sample from population 1 (i.e., fixed in the
    sample from population 2). \cr
    \code{"F(S_sl)"} \tab 
    : The mean derived allele frequency over \eqn{n} haplotypes for
    the \eqn{S_sl} segregating sites with the derived allele shared
    AND at low frequency (i.e., \eqn{p_s \le 10\%}). \cr 
    \code{"F(S_sh)"} \tab 
    : The mean derived allele frequency over \eqn{n} haplotypes for
    the \eqn{S_sh} segregating sites with the derived allele shared
    AND at high frequency (i.e., \eqn{p_s > 10\%}).\cr
    \code{"F_st"} \tab 
    : The level of differentiation between the two population
    samples. \cr 
    \tab \eqn{F_st=1-H_w/H_b}, where \eqn{H_w=H_w1+H_w2} (Wright 1931;
    Hudson et al. 1992). \cr 
    \code{"H_b"} \tab 
    : The mean pairwise differences between both population
    samples. \cr 
    \code{"H_w1"} \tab 
    : The mean pairwise differences within the sample from population
    1. \cr 
    \code{"H_w2"} \tab 
    : The mean pairwise differences within the sample from population
    2. \cr 
    \code{"S_nn"} \tab 
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
    \eqn{i} that are from the same population sample as haplotype
    \eqn{i}. \cr
    \code{"pi"} \tab
    : \eqn{\pi} for the data set for
    a genomic region or locus. \cr
    \tab \eqn{\pi=\Sigma[2*p_s*(1-p_s)*n_s/(n_s-1)]}, \eqn{s}
    \eqn{\in} \eqn{[1,S]}. \cr 
    \code{"pi_1"} \tab 
    : \eqn{\pi} for the sample from population 1. \cr
    \tab \eqn{\pi_1=\Sigma[2*q_s*(1-q_s)*n_1s/(n_1s-1)]}, \eqn{s}
    \eqn{\in} \eqn{[1,S_1]}. \cr 
    \code{"pi_2"} \tab 
    : \eqn{\pi} for the sample from population 2. \cr
    \tab \eqn{\pi_2=\Sigma[2*r_s*(1-r_s)*n_2s/(n_2s-1)]}, \eqn{s}
    \eqn{\in} \eqn{[1,S_2]}. \cr
    \code{"theta_W"} \tab 
    : \eqn{\theta_W} for the data set for
    a genomic region or locus. \cr
    \tab \eqn{\theta_W=\Sigma[1/(1/n_s)]}, \eqn{s} \eqn{\in}
    \eqn{[1,S]}. \cr 
    \code{"theta_W1"} \tab 
    : \eqn{\theta_W} for the sample from population 1. \cr
    \tab \eqn{\theta_W1=\Sigma[1/(1/n_1s)]}, \eqn{s} \eqn{\in}
    \eqn{[1,S_1]}. \cr 
    \code{"theta_W2"} \tab 
    : \eqn{\theta_W} for the sample from population 2. \cr
    \tab \eqn{\theta_W2=\Sigma[1/(1/n_2s)]}, \eqn{s} \eqn{\in}
    \eqn{[1,S_2]}. \cr
    \code{"theta_H"} \tab 
    : \eqn{\theta_H} for the data set for
    a genomic region or locus. \cr
    \tab \eqn{\theta_H=\Sigma[2*(p_s^2)/[n_s(n_s-1)]]}, \eqn{s}
    \eqn{\in} \eqn{[1,S]}. \cr 
    \code{"theta_H1"} \tab 
    : \eqn{\theta_H} for the sample from population 1. \cr
    \tab \eqn{\theta_H1=\Sigma[2*(q_s^2)/[n_1s(n_1s-1)]]}, \eqn{s}
    \eqn{\in} \eqn{[1,S_1]}. \cr 
    \code{"theta_H2"} \tab 
    : \eqn{\theta_H} for the sample from population 2. \cr
    \tab \eqn{\theta_H2=\Sigma[2*(r_s^2)/[n_2s(n_2s-1)]]}, \eqn{s}
    \eqn{\in} \eqn{[1,S_2]}. \cr
    \code{"D"} \tab 
    : Tajima's \eqn{D} for the data set for
    a genomic region or locus. \cr
    \tab \eqn{D=(\pi-\theta_W)/\sqrt(var(\pi-\theta_W))}. \cr 
    \code{"D_1"} \tab 
    : Tajima's \eqn{D} for the sample from population 1. \cr
    \tab \eqn{D_1=(\pi_1-\theta_W1)/\sqrt(var(\pi_1-\theta_W1))}. \cr 
    \code{"D_2"} \tab 
    : Tajima's \eqn{D} for the sample from population 2. \cr 
    \tab \eqn{D_2=(\pi_2-\theta_W2)/\sqrt(var(\pi_2-\theta_W2))}. \cr
    \code{"H"} \tab 
    : Fay & Wu's \eqn{H} for the data set for
    a genomic region or locus. \cr
    \tab \eqn{H=\pi-\theta_H=-2*p_s(2*p_s-1)n_s/(n-1)}, \eqn{s}
    \eqn{\in} \eqn{[1,S]}. \cr 
    \code{"H_1"} \tab 
    : Fay & Wu's \eqn{H} for the sample from population 1. \cr
    \tab \eqn{H_1=\pi_1-\theta_H1=-2*q_s(2*q_s-1)n_1s/(n_1s-1)},
    \eqn{s} \eqn{\in} \eqn{[1,S_1]}. \cr 
    \code{"H_2"} \tab 
    : Fay & Wu's \eqn{H} for the sample from population 2. \cr 
    \tab \eqn{H_2=\pi_2-\theta_H2=-2*r_s(2*r_s-1)n_2s/(n_2s-1)},
    \eqn{s} \eqn{\in} \eqn{[1,S_2]}. \cr
    \code{"D_star"} \tab 
    : Fu & Li's \eqn{D*} for the data set for
    a genomic region or locus. \cr
    \tab \eqn{D*=(S-S_o)/\sqrt(var(S-S_o))}. \cr
    \code{"D_star1"} \tab 
    : Fu & Li's \eqn{D*} for the sample from population 1. \cr
    \tab \eqn{D_1*=(S_1-S_o1)/\sqrt(var(S_1-S_o1))}. \cr
    \code{"D_star2"} \tab 
    : Fu & Li's \eqn{D*} for the sample from population 2. \cr
    \tab \eqn{D_2*=(S_2-S_o2)/\sqrt(var(S_2-S_o2))}.
  }Where: \tabular{ll}{
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
    \eqn{s} in the data set for
    a genomic region or locus, \eqn{s} \eqn{\in} \eqn{[1,S]}. \cr 
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
  The vector of summary statistic names \code{\link{list_stats}} is
  within the the functions \code{\link{estimate_IMc}} and
  \code{\link{simulate_data}} and is used as input for
  the function \code{\link{get_liststats}}.\cr
  \code{\link{list_stats}} is also used in the examples of
  \code{\link{get_stats}} and \code{\link{output_stats}}.
}

\examples{
### Create the data file.
data(simulation_files) # download the data
list_stats
}

\keyword{datasets}
