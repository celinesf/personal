\name{estimate_IMc}
\alias{estimate_IMc}

\title{
  Function to estimate the parameters of extensions of the
  Isolation-with-migration model. 
}
\description{
  The function \code{\link{estimate_IMc}} estimates the parameters of
  extensions of the Isolation-with-migration model. \cr 
  \code{\link{estimate_IMc}} uses the values and descriptions of
  the prior distributions for the fixed and estimated parameters,
  respectively, specified in the file \code{paramfile}, the data
  specified in the files \code{"}\code{datafile}\eqn{r}\code{"}, \eqn{r}
  \eqn{\in} \eqn{[1,R]}, as well as the information on
  the \eqn{R} independent genomic regions specified in the file
  \code{regfile}. \cr 
  If there are multiple loci with gaps in one or more genomic regions
  considered (as specified by \eqn{Y>1} in the file \code{regfile}), the
  function \code{\link{estimate_IMc}} also requires the information on
  the multiple loci to be specified in the file \code{locifile}. \cr 
  The function \code{\link{estimate_IMc}} writes the values at the 
  summary statistics listed in \code{summaries} of the simulated
  polymorphisms for the independent genomic 
  regions considered (and their loci is applicable) in the file named
  \code{"datafile-est-S"} (\code{"datafile-loci-est-S"}).\cr
  \code{\link{estimate_IMc}} also prints the string of characters
  with the parameter values and other information required to run the C
  program \code{estlike} into text files starting with  \code{"datafile-lik-"}.
}

\usage{
estimate_IMc(paramfile="estimation.par", datafile="datafiles/reg",
regfile="datafiles/info.reg", locifile="datafiles/info.loc", type="anc",
summaries="S_1 S_2 S_sl S_sh F_st D_1 D_2 D_star1 D_star2")
}

\arguments{
  \item{paramfile}{
    The name of the file with the values and descriptions of the
    prior distributions for the fixed and estimated parameters,
    respectively, required to estimate the parameters of extensions of the
    Isolation-with-migration model with the \R function
    \code{\link{estimate_IMc}}.\cr 
    By default the file name is \code{"estimation.par"}. \cr
    Comments start with \code{"\#"} and are ignored by the function
    \code{\link{estimate_IMc}}. The keywords and the values are separated
    by space or tabulations. \cr 
    \item{\strong{Specifying estimated parameters: }}{
      Any of the demographic parameters listed below can be either fixed
      (i.e., only one value specified after the keyword) or estimated. \cr 
      If a parameter is to be estimated, the information on the uniform prior
      distribution is specified with three values, two floating numbers and an
      integer as follows: \cr 
      \code{"}\emph{\code{param_to_estimate }}\eqn{\Theta_l}\code{ }\eqn{\Theta_u}\code{ }\eqn{\Theta_n}\code{"}. 
      \tabular{ll}{
        \eqn{\Theta_l}\tab 
        : The lower limit on the prior distribution
        range. \eqn{\Theta_l}\eqn{\ge}\eqn{0}. \cr 
        \eqn{\Theta_u} \tab 
        : The upper limit on the prior distribution range.
        \eqn{\Theta_l}\eqn{\le}\eqn{\Theta_u}. \cr 
        \eqn{\Theta_n} \tab 
        : The number of values to consider along the prior distribution
        range (including the bounderies). \eqn{\Theta_n>0} \cr 
      } ATTENTION: \code{\link{estimate_IMc}} requires that
      \eqn{0\le}\eqn{\Theta_l}\eqn{\le}\eqn{\Theta_u} and
      \eqn{\Theta_n>0}. \cr 
      If only one value is specified (i.e., \eqn{\Theta_l>0} and
      \eqn{\Theta_u} not specified or \eqn{\Theta_u=\Theta_l}), the
      parameter is considered fixed to \eqn{\Theta_l}.
    }
    The information on the fixed and/or estimated parameters for the
    IM-like model and optimization of the estimation method are specified
    with the eleven following keywords: 
    \item{\strong{REQUIRED PARAMETERS: }}{
      \item{\emph{\code{theta_1}}}{
        ATTENTION: \code{\link{estimate_IMc}} requires that
        either the parameter value  of the population mutation rate per bp per
        generation for population 1 (i.e. \eqn{\theta_1>0}) or the prior
        distribution for \eqn{\theta_1} (see above for the restrictions)
        be specified: \cr 
        \code{"}\emph{\code{theta_1}}\code{ }\eqn{\theta_1}\code{"} or
        \code{"}\emph{\code{theta_1}}\code{ }\eqn{\theta_1l}\code{ }\eqn{\theta_1u}\code{ }\eqn{\theta_1n}\code{"}. 
        \tabular{ll}{
          ---> Where: \tab \cr
          \eqn{\theta_1=4N_1*\mu} \tab
          : The population mutation rate per bp per generation for
          population 1 (required). \cr 
          \eqn{N_1} \tab 
          : The effective population size in population 1 (the reference
          population so by default and unless specified, \eqn{N_1=N_2=N_A}). \cr
          \eqn{\mu} \tab
          : The genomic generational mutation rate per bp. 
        }NOTE: If \eqn{\theta_1l=0}, \code{\link{estimate_IMc}}
        will generate \eqn{\theta_1n} possible values of \eqn{\theta_1} \eqn{\in}
        \eqn{[(\theta_1u-\theta_1l)/\theta_1n,\theta_1u]}. 
      }% end theta
      \item{\emph{\code{M_present}}}{
        Use \code{"}\emph{\code{M_present}}\code{ }\eqn{M_p}\code{"} or 
        \code{"}\emph{\code{M_present}}\code{ }\eqn{M_pl}\code{ }\eqn{M_pu}\code{ }\eqn{M_pn}\code{"}\cr
        to specify \eqn{M_p=4N_1*m_p}, the number of migrants exchanged
        each generation by the TWO populations at present. \cr
        By default there is no gene flow. \cr
        \eqn{M_p} is ignored unless there is population structure (as
        specified with \eqn{n_1>0} and \eqn{n_2>0} for the genomic region
        considered in the file \code{regfile}). \cr
        \eqn{M_p} has different definitions depending of the model
        considered: \cr 
        * \eqn{M_p} is the symmetrical rate of gene flow between
        TWO populations in an island model (i.e., \eqn{0<t<\inf}) as specified
        with: \cr 
        -- At least one of the independent genomic regions has
        \eqn{n_1>0} and \eqn{n_2>0} as specified in the file
        \code{regfile}. \cr 
        -- Either \eqn{T_s=0} when fixed and specified with the keyword
        \emph{\code{T_split}} or unspecified. \cr
        ATTENTION: \code{\link{estimate_IMc}} requires that
        \eqn{M_p>0} when fixed or
        \eqn{0\le}\eqn{M_pl}\eqn{\le}\eqn{M_pu} and \eqn{M_pn>0} when
        estimated be specified in case of an island model. \cr 
        * \eqn{M_p} is the constant symmetrical rate of gene flow
        since the split until present (i.e., \eqn{0<t<T_s}) if there is a
        population split (as specified with the keyword
        \emph{\code{T_split}} with \eqn{T_s>0} when fixed or
        \eqn{0\le}\eqn{T_sl}\eqn{\le}\eqn{T_su} and \eqn{T_sn>0} when
        estimated). \cr
        * \eqn{M_p} is the constant symmetrical rate of gene flow since
        the time of gene flow rate change until present (i.e., \eqn{0<t<T_c})
        if a time at which the gene flow rate changed is specified with: \cr 
        -- \eqn{0<T_c<T_s} (as specified with the keyword
        \emph{\code{T_change}} with \eqn{0<\epsilon<1} when fixed or
        \eqn{0\le}\eqn{\epsilon_l}\eqn{\le}\eqn{\epsilon_u}\eqn{\le1} and \eqn{\epsilon_n>0}
        when estimated). \cr
        -- \eqn{M_c \not=M_p} (as specified with the keywords
        \emph{\code{M_present}} and \emph{\code{M_change}} with
        \eqn{0\le}\eqn{M_c \not=M_p} when fixed).
        \tabular{l}{
          ---> Where, \eqn{m_p}
          : The generational fraction of migrant individuals at
          present. 
        }  
      }% end M_p
    }% end req arg
    \item{\strong{OPTIONAL DEMOGRAPHIC PARAMETERS: }}{
      All the parameters listed below are optional and can be either
      fixed or estimated. 
      \item{\emph{\code{theta_2}}}{
        Use \code{"}\emph{\code{theta_2}}\code{ }\eqn{\theta_2}\code{"} or
        \code{"}\emph{\code{theta_2}}\code{ }\eqn{\theta_2l}\code{ }\eqn{\theta_2u}\code{ }\eqn{\theta_2n}\code{"}\cr
        to specify \eqn{\theta_2=4N_2*\mu}, the population mutation rate
        per bp per generation for population 2. \cr 
        \eqn{\theta_2} is ignored unless there is population structure
        (as specified with \eqn{n_1>0} and \eqn{n_2>0} for the genomic region
        considered in the file \code{regfile}). 
        \tabular{l}{
          ---> Where, \eqn{N_2} 
          : The effective population size in population 2. 
        }
        ATTENTION: \code{\link{estimate_IMc}} requires that
        \eqn{\theta_2>0} when specified and fixed. \cr
        NOTE: If \eqn{\theta_2l=0}, \code{\link{estimate_IMc}}
        will sample a value of \eqn{\theta_2} among the \eqn{\theta_2n}
        possible values \eqn{\in}
        \eqn{[(\theta_2u-\theta_2l)/\theta_2n,\theta_2u]}.\cr
        If \eqn{\theta_2} is not specified the \eqn{\theta_2=\theta_1}.
      }% end theta2 
      \item{\emph{\code{theta_A}}}{
        Use \code{"}\emph{\code{theta_A}}\code{ }\eqn{\theta_A}\code{"} or
        \code{"}\emph{\code{theta_A}}\code{ }\eqn{\theta_Al}\code{ }\eqn{\theta_Au}\code{ }\eqn{\theta_An}\code{"}\cr
        to specify \eqn{\theta_A=4N_A*\mu}, the ancestral population
        mutation rate per bp per generation. \cr 
        \eqn{\theta_A} is ignored unless there is a population split (as
        specified with the keyword \emph{\code{T_split}} with \eqn{T_s>0} or
        \eqn{0\le}\eqn{T_sl}\eqn{\le}\eqn{T_su} and \eqn{T_sn>0} when
        estimated).
        \tabular{l}{
          ---> Where, \eqn{N_A} 
          : The ancestral effective population size. 
        }
        ATTENTION: \code{\link{estimate_IMc}} requires that
        \eqn{\theta_A>0} when specified and fixed.  \cr
        NOTE: If \eqn{\theta_Al=0}, \code{\link{estimate_IMc}}
        will sample a value of \eqn{\theta_A} among the \eqn{\theta_An}
        possible values \eqn{\in} \eqn{[(\theta_Au-\theta_Al)/\theta_An,\theta_Au]}.\cr
        If \eqn{\theta_A} is not specified the \eqn{\theta_A=\theta_1}.
      }% end thetaA
      \item{\emph{\code{T_split}}}{
        Use \code{"}\emph{\code{T_split}}\code{ }\eqn{T_s}\code{"} or 
        \code{"}\emph{\code{T_split}}\code{ }\eqn{T_sl}\code{ }\eqn{T_su}\code{ }\eqn{T_sn}\code{"}\cr
        to specify \eqn{T_s}, the split time in unit of \eqn{4N_1}
        generations between the TWO populations. \cr 
        By default \eqn{T_s=0} to specify a model without population split. \cr
        \eqn{T_s} is ignored unless there is population structure (as
        specified with \eqn{n_1>0} and \eqn{n_2>0} for the genomic region
        considered in the file \code{regfile}).
      }% end T_split
      \item{\emph{\code{T_change}}}{
        Use \code{"}\emph{\code{T_change}}\code{ }\eqn{\epsilon}\code{"} or 
        \code{"}\emph{\code{T_change}}\code{ }\eqn{\epsilon_l}\code{ }\eqn{\epsilon_u}\code{ }\eqn{\epsilon_n}\code{"}\cr
        to specify the ratio \eqn{\epsilon=T_c/T_s}. \cr
        By default there is no change of gene flow rate since the
        split. \cr 
        * \eqn{\epsilon} is ignored unless: \cr
        -- There is a population split (as specified with the keyword
        \emph{\code{T_split}} with \eqn{T_s>0} when fixed or
        \eqn{0\le}\eqn{T_sl}\eqn{\le}\eqn{T_su} and \eqn{T_sn>0} when
        estimated). \cr 
        -- \eqn{0<T_c<T_s} (as specified with \eqn{0<\epsilon<1} when fixed or
        \eqn{0\le}\eqn{\epsilon_l}\eqn{\le}\eqn{\epsilon_u}\eqn{\le1} and \eqn{\epsilon_n>0}
        when estimated). \cr 
        -- \eqn{M_c \not=M_p} (as specified with the keywords
        \emph{\code{M_present}} and \emph{\code{M_change}} with
        \eqn{0\le}\eqn{M_c \not=M_p} when fixed). 
        \tabular{l}{
          ---> Where, \eqn{T_c} 
          : The time at which the rate of gene flow changed between the
          two populations in unit of \eqn{4N_1} generations.
        } ATTENTION: \code{\link{estimate_IMc}} requires that
        \eqn{0}\eqn{\le}\eqn{T_c}\eqn{\le}\eqn{T_s} when specified and fixed
        with the keywords \emph{\code{T_split}} and
        \emph{\code{T_change}}.\cr
        NOTE: If \eqn{\epsilon_u=1}, \code{\link{estimate_IMc}}
        will sample a value of \eqn{\epsilon} among the \eqn{\epsilon_n}
        possible values \eqn{\in} \eqn{[\epsilon_l),\epsilon_u-(\epsilon_u-\epsilon_l)/\epsilon_n,]}.
      }% end T_change
      \item{\emph{\code{M_change}}}{
        Use \code{"}\emph{\code{M_change}}\code{ }\eqn{M_c}\code{"} or 
        \code{"}\emph{\code{M_change}}\code{ }\eqn{M_cl}\code{ }\eqn{M_cu}\code{ }\eqn{M_cn}\code{"}\cr
        to specify \eqn{M_c=4N_1*m_c}, the number of migrants exchanged
        each generation by the TWO populations since the split until the time
        of gene flow rate change (i.e., \eqn{T_c<t<T_s}). \cr 
        By default there is no change of gene flow rate since the
        split. \cr 
        * \eqn{M_c} is ignored unless: \cr
        -- There is a population split (as specified with the keyword
        \emph{\code{T_split}} with \eqn{T_s>0} when fixed or
        \eqn{0\le}\eqn{T_sl}\eqn{\le}\eqn{T_su} and \eqn{T_sn>0} when
        estimated). \cr 
        -- \eqn{0<T_c<T_s} (as specified with the keyword
        \emph{\code{T_change}} with \eqn{0<\epsilon<1} when fixed or
        \eqn{0\le}\eqn{\epsilon_l}\eqn{\le}\eqn{\epsilon_u}\eqn{\le1} and \eqn{\epsilon_n>0}
        when estimated). \cr 
        -- \eqn{M_c \not=M_p} (as specified with the keywords
        \emph{\code{M_present}} and \emph{\code{M_change}} with
        \eqn{0\le}\eqn{M_c \not=M_p} when fixed).
        \tabular{l}{
          ---> Where, \eqn{m_c}
          : The generational fraction of migrant individuals between
          \eqn{T_c<t<T_s}.
        }
      }% end M_c
    }% end Optional demo
    \item{\strong{OPTIONAL NUISANCE PARAMETER: }}{
      \item{\emph{\code{rho}}}{
        Use this keyword to specify the parameter on the intra-region
        recombination rate.\cr
        By default there is no recombination. \cr
        \emph{\code{rho}} is directly followed by either a floating
        number, \eqn{\rho=4N_1*c}, or the values
        \eqn{{\code{"-2"},\code{"-1"},\code{"1"},\code{"2"}}}. In all
        cases, the function \code{\link{estimate_IMc}} obtains the
        region-specific population recombination rate
        \eqn{\rho_r=\beta*4N_1c*(Z-1)} for the genomic region \eqn{r}
        \eqn{\in} \eqn{[1,R]} (as specified with the keyword
        \emph{\code{nregions}}):
        \tabular{lllll}{
          \code{"}\emph{\code{rho}} \tab \eqn{\rho}\code{"} \tab \tab \tab
          : This specifies that the genomic average population
          intra-region recombination rate per bp per generation is fixed to the
          value \eqn{\rho=4N_1*c}. \cr 
          \tab \tab \tab \tab In this case, \code{\link{estimate_IMc}} calculates
          \eqn{\rho_r=\rho*w*(Z-1)} fixed across estimation steps for the genomic region \eqn{r} (here
          \eqn{w=\beta}. \cr 
          \code{"}\emph{\code{rho}} \tab \code{1"} \tab \tab \tab 
          : This specifies that an estimate of the region-specific population
          recombination rate per bp, \eqn{\rho_o=4N_1*c_o}, is KNOWN from
          linkage disequilibrium analysis and specified with the parameter
          \eqn{w=\beta*\rho_o} for each recombining region (in the file
          \code{regfile}).\cr 
          \tab \tab \tab \tab In this case, \code{\link{estimate_IMc}}
          calculates \eqn{\rho_r=w*(Z-1)} fixed across estimation steps for the genomic region
          \eqn{r}. \cr 
          \code{"}\emph{\code{rho}} \tab \code{2} \tab \eqn{\mu}\code{"} \tab \tab
          : This specifies that an estimate of the region-specific
          recombination rate per bp, \eqn{c_o}, is KNOWN from pedigree
          analysis and specified with the parameter \eqn{w=\beta*c_o}
          for each recombining region (in the file \code{regfile}).\cr
          \tab \tab \tab \tab This also specifies \eqn{\mu}, an independent estimate of the
          genomic generational mutation rate per bp.\cr 
          \tab \tab \tab \tab In this case, \code{\link{estimate_IMc}}
          calculates \eqn{\rho_r=w*(Z-1)\theta_1/\mu} every estimation
          steps for the genomic region \eqn{r}. \cr 
          \code{"}\emph{\code{rho}} \tab \code{-1} \tab \eqn{\lambda}\code{"} \tab \tab 
          : This specifies that the intra-region recombination rate is
          UNKNOWN and the ratio of recombination over mutation
          rate for the genomic region \eqn{r}, \eqn{\alpha=c_r/\mu}, is drawn
          from an exponential distribution with mean \eqn{1/}\eqn{\lambda}. \cr 
          \tab \tab \tab \tab In this case, \code{\link{estimate_IMc}}
          calculates \eqn{\rho_r=w*\alpha*(Z-1)*theta_1} every estimation
          steps for the genomic region \eqn{r}. \cr 
          \code{"}\emph{\code{rho}} \tab \code{-2} \tab \eqn{\nu} \tab \eqn{\sigma}\code{"} \tab
          : This specifies that the intra-region recombination rate is
          UNKNOWN and the ratio of recombination over mutation
          rate for the genomic region \eqn{r}, \eqn{\alpha=c_r/\mu}, is drawn
          from an normal distribution with mean
          \eqn{\nu} and standard deviation \eqn{\sigma}. \cr
          \tab \tab \tab \tab In this case, \code{\link{estimate_IMc}}
          calculates \eqn{\rho_r=w*\alpha*(Z-1)*theta_1} every estimation
          steps for the genomic region \eqn{r}.
        }
        Where:
        \tabular{ll}{ 
          \eqn{\rho=4N_1*c}\tab
          : The genomic average population intra-region recombination rate per bp
          per generation. \cr 
          \eqn{\rho_r=\beta*(Z}\eqn{-}\eqn{1)*4N_1c}\tab 
          : The region-specific recombination rate per generation for the
          genomic region considered. \cr 
          \eqn{\rho_o=4N_1*c_o} \tab 
          : The estimate of the region-specific population recombination
          rate per bp per generation for the genomic region considered from
          linkage disequilibrium analysis. \cr
          \eqn{c} \tab 
          : The genomic generational cross-over rate per bp. \cr
          \eqn{c_o} \tab 
          : The estimate of the region-specific cross-over rate per bp
          per generation for the genomic region considered from pedigree
          analysis. \cr 
          \eqn{\alpha=c_r/\mu} \tab :Drawn from a prior distribution. \cr
          \eqn{c_r} \tab
          : The generational region-specific cross-over rate per bp for the
          genomic region \eqn{r}. \cr 
          \eqn{w} \tab 
          : The recombination scalar for the genomic region considered
          specified in the file \code{regfile}.\cr
          \eqn{\beta} \tab
          : The ratio of the region-specific population recombination
          rate per bp over \eqn{\rho=4N_1*c} for the genomic region
          considered.\cr
          \tab \eqn{\beta=}\code{"1"} (\code{"0.5"} in
          \emph{Drosophila}) for autosomal region, \code{"0.5"} for X- and
          \code{"0"} for Y- and mtDNA-linked region. \cr
          \eqn{Z} \tab
          : The size in bp of the genomic region considered. \cr
          \tab \eqn{z_s} and \eqn{z_e} (such as \eqn{Z=z_e-z_s}) are
          specified in the file \code{regfile}. 
        }
      }% end rho
    }% end nuisance
    \item{\strong{OPTIONAL OTHER PARAMETERS: }}{
      \item{\emph{\code{nregions}}}{
        Use \code{"}\emph{\code{nregions}}\code{ }\eqn{R}\code{"} to specify
        the number of independent genomic regions considered. \cr
        By default, \eqn{R=1}.\cr
        ATTENTION: \code{\link{estimate_IMc}} requires that \eqn{R>0}
        when specified with the keyword \emph{\code{nregions}}. \cr 
      }
      \item{\emph{\code{howmany}}}{
        Use \code{"}\emph{\code{howmany }}\eqn{H}\code{"} to specify
        \eqn{H}, the number of data sets to simulate per set of parameters
        (i.e., per grid point) to estimate the likelihood of the data given
        the set of parameters of the extension of the isolation-with-migration
        model.\cr 
        By default, \eqn{H=1000}.\cr
        ATTENTION: \code{\link{estimate_IMc}} requires that
        \eqn{H>0} when specified with the keyword \emph{\code{howmany}}. \cr 
      }% end howmany
      \item{\emph{\code{parallel}}}{
        Use \code{"}\emph{\code{parallel }}\eqn{J}\code{"} to specify
        \eqn{J}, the number of jobs to run in parallel to perform the
        estimation of the posterior distribution of the parameters of the
        extension of the isolation-with-migration model.\cr 
        By default, \eqn{J=1}, i.e., no parallelization.\cr
        ATTENTION: \code{\link{estimate_IMc}} requires that
        \eqn{J>0} when specified with the keyword
        \emph{\code{parallel}}. \cr 
      }% end parallel
    }% end extra
  }
  \item{datafile}{
    The argument \code{datafile} specifies the location and the string
    of characters that initiates the name of the text files with the data
    for each independent genomic region. \cr 
    By default the file name is \code{"reg"}
    located in the directory \code{"datafiles/"}. \cr
    In this case, the function \code{\link{estimate_IMc}} assumes that
    the data for the \eqn{R} independent genomic regions are in text files
    named \code{"}\code{datafile}\eqn{r}\code{"}, where \eqn{r} is the
    number of the genomic region considered, \eqn{r} \eqn{\in} \eqn{[1,R]}
    and \eqn{R} is specified with the keyword \emph{\code{nregions}} in
    the file \code{paramfile}. \cr 
    Each file should contain the data for the genomic region \eqn{r} in
    the usual \code{ms} format (with some variations from Hudson's
    program):
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
        \eqn{r}. \eqn{n=n_1+n_2}. \cr 
        \eqn{n_1} \tab
        : The sample size from population 1 for the genomic region
        \eqn{r}. \cr 
        \eqn{n_2} \tab
        : The sample size from population 2 for the genomic region
        \eqn{r}.
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
      \eqn{r} starts below.
    }
    \item{Line 5}{
      \code{"}\code{segsites: } \eqn{S}\code{"}, the number of
      segregating sites found in the data set for the genomic region
      \eqn{r}.
    }
    \item{Line 6}{
      \code{"}\code{position: }\code{"} followed by the \eqn{S}
      positions in base pairs of the \eqn{S} segregating sites found in the
      data set along the genomic region \eqn{r}. \cr 
      The position of the segregating site \eqn{s}
      \eqn{\in} \eqn{[1,S]} is \eqn{p_s} \eqn{\in} \eqn{[0,Z]}. \cr
      \eqn{Z=z_e-z_s} is the size of the genomic region
      \eqn{r} as specified in the file \code{regfile}. 
    }
    \item{Lines 7 to 7\emph{+n}}{
      \eqn{n=n_1+n_2} lines follow corresponding to the \eqn{n}
      haplotypes sampled for the genomic 
      region \eqn{r}. \cr 
      \eqn{n_1} and \eqn{n_2} are the sample sizes in the
      two populations for the genomic region \eqn{r} specified in the file
      \code{regfile}. \cr 
      Each line \eqn{j=i+7} \eqn{\in} \eqn{[7,7+n]} of
      haplotype is a string of \eqn{S} characters corresponding to the
      alleles found at the \eqn{S} segregating sites in the data set
      for the genomic region \eqn{r}. \cr 
      i.e., \code{"0"}, \code{"1"} or \code{"?"}, for
      ancestral (when known i.e., when the argument
      \code{type="anc"}), derived and missing allele, respectively. 
    }
  }
  \item{regfile}{
    The location and name of the file with the information on the
    independent genomic regions considered for the estimation. \cr 
    By default the file name is \code{"info.reg"}
    located in the directory \code{"datafiles/"}. \cr 
    The first line of the file is ignored. The information on the
    independent genomic regions starts on the SECOND line. \cr 
    Each genomic region is described by ten values: 
    \tabular{ll}{
      \eqn{r} \tab 
      : The genomic region number, \eqn{r} \eqn{\in} \eqn{[1,R]}. \cr 
      \tab \eqn{R} is the number of independent genomic regions
      considered as specified with the keyword \emph{\code{nregions}} in the
      file \code{paramfile}. \cr 
      \tab ATTENTION: The independent genomic region numbers, \eqn{r},
      need to be in order, i.e., region 1 starts on line 2, region 2 on
      line 3 \dots region \eqn{R} on line \eqn{R+1}. \cr 
      \emph{\code{Region name}} \tab 
      : The name of the genomic region \eqn{r}, which should contain at
      least ONE non-numerical character. \cr 
      \eqn{x} \tab 
      : The inheritance scalar for the genomic region \eqn{r} (i.e.,
      \code{"1"} for autosomal region, \code{"0.75"} for X- and \code{"0.5"}
      for Y- and mtDNA-linked region). \cr 
      \eqn{v} \tab 
      : The mutation rate scalar for the genomic region \eqn{r} (which
      can be estimated e.g., from divergence data). \cr
      \eqn{w} \tab 
      : The recombination scalar for the genomic region \eqn{r}.\cr
      \tab -- Usually \eqn{w=\beta}, the ratio of the locus-specific
      population recombination rate per bp over \eqn{\rho=4N_1*c}. \cr 
      \tab -- If an estimate of the region-specific population recombination
      rate per bp is available for each region from linkage
      disequilibrium analysis, \eqn{\rho_o=4N_1*c_o}, set
      \eqn{w=\beta*\rho_o} to incorporate this knowledge in the
      estimation (with \code{"}\emph{\code{rho}} \code{1"} in the
      file \code{paramfile}). \cr
      \tab In this case, \eqn{w} is the scaled
      sex-averaged region-specific population recombination rate per
      bp, i.e., for an X-linked locus \eqn{c_o} is the female
      recombination rate and \eqn{\beta=0.5} so that
      \eqn{\beta*\rho_o=2N_1*c_o}. \cr 
      \tab -- If an estimate of the region-specific recombination rate per bp
      is available for each region from pedigree analysis, \eqn{c_o},
      set \eqn{w=\beta*c_o} to incorporate this knowledge in the estimation
      (with \code{"}\emph{\code{rho}} \code{2 }\eqn{\mu}\code{"} in the file
      \code{paramfile}).\cr
      \tab In this case, \eqn{w} is the scaled sex-averaged
      region-specific recombination rate per bp, i.e., for an X-linked locus,
      \eqn{c_o} is the estimated female recombination rate so the
      scaled sex-averaged recombination rate is \eqn{\beta*c_o=0.5c}. \cr 
      \eqn{n_1} \tab 
      : The sample size from population 1 for the genomic region
      \eqn{r}. \cr 
      \eqn{n_2} \tab 
      : The sample size from population 2 for the genomic region
      \eqn{r}. \cr 
      \eqn{z_s} \tab 
      : The start position of the genomic region \eqn{r} in bp. \cr
      \eqn{z_e} \tab 
      : The end position of the genomic region \eqn{r} in bp. \cr
      \eqn{Y} \tab 
      : The number of loci spanning the genomic region \eqn{r}. \cr
      \tab \eqn{Y=1} means that the whole genomic region \eqn{r} is
      considered as a single locus. \cr 
      \tab \eqn{Y>1} means that multiple loci are considered within the
      genomic region \eqn{r} (the information on the loci should be given in
      the file \code{locifile}). 
    }
  }
  \item{locifile}{
    The location and name of the file with the information on the loci
    for the multi-locus genomic regions (as specified with \eqn{Y>1} in
    the file \code{regfile}). \cr 
    By default the file name is \code{"info.loc"}
    located in the directory \code{"datafiles/"}. \cr
    The first line of the file is ignored. The information on the loci
    starts on the SECOND line. \cr 
    Each locus is described by six values: 
    \tabular{ll}{
      \eqn{r} \tab 
      : The multi-locus genomic region number, \eqn{r} \eqn{\in}
      \eqn{[1,R]}, that this locus is part of. \cr 
      \tab \eqn{R} is the number of independent genomic regions
      considered as specified with the keyword \emph{\code{nregions}} in the
      file \code{paramfile}. \cr 
      \tab ATTENTION: The multi-locus genomic region numbers, \eqn{r},
      need to be in order, i.e., the information for the loci for the first
      multi-locus genomic region (\eqn{a}) are from line 2 to \eqn{Y_a+1},
      for the second multi-locus genomic region (\eqn{b}), the loci
      information starts on line \eqn{Y_a+2} to \eqn{Y_a+Y_b+1} \dots \eqn{a}
      and \eqn{b} \eqn{\in} \eqn{[1,R]}. \cr 
      \eqn{y} \tab 
      : The locus number, \eqn{y} \eqn{\in} \eqn{[1,Y]}. \cr
      \tab \eqn{Y} is the total number of loci spanning the multi-locus
      genomic region \eqn{r} as specified in the file \code{regfile}. \cr 
      \tab ATTENTION: The loci numbers, \eqn{y}, need to be in order,
      i.e., information for locus 1 of the first multi-locus genomic region
      \eqn{r} is on line 2, locus 2 on line 3 \dots locus \eqn{Y} on line
      \eqn{Y+1}. \cr 
      \eqn{n_1y} \tab 
      : The sample size from population 1 for the locus \eqn{y}. \cr
      \tab \eqn{n_1y}\eqn{\le}\eqn{n_1}, where \eqn{n_1} for the multi-locus
      genomic region \eqn{r} is specified in the file \code{regfile}. \cr 
      \eqn{n_2y} \tab 
      : The sample size from population 2 for the locus \eqn{y}. \cr
      \tab \eqn{n_2y}\eqn{\le}\eqn{n_2}, where \eqn{n_2} for the multi-locus
      genomic region \eqn{r} is specified in the file \code{regfile}. \cr 
      \eqn{z_sy} \tab 
      : The start position of the locus \eqn{y} in bp. \cr
      \tab \eqn{z_s1}\eqn{\ge}\eqn{z_s}, where \eqn{z_s} is the start position of
      the multi-locus genomic region \eqn{r} specified in the file
      \code{regfile}. \cr
      \eqn{z_ey} \tab 
      : The end position of the locus \eqn{y} in bp. \cr
      \tab \eqn{z_eY}\eqn{\le}\eqn{z_e} where \eqn{z_e} is the end position of
      the multi-locus genomic region \eqn{r} specified in the file
      \code{regfile}. 
    }
  }% end loci_file
  \item{type}{
    Use this argument to specify whether the ancestral alleles of the
    \eqn{R} independent genomic regions considered is known when the
    function \code{\link{estimate_IMc}} calculates the summary statistics
    of the data. \cr 
    \code{type} can take the two following values: 
    \item{\code{type="anc"}}{
      By default, \code{\link{estimate_IMc}} calculates the
      summary statistics for the data set of the genomic region \eqn{r}
      considering the character \code{"0"} at segregating site \eqn{s} as
      the ancestral allele (as specified in lines 7 to 7\eqn{+n} in the
      file \code{"}\code{datafile}\eqn{r}\code{"}), \eqn{s} \eqn{\in} \eqn{[1,S]} (as specified in
      line 5 in the file \code{"}\code{datafile}\eqn{r}\code{"}). \cr 
      \tab In turn, \code{\link{estimate_IMc}} considers the
      character \code{"1"} as the derived allele at segregating site
      \eqn{s}.
    }
    \item{\code{type="maf"}}{
      \code{\link{estimate_IMc}} considers the allele with the
      minimum frequency at segregating site \eqn{s} as the derived allele
      (regardless of the character specified in lines 7 to 7\eqn{+n} in the
      file \code{"}\code{datafile}\eqn{r}\code{"}).
    }
  }% end type

  \item{summaries}{
    The string of characters with the list of the summary statistics
    calculated for the \eqn{R} genomic regions and the \eqn{\Sigma Y_r} 
    loci, \eqn{r} \eqn{\in} \eqn{[1,R]} for the multi-loci genomic regions
    that \code{\link{estimate_IMc}} will use for the estimation.\cr
    By default \code{summaries="S_1 S_2 S_sl S_sh F_st D_1 D_2 D_star1 D_star2"}.
    The statistics are separated by spaces, can be either upper or
    lower case and  in any order. \eqn{47} summary statistics "words" are
    allowed and reconized (i.e., the ones output by the function
    \code{\link{stats_pop}}, see \code{\link{list_stats}}):
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
  }% end summaries
}% end of arguments

\details{
  The function \code{\link{estimate_IMc}} uses the function
  \code{\link{get_ms_output}} to transform the data in the
  \code{"}\code{datafile}\eqn{r}\code{"}s in a data frame in the same format as output by the
  \R function \code{\link{ms}}. \cr 
  The function calls the function \code{\link{stats_pop}} to calculate
  the summary statistics specified in the argument \code{summaries} for
  the  genomic region \eqn{r} and for the \eqn{Y} loci if the genomic
  region \eqn{r} has multiple loci (as specified with \eqn{Y>1} in the
  file \code{regfile}). \cr  
  \code{\link{stats_pop}} takes into account missing data in the
  calculations and whether or not the ancestral allele is known (as
  specified with the argument \code{type}). \cr
  \code{\link{estimate_IMc}} calls the function
  \code{\link{output_stats}} to output the summary statistics into
  files, as well as the function \code{\link{output_cmd}} to output the
  parameters values for each estimation steps considered.  
}

\value{
  \item{Summaries}{
    The function  \code{\link{estimate_IMc}} writes the summary statistics
    considered for the \eqn{R} independent genomic regions and the
    \eqn{\Sigma Y_r} uni-locus genomic
    regions and multiple locis, into two files named \code{"datafile-est-S"} and
    \code{"datafile-loci-est-S"}, respectively. \cr
    The statistics for the genomic regions without multiple loci are 
    reported in both files (i.e., if the \eqn{R} genomic regions are
    uni-locus, \code{"datafile-est-S"} and
    \code{"datafile-loci-est-S"} are indentical).
    The format of the files is as follows:
    \item{Line 1}{ The word "\code{type}".}
    \item{Line 2}{ The sentence describing the type of summary statistics
      calculated as specified with the argument \code{type}:
      \tabular{ll}{
        \code{"Known ancestral allele"}\tab
        : Specifies that  \code{\link{estimate_IMc}} calculated the summary statistics
        for the data sets considering the character \code{"0"} at each
        segregating site (as specified in the files \code{"}\code{datafile}\eqn{r}\code{"}) as the
        ancestral allele. \cr     
        \tab In turn,  \code{\link{estimate_IMc}} considered the character
        \code{"1"} as the derived allele at each segregating site.\cr
        \code{"Unknown ancestral allele"} \tab     
        : Specifies that  \code{\link{estimate_IMc}} considered the allele with
        the minimum frequency at each segregating site as the derived allele
        (regardless of the character specified in the files \code{"}\code{datafile}\eqn{r}\code{"}).
      }
    }
    \item{Line 3}{ The header for the summary statistics
      separated by tabulations.
      \code{"nregions nloci ..."} followed by the list of
      summary statistics specified in the argument \code{summaries} in
      order.
    }
    \item{Line  4 to 4\emph{+X}}{ The values for the summary statistics
      considered  calculated by the function \code{\link{stats_pop}} separated
      by tabulations
      for \eqn{X=R} independent genomic regions in
      \code{"datafile-est-S"} and  \eqn{X=\Sigma Y_r} uni-locus genomic
      regions and linked loci in  \code{"datafile-loci-est-S"}.
    }
     e.g.:\cr
    \code{$type}\cr
    \code{Known ancestral allele}\cr
     \tabular{lllllllllll}{
        \code{nregion}\tab \code{nloci}\tab \code{S_1}\tab \code{S_2}\tab \code{S_sl}\tab \code{S_sh}\tab \code{F_st}\tab \code{D_1}\tab \code{D_2}\tab \code{D_star1}\tab \code{D_star2}\cr
        \code{1}\tab \code{1}\tab \code{1}\tab \code{3}\tab \code{0}\tab \code{0}\tab \code{0.29}\tab \code{1.63}\tab \code{1.09}\tab \code{1.10}\tab \code{0.64}\cr
        \code{2}\tab \code{2}\tab \code{1}\tab \code{1}\tab \code{0}\tab \code{1}\tab \code{0.00}\tab \code{1.89}\tab \code{0.59}\tab \code{1.44}\tab \code{0.12}\cr
        \code{3}\tab \code{1}\tab \code{0}\tab \code{0}\tab \code{0}\tab \code{6}\tab \code{0.50}\tab \code{NA}\tab \code{NA}\tab \code{NA}\tab \code{-1.49}\cr
        \code{3}\tab \code{2}\tab \code{0}\tab \code{0}\tab \code{0}\tab \code{16}\tab \code{0.00}\tab \code{NA}\tab \code{NA}\tab \code{NA}\tab \code{3.93}\cr
        \code{4}\tab \code{1}\tab \code{1}\tab \code{0}\tab \code{0}\tab \code{0}\tab \code{0.00}\tab \code{-0.61}\tab \code{NA}\tab \code{-0.91}\tab \code{NA}\cr
        \code{4}\tab \code{2}\tab \code{0}\tab \code{0}\tab \code{0}\tab \code{0}\tab \code{NA}\tab \code{0.00}\tab \code{0.00}\tab \code{0.00}\tab \code{0.00}\cr
        \code{4}\tab \code{3}\tab \code{0}\tab \code{0}\tab \code{0}\tab \code{0}\tab \code{NA}\tab \code{0.00}\tab \code{0.00}\tab \code{0.00}\tab \code{0.00}\cr
        \code{4}\tab \code{4}\tab \code{0}\tab \code{0}\tab \code{0}\tab \code{3}\tab \code{-0.67}\tab \code{NA}\tab \code{NA}\tab \code{-1.15}\tab \code{NA}\cr
      }
  }% end summary files
  \item{Parameters}{
    The function  \code{\link{estimate_IMc}} writes the 
    parameter values for each of the estimation steps considered into
    \eqn{J} text files named  \code{"datafile-lik-}\eqn{j}\code{"}, \eqn{j}
    \eqn{\in} \eqn{[0,J-1]} (\eqn{J} is provided with the keyword
  \emph{\code{parallel}} in the file \code{paramfile}).\cr
    Each file contains \eqn{\Pi(\Theta_in}\eqn{)/J} sets of parameters,
  where \eqn{\Theta_i}
  is any of the estimated parameters out of seven parameters of the model.\cr
    The file numbered \eqn{J-1} may have less lines since some sets of
    parameters are ignored (see below).\cr
    Each line is as follows:
    \code{"-h }\eqn{h}\eqn{ -p }\eqn{theta_1}\code{ }\eqn{theta_2}\code{ }\eqn{theta_A}\code{ }\eqn{T_s}\code{ }\eqn{T_c}\code{ }\eqn{M_p}\code{ }\eqn{M_c}\code{"}\cr 
    Each line contains a unic set of parameter values.\cr All the combination
    of values for the estimate parameters are recorded: i.e. the
    combination of the \eqn{\Theta_in} possible values \eqn{\in}
    \eqn{[\Theta_il,\Theta_iu]}, with the following expections:\cr
    --- If \eqn{\theta_il=0} , the
    \eqn{\theta_in}  possible values are \eqn{\theta_i} \eqn{\in}
    \eqn{[(\theta_iu-\theta_il)/\theta_in,\theta_iu]}, \eqn{i} \eqn{\in} \eqn{{1,2,A}}.\cr
   If \eqn{\theta_i} is not specified the \eqn{\theta_i=\theta_1}, \eqn{i} \eqn{\in} \eqn{{2,A}}.
    The sets of parameters for which  \eqn{T_s=M_p=0} are not
    recorded.\cr
   --- For the steps with \eqn{T_s=0} only the steps with \eqn{\epsilon=M_c=0}
    are recorded.\cr
   --- For the steps with \eqn{\epsilon=0} only the steps with \eqn{M_c=0}
    are recorded.\cr
    --- For the steps with \eqn{M_c=M_p>0} only the steps with \eqn{\epsilon=0}
    are recorded.\cr
     e.g. of two lines:\cr
    \code{-h 1 -r 0.0005 -p 0.0001 0.0001 0.00015 0.222222 0 0 0}\cr
    \code{-h 1 -r 0.0005 -p 0.0001 0.0001 0.00015 0.222222 0.8 0 2.22222}
  }% end of files on parameters
}

\note{
  \item{ ATTENTION: It is the user's responsibility to mind the following
    restrictions:}{
    -> The user needs to be comfortable with the use of the C program
    \code{ms} and needs to be familiar with the documentation at
    \url{http://pps-spud.uchicago.edu/~cbecquet/msdoc.pdf}.
  }
  \item{In the file \code{paramfile}: }{
    -> \code{\link{estimate_IMc}} requires that either the parameter value
    (i.e. \eqn{\theta_1>0}) or the prior distribution for \eqn{\theta_1}
    (see above for the restrictions) be specified with the keyword
    \emph{\code{theta_1}}. \cr 
    -> \code{\link{estimate_IMc}} requires that \eqn{\theta_i>0} for any
    \eqn{i} \eqn{\in} \eqn{{1,2,A}} when these parameters are fixed (as specified
    with the keywords \emph{\code{theta_i}}). \cr
    -> \code{\link{estimate_IMc}} requires that \eqn{M_p>0} when fixed
    or \eqn{0\le}\eqn{M_pl}\eqn{\le}\eqn{M_pu} and \eqn{M_pn>0} when
    estimated be specified in case of an island model. \cr 
    -> \code{\link{estimate_IMc}} requires that
    \eqn{0}\eqn{\le}\eqn{T_c}\eqn{\le}\eqn{T_s} when specified and fixed
    with the keywords \emph{\code{T_split}} and
    \emph{\code{T_change}}.\cr
    -> \code{\link{estimate_IMc}} requires that \eqn{R>0} when specified
    with the keyword \emph{\code{nregions}}. \cr
    -> \code{\link{estimate_IMc}} requires that \eqn{H>0} when specified
    with the keyword \emph{\code{howmany}}. \cr 
    -> \code{\link{estimate_IMc}} requires that \eqn{J>0} when specified
    with the keyword \emph{\code{parallel}}. \cr
  }
  \item{ In the file \code{"}\code{datafile}\eqn{r}\code{"}:}{
    -> The \code{ms}-like command line (i.e., line 1 of the file
    \code{"}\code{datafile}\eqn{r}\code{"}) needs to contain at minimum these SEVEN "words"
    (The first word can be any string of characters without white space):
    \cr
    \code{"./ms }\eqn{n}\code{ 1 -I 2 }\eqn{n_1}\code{ }\eqn{n_2}\code{"}. \cr 
  }
  \item{In the file \code{regfile}: }{
    -> The information on the first genomic region starts on the SECOND
    line of the file. \cr 
    -> The independent genomic region numbers, \eqn{r}, need to be in
    order, i.e., region 1 starts on line 2, region 2 on line 3 \dots region
    \eqn{R} on line \eqn{R+1}. \cr 
    -> Reasonable values need to be specified for each genomic region. 
  }
  \item{In the file \code{locifile}:}{
    -> The information on the first locus of the first multi-locus
    genomic region starts on the SECOND line of the file. \cr 
    -> The loci numbers, \eqn{y}, need to be in order, i.e., information
    for locus 1 of the first multi-locus genomic region is on line 2,
    locus 2 on line 3 \dots locus \eqn{Y} on line \eqn{Y+1}. \cr 
    -> The multi-locus genomic region numbers, \eqn{r}, need to be in
    order, i.e., the information for the loci for the first multi-locus
    genomic region (\eqn{a}) are from line 2 to \eqn{Y_a+1}, for the
    second multi-locus genomic region (\eqn{b}), the loci information
    starts on line \eqn{Y_a+2} to \eqn{Y_a+Y_b+1} \dots \eqn{a} and \eqn{b}
    are \eqn{\in} \eqn{[1,R]}. \cr 
    -> Reasonable values need to be specified for each locus. 
  }
}

\author{
  Celine Becquet - \email{celine.becquet@gmail.com}. 
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

\seealso{
  The function \code{\link{estimate_IMc}} calls the following functions:
  \cr
  \code{\link{check_param}}, \code{\link{check_reg}},
  \code{\link{cmd_line}}, \code{\link{error_message}}, \code{\link{get_liststats}},
  \code{\link{get_loci}},
  \code{\link{get_ms_output}}, \code{\link{get_stats}}, \code{\link{ms}},
  \code{\link{order_param}},   \code{\link{output_cmd}},
  \code{\link{output_stats}} and \code{\link{stats_pop}}. \cr 
  Data file examples required as input for the function
  \code{\link{estimate_IMc}}: \cr 
  \code{\link{param_est}}, \code{\link{data_file}}, \code{\link{info_region}} and
  \code{\link{info_loci}}.\cr
  \code{\link{estimate_IMc}} uses as input data files in the same format
  as output by the functions \code{\link{simulate_data}} and
  \code{\link{write_ms_output}}. \cr
  Lists of definitions of the symbols and parameters mentioned in
  this file are found in \code{\link{Rmspack-package}}.
}

\examples{
### Write the files of information in the local directory. 
data(simulation_files) # download the data

write.table(file="sim", x=param_sim, row.name=FALSE, col.names=FALSE, quote=FALSE)
# Creates the file "sim" containing (with comments):
read.csv("sim",header =FALSE,sep="", comment="#")

write.table(file="est", x=param_est, row.name=FALSE, col.names=FALSE, quote=FALSE)
# Creates the file "est" containing (with comments):
read.csv("est",header =FALSE,sep="", comment="#")

write.table(file="regions", x=info_region, row.name=FALSE, col.names=FALSE, quote=FALSE)
# Creates the file "regions" containing:
read.csv("regions",header =FALSE,sep="", comment="#")

write.table(file="loci", x=info_loci, row.name=FALSE, col.names=FALSE, quote=FALSE)
# Creates the file "loci" containing:
read.csv("loci",header =FALSE,sep="", comment="#")

## Create 4 files named regr (r in [1,4]) with simulated data and files
## of summaries reg-S and reg-loci-S. 
simulate_data(paramfile="sim", regfile="regions", locifile="loci");

## Estimate parameters & create files of summaries "reg-est-S" and "reg-loci-est-S" and of cmd lines "reg-lik-0" (containing 281930 lines - the last five shown).
estimate_IMc(paramfile="est", datafile="reg", regfile="regions", locifile="loci");
read.csv("reg-est-S",header =FALSE,sep="", comment="#")
read.csv("reg-loci-est-S",header =FALSE,sep="", comment="#")
read.csv("reg-lik-0",header =FALSE,sep="", comment="#", skip=181895)

# Clean up the directory.
unlink(c("est", "sim", "regions", "loci", "reg1", "reg2", "reg3", "reg4","reg-S","reg-loci-S","reg-loci-est-S","reg-est-S","reg-lik-0"))
}

\keyword{IO}
\keyword{programming}
\keyword{methods}