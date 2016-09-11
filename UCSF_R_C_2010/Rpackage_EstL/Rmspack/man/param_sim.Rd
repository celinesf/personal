\name{param_sim}
\alias{param_sim}
\docType{data}

\title{
  Example of a file with parameter values to simulate data. 
}

\description{
  The file \code{\link{param_sim}} provides an example of the parameter file
  \code{paramfile}, which contains the parameter values required to
  simulate data with the \R function \code{\link{simulate_data}}. 
}

\usage{
data(param_sim)
}

\format{
  The file \code{\link{param_sim}} provides an example of the parameter file
  \code{paramfile}, which contains the parameter values required to
  simulate data with the \R function \code{\link{simulate_data}}. \cr 
  Comments start with \code{"\#"} but are absent in
  \code{\link{param_est}}. The keywords and the values are separated by
  space or tabulations. \cr 
  \item{\strong{Specifying variable parameters: }}{
    Any of the demographic parameters listed below can be either fixed
    (i.e., only one value specified after the keyword) or variable accross independent genomic regions. \cr 
    If a parameter is to be variable accross independent genomic regions, the information on the uniform prior
    distribution is specify with three values, two floating numbers and an
    integer as follows: \cr 
    \code{"}\emph{\code{variable_param }}\eqn{\Theta_l}\code{ }\eqn{\Theta_u}\code{ }\eqn{\Theta_n}\code{"}. 
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
    } ATTENTION: It is required that
    \eqn{0\le}\eqn{\Theta_l}\eqn{\le}\eqn{\Theta_u} and
    \eqn{\Theta_n>0}. \cr 
    If only one value is specified (i.e., \eqn{\Theta_l>0} and
    \eqn{\Theta_u} not specified or \eqn{\Theta_u=\Theta_l}), the
    parameter is considered fixed to \eqn{\Theta_l}.
  }
   The information on the fixed and/or variable accross independent
    genomic regions parameters for the IM-like model  are specified
    with the nine following keywords: 
  \item{\strong{REQUIRED PARAMETERS: }}{
    \item{\emph{\code{theta_1}}}{
      ATTENTION: It is required that either the parameter value  of the population mutation rate per bp per
      generation for population 1
      (i.e. \eqn{\theta_1>0}) or the prior distribution for \eqn{\theta_1}
      (see above for the restrictions) be specified: \cr 
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
      } NOTE: If \eqn{\theta_1l=0}, \code{\link{simulate_data}}
        will sample a value of \eqn{\theta_1} among the \eqn{\theta_1n}
        possible values \eqn{\in} \eqn{[(\theta_1u-\theta_1l)/\theta_1n,\theta_1u]}. 
    }% end theta 
    \item{\emph{\code{M_present}}}{
      Use \code{"}\emph{\code{M_present}}\code{ }\eqn{M_p}\code{"} or 
      \code{"}\emph{\code{M_present}}\code{ }\eqn{M_pl}\code{ }\eqn{M_pu}\code{ }\eqn{M_pn}\code{"}\cr
      to specify \eqn{M_p=4N_1*m_p}, the number of migrants exchanged
      each generation by the TWO populations at present. \cr 
      By default there is no gene flow. \cr
      \eqn{M_p} is ignored unless there is population structure (as
      specified with \eqn{n_1>0} and \eqn{n_2>0} for the genomic region
      considered in the file \code{regfile}, see \code{\link{info_region}}). \cr
      \eqn{M_p} has different definitions depending of the model considered: \cr 
      * \eqn{M_p} is the symmetrical rate of gene flow between
      TWO populations in an island model (i.e., \eqn{0<t<\inf}) as specified
      with: \cr 
      -- At least one of the independent genomic regions has
      \eqn{n_1>0} and \eqn{n_2>0} as specified in the file
      \code{regfile} (see \code{\link{info_region}}). \cr 
      -- Either \eqn{T_s=0} when fixed and specified with the keyword
      \emph{\code{T_split}} or unspecified. \cr
      ATTENTION: It is required that \eqn{M_p>0} when fixed or
      \eqn{0\le}\eqn{M_pl}\eqn{\le}\eqn{M_pu} and \eqn{M_pn>0} when
      variable accross independent genomic regions be specified in case
      of an island model.
      In this case and even if
  \eqn{M_pl=0}, \code{\link{simulate_data}} will sample a value
  \eqn{M_p>0} among the \eqn{M_pn} 
        possible values \eqn{\in} \eqn{[(M_pu-M_pl)/M_pn,M_pu]}. \cr
      * \eqn{M_p} is the constant symmetrical rate of gene flow
      since the split until present (i.e., \eqn{0<t<T_s}) if
      there is a population split (as specified with the keyword
      \emph{\code{T_split}} with \eqn{T_s>0} when fixed or
      \eqn{0\le}\eqn{T_sl}\eqn{\le}\eqn{T_su} and \eqn{T_sn>0} when
      variable accross independent genomic regions). \cr 
      * \eqn{M_p} is the constant symmetrical rate of gene flow since
      the time of gene flow rate change until present (i.e., \eqn{0<t<T_c})
      if a time at which the gene flow rate changed is specified with: \cr 
      -- \eqn{0<T_c<T_s} (as specified with the keyword
      \emph{\code{T_change}} with \eqn{0<\epsilon<1} when fixed or
      \eqn{0\le}\eqn{\epsilon_l}\eqn{\le}\eqn{\epsilon_u}\eqn{\le1} and \eqn{\epsilon_n>0}
      when variable accross independent genomic regions). \cr 
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
    fixed or variable accross independent genomic regions. 
    \item{\emph{\code{theta_2}}}{
      Use \code{"}\emph{\code{theta_2}}\code{ }\eqn{\theta_2}\code{"} or
      \code{"}\emph{\code{theta_2}}\code{ }\eqn{\theta_2l}\code{ }\eqn{\theta_2u}\code{ }\eqn{\theta_2n}\code{"}\cr
      to specify \eqn{\theta_2=4N_2*\mu}, the population mutation rate
      per bp per generation for population 2. \cr 
      \eqn{\theta_2} is ignored unless there is population structure
      (as specified with \eqn{n_1>0} and \eqn{n_2>0} for the genomic region
      considered in the file \code{regfile}, see
      \code{\link{info_region}}). 
      \tabular{l}{
        ---> Where, \eqn{N_2} 
        : The effective population size in population 2. 
      } ATTENTION: It is required that \eqn{\theta_2>0} when specified and fixed.\cr
        NOTE: If \eqn{\theta_2l=0}, \code{\link{simulate_data}}
        will sample a value of \eqn{\theta_2} among the \eqn{\theta_2n}
        possible values \eqn{\in} \eqn{[(\theta_2u-\theta_2l)/\theta_2n,\theta_2u]}.\cr
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
      variable accross independent genomic regions). 
      \tabular{l}{
        ---> Where, \eqn{N_A} 
        : The ancestral effective population size. 
      } ATTENTION: It is required that \eqn{\theta_A>0} when specified
      and fixed.\cr
        NOTE: If \eqn{\theta_Al=0}, \code{\link{simulate_data}}
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
      considered in the file \code{regfile}, see \code{\link{info_region}}).\cr
        NOTE, in the case were a value \eqn{M_p=0} has been sampled,
  \code{\link{simulate_data}} will sample a value \eqn{T_s>0} among the \eqn{T_sn}
        possible values \eqn{\in} \eqn{[(T_su-T_sl)/T_sn,T_su]}.\cr  
    }% end T_split
    \item{\emph{\code{T_change}}}{
      Use \code{"}\emph{\code{T_change}}\code{ }\eqn{\epsilon}\code{"} or 
      \code{"}\emph{\code{T_change}}\code{ }\eqn{\epsilon_l}\code{ }\eqn{\epsilon_u}\code{ }\eqn{\epsilon_n}\code{"}\cr
      to specify the ratio \eqn{\epsilon=T_c/T_s}.\cr 
      By default there is no change of gene flow rate since the
      split. \cr 
      * \eqn{\epsilon} is ignored unless: \cr
      -- There is a population split (as specified with the keyword
      \emph{\code{T_split}} with \eqn{T_s>0} when fixed or
      \eqn{0\le}\eqn{T_sl}\eqn{\le}\eqn{T_su} and \eqn{T_sn>0} when
      variable accross independent genomic regions). \cr 
      -- \eqn{0<T_c<T_s} (as specified with \eqn{0<\epsilon<1} when fixed or
      \eqn{0\le}\eqn{\epsilon_l}\eqn{\le}\eqn{\epsilon_u}\eqn{\le1} and \eqn{\epsilon_n>0}
      when variable accross independent genomic regions). \cr
      -- \eqn{M_c \not=M_p} (as specified with the keywords
      \emph{\code{M_present}} and \emph{\code{M_change}} with
      \eqn{0\le}\eqn{M_c \not=M_p} when fixed). 
      \tabular{l}{
        ---> Where, \eqn{T_c}
        : The time at which the rate of gene flow changed between the
        two populations in unit of \eqn{4N_1} generations.
      } ATTENTION: It is required that
      \eqn{0}\eqn{\le}\eqn{T_c}\eqn{\le}\eqn{T_s} when specified and fixed
      with the keywords \emph{\code{T_split}} and
      \emph{\code{T_change}}. \cr
         NOTE: If \eqn{\epsilon_u=1}, \code{\link{simulate_data}}
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
      variable accross independent genomic regions). \cr 
      -- \eqn{0<T_c<T_s} (as specified with the keyword
      \emph{\code{T_change}} with \eqn{0<\epsilon<1} when fixed or
      \eqn{0\le}\eqn{\epsilon_l}\eqn{\le}\eqn{\epsilon_u}\eqn{\le1} and \eqn{\epsilon_n>0}
      when variable accross independent genomic regions). \cr 
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
      Use the keyword \code{"}\emph{\code{rho}}\code{"} to
      specify the intra-region recombination rate. \cr 
      By default there is no recombination. \cr 
      \emph{\code{rho}} is directly followed by either a floating
      number, \eqn{\rho=4N_1*c}, or the values
      \eqn{{\code{"-2"},\code{"-1"},\code{"1"},\code{"2"}}}. In all cases,
      the function \code{\link{simulate_data}} obtains the region-specific
      population recombination rate \eqn{\rho_r=\beta*4N_1c*(Z-1)} for the
      genomic region \eqn{r} \eqn{\in} \eqn{[1,R]} (as specified with the
      keyword \emph{\code{nregions}}):
      \tabular{lllll}{
        \code{"}\emph{\code{rho}} \tab \eqn{\rho}\code{"} \tab \tab \tab
        : This specifies that the genomic average population
        intra-region recombination rate per bp per generation is fixed to the
        value \eqn{\rho=4N_1*c}. \cr 
        \tab \tab \tab \tab In this case, \code{\link{simulate_data}} calculates
        \eqn{\rho_r=\rho*w*(Z-1)} for the genomic region \eqn{r} (here
        \eqn{w=\beta}). \cr 
        \code{"}\emph{\code{rho}} \tab \code{1"} \tab \tab \tab 
        : This specifies that an estimate of the region-specific population
        recombination rate per bp, \eqn{\rho_o=4N_1*c_o}, is KNOWN from
        linkage disequilibrium analysis and specified with the parameter
        \eqn{w=\beta*\rho_o} for each recombining region (in the file
        \code{regfile}, see \code{\link{info_region}}).\cr 
        \tab \tab \tab \tab In this case, \code{\link{simulate_data}}
        calculates \eqn{\rho_r=w*(Z-1)} for the genomic region \eqn{r}. \cr 
        \code{"}\emph{\code{rho}} \tab \code{2} \tab \eqn{\mu}\code{"} \tab \tab
        : This specifies that an estimate of the region-specific
        recombination rate per bp, \eqn{c_o}, is KNOWN from pedigree
        analysis and specified with the parameter \eqn{w=\beta*c_o}
        for each recombining region (in the file \code{regfile}, see
        \code{\link{info_region}}). \cr
        \tab \tab \tab \tab This also specifies \eqn{\mu}, an independent estimate of the
        genomic generational mutation rate per bp.\cr 
        \tab \tab \tab \tab In this case, \code{\link{simulate_data}}
        calculates \eqn{\rho_r=w*(Z-1)\theta_1/\mu} for the genomic
        region \eqn{r}. \cr 
        \code{"}\emph{\code{rho}} \tab \code{-1} \tab \eqn{\lambda}\code{"} \tab \tab 
        : This specifies that the intra-region recombination rate is
        UNKNOWN and the ratio of recombination over mutation
        rate for the genomic region \eqn{r}, \eqn{\alpha=c_r/\mu}, is drawn
        from an exponential distribution with mean \eqn{1/}\eqn{\lambda}. \cr 
        \tab \tab \tab \tab In this case, \code{\link{simulate_data}}
        calculates \eqn{\rho_r=w*\alpha*(Z-1)*theta_1} for the genomic
        region \eqn{r}. \cr   
        \code{"}\emph{\code{rho}} \tab \code{-2} \tab \eqn{\nu} \tab \eqn{\sigma}\code{"} \tab
        : This specifies that the intra-region recombination rate is
        UNKNOWN and the ratio of recombination over mutation
        rate for the genomic region \eqn{r}, \eqn{\alpha=c_r/\mu}, is drawn
        from an normal distribution with mean
        \eqn{\nu} and standard deviation \eqn{\sigma}. \cr
        \tab \tab \tab \tab In this case, \code{\link{simulate_data}}
        calculates \eqn{\rho_r=w*\alpha*(Z-1)*theta_1} for the genomic
        region \eqn{r}.
      }
      Where:
      \tabular{ll}{
        \eqn{\rho=4N_1*c}\tab
        : The genomic average population intra-region recombination rate per bp per
        generation. \cr
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
        \eqn{\alpha=c_r/\mu} \tab : Drawn from a prior distribution. \cr
        \eqn{c_r} \tab
        : The generational region-specific cross-over rate per bp for the
        genomic region \eqn{r}. \cr   
        \eqn{w} \tab 
        : The recombination scalar for the genomic region considered
        specified in the file \code{regfile} (see
        \code{\link{info_region}}).\cr
        \eqn{\beta} \tab
        : The ratio of the region-specific population recombination
        rate per bp over \eqn{\rho=4N_1*c} for the genomic region
        considered.\cr
        \tab \eqn{\beta=}\code{"1"} (\code{"0.5"} in
        \emph{Drosophila}) for autosomal region, \code{"0.5"} for X- and
        \code{"0"} for Y- and mtDNA-linked region.\cr
        \eqn{Z} \tab
        : The size in bp of the genomic region considered. \cr
        \tab \eqn{z_s} and \eqn{z_e} (such as \eqn{Z=z_e-z_s}) are
        specified in the file \code{regfile} (see \code{\link{info_region}}). 
      }
    }
  }
  \item{\strong{OPTIONAL OTHER PARAMETERS: }}{
    \item{\emph{\code{nregions}}}{
      \item{\emph{\code{nregions}}}{
        Use \code{"}\emph{\code{nregions}}\code{ }\eqn{R}\code{"} to specify
        the number of independent genomic regions to simulate. \cr
        By default, \eqn{R=1}.\cr
        ATTENTION: It is required that \eqn{R>0} when specified with the
        keyword \emph{\code{nregions}}. 
      }
    }
  }
}% end format

\note{
  \item{ATTENTION: It is the user's responsibility to mind the following
    restrictions:}{ 
    -> It is required that \eqn{\theta_1>0} be specified
    with the keyword \emph{\code{theta}}. \cr 
    -> It is required that \eqn{\theta_i>0} for any \eqn{i} \eqn{\in}
    \eqn{{1,2,A}} when specified with the keyword \emph{\code{theta}}. \cr 
    -> It is required that \eqn{M_p>0} be specified in case of an island model. \cr 
    -> It is required that \eqn{0}\eqn{\le}\eqn{T_c}\eqn{\le}\eqn{T_s}
    when specified with the keyword \emph{\code{time}}.\cr 
    -> It is required that \eqn{R>0} when specified with the keyword
    \emph{\code{nregions}}. 
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
  Files in the format of \code{\link{param_sim}} are used as input for
  the function \code{\link{simulate_data}}. \cr 
  The function \code{\link{check_param}} checks the format of the
  files like \code{\link{param_sim}}. \cr 
  \code{\link{info_region}} and \code{\link{info_loci}} are other data
  file examples used as input by \code{\link{simulate_data}}. \cr 
  \code{\link{param_est}} is similar to \code{\link{param_sim}} but used
  by the \R function \code{\link{estimate_IMc}}. \cr
  Lists of definitions of the symbols and parameters mentioned in
  this file are found in \code{\link{Rmspack-package}}.
}

\examples{
### Create the data file.
data(simulation_files) # download the data

## Write the file of information on the parameters in the local directory.
write.table(file="sim", x=param_sim, row.name=FALSE, col.names=FALSE,
quote=FALSE)
# Creates the file "sim" containing (with comments): 
read.csv("sim",header =FALSE,sep="", comment="#")

# Clean up the directory.
unlink("sim")
}

\keyword{datasets}
