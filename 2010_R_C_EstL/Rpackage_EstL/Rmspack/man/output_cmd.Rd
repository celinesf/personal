\name{output_cmd}
\alias{output_cmd}

\title{
  Function to write the parameter values for each
  estimation step and the information required to run
  parallel jobs of the estimation program into text files.
}

\description{
  The function  \code{\link{output_cmd}} prints the parameter values for
  each sets of parameters into text files named
  \code{"datafile-lik-}\eqn{j}\code{"} (\eqn{j} \eqn{\in} \eqn{[0,J)}). \cr
  It also prints into the text file  \code{"datafile-estlikcmd"} the command line with other information required to run
  the estimation program.\cr
  This function is called within the \R function
  \code{\link{estimate_IMc}}. 
}
\usage{
output_cmd(datafile, param)
}

\arguments{
  \item{datafile}{
    The argument \code{datafile} specifies the location and the string
    of characters that initiates the name of the text files output by
    the function \code{\link{output_cmd}}. 
  }
  \item{param}{
    The matrix \eqn{11} \eqn{\times} \eqn{3} of either the values or
    descriptions of the prior distributions for the fixed and  estimated
    parameters, respectively.\cr
    The matrix  \code{param} is lightly modified from the matrix
    \code{$param} output by the function \code{\link{check_param}}: The
    \code{"NA"} values should replaced with \code{"-1"}. In
    \code{param[1:7,3]}, \code{"NA"}  should become \code{"1"} instead. 
    \item{\strong{Specifying variable (or estimated) parameters: }}{
      Any of the demographic parameters listed below can be either fixed
      (i.e., only one value specified after the keyword) or variable (or estimated). \cr 
      If a parameter \eqn{i} \eqn{\in} \eqn{[1:8]} is to be variable (or estimated),
      the information on the uniform prior distribution is specify with
      three values, two floating numbers and an integer as follows: \cr
      \code{param[i,]=c(}\eqn{\Theta_l}\code{, }\eqn{\Theta_u}\code{, }\eqn{\Theta_n}\code{)}. 
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
      } ATTENTION:
      \eqn{0\le}\eqn{\Theta_l}\eqn{\le}\eqn{\Theta_u} and
      \eqn{\Theta_n>0}. \cr 
      If only one value is specified (i.e., \eqn{\Theta_l>0} and
      \eqn{\Theta_u} not specified or \eqn{\Theta_u=\Theta_l}), the
      parameter is considered fixed to \eqn{\Theta_l}. 
    }
    The twelve rows correspond to the twelve following parameters and keywords:
    \item{\strong{REQUIRED PARAMETERS: }}{
      \item{\code{param[1,]~}\emph{\code{theta_1}}}{
        Either the parameter value of the population mutation rate per bp per
        generation for population 1 
        (i.e. \code{param[1,1]=}\eqn{\theta_1>0}) or the prior distribution
        for \eqn{theta_1} (see above for the restrictions) should be specified:\cr
        \code{param[1,1]}\eqn{\theta_1} or
        \code{param[1,]=c(}\eqn{\theta_1l}\code{, }\eqn{\theta_1u}\code{, }\eqn{\theta_1n}\code{)}. \cr
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
        } 
        The information on the parameter \eqn{\theta_1} can be
        specified with the keyword \emph{\code{theta_1}} in the file
        \code{paramfile} (see \code{\link{param_est}}).\cr
        NOTE: If \eqn{\theta_1l=0}, \code{\link{output_cmd}}
        will generate \eqn{\theta_1n} possible values of \eqn{\theta_1} \eqn{\in}
        \eqn{[(\theta_1u-\theta_1l)/\theta_1n,\theta_1u]}. 
      }% end theta
      \item{\code{param[2,]~}\emph{\code{M_present}}}{
        \code{param[2,1]=}\eqn{M_p} or
        \code{param[2,]=c(}\eqn{M_pl}\code{, }\eqn{M_pu}\code{, }\eqn{M_pn}\code{)}\cr
        specifies \eqn{M_p=4N_1*m_p}, the number of migrants exchanged
        each generation by the TWO populations at present. \cr 
        \eqn{M_p} has different definitions depending of the model considered: \cr 
        * \eqn{M_p} is the symmetrical rate of gene flow between
        TWO populations in an island model (i.e., \eqn{0<t<\inf}) as specified
        with: \cr 
        -- At least one of the independent genomic regions has
        \eqn{n_1>0} and \eqn{n_2>0} as specified in the file
        \code{regfile} (see \code{\link{info_region}}). \cr
        -- Either \code{param[5,1]=}\eqn{T_s=0} or \code{param[5,1]=0}. \cr
        \eqn{M_p>0} when fixed or
        \eqn{0\le}\eqn{M_pl}\eqn{\le}\eqn{M_pu} and \eqn{M_pn>0} when
        variable (or estimated) be specified in \code{param[2,]} in case of an island
        model. \cr 
        * \eqn{M_p} is the constant symmetrical rate of gene flow
        since the split until present (i.e., \eqn{0<t<T_s}) if
        there is a population split (as specified with \eqn{T_s>0} when fixed or
        \eqn{0\le}\eqn{T_sl}\eqn{\le}\eqn{T_su} and \eqn{T_sn>0} when
        variable (or estimated) in \code{param[5,]}). \cr 
        * \eqn{M_p} is the constant symmetrical rate of gene flow since
        the time of gene flow rate change until present (i.e., \eqn{0<t<T_c})
        if a time at which the gene flow rate changed is specified with: \cr 
        -- \eqn{0<T_c<T_s} (as specified with \eqn{0<\epsilon<1} when fixed or
        \eqn{0\le}\eqn{\epsilon_l}\eqn{\le}\eqn{\epsilon_u}\eqn{\le1} and \eqn{\epsilon_n>0}
        when variable (or estimated) in \code{param[6,]}). \cr 
        -- \eqn{M_c \not=M_p} (as specified with \eqn{0\le}\eqn{M_c
          \not=M_p} when fixed in \code{param[2,1]} and \code{param[7,1]}).\cr 
        \tabular{l}{
          ---> Where, \eqn{m_p} 
          : The generational fraction of migrant individuals at
          present. 
        }
        The information on the parameter \eqn{M_p} can be
        specified with the keyword \emph{\code{M_present}} in the file
        \code{paramfile} (see \code{\link{param_est}}).
      }% end M_p
    }% end req arg
    \item{\strong{OPTIONAL DEMOGRAPHIC PARAMETERS: }}{
      All the parameters listed below are optional and can be either
      fixed or variable (or estimated). 
      \item{\code{param[3,]~}\emph{\code{theta_2}}}{
        \code{param[3,1]=}\eqn{\theta_2} or
        \code{param[3,]=c(}\eqn{\theta_2l}\code{, }\eqn{\theta_2u}\code{, }\eqn{\theta_2n}\code{)}\cr
        specifies \eqn{\theta_2=4N_2*\mu}, the population mutation rate
        per bp per generation for population 2. \cr
        \tabular{l}{
          ---> Where, \eqn{N_2} 
          : The effective population size in population 2. 
        }
        \eqn{theta_2>0} when specified and fixed in \code{param[3,1]}. \cr 
        The information on the parameter \eqn{\theta_2} can be
        specified with the keyword \emph{\code{theta_2}} in the file
        \code{paramfile} (see \code{\link{param_est}}).\cr
        NOTE: If \eqn{\theta_2l=0}, \code{\link{output_cmd}}
        will sample a value of \eqn{\theta_2} among the \eqn{\theta_2n}
        possible values \eqn{\in}
        \eqn{[(\theta_2u-\theta_2l)/\theta_2n,\theta_2u]}.\cr
        If \eqn{\theta_2} is not specified the \eqn{\theta_2=\theta_1}.
      }% end theta2 
      \item{\code{param[4,]~}\emph{\code{theta_A}}}{
        \code{param[4,1]=}\eqn{\theta_A} or
        \code{param[4,]=c(}\eqn{\theta_Al}\code{, }\eqn{\theta_Au}\code{, }\eqn{\theta_An}\code{)}\cr
        specifies \eqn{\theta_A=4N_A*\mu}, the ancestral population
        mutation rate per bp per generation. 
        \tabular{l}{
          ---> Where, \eqn{N_A}
          : The ancestral effective population size. 
        }
        \eqn{theta_A>0} when specified and fixed in \code{param[4,1]}. \cr
        The information on the parameter \eqn{\theta_A} can be
        specified with the keyword \emph{\code{theta_A}} in the file
        \code{paramfile} (see \code{\link{param_est}}).\cr
        NOTE: If \eqn{\theta_Al=0}, \code{\link{output_cmd}}
        will sample a value of \eqn{\theta_A} among the \eqn{\theta_An}
        possible values \eqn{\in} \eqn{[(\theta_Au-\theta_Al)/\theta_An,\theta_Au]}.\cr
        If \eqn{\theta_A} is not specified the \eqn{\theta_A=\theta_1}.
      }% end thetaA
      \item{\code{param[5,]~}\emph{\code{T_split}}}{
        \code{param[5,1]=}\eqn{T_s} or
        \code{param[5,]=c(}\eqn{T_sl}\code{, }\eqn{T_su}\code{, }\eqn{T_sn}\code{)}\cr
        specifies \eqn{T_s}, the split time in unit of \eqn{4N_1}
        generations between the TWO populations. \cr
        The information on the parameter \eqn{T_s} can be specified with the keyword
        \emph{\code{T_split}} in the file \code{paramfile} (see
        \code{\link{param_est}}).
      }% end T_split
      \item{\code{param[6,]~}\emph{\code{T_change}}}{
        \code{param[6,1]=}\eqn{\epsilon} or
        \code{param[6,]=c(}\eqn{\epsilon_l}\code{, }\eqn{\epsilon_u}\code{, }\eqn{\epsilon_n}\code{)}\cr
        specifies the ratio \eqn{\epsilon=T_c/T_s}.
        \tabular{l}{
          ---> Where, \eqn{T_c}
          : The time at which the rate of gene flow changed between the
          two populations in unit of \eqn{4N_1} generations.
        }
        The information on the parameter \eqn{\epsilon} can be
        specified with the keyword \emph{\code{T_change}} in the file
        \code{paramfile} (see \code{\link{param_est}}).\cr
        \eqn{0}\eqn{\le}\eqn{T_c}\eqn{\le}\eqn{T_s} when specified and fixed
        in \code{param[5,1]} and \code{param[6,1]}). \cr
        NOTE: If \eqn{\epsilon_u=1}, \code{\link{output_cmd}}
        will sample a value of \eqn{\epsilon} among the \eqn{\epsilon_n}
        possible values \eqn{\in} \eqn{[\epsilon_l),\epsilon_u-(\epsilon_u-\epsilon_l)/\epsilon_n,]}.
      }% end T_change
      \item{\code{param[7,]~}\emph{\code{M_change}}}{
        \code{param[7,1]=}\eqn{M_c} or
        \code{param[7,]=c(}\eqn{M_cl}\code{, }\eqn{M_cu}\code{, }\eqn{M_cn}\code{)}\cr
        specifies \eqn{M_c=4N_1*m_c}, the number of migrants exchanged
        each generation by the TWO populations since the split until the time
        of gene flow rate change (i.e., \eqn{T_c<t<T_s}).
        \tabular{l}{
          ---> Where, \eqn{m_c} 
          : The generational fraction of migrant individuals between
          \eqn{T_c<t<T_s}.
        }
        The information on the parameter \eqn{M_c} can be
        specified with the keyword \emph{\code{M_change}} in the file
        \code{paramfile} (see \code{\link{param_est}}).
      }% end M_c
    }% end Optional demo
    \item{\strong{OPTIONAL NUISANCE PARAMETER: }}{
      \item{\code{param[8,]~}\emph{\code{rho}}}{
        This specifies the parameter on the intra-region recombination
        rate.\cr
        The region-specific recombination rate per generation is calculated
        with \eqn{\rho_r=\beta*4N_1c*(Z-1)} for the genomic region \eqn{r}
        \eqn{\in} \eqn{[1,R]} (as specified in \code{param[9,1]}). \cr 
        \code{param[8,]} can have the five following forms:
        \tabular{lllll}{
          \code{param[8,1]=}\tab \eqn{\rho} \tab \tab \tab
          : This specifies that the genomic average population
          intra-region recombination rate per bp per generation is fixed to the
          value \eqn{\rho=4N_1*c}. \cr
          \tab \tab \tab \tab In this case, \eqn{\rho_r=\rho*w*(Z-1)}
          for the genomic region \eqn{r} (here
          \eqn{w=\beta}). \eqn{\rho_r}  is fixed across estimation steps in the
          function \code{\link{output_cmd}}. \cr  
          \code{param[8,1]=} \tab \code{1} \tab \tab \tab
          : This specifies that an estimate of the region-specific population
          recombination rate per bp, \eqn{\rho_o=4N_1*c_o}, is KNOWN from
          linkage disequilibrium analysis and specified with the parameter
          \eqn{w=\beta*\rho_o} for each recombining region (in the file
          \code{regfile}, see \code{\link{info_region}}).\cr
          \tab \tab \tab \tab In this case, \eqn{\rho_r=w*(Z-1)} for the
          genomic region \eqn{r}. \eqn{\rho_r}  is fixed across estimation steps
          in the function \code{\link{output_cmd}}. \cr 
          \code{param[8,]=c(}\tab \code{2}\code{,}\tab \eqn{\mu}\code{,}\tab \code{NA)}\tab
          : This specifies that an estimate of the region-specific
          recombination rate per bp, \eqn{c_o}, is KNOWN from pedigree
          analysis and specified with the parameter \eqn{w=\beta*c_o}
          for each recombining region (in the file \code{regfile}, see
          \code{\link{info_region}}). \cr
          \tab \tab \tab \tab This also specifies \eqn{\mu}, an
          independent estimate of the genomic generational mutation
          rate per bp.\cr 
          \tab \tab \tab \tab In this case,
          \eqn{\rho_r=w*(Z-1)\theta_1/\mu}  for the genomic region
          \eqn{r}. \eqn{\rho_r}  varies across estimation steps if \eqn{theta_1}
          is estimated in the function \code{\link{output_cmd}}. \cr 
          \code{param[8,]=c(}\tab \code{-}\code{1}\code{,}\tab \eqn{1/}\eqn{\lambda}\code{,}\tab \code{NA)} \tab
          : This specifies that the intra-region recombination rate is
          UNKNOWN and the ratio of recombination over mutation
          rate for the genomic region \eqn{r}, \eqn{\alpha=c_r/\mu}, is drawn
          from an exponential distribution with mean \eqn{1/}\eqn{\lambda}. \cr 
          \tab \tab \tab \tab In this case,
          \eqn{\rho_r=w*\alpha*(Z-1)*theta_1} for the genomic region
          \eqn{r}. \eqn{\rho_r}  varies across estimation steps if \eqn{theta_1}
          is estimated in the function \code{\link{output_cmd}}. \cr
          \code{param[8,]=c(}\tab \code{-}\code{2}\code{,}\tab \eqn{\nu}\code{,}\tab \eqn{\sigma}\code{)} \tab
          : This specifies that the intra-region recombination rate is
          UNKNOWN and the ratio of recombination over mutation
          rate for the genomic region \eqn{r}, \eqn{\alpha=c_r/\mu}, is drawn
          from an normal distribution with mean
          \eqn{\nu} and standard deviation \eqn{\sigma}. \cr
          \tab \tab \tab \tab In this case,
          \eqn{\rho_r=w*\alpha*(Z-1)*theta_1} for the genomic region
          \eqn{r}. \eqn{\rho_r}  varies across estimation steps if \eqn{theta_1}
          is estimated in the function \code{\link{output_cmd}}.
        } Where:
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
          \code{"0"} for Y- and mtDNA-linked region. \cr
          \eqn{Z} \tab
          : The size in bp of the genomic region considered. \cr
          \tab \eqn{z_s} and \eqn{z_e} (such as \eqn{Z=z_e-z_s}) are
          specified in the file \code{regfile} (see
          \code{\link{info_region}}). 
        }
        The information on the parameter \eqn{\rho} can be
        specified with the keyword \emph{\code{rho}} in the file
        \code{paramfile} (see  \code{\link{param_est}}).
      }% end rho
    }% end nuisance
    \item{\strong{OPTIONAL OTHER PARAMETERS: }}{
      \item{\code{param[9,1]~}\emph{\code{nregions}}}{
        \code{param[9,1]}\eqn{R} is the number of independent
        genomic regions considered.\cr
        \eqn{R} can be specified with the keyword
        \emph{\code{nregions}} in the file \code{paramfile} (see
        \code{\link{param_est}}). \cr 
        \code{param[9,1]=}\eqn{R>0} when specified. \cr 
      }
    }
    \item{\strong{OPTIONAL PARAMETERS SPECIFIC TO \code{\link{estimate_IMc}}: }}{
      \item{\code{param[10,1]~}\emph{\code{howmany}}}{
        \code{param[10,1]=}\eqn{H} specifies \eqn{H}, the number of
        data sets to simulate per set of parameters (i.e., per grid point) to
        estimate the likelihood of the data given the set of parameters of the
        extension of the isolation-with-migration model.\cr 
        By default, \code{param[10,1]=}\eqn{1000}.\cr
        \code{param[10,1]=}\eqn{H>0} when specified. \cr
        The information on \eqn{J} can be specified with the keyword
        \emph{\code{howmany}} in the file \code{paramfile} (see
        \code{\link{param_est}}). \cr 
      }% end howmany
      \item{\code{param[11,1]~}\emph{\code{parallel}}}{
        \code{param[11,1]=}\eqn{J} specifies \eqn{J}, the number of jobs
        to run in parallel to perform the estimation of the posterior
        distribution of the parameters of the extension of the
        isolation-with-migration model.\cr 
        By default, \code{param[11,1]=}\eqn{1}, i.e., no parallelization.\cr
        \code{param[11,1]=}\eqn{J>0} when specified.\cr
        The information on \eqn{J} can be specified with the keyword
        \emph{\code{parallel}} in the file \code{paramfile} (see
        \code{\link{param_est}}). \cr 
      }% end parallel
    }% end extra 
  }
}

\value{
  \item{\code{"datafile-lik-}\eqn{j}\code{"}}{
    The function  \code{\link{output_cmd}} writes the 
    parameter values for each of the estimation steps considered divided in
    \eqn{J} text files named \code{"datafile-lik-}\eqn{j}\code{"}, \eqn{j}
    \eqn{\in} \eqn{[0,J-1]} (\eqn{J} is provided in \code{param[11,1]}).\cr
    Each file contains \eqn{\Pi(}\code{param[i,3]}\eqn{)/J} lines, \eqn{i}
    \eqn{\in} \eqn{[1,7]}, but the file \eqn{J-1} may have less lines
    since some sets of parameters are ignored (see below).\cr
    Each line is as follows:
    \code{"}\eqn{theta_1}\code{ }\eqn{theta_2}\code{ }\eqn{theta_A}\code{ }\eqn{T_s}\code{ }\eqn{T_c}\code{ }\eqn{M_p}\code{ }\eqn{M_c}\code{"}\cr 
    Each line contains a unic set of parameter values.\cr All the combination
    of values for the estimate parameters are recorded: i.e. the
    combination of the \code{param[i,3]} possible values \eqn{\in}
    \code{[param[i,1],param[i,2]]}, \eqn{i}
    \eqn{\in} \eqn{[1,7]}, with the following expections:\cr
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
    e.g. of the first two lines of \code{"reg-lik-0"}:\cr
    \code{0.0001 0.0001 0.00015 0.222222 0 0 0}\cr
    \code{0.0001 0.0001 0.00015 0.222222 0.8 0 2.22222}
  }
  \item{\code{"datafile-estlikcmd"}}{
    The function  \code{\link{output_cmd}} writes the information required
    for the estimation program to be run \eqn{J} times on the \eqn{J} text
    files named \code{"datafile-lik-}\eqn{j}\code{"}, into the text file
    name \code{"datafile-estlikcmd"}.\cr
    This file contains \eqn{J} lines, as follows:
    \code{"./estlik -a }\eqn{h}\code{ -r }\eqn{rho}\code{ -p datafile-lik-}\eqn{j}\code{"}.\cr
    e.g. for \code{"reg-estlikecmd"}:\cr
    \code{./estlik -a 1 -r .005 -p reg-lik-0}\cr
  }
}

\author{
  Celine Becquet - \email{celine.becquet@gmail.com}. 
}

\seealso{
  The function \code{\link{estimate_IMc}} calls the function
  \code{\link{output_cmd}}.\cr
  The function \code{\link{output_cmd}} uses as input a modified data frame
  output by \code{\link{check_param}}. \cr
}

\examples{
### Write the file of information on the parameters in the local directory.
data(simulation_files) # download the data

write.table(file="est", x=param_est, row.name=FALSE, col.names=FALSE, quote=FALSE)
# Creates the file "est" containing (with comments):
read.csv("est",header =FALSE,sep="", comment="#")

### Create the input for the function.
listparam=c("theta_1", "M_present", "theta_2", "theta_A", "T_split", "T_change", "M_change", "rho", "nregions", "howmany", "parallel")
param=scan("est", comment.char="#", what=c("c", "n", "n", "n"), allowEscapes=FALSE, fill=TRUE, sep="\n", strip.white=TRUE, quiet=TRUE) 
vparam=order_param(param, listparam) 
vparam=check_param(param=vparam, paramfile="est", listparam=listparam)$param
vparam[is.na(vparam[])]=-1
vparam[which(vparam[1:7,3]==-1),3]=1

## case for a sigle estlike run.
output_cmd(datafile="reg", param=vparam)
# Creates the file "reg-0"  containing 181900 lines (the last five shown):
read.csv("reg-0",header =FALSE,sep="", comment="#", skip=181895)
read.csv("reg-estlikcmd",header =FALSE,sep="")

## Case for 2 parellel estlike runs.
vparam[11,1]=2
output_cmd(datafile="reg", param=vparam)
# Creates the files "reg-0" and "reg-1" containing 100001 and 81899 lines (the last five shown):
read.csv("reg-0",header =FALSE,sep="", comment="#", skip=99996)
read.csv("reg-1",header =FALSE,sep="", comment="#", skip=81894)
read.csv("reg-estlikcmd",header =FALSE,sep="")

# Clean up the directory.
unlink(c("est","reg-0","reg-1"))
}

\keyword{IO}
\keyword{print}

