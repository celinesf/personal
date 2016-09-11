

#' Function to check the format of the file of parameter used to simulate data
#' under or estimate parameters of extensions of the Isolation-with-migration
#' model.
#' The function \code{\link{check_param}} checks that the information on the
#' values and descriptions of the prior distributions for the fixed and
#' variable (or estimated) parameters, respectively, specified in the file
#' \code{paramfile} (see \code{\link{param_est}} and \code{\link{param_sim}})
#' are in the correct format to simulate data under (or to estimate parameters
#' of) extensions of the Isolation-with-migration model with the function
#' \code{\link{simulate_data}} (or \code{\link{estimate_IMc}}). \cr This
#' function is called within the functions \code{\link{simulate_data}} and
#' \code{\link{estimate_IMc}}.
#' 
#' 
#' @param param The matrix \eqn{9} (or \eqn{11}) \eqn{\times} \eqn{3} of either
#'   the values or descriptions of the prior distributions for the fixed and
#'   variable (or estimated) parameters, respectively, required to simulate
#'   data under (or to estimate parameters of) extensions of the
#'   Isolation-with-migration model with the functions
#'   \code{\link{simulate_data}} and \code{\link{estimate_IMc}}.
#'   \item\strong{Specifying variable (or estimated) parameters: } Any of the
#'   demographic parameters listed below can be either fixed (i.e., only one
#'   value specified after the keyword) or variable (or estimated). \cr If a
#'   parameter \eqn{i} \eqn{\in} \eqn{[1:7]} is to be variable (or estimated),
#'   the information on the uniform prior distribution is specify with three
#'   values, two floating numbers and an integer as follows: \cr
#'   \code{param[i,]=c(}\eqn{\Theta_l}\code{, }\eqn{\Theta_u}\code{,
#'   }\eqn{\Theta_n}\code{)}.  \tabular{ll}{ \eqn{\Theta_l}\tab : The lower
#'   limit on the prior distribution range. \eqn{\Theta_l}\eqn{\ge}\eqn{0}. \cr
#'   \eqn{\Theta_u} \tab : The upper limit on the prior distribution range.
#'   \eqn{\Theta_l}\eqn{\le}\eqn{\Theta_u}. \cr \eqn{\Theta_n} \tab : The
#'   number of values to consider along the prior distribution range.
#'   \eqn{\Theta_n>0} \cr } ATTENTION: \code{\link{check_param}} requires that
#'   \eqn{0\le}\eqn{\Theta_l}\eqn{\le}\eqn{\Theta_u} and \eqn{\Theta_n>0}. \cr
#'   If only one value is specified (i.e., \eqn{\Theta_l>0} and \eqn{\Theta_u}
#'   not specified or \eqn{\Theta_u=\Theta_l}), the parameter is considered
#'   fixed to \eqn{\Theta_l}.
#' 
#' The twelve rows correspond to the twelve following parameters and keywords:
#'   \item\strong{REQUIRED PARAMETERS: }
#'   \item\code{param[1,]~}\emph{\code{theta_1}} ATTENTION:
#'   \code{\link{check_param}} requires that either the parameter value of the
#'   population mutation rate per bp per generation for population 1 (i.e.
#'   \code{param[1,1]=}\eqn{\theta_1>0}) or the prior distribution for
#'   \eqn{theta_1} (see above for the restrictions) be specified:\cr
#'   \code{param[1,1]}\eqn{\theta_1} or
#'   \code{param[1,]=c(}\eqn{\theta_1l}\code{, }\eqn{\theta_1u}\code{,
#'   }\eqn{\theta_1n}\code{)}. \cr \tabular{ll}{ ---> Where: \tab \cr
#'   \eqn{\theta_1=4N_1*\mu} \tab : The population mutation rate per bp per
#'   generation for population 1 (required). \cr \eqn{N_1} \tab : The effective
#'   population size in population 1 (the reference population so by default
#'   and unless specified, \eqn{N_1=N_2=N_A}). \cr \eqn{\mu} \tab : The genomic
#'   generational mutation rate per bp.  } The information on the parameter
#'   \eqn{\theta_1} can be specified with the keyword \emph{\code{theta_1}} in
#'   the file \code{paramfile} (see \code{\link{param_sim}} and
#'   \code{\link{param_est}}).  % end theta
#'   \item\code{param[2,]~}\emph{\code{M_present}} \code{param[2,1]=}\eqn{M_p}
#'   or \code{param[2,]=c(}\eqn{M_pl}\code{, }\eqn{M_pu}\code{,
#'   }\eqn{M_pn}\code{)}\cr specifies \eqn{M_p=4N_1*m_p}, the number of
#'   migrants exchanged each generation by the TWO populations at present. \cr
#'   \eqn{M_p} has different definitions depending of the model considered: \cr
#'   * \eqn{M_p} is the symmetrical rate of gene flow between TWO populations
#'   in an island model (i.e., \eqn{0<t<\inf}) as specified with: \cr -- At
#'   least one of the independent genomic regions has \eqn{n_1>0} and
#'   \eqn{n_2>0} as specified in the file \code{regfile} (see
#'   \code{\link{info_region}}). \cr -- Either \code{param[5,1]=}\eqn{T_s=0} or
#'   \code{param[5,1]=NA}. \cr ATTENTION: \code{\link{check_param}} requires
#'   that \eqn{M_p>0} when fixed or \eqn{0\le}\eqn{M_pl}\eqn{\le}\eqn{M_pu} and
#'   \eqn{M_pn>0} when variable (or estimated) be specified in \code{param[2,]}
#'   in case of an island model. \cr * \eqn{M_p} is the constant symmetrical
#'   rate of gene flow since the split until present (i.e., \eqn{0<t<T_s}) if
#'   there is a population split (as specified with \eqn{T_s>0} when fixed or
#'   \eqn{0\le}\eqn{T_sl}\eqn{\le}\eqn{T_su} and \eqn{T_sn>0} when variable (or
#'   estimated) in \code{param[5,]}). \cr * \eqn{M_p} is the constant
#'   symmetrical rate of gene flow since the time of gene flow rate change
#'   until present (i.e., \eqn{0<t<T_c}) if a time at which the gene flow rate
#'   changed is specified with: \cr -- \eqn{0<T_c<T_s} (as specified with
#'   \eqn{0<\epsilon<1} when fixed or
#'   \eqn{0\le}\eqn{\epsilon_l}\eqn{\le}\eqn{\epsilon_u}\eqn{\le1} and
#'   \eqn{\epsilon_n>0} when variable (or estimated) in \code{param[6,]}). \cr
#'   -- \eqn{M_c \not=M_p} (as specified with \eqn{0\le}\eqn{M_c \not=M_p} when
#'   fixed in \code{param[2,1]} and \code{param[7,1]}).\cr \tabular{l}{ --->
#'   Where, \eqn{m_p} : The generational fraction of migrant individuals at
#'   present.  } The information on the parameter \eqn{M_p} can be specified
#'   with the keyword \emph{\code{M_present}} in the file \code{paramfile} (see
#'   \code{\link{param_sim}} and \code{\link{param_est}}).  % end M_p % end req
#'   arg \item\strong{OPTIONAL DEMOGRAPHIC PARAMETERS: } All the parameters
#'   listed below are optional and can be either fixed or variable (or
#'   estimated).  \item\code{param[3,]~}\emph{\code{theta_2}}
#'   \code{param[3,1]=}\eqn{\theta_2} or
#'   \code{param[3,]=c(}\eqn{\theta_2l}\code{, }\eqn{\theta_2u}\code{,
#'   }\eqn{\theta_2n}\code{)}\cr specifies \eqn{\theta_2=4N_2*\mu}, the
#'   population mutation rate per bp per generation for population 2. \cr
#'   \tabular{l}{ ---> Where, \eqn{N_2} : The effective population size in
#'   population 2.  } ATTENTION: \code{\link{check_param}} requires that
#'   \eqn{theta_2>0} when specified and fixed in \code{param[3,1]}. \cr The
#'   information on the parameter \eqn{\theta_2} can be specified with the
#'   keyword \emph{\code{theta_2}} in the file \code{paramfile} (see
#'   \code{\link{param_sim}} and \code{\link{param_est}}).  % end theta2
#'   \item\code{param[4,]~}\emph{\code{theta_A}}
#'   \code{param[4,1]=}\eqn{\theta_A} or
#'   \code{param[4,]=c(}\eqn{\theta_Al}\code{, }\eqn{\theta_Au}\code{,
#'   }\eqn{\theta_An}\code{)}\cr specifies \eqn{\theta_A=4N_A*\mu}, the
#'   ancestral population mutation rate per bp per generation.  \tabular{l}{
#'   ---> Where, \eqn{N_A} : The ancestral effective population size.  }
#'   ATTENTION: \code{\link{check_param}} requires that \eqn{theta_A>0} when
#'   specified and fixed in \code{param[4,1]}. \cr The information on the
#'   parameter \eqn{\theta_A} can be specified with the keyword
#'   \emph{\code{theta_A}} in the file \code{paramfile} (see
#'   \code{\link{param_sim}} and \code{\link{param_est}}).  % end thetaA
#'   \item\code{param[5,]~}\emph{\code{T_split}} \code{param[5,1]=}\eqn{T_s} or
#'   \code{param[5,]=c(}\eqn{T_sl}\code{, }\eqn{T_su}\code{,
#'   }\eqn{T_sn}\code{)}\cr specifies \eqn{T_s}, the split time in unit of
#'   \eqn{4N_1} generations between the TWO populations. \cr The information on
#'   the parameter \eqn{T_s} can be specified with the keyword
#'   \emph{\code{T_split}} in the file \code{paramfile} (see
#'   \code{\link{param_sim}} and \code{\link{param_est}}).  % end T_split
#'   \item\code{param[6,]~}\emph{\code{T_change}}
#'   \code{param[6,1]=}\eqn{\epsilon} or
#'   \code{param[6,]=c(}\eqn{\epsilon_l}\code{, }\eqn{\epsilon_u}\code{,
#'   }\eqn{\epsilon_n}\code{)}\cr specifies the ratio \eqn{\epsilon=T_c/T_s}.
#'   \tabular{l}{ ---> Where, \eqn{T_c} : The time at which the rate of gene
#'   flow changed between the two populations in unit of \eqn{4N_1}
#'   generations.  } The information on the parameter \eqn{\epsilon} can be
#'   specified with the keyword \emph{\code{T_change}} in the file
#'   \code{paramfile} (see \code{\link{param_sim}} and
#'   \code{\link{param_est}}).\cr ATTENTION: \code{\link{check_param}} requires
#'   that \eqn{0}\eqn{\le}\eqn{T_c}\eqn{\le}\eqn{T_s} when specified and fixed
#'   in \code{param[5,1]} and \code{param[6,1]}).  % end T_change
#'   \item\code{param[7,]~}\emph{\code{M_change}} \code{param[7,1]=}\eqn{M_c}
#'   or \code{param[7,]=c(}\eqn{M_cl}\code{, }\eqn{M_cu}\code{,
#'   }\eqn{M_cn}\code{)}\cr specifies \eqn{M_c=4N_1*m_c}, the number of
#'   migrants exchanged each generation by the TWO populations since the split
#'   until the time of gene flow rate change (i.e., \eqn{T_c<t<T_s}).
#'   \tabular{l}{ ---> Where, \eqn{m_c} : The generational fraction of migrant
#'   individuals between \eqn{T_c<t<T_s}.  } The information on the parameter
#'   \eqn{M_c} can be specified with the keyword \emph{\code{M_change}} in the
#'   file \code{paramfile} (see \code{\link{param_sim}} and
#'   \code{\link{param_est}}).  % end M_c % end Optional demo
#'   \item\strong{OPTIONAL NUISANCE PARAMETER: }
#'   \item\code{param[8,]~}\emph{\code{rho}} This specifies the parameter on
#'   the intra-region recombination rate.\cr The region-specific recombination
#'   rate per generation is calculated with \eqn{\rho_r=\beta*4N_1c*(Z-1)} for
#'   the genomic region \eqn{r} \eqn{\in} \eqn{[1,R]} (as specified in
#'   \code{param[9,1]}). \cr \code{param[8,]} can have the five following
#'   forms: \tabular{lllll}{ \code{param[8,1]=}\tab \eqn{\rho} \tab \tab \tab :
#'   This specifies that the genomic average population intra-region
#'   recombination rate per bp per generation is fixed to the value
#'   \eqn{\rho=4N_1*c}. \cr \tab \tab \tab \tab In this case,
#'   \eqn{\rho_r=\rho*w*(Z-1)} for the genomic region \eqn{r} (here
#'   \eqn{w=\beta}). \eqn{\rho_r} is fixed across estimation steps in the
#'   function \code{\link{estimate_IMc}}. \cr \code{param[8,1]=} \tab \code{1}
#'   \tab \tab \tab : This specifies that an estimate of the region-specific
#'   population recombination rate per bp, \eqn{\rho_o=4N_1*c_o}, is KNOWN from
#'   linkage disequilibrium analysis and specified with the parameter
#'   \eqn{w=\beta*\rho_o} for each recombining region (in the file
#'   \code{regfile}, see \code{\link{info_region}}).\cr \tab \tab \tab \tab In
#'   this case, \eqn{\rho_r=w*(Z-1)} for the genomic region \eqn{r}.
#'   \eqn{\rho_r} is fixed across estimation steps in the function
#'   \code{\link{estimate_IMc}}. \cr \code{param[8,]=c(}\tab
#'   \code{2}\code{,}\tab \eqn{\mu}\code{,}\tab \code{NA)}\tab : This specifies
#'   that an estimate of the region-specific recombination rate per bp,
#'   \eqn{c_o}, is KNOWN from pedigree analysis and specified with the
#'   parameter \eqn{w=\beta*c_o} for each recombining region (in the file
#'   \code{regfile}, see \code{\link{info_region}}). \cr \tab \tab \tab \tab
#'   This also specifies \eqn{\mu}, an independent estimate of the genomic
#'   generational mutation rate per bp.\cr \tab \tab \tab \tab In this case,
#'   \eqn{\rho_r=w*(Z-1)\theta_1/\mu} for the genomic region \eqn{r}.
#'   \eqn{\rho_r} varies across estimation steps if \eqn{theta_1} is estimated
#'   in the function \code{\link{estimate_IMc}}. \cr \code{param[8,]=c(}\tab
#'   \code{-}\code{1}\code{,}\tab \eqn{1/}\eqn{\lambda}\code{,}\tab \code{NA)}
#'   \tab : This specifies that the intra-region recombination rate is UNKNOWN
#'   and the ratio of recombination over mutation rate for the genomic region
#'   \eqn{r}, \eqn{\alpha=c_r/\mu}, is drawn from an exponential distribution
#'   with mean \eqn{1/}\eqn{\lambda}. \cr \tab \tab \tab \tab In this case,
#'   \eqn{\rho_r=w*\alpha*(Z-1)*theta_1} for the genomic region \eqn{r}.
#'   \eqn{\rho_r} varies across estimation steps if \eqn{theta_1} is estimated
#'   in the function \code{\link{estimate_IMc}}. \cr \code{param[8,]=c(}\tab
#'   \code{-}\code{2}\code{,}\tab \eqn{\nu}\code{,}\tab \eqn{\sigma}\code{)}
#'   \tab : This specifies that the intra-region recombination rate is UNKNOWN
#'   and the ratio of recombination over mutation rate for the genomic region
#'   \eqn{r}, \eqn{\alpha=c_r/\mu}, is drawn from an normal distribution with
#'   mean \eqn{\nu} and standard deviation \eqn{\sigma}. \cr \tab \tab \tab
#'   \tab In this case, \eqn{\rho_r=w*\alpha*(Z-1)*theta_1} for the genomic
#'   region \eqn{r}. \eqn{\rho_r} varies across estimation steps if
#'   \eqn{theta_1} is estimated in the function \code{\link{estimate_IMc}}.  }
#'   Where: \tabular{ll}{ \eqn{\rho=4N_1*c}\tab : The genomic average
#'   population intra-region recombination rate per bp per generation. \cr
#'   \eqn{\rho_r=\beta*(Z}\eqn{-}\eqn{1)*4N_1c}\tab : The region-specific
#'   recombination rate per generation for the genomic region considered. \cr
#'   \eqn{\rho_o=4N_1*c_o} \tab : The estimate of the region-specific
#'   population recombination rate per bp per generation for the genomic region
#'   considered from linkage disequilibrium analysis. \cr \eqn{c} \tab : The
#'   genomic generational cross-over rate per bp. \cr \eqn{c_o} \tab : The
#'   estimate of the region-specific cross-over rate per bp per generation for
#'   the genomic region considered from pedigree analysis. \cr
#'   \eqn{\alpha=c_r/\mu} \tab : Drawn from a prior distribution. \cr \eqn{c_r}
#'   \tab : The generational region-specific cross-over rate per bp for the
#'   genomic region \eqn{r}. \cr \eqn{w} \tab : The recombination scalar for
#'   the genomic region considered specified in the file \code{regfile} (see
#'   \code{\link{info_region}}).\cr \eqn{\beta} \tab : The ratio of the
#'   region-specific population recombination rate per bp over
#'   \eqn{\rho=4N_1*c} for the genomic region considered.\cr \tab
#'   \eqn{\beta=}\code{"1"} (\code{"0.5"} in \emph{Drosophila}) for autosomal
#'   region, \code{"0.5"} for X- and \code{"0"} for Y- and mtDNA-linked region.
#'   \cr \eqn{Z} \tab : The size in bp of the genomic region considered. \cr
#'   \tab \eqn{z_s} and \eqn{z_e} (such as \eqn{Z=z_e-z_s}) are specified in
#'   the file \code{regfile} (see \code{\link{info_region}}).  } The
#'   information on the parameter \eqn{\rho} can be specified with the keyword
#'   \emph{\code{rho}} in the file \code{paramfile} (see
#'   \code{\link{param_sim}} and \code{\link{param_est}}).  % end rho % end
#'   nuisance \item\strong{OPTIONAL OTHER PARAMETERS: }
#'   \item\code{param[9,1]~}\emph{\code{nregions}} \code{param[9,1]}\eqn{R} is
#'   the number of independent genomic regions considered.\cr \eqn{R} can be
#'   specified with the keyword \emph{\code{nregions}} in the file
#'   \code{paramfile} (see \code{\link{param_sim}} and
#'   \code{\link{param_est}}). \cr ATTENTION: \code{\link{check_param}}
#'   requires that \code{param[9,1]=}\eqn{R>0} be specified. \cr
#' 
#' \item\strong{OPTIONAL PARAMETERS SPECIFIC TO \code{\link{estimate_IMc}}: }
#'   \item\code{param[10,1]~}\emph{\code{howmany}} \code{param[10,1]=}\eqn{H}
#'   specifies \eqn{H}, the number of data sets to simulate per set of
#'   parameters (i.e., per grid point) to estimate the likelihood of the data
#'   given the set of parameters of the extension of the
#'   isolation-with-migration model.\cr By default,
#'   \code{param[10,1]=}\eqn{1000}.\cr ATTENTION: \code{\link{check_param}}
#'   requires that \code{param[10,1]=}\eqn{H>0} when specified. \cr The
#'   information on \eqn{J} can be specified with the keyword
#'   \emph{\code{howmany}} in the file \code{paramfile} (see
#'   \code{\link{param_sim}} and \code{\link{param_est}}). \cr % end howmany
#'   \item\code{param[11,1]~}\emph{\code{parallel}} \code{param[11,1]=}\eqn{J}
#'   specifies \eqn{J}, the number of jobs to run in parallel to perform the
#'   estimation of the posterior distribution of the parameters of the
#'   extension of the isolation-with-migration model.\cr By default,
#'   \code{param[11,1]=}\eqn{1}, i.e., no parallelization.\cr ATTENTION:
#'   \code{\link{check_param}} requires that \code{param[11,1]=}\eqn{J>0} when
#'   specified.\cr The information on \eqn{J} can be specified with the keyword
#'   \emph{\code{parallel}} in the file \code{paramfile} (see
#'   \code{\link{param_sim}} and \code{\link{param_est}}). \cr % end parallel %
#'   end extra
#' @param paramfile The name of the file with the values and descriptions of
#'   the prior distributions for the fixed and variable (or estimated)
#'   parameters, respectively, required to estimate the parameters of
#'   extensions of the Isolation-with-migration model with the function
#'   \code{\link{estimate_IMc}}.\cr By default the file name is
#'   \code{"estimation.par"}.
#' @param listparam The list of possible keywords/parameters in the required
#'   order: \cr For function \code{\link{simulate_data}}:\cr
#'   \code{c("}\emph{\code{theta_1}}\code{", "}\emph{\code{M_present}}\code{",
#'   "}\emph{\code{theta_2}}\code{", "}\emph{\code{theta_A}}\code{",
#'   "}\emph{\code{T_split}}\code{", "}\emph{\code{T_change}}\code{",
#'   "}\emph{\code{M_change}}\code{", "}\emph{\code{rho}}\code{",
#'   "}\emph{\code{nregions}}\code{")}.\cr For function
#'   \code{\link{estimate_IMc}}:\cr \code{c("}\emph{\code{theta_1}}\code{",
#'   "}\emph{\code{M_present}}\code{", "}\emph{\code{theta_2}}\code{",
#'   "}\emph{\code{theta_A}}\code{", "}\emph{\code{T_split}}\code{",
#'   "}\emph{\code{T_change}}\code{", "}\emph{\code{M_change}}\code{",
#'   "}\emph{\code{rho}}\code{", "}\emph{\code{nregions}}\code{",
#'   "}\emph{\code{howmany}}\code{", "}\emph{\code{parallel}}\code{")}.
#' @return The function \code{\link{check_param}} outputs error and warning
#'   messages regarding the format/information in the matrix of values
#'   \code{param} as well as the following data frame: \item{list("$ok")}{
#'   \code{$ok} takes the value \code{"1"} if the format of the input file
#'   \code{paramfile} for the function \code{\link{simulate_data}} (or
#'   \code{\link{estimate_IMc}}) has the format/information required. \cr
#'   Otherwise, \code{$ok} takes the value \code{"0"} to inform the function
#'   \code{\link{simulate_data}} (or \code{\link{estimate_IMc}}) that it should
#'   stop. \cr In this later, case,\code{\link{check_param}} outputs an error
#'   message explaining what part of the format/information is incorrect.  }
#'   \item{list("$param")}{ The matrix \eqn{9} (or \eqn{11}) \eqn{\times}
#'   \eqn{3} of the information for the parameters of the model that will be
#'   simulated provided as input to the function \code{\link{check_param}}, but
#'   updated to take into account parameters that will be ignored by the
#'   function \code{\link{simulate_data}} (or \code{\link{estimate_IMc}}).\cr
#'   When the matrix of values is updated, \code{\link{check_param}} outputs a
#'   "WARNING" message.\cr } e.g. for the function \code{\link{estimate_IMc}}:
#'   \code{[1] "PROBLEM: Error message on the parameter values and prior
#'   distributions for the }\R\code{ function estimate_IMc."}\cr \code{$ok} \cr
#'   \code{[1] 0}\cr \code{$param}\cr \tabular{llll}{ \tab \code{[,1]} \tab
#'   \code{[,2]} \tab \code{[,3]} \cr \code{[1,]} \tab \code{0.0001} \tab
#'   \code{1e-03} \tab \code{10}\cr \code{[2,]} \tab \code{0.0000} \tab
#'   \code{1e+01} \tab \code{20}\cr \code{[3,]} \tab \code{ 1.5000} \tab
#'   \code{1e-03} \tab \code{10}\cr \code{[4,]} \tab \code{0.0001} \tab
#'   \code{1e-03} \tab \code{10}\cr \code{[5,]} \tab \code{0.0000} \tab
#'   \code{2e+00} \tab \code{10}\cr \code{[6,]} \tab \code{ 0.0000} \tab
#'   \code{1e+00} \tab \code{10}\cr \code{[7,]} \tab \code{0.0000} \tab
#'   \code{2e+01} \tab \code{20}\cr \code{[8,]} \tab \code{0.0005} \tab
#'   \code{NA} \tab \code{NA}\cr \code{[9,]} \tab \code{4.0000} \tab \code{NA}
#'   \tab \code{NA}\cr \code{[10,]} \tab \code{NA} \tab \code{NA} \tab \code{
#'   NA}\cr \code{[11,]} \tab \code{NA} \tab \code{NA} \tab \code{ NA} }
#' @note \itemATTENTION: It is the user's responsibility to mind the following
#'   restrictions: -> \code{\link{check_param}} requires that
#'   \code{param[1,1]=}\eqn{R>0} be specified. \cr -> \code{\link{check_param}}
#'   requires that either the parameter value (i.e.
#'   \code{param[1,1]=}\eqn{\theta_1>0}) or the prior distribution for
#'   \eqn{theta_1} (see above for the restrictions) be specified in
#'   \code{param[1,]}. \cr -> \code{\link{check_param}} requires that
#'   \eqn{M_p>0} when fixed or \eqn{0\le}\eqn{M_pl}\eqn{\le}\eqn{M_pu} and
#'   \eqn{M_pn>0} when variable (or estimated) be specified in \code{param[2,]}
#'   in case of an island model. \cr -> \code{\link{check_param}} requires that
#'   \code{param[j,i]=}\eqn{theta_i}\eqn{\ge}\eqn{0} for any \eqn{i} \eqn{\in}
#'   \eqn{{2,A}} and \eqn{j} \eqn{\in} \eqn{{2,3}}. \cr ->
#'   \code{\link{check_param}} requires that
#'   \eqn{0}\eqn{\le}\eqn{T_c}\eqn{\le}\eqn{T_s} when specified and fixed in
#'   \code{param[5,1]} and \code{param[6,1]}.\cr -> \code{\link{check_param}}
#'   requires that \code{param[10,1]=}\eqn{H>0} when specified. \cr ->
#'   \code{\link{check_param}} requires that \code{param[11,1]=}\eqn{J>0} when
#'   specified.\cr
#' @author Celine Becquet - \email{celine.becquet@@gmail.com}.
#' @seealso The functions \code{\link{simulate_data}} and
#'   \code{\link{estimate_IMc}} calls the function \code{\link{check_param}},
#'   which in turn calls the function \code{\link{error_message}}.\cr The
#'   function \code{\link{check_param}} checks the format of the files like
#'   \code{\link{param_sim}} and \code{\link{param_est}}.\cr Other functions to
#'   check the format of input files: \cr \code{\link{get_rho}} and
#'   \code{\link{check_reg}}. \cr Lists of definitions of the symbols and
#'   parameters mentioned in this file are found in
#'   \code{\link{Rmspack-package}}.
#' @keywords error print
#' @examples
#' 
#' ### Write the file of information on the parameters in the local directory.
#' data(simulation_files) # download the data
#' 
#' write.table(file="est", x=param_est, row.name=FALSE, col.names=FALSE, quote=FALSE)
#' # Creates the file "est" containing (with comments):
#' read.csv("est",header =FALSE,sep="", comment="#")
#' 
#' ### Create the inputs for the function.
#' listparam=c("theta_1", "M_present", "theta_2", "theta_A", "T_split", "T_change", "M_change", "rho", "nregions", "howmany", "parallel")
#' param=scan("est", comment.char="#", what=c("c", "n", "n", "n"), allowEscapes=FALSE, fill=TRUE, sep="\n", strip.white=TRUE, quiet=TRUE) 
#' vparam=order_param(param, listparam) 
#' vparam
#' 
#' ## Case with no errors.
#' check_param(param=vparam, paramfile="est", listparam=listparam)
#' 
#' ## Case with warning
#' vparam[3, 1]=.001 # theta_2l==theta_2u
#' check_param(param=vparam, paramfile="est", listparam=listparam)
#' 
#' 
#' ## Case with error in defining the prior of theta_2.
#' vparam[3, 1]=.002 # theta_2l>theta_2u
#' check_param(param=vparam, paramfile="est", listparam=listparam)
#' 
#' # Clean up the directory.
#' unlink("est")
#' 
NULL





#' Function to check the format of the files with the information on the
#' genomic regions and multiple loci.
#' The function \code{\link{check_reg}} checks that the information on the
#' independent genomic regions described in the file \code{regfile} and the
#' information on the loci for the multi-locus genomic regions in the file
#' \code{locifile} are in the correct format to simulate data with the function
#' \code{\link{simulate_data}} or estimate models of Isolation-with-migration
#' and possible extensions with the function \code{\link{estimate_IMc}}. \cr
#' This function is called by the functions \code{\link{simulate_data}} and
#' \code{\link{estimate_IMc}}.
#' 
#' 
#' @param nregions The number of independent genomic regions considered,
#'   \eqn{R}. \cr \eqn{R} is specified with the keyword \emph{\code{nregions}}
#'   in the file \code{paramfile} (see \code{\link{param_sim}} and
#'   \code{\link{param_est}}). \cr
#' @param info_region The matrix \eqn{R} \eqn{\times} \eqn{10} of information
#'   on the \eqn{R} independent genomic regions. \cr Each genomic region is
#'   described by ten values: \tabular{ll}{ \eqn{r} \tab : The genomic region
#'   number, \eqn{r} \eqn{\in} \eqn{[1,R]}. \cr \tab \eqn{R} is the number of
#'   independent genomic regions specified with the argument
#'   \emph{\code{nregions}}. \cr \tab ATTENTION: The independent genomic region
#'   numbers, \eqn{r}, need to be in order, i.e., region 1 starts on line 2,
#'   region 2 on line 3 \dots{} region \eqn{R} on line \eqn{R+1}. \cr
#'   \emph{\code{Region name}} \tab : The name of the genomic region \eqn{r},
#'   which should contain at least ONE non-numerical character. \cr \eqn{x}
#'   \tab : The inheritance scalar for the genomic region \eqn{r} (i.e.,
#'   \code{"1"} for autosomal region, \code{"0.75"} for X- and \code{"0.5"} for
#'   Y- and mtDNA-linked region). \cr \eqn{v} \tab : The mutation rate scalar
#'   for the genomic region \eqn{r} (which can be estimated e.g., from
#'   divergence data). \cr \eqn{w} \tab : The recombination scalar for the
#'   genomic region \eqn{r}.\cr \tab -- Usually \eqn{w=\beta}, the ratio of the
#'   locus-specific population recombination rate per bp over
#'   \eqn{\rho=4N_1*c}. \cr \tab -- If an estimate of the region-specific
#'   population recombination rate per bp is available for each region from
#'   linkage disequilibrium analysis, \eqn{\rho_o=4N_1*c_o}, set
#'   \eqn{w=\beta*\rho_o} to incorporate this knowledge in the simulation or
#'   estimation (with \code{"}\emph{\code{rho}} \code{1"} in the file
#'   \code{paramfile}, see \code{\link{param_sim}} and
#'   \code{\link{param_est}}).\cr \tab In this case, \eqn{w} is the scaled
#'   sex-averaged region-specific population recombination rate per bp, i.e.,
#'   for an X-linked locus \eqn{c_o} is the female recombination rate and
#'   \eqn{\beta=0.5} so that \eqn{\beta*\rho_o=2N_1*c_o}. \cr \tab -- If an
#'   estimate of the region-specific recombination rate per bp is available for
#'   each region from pedigree analysis, \eqn{c_o}, set \eqn{w=\beta*c_o} to
#'   incorporate this knowledge in the simulation or estimation (with
#'   \code{"}\emph{\code{rho}} \code{2"} in the file \code{paramfile}, see
#'   \code{\link{param_sim}} and \code{\link{param_est}}). \cr \tab In this
#'   case, \eqn{w} is the scaled sex-averaged region-specific recombination
#'   rate per bp, i.e., for an X-linked locus, \eqn{c_o} is the estimated
#'   female recombination rate so the scaled sex-averaged recombination rate is
#'   \eqn{\beta*c_o=0.5c}. \cr \eqn{n_1} \tab : The sample size from population
#'   1 for the genomic region \eqn{r}. \cr \eqn{n_2} \tab : The sample size
#'   from population 2 for the genomic region \eqn{r}. \cr \eqn{z_s} \tab : The
#'   start position of the genomic region \eqn{r} in bp. \cr \eqn{z_e} \tab :
#'   The end position of the genomic region \eqn{r} in bp. \cr \eqn{Y} \tab :
#'   The number of loci spanning the genomic region \eqn{r}. \cr }Where:
#'   \tabular{ll}{ \eqn{\rho=4N_1*c}\tab : The genomic average population
#'   intra-region recombination rate per bp per generation. \cr
#'   \eqn{\rho_o=4N_1*c_o} \tab : The estimate of the region-specific
#'   population recombination rate per bp per generation for the genomic region
#'   considered from linkage disequilibrium analysis. \cr \eqn{N_1} \tab : The
#'   effective population size in population 1 (the reference population). \cr
#'   \eqn{c} \tab : The genomic generational cross-over rate per bp. \cr
#'   \eqn{c_o} \tab : The estimate of the region-specific cross-over rate per
#'   bp per generation for the genomic region considered from pedigree
#'   analysis.\cr \eqn{\beta} \tab : The ratio of the region-specific
#'   population recombination rate per bp over \eqn{\rho=4N_1*c} for the
#'   genomic region considered.\cr \tab \eqn{\beta=}\code{"1"} (\code{"0.5"} in
#'   \emph{Drosophila}) for autosomal region, \code{"0.5"} for X- and
#'   \code{"0"} for Y- and mtDNA-linked region.  } See
#'   \code{\link{info_region}} for further details.
#' @param info_loci The matrix \eqn{\Sigma{Y_r}} \eqn{\times} \eqn{6} of
#'   information on the loci for the multi-locus genomic regions (as specified
#'   with \eqn{Y_r>1} for the multi-locus region \eqn{r} in
#'   \code{info_region[}\eqn{r}\code{]$V10}). \cr This matrix is empty unless
#'   \eqn{\Sigma{Y_r}>0}. \cr Each locus is described by six values:
#'   \tabular{ll}{ \eqn{r} \tab : The multi-locus genomic region number,
#'   \eqn{r} \eqn{\in} \eqn{[1,R]}, that this locus is part of. \cr \tab
#'   ATTENTION: The multi-locus genomic region numbers, \eqn{r}, need to be in
#'   order, i.e., the information for the loci for the first multi-locus
#'   genomic region (\eqn{a}) are from line 2 to \eqn{Y_a+1}, for the second
#'   multi-locus genomic region (\eqn{b}), the loci information starts on line
#'   \eqn{Y_a+2} to \eqn{Y_a+Y_b+1} \dots{} \eqn{a} and \eqn{b} \eqn{\in}
#'   \eqn{[1,R]}. \cr \eqn{y} \tab : The locus number, \eqn{y} \eqn{\in}
#'   \eqn{[1,Y]}. \cr \tab \eqn{Y} is the total number of loci spanning the
#'   multi-locus genomic region \eqn{r} as specified in the matrix
#'   \code{info_region}. \cr \tab ATTENTION: The loci numbers, \eqn{y}, need to
#'   be in order, i.e., information for locus 1 of the first multi-locus
#'   genomic region \eqn{r} is on line 2, locus 2 on line 3 \dots{} locus
#'   \eqn{Y} on line \eqn{Y+1}. \cr \eqn{n_1y} \tab : The sample size from
#'   population 1 for the locus \eqn{y}.\cr \tab \eqn{n_1y}\eqn{\le}\eqn{n_1},
#'   where \eqn{n_1} for the multi-locus genomic region \eqn{r} is specified in
#'   the matrix \code{info_region}. \cr \eqn{n_2y} \tab : The sample size from
#'   population 2 for the locus \eqn{y}.\cr \tab \eqn{n_2y}\eqn{\le}\eqn{n_2},
#'   where \eqn{n_2} for the multi-locus genomic region \eqn{r} is specified in
#'   the matrix \code{info_region}. \cr \eqn{z_sy} \tab : The start position of
#'   the locus \eqn{y} in bp.\cr \tab \eqn{z_s1}\eqn{\ge}\eqn{z_s}, where
#'   \eqn{z_s} is the start position of the multi-locus genomic region \eqn{r}
#'   specified in the matrix \code{info_region}. \cr \eqn{z_ey} \tab : The end
#'   position of the locus \eqn{y} in bp.\cr \tab \eqn{z_eY}\eqn{\le}\eqn{z_e}
#'   where \eqn{z_e} is the end position of the multi-locus genomic region
#'   \eqn{r} specified in the matrix \code{info_region}.  } See
#'   \code{\link{info_loci}} for further details.
#' @param locifile The name of the file with the information on the loci for
#'   the multi-locus genomic regions (as specified with \eqn{Y>1} in the file
#'   \code{regfile}). \cr By default the file name is \code{"info.loc"}.
#' @return The function \code{\link{check_reg}} outputs \code{"1"} if the
#'   format of the matrices \code{info_region} (and \code{info_loci}) (from the
#'   input file \code{regfile} (and \code{locifile})) for the calling function
#'   (either \code{\link{simulate_data}} or \code{\link{estimate_IMc}}) have
#'   the format/information required. \cr Otherwise, \code{\link{check_reg}}
#'   informs the calling function that it should stop and outputs an error
#'   message explaining what part of the format/information is incorrect: \cr
#'   \tabular{ll}{ e.g.:\tab \code{[1] "Error message on independent
#'   (potentially multi-locus) genomic regions."}\cr \tab \code{[1] 0} }
#' @note \itemATTENTION: It is the user's responsibility to mind the following
#'   restrictions: \itemIn the matrix \code{info_region}: -> The independent
#'   genomic region numbers, \eqn{r}, need to be in order, i.e., region 1
#'   starts on line 2, region 2 on line 3 \dots{} region \eqn{R} on line
#'   \eqn{R+1}. \cr -> Reasonable values need to be specified for each genomic
#'   region.
#' 
#' \itemIn the matrix \code{info_loci}: -> The loci numbers, \eqn{y}, need to
#'   be in order, i.e., information for locus 1 of the first multi-locus
#'   genomic region is on line 2, locus 2 on line 3 \dots{} locus \eqn{Y} on
#'   line \eqn{Y+1}. \cr -> The multi-locus genomic region numbers, \eqn{r},
#'   need to be in order, i.e., the information for the loci for the first
#'   multi-locus genomic region (\eqn{a}) are from line 2 to \eqn{Y_a+1}, for
#'   the second multi-locus genomic region (\eqn{b}), the loci information
#'   starts on line \eqn{Y_a+2} to \eqn{Y_a+Y_b+1} \dots{} \eqn{a} and \eqn{b}
#'   are \eqn{\in} \eqn{[1,R]}. \cr -> Reasonable values need to be specified
#'   for each locus.
#' @author Celine Becquet - \email{celine.becquet@@gmail.com}.
#' @seealso The functions \code{\link{simulate_data}} and
#'   \code{\link{estimate_IMc}} call the function \code{\link{check_reg}},
#'   which in turn calls the function \code{\link{error_message}}. \cr Other
#'   functions to check the format of input files: \cr \code{\link{get_rho}}
#'   and \code{\link{check_param}}. \cr Lists of definitions of the symbols and
#'   parameters mentioned in this file are found in
#'   \code{\link{Rmspack-package}}.
#' @keywords error print
#' @examples
#' 
#' ### Write the files of information in the local directory.
#' data(simulation_files) # download the data
#' 
#' write.table(file="regions", x=info_region, row.name=FALSE, col.names=FALSE, quote=FALSE)
#' # Creates the file "regions" containing: 
#' read.csv("regions",header =FALSE,sep="", comment="#")
#' 
#' write.table(file="loci", x=info_loci, row.name=FALSE, col.names=FALSE, quote=FALSE)
#' # Creates the file "loci" containing: 
#' read.csv("loci",header =FALSE,sep="", comment="#")
#' 
#' ### Create the inputs for the function.
#' reg=read.table("regions", skip=1, fill=TRUE)
#' loci=read.table("loci", skip=1, fill=TRUE)
#' 
#' ## Case with no errors.
#' check_reg(nregions=4, info_region=reg, info_loci=loci, locifile="loci")
#' 
#' ## Case with error in the # of regions
#' reg[1, 1]=2
#' check_reg(nregions=4, info_region=reg, info_loci=loci, locifile="loci")
#' 
#' # Clean up the directory.
#' unlink(c("regions", "loci"))
#' 
NULL



