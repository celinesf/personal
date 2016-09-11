\name{Rmspack-package}
\alias{Rmspack-package}
\alias{Rmspack}
\docType{package}

\title{
  A package to estimate the parameters of extensions of the
  Isolation-with-migration model. 
}

\section{Symbols}{
  \item{The following tables list the definitions of the symbols mentioned in
    this file.}{
    \item{\bold{\acronym{Population mutation rate parameters}}}{
      \tabular{ll}{
        \eqn{\theta_1*x*v*Z} \tab
        : The region-specific population mutation rate per generation for
        population 1. \cr
        \eqn{\theta_2*x*v*Z} \tab
        : The region-specific population mutation rate per generation for
        population 2. \cr
        \eqn{\theta_A*x*v*Z} \tab
        : The region-specific ancestral population mutation rate per
        generation. \cr 
        \eqn{\theta_1=4N_1*\mu} \tab
        : The population mutation rate per bp per generation for
        population 1. \cr 
        \eqn{\theta_2=4N_2*\mu} \tab
        : The population mutation rate per bp per generation for
        population 2. \cr 
        \eqn{\theta_A=4N_A*\mu} \tab
        : The ancestral population mutation rate per bp per
        generation. \cr 
        \eqn{\mu} \tab
        : The genomic generational mutation rate per bp. \cr
        \eqn{x} \tab 
        : The inheritance scalar for the genomic region considered (i.e.,
        \code{"1"} for autosomal region, \code{"0.75"} for X- and \code{"0.5"}
        for Y- and mtDNA-linked region). \cr 
        \eqn{v} \tab 
        : The mutation rate scalar for the genomic region considered (which
        can be estimated e.g., from divergence data). \cr 
      }
    }
    \item{\bold{\acronym{Population recombination rate parameters}}}{
      \tabular{ll}{
        \eqn{\rho=4N_1*c} \tab
        : The genomic average population intra-region recombination rate
        per bp per generation. \cr 
        \eqn{c} \tab
        : The generational cross-over rate per bp. \cr
        \eqn{w} \tab
        : The recombination scalar for the genomic region
        considered. Usually \eqn{w=\beta}. \cr
        \eqn{\beta} \tab
        : The ratio of the region-specific population recombination
        rate per bp over \eqn{\rho=4N_1*c} for the genomic region
        considered.\cr
        \tab i.e., \code{"1"} (\code{"0.5"} in \emph{Drosophila}) for
        autosomal region, \code{"0.5"} for X- and \code{"0"} for Y- and
        mtDNA-linked region.\cr 
        \eqn{\rho_r=\rho*\beta*(Z}\eqn{-}\eqn{1)}\tab
        : The region-specific recombination rate per generation for the
        genomic region \eqn{r}. \cr  
        \eqn{\rho_o=4N_1*c_o} \tab
        : The estimate of the region-specific population recombination
        rate per bp per generation for the genomic region considered from
        linkage disequilibrium analysis. \cr
        \eqn{c_o} \tab 
        : The estimate of the region-specific cross-over rate per bp
        per generation for the genomic region considered from pedigree
        analysis. \cr
        \eqn{\alpha=c_r/\mu} \tab Drawn from a distribution. \cr
        \eqn{c_r} \tab
        : The generational region-specific cross-over rate per bp for the
        genomic region \eqn{r}. \cr 
      }
    }
    \item{\bold{\acronym{Parameters of time of events}}}{
      \tabular{ll}{
        \eqn{T_s} \tab
        : The split time in unit of \eqn{4N_1} generations between the
        TWO populations. \cr
        \eqn{T_c} \tab
        : The time at which the rate of gene flow changed between the
        two populations in unit of \eqn{4N_1} generations. \cr
        \eqn{\epsilon=T_c/T_s} \tab 
      }
    }
    \item{\bold{\acronym{Parameters of gene flow rates}}}{
      \tabular{ll}{
        \eqn{M_p=4N_1*m_p} \tab 
        : The number of migrants exchanged each generation by the TWO
        populations at present. \cr
        \tab \eqn{M_p} has different definitions depending of the model
        considered: \cr 
        \tab * \eqn{M_p} is the symmetrical rate of gene flow between
        TWO populations in an island model (i.e., \eqn{0<t<\inf}) as specified
        with: \cr 
        \tab -- At least one of the independent genomic regions has
        \eqn{n_1>0} and \eqn{n_2>0}. \cr 
        \tab -- Either \eqn{T_s=0} when specified or unspecified. \cr
        \tab * \eqn{M_p} is the constant symmetrical rate of gene flow
        since the split until present (i.e., \eqn{0<t<T_s}) if
        there is a population split (as specified with \eqn{T_s>0}). \cr 
        \tab * \eqn{M_p} is the constant symmetrical rate of gene flow
        since the time of gene flow rate change until present (i.e.,
        \eqn{0<t<T_c}) if a time at which the gene flow rate changed is
        specified with: \cr 
        \tab -- \eqn{0<T_c<T_s}. \cr 
        \tab -- \eqn{M_c \not=M_p}. \cr 
        \eqn{M_c=4N_1*m_c} \tab 
        : The number of migrants exchanged each generation by the TWO
        populations since the split until the time of gene flow rate
        change (i.e., \eqn{T_c<t<T_s}). \cr 
        \eqn{m_p} \tab 
        : The generational fraction of migrant individuals at present. \cr 
        \eqn{m_c} \tab
        : The generational fraction of migrant individuals between
        \eqn{T_c<t<T_s}. 
      }
    }
    \item{\bold{\acronym{List of other symbols}}}{
      \tabular{ll}{
        \eqn{n_1} \tab
        : The sample size from population 1 for the genomic region
        considered. \cr 
        \eqn{n_2} \tab
        : The sample size from population 2 for the genomic region
        considered.\cr 
        \eqn{N_1} \tab
        : The effective population size in population 1 (the reference
        population so by default and unless specified, \eqn{N_1=N_2=N_A}). \cr 
        \eqn{N_2} \tab
        : The effective population size in population 2. \cr 
        \eqn{N_A} \tab
        : The ancestral effective population size. \cr
        \eqn{R} \tab
        : The number of independent genomic regions considered. \cr
        \eqn{Z} \tab
        : The size in bp of the genomic region considered. \cr
        \eqn{Y} \tab
        : The number of loci spanning the genomic region considered. 
      }
    }
  }
}% end symbols

\description{
  \item{The package \code{\link{Rmspack}} is designed to estimate parameters
    of extensions of the Isolation-with-migration model. More
    specifically, it can estimates the following parameters with the
    function \code{\link{estimate_IMc}}:}{
    -- The population mutation rates per bp per generation for the population
    1, 2 and the ancestral population (\eqn{\theta_1}, \eqn{\theta_2} and
    \eqn{\theta_A}, respectively).\cr
    -- The time of events that ocured during the demographic history of the
    populations, \eqn{T_s} and \eqn{\epsilon=T_c/T_s}.\cr
    -- The rates of gene flow during different time periods, \eqn{M_p} and \eqn{M_c}.\cr
    See the lists of definitions of parameters in the section
    \code{Symbols}.
  }\cr
  \item{\code{\link{Rmspack}} also helps run of Dick Hudson C program
    \code{ms}. It allows the following:}{
    -- To generate command lines to run the program of Dick Hudson
    \code{ms} or the \R function \code{\link{ms}} (\R function
    \code{\link{cmd_line}}).\cr 
    -- To transform \code{ms}-like data from a text file into the data
    frame output by the \R function \code{\link{ms}} (\R function
    \code{\link{get_ms_output}}).\cr 
    -- To calculate sample statistics and population statistics from the
    data frame output by the \R function \code{\link{ms}} (\R function
    and \code{\link{stats_pop}}).\cr  
    -- To simulate data sets of independent regions in the format required
    by the estimation function \code{\link{estimate_IMc}} (\R function
    \code{\link{simulate_data}}).\cr 
  }
}

\details{
  \tabular{ll}{
    Package: \tab Rmspack\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2009-09-18\cr
    License: \tab Artistic License, v2.0\cr
    LazyLoad: \tab no
  }
  \item{The ready to use functions in the package are:}{
    -- \code{estimate_IMc(paramfile = "estimation.par", datafile = "datafiles/reg", regfile = "datafiles/info.reg", locifile = "datafiles/info.loc")}\cr
    -- \code{cmd_line(nsam,theta,howmany=1,rho=0,struct=1,seeds="NA",time=0,mig=0,extra="NA")}\cr
    -- \code{ms(cmdline)}\cr
    -- \code{stats_pop(msframe, type = "anc")}\cr
    -- \code{simulate_data(argfile = "simulation.arg", regfile = "info.reg", locifile = "info.loc")}\cr
    -- \code{get_ms_output(msoutput = "msout")}\cr
    -- \code{write_ms_output(msframe, msoutput="msout")}
  }
}
\author{Celine Becquet - \email{celine.becquet@gmail.com}.
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
  \item{All the functions in \code{\link{Rmspack}} are listed below:}{
    \item{Function \code{\link{check_param}}}{
      The function \code{\link{check_param}} checks that the information
      on the values and descriptions of the prior distributions for the
      fixed and variable (or estimated) parameters, respectively, specified
      in the file \code{paramfile} (see \code{\link{param_est}} and
      \code{\link{param_sim}}) are in the correct format to simulate data
      under (or to estimate parameters of) extensions of the
      Isolation-with-migration model with the function
      \code{\link{simulate_data}} (or  \code{\link{estimate_IMc}}). \cr  
      This function is called within the \R functions
      \code{\link{simulate_data}} and \code{\link{estimate_IMc}}.
    }
    \item{Function \code{\link{check_reg}}}{
      The function \code{\link{check_reg}} checks that the information on the
      independent genomic regions described in the file \code{regfile} and
      the information on the loci for the multi-locus genomic regions in the
      file \code{locifile} are in the correct format to simulate data with
      the function \code{\link{simulate_data}} or estimate models of
      Isolation-with-migration and possible extensions with the function
      \code{\link{estimate_IMc}}. \cr 
      This function is called by the \R functions \code{\link{simulate_data}}
      and \code{\link{estimate_IMc}}. 
    }
    \item{Function \code{\link{cmd_line}}}{
      The function \code{\link{cmd_line}} takes values for the number of
      samples, the region-specific population mutation and recombination
      rates, the parameters of the demographic model to simulate... and
      generates the \code{ms}-like command line with the tags and values
      required as input by the \R function \code{\link{ms}}. \cr 
      This function is called by the \R functions \code{\link{simulate_data}}
      and \code{\link{estimate_IMc}}. 
    }
    \item{Function \code{\link{error_message}}}{
      The function \code{\link{error_message}} pastes together and prints
      strings of characters, value of parameters and other relevant
      information into an error message. \cr 
      This function is called by the \R functions checking the format of data
      files.
    }
    \item{Function \code{\link{estimate_IMc}}}{
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
      The function \code{\link{estimate_IMc}} write values at the 
      summary statistics listed in \code{summaries} of the simulated polymorphisms for the independent genomic
      regions considered (and their loci is applicable) in the file named
      \code{"datafile-est-S"} (\code{"datafile-loci-est-S"}). \cr
  \code{\link{estimate_IMc}} also prints the string of characters
  with the parameter values and other information required to run the C
  program \code{estlike} into text files starting with  \code{"datafile-lik-"}.
    }
    \item{Function \code{\link{get_liststats}}}{
      The function  \code{\link{get_liststats}} compares the vectors of
      summary statistics \code{liststats} and \code{summaries} to recover
      the positions in the vector \code{liststats} of the summary statistics
      listed in  vector \code{summaries}.\cr
      \code{\link{get_liststats}} outputs a data frame
      listing the position in of the summary statistics considered in the 
      data frame of the summary statistics output by the function
      \code{\link{stats_pop}}.\cr  
      This function is called within the \R functions
      \code{\link{simulate_data}} and \code{\link{estimate_IMc}}.   
    }
    \item{Function \code{\link{get_loci}}}{
      The function  \code{\link{get_loci}} extracts the data for the locus
      \code{nlocus} from a data frame \code{msframe} of the kind output by
      the \R function 
      \code{\link{ms}} for a multi-locus region.\cr
      This function is called within the \R functions
      \code{\link{simulate_data}} and \code{\link{estimate_IMc}}. 
    }
    \item{Function \code{\link{get_ms_output}}}{
      The function \code{\link{get_ms_output}} transforms text file in the
      format as output by Hudson's C program \code{ms} into the data frame
      output by the \R function \code{\link{ms}}. \cr 
      This function is called by the \R function
      \code{\link{estimate_IMc}}. 
    }
    \item{Function \code{\link{get_nparam}}}{
      The function \code{\link{get_nparam}} outputs the position in the list
      \code{listparam} of the keyword specified in the argument
      \code{"keyword"}. \cr 
      This function is called by the \R function
      \code{\link{order_param}}.
    }
    \item{Function \code{\link{get_rho}}}{
      The function \code{\link{get_rho}} checks that the recombination
      scalar \eqn{w} specified for the recombining genomic region \eqn{r}
      in the file \code{regfile} (see \code{\link{info_region}})
      is in the correct format to simulate data with the function
      \code{\link{simulate_data}}.\cr 
      If so \code{\link{get_rho}} outputs the region-specific
      population recombination rate, \eqn{\rho_r=\beta*4N_1c*(Z-1)} for
      the recombining region \eqn{r}. \cr
      This function is called by the \R function \code{\link{simulate_data}}.
    }
    \item{Function \code{\link{get_stats}}}{
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
    \item{Function \code{\link{ms}}}{
      The function \code{\link{ms}} calls Hudson's C program \code{ms} to
      simulate independent data sets for a genomic region of interest.
      \code{\link{ms}} uses as input an \code{ms}-like command line defined
      by the user. \cr 
      This function is called by the \R functions \code{\link{simulate_data}}
      and \code{\link{estimate_IMc}}.
    }
    \item{Function \code{\link{order_param}}}{
      The function \code{\link{order_param}} extracts the numerical
      information from an unordered vector of strings of characters,
      \code{param}, into a matrix following the order of
      keywords/arguments/parameters specified by \code{listparam}. \cr 
      This function is called by the \R functions \code{\link{simulate_data}}
      and \code{\link{estimate_IMc}}.
    }
    \item{Function \code{\link{out_ms_loci}}}{
      The function \code{\link{out_ms_loci}} updates the data frame returned
      by the \R function \code{\link{ms}} to account for missing data and gaps
      when there are multiple loci in a specific genomic region. \cr 
      This function is called by the \R function
      \code{\link{output_ms}}.
    }
    \item{Function \code{\link{output_ms}}}{
      The function \code{\link{output_ms}} outputs and writes the data frame
      \code{msframe} output by the \R function \code{\link{ms}} updated for
      multiple-loci if applicable into an
      \code{ms}-like text file named \code{"}\code{datafile}\eqn{r}\code{"},
      where \eqn{r} is the number of the genomic region considered specified
      by the argument \code{nregion}. \cr 
      This function is called by the \R function \code{\link{simulate_data}}.
    }
    \item{Function \code{\link{output_cmd}}}{
      The function  \code{\link{output_cmd}} prints the string of characters
      with the parameter values and other information required to run the C
      program \code{estlike} into text files starting with  \code{"datafile-"}. \cr
      This function is called within the \R function
      \code{\link{estimate_IMc}}. 
      }
    \item{Function \code{\link{output_stats}}}{
      The function  \code{\link{output_stats}} prints the considered summary
      statistics (specified in the vector \code{liststats}) of the
      independent genomic regions and/or linked loci  
      the data frame \code{sumstats} \code{\link{get_stats}}. \cr
      This function is called within the \R functions
      \code{\link{simulate_data}} and \code{\link{estimate_IMc}}. 
    }
    \item{Function \code{\link{simulate_data}}}{
      The function \code{\link{simulate_data}} simulates data in the format
      required by the \R function
      \code{\link{estimate_IMc}}. \code{\link{simulate_data}} uses the \R
      function \code{\link{ms}} with the set of parameters specified in the
      file \code{paramfile} as well as the information on independent
      genomic regions specified in the file \code{regfile}. \cr 
      If there are multiple loci with gaps in one or more regions considered
      (as specified by \eqn{Y>0} in the file \code{regfile}), the function
      \code{\link{simulate_data}} also requires the information on the
      multiple loci to be specified in the file \code{locifile}. \cr
      The parameters of the model specified in \code{paramfile} can be
      either fixed or variable.
      \code{\link{simulate_data}} writes the simulated data in \code{ms}-like text files
      named \code{"datafile}\eqn{r}\code{"}, where \eqn{r} \eqn{\in}
      \eqn{[1,R]} is the genomic region number. It also writes the \eqn{43}
      summary statistics of the simulated polymorphisms calculated by
      \code{\link{stats_pop}}
      for the independent genomic
      regions considered (and their loci is applicable) in the file named
      \code{"datafile-S"} (\code{"datafile-loci-S"}).
    }
    \item{Function \code{\link{stats_pop}}}{
      The function \code{\link{stats_pop}} calculates various summary
      statistics for the data frame \code{msframe} of the same format as
      output by the \R function \code{\link{ms}} for a model with TWO
      populations.\cr
      This function is called by the \R functions
      \code{\link{simulate_data}} and \code{\link{estimate_IMc}}.
    }
    \item{Function \code{\link{write_ms_output}}}{
      The function \code{\link{write_ms_output}} writes the data frame
      \code{msframe} output by the \R function \code{\link{ms}} into an
      \code{ms}-like text file named \code{msoutput}. \cr\cr
    }
  }% end of functions
  \item{\code{\link{Rmspack}} also contains example text files of the data
    files used as argument by the functions \code{\link{estimate_IMc}},
    \code{\link{simulate_data}} and \code{\link{get_ms_output}}:}{ 
    \item{File of data \code{\link{data_file}}}{
      The file \code{\link{data_file}} provides an example of the file
      \code{datafile}, which contains the data for a genomic region in the
      format as output by Hudson's C program \code{ms} required to estimate
      parameters of extensions of the Isolation-with-migration model with
      the \R function \code{\link{estimate_IMc}} or to transform into the data
      frame output by the \R function \code{\link{ms}} with the \R function
      \code{\link{get_ms_output}}.\cr 
      Text files in the format of \code{\link{data_file}} are output by the \R
      functions:\cr
      \code{\link{output_ms}}, \code{\link{simulate_data}} and
      \code{\link{write_ms_output}}. 
    }
    \item{File of information \code{\link{info_loci}}}{
      The file \code{\link{info_loci}} provides an example of the file
      \code{locifile}, which contains information on the loci for the
      multi-locus genomic regions (as specified with \eqn{Y>1} in the file
      \code{regfile}, see \code{\link{info_region}}) required to simulate
      data with the \R function \code{\link{simulate_data}} or estimate
      parameters of extensions of the Isolation-with-migration model with the
      \R function \code{\link{estimate_IMc}}. 
    }
    \item{File of information \code{\link{info_region}}}{
      The file \code{\link{info_region}} provides an example of the file
      \code{regfile}, which contains information on the independent genomic
      regions required to simulate data with the \R function
      \code{\link{simulate_data}} or estimate parameters of extensions of the
      Isolation-with-migration model with the \R function
      \code{\link{estimate_IMc}}. 
    }
    \item{File of information \code{\link{param_est}}}{
      The file \code{\link{param_est}} provides an example of the
      parameter file \code{paramfile}. This file contains the values and
      descriptions of the prior distributions for the fixed and estimated parameters,
      respectively, required to estimate the parameters of extensions of the
      Isolation-with-migration model with the \R function
      \code{\link{estimate_IMc}}.
    }
    \item{File of information \code{\link{param_sim}}}{ 
      The file \code{\link{param_sim}} provides an example of the parameter file
      \code{paramfile}, which contains the parameter values required to
      simulate data with the \R function \code{\link{simulate_data}}.
    }
  }% end data files
}

\examples{
##########
######### RUN ms ##########
np=c(2, 5, 6) # npop n1 n2.
ns=np[2]+np[3] # total number of samples.
h=2 # independent gene genealogies.
s=c("1", "2", "3") # Seed numbers.

th=c(5, 10) # 4Ne*mu*Z for pop1 and pop2, here theta2/Theta1=2.
rh=c(1, 1000) # 4Ne*c*Z and Z, recombination rate and size of locus in bp.

T=c(1, 0.5) # Time of divergence and time of migration rate change in unit of 4Ne generations. 
M=c(5, 0) # Migration rate until present and eM, the migration rate between T_div and T_change, forward in time.

## Generate a command for an IM model with a change in gene flow rate 0.5/4N_1 generations ago.
cmdline=cmd_line(nsam=ns, howmany=h, theta=th, rho=rh, structure=np, seeds=s, time=T, migration=M) 
cmdline

## Run ms.
msframe=ms(cmdline=cmdline)
msframe

##########
######### Calculate population statistics ########
stats_pop(msframe=msframe)

## unknown ancestral allele
stats_pop(msframe=msframe, type="maf")

##########
######### Simulate data in the format for estimation ########
data(simulation_files) # download the data

### Write the information on text files in the local directory
write.table(file="sim", x=param_sim, row.name=FALSE, col.names=FALSE, quote=FALSE )
# Creates the file "sim" (with comments):
read.csv("sim",header =FALSE,sep="", comment="#")

write.table(file="regions", x=info_region, row.name=FALSE, col.names=FALSE, quote=FALSE )
# Creates the file "regions" containing:
read.csv("regions",header =FALSE,sep="", comment="#")

write.table(file="loci", x=info_loci, row.name=FALSE, col.names=FALSE, quote=FALSE )
# Creates the file "loci" containing:
read.csv("loci",header =FALSE,sep="", comment="#")

## Create 4 files named regr (r in [1,4]) with simulated data and files
## of summaries reg-S and reg-loci-S.
simulate_data(paramfile="sim", regfile="regions", locifile="loci");

## download the sampled regions with function get_ms_output
get_ms_output("reg1")
get_ms_output("reg2")
get_ms_output("reg3")
get_ms_output("reg4")

##########
######### Estimate parameters of IMc ########
write.table(file="est", x=param_est, row.name=FALSE, col.names=FALSE, quote=FALSE )
# Creates the file "est" (with comments):
read.csv("est",header =FALSE,sep="", comment="#")

## Estimate parameters & create files of summaries reg-est-S and reg-loci-est-S.
estimate_IMc(paramfile="est", datafile="reg", regfile="regions", locifile="loci");

# Clean up the directory.
unlink(c("sim", "est", "regions", "loci", "reg1", "reg2", "reg3", "reg4","reg-S","reg-loci-S","reg-loci-est-S","reg-est-S","reg-lik-0"))
}

\keyword{package}
