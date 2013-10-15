/*! \file msonly.h
  \brief Functions in msC.c and msR.c only. 
*/

double segfac ;//!< ms specific when S fixed
int count;//!< # independent samples
int ntbs;//!< # Used when reading cmd line from file

/***/
/*! \brief cmd line Word Check 

  -  Checks that arg is not "-" before one expects a flag 
  - Called in getpars() in msC.c and msR.c
*/
void argcheck( 
              int arg, //!< arg/word #
              int argc, //!< tot # of words in cmd line
              char *argv[] //!< cmd line
               );


/***/
/*! \brief Memory reallocation for the list of haplotypes 

  Used when S>maxsites
*/
int biggerlist(
               int nsam,  //!< number of chromosomes sampled
               char **list //!< list of haplotypes
               );

/***/
/*! \brief Memory allocation for haplotype list.

  Called by main() in msC.c and msMain() in msR.c
 */
char ** cmatrix(
                int nsam ,//!< # of chromosomes
                int len//!< #number of seg sites
);

/*! \brief ARG and sample genrator */
int gensam( 
           char **list, //!< haplotype data
           double *pprobss, //!< arrays of prob of observing S_fixed
           double *ptmrca,  //!< tMRCA when -L used
           double *pttot //!< Sum of time of branches when -L used
) ;


/*! \brief Parameter recovery from cmd line */
void getpars(
             int argc,//!< # of word in cmd line
             char *argv[], //!< cmd line
             int *phowmany //!< # of independent sample to generate
);

/*!\brief List usage of options and tags.*/
	
/*! 
  usage: ./msC nsam howmany 
  - -t theta   (this option and/or the next must be used. Theta = 4*N0*u)
  - -s segsites   ( fixed number of segregating sites)
  - -T          (Output gene tree.)
  - -F minfreq     Output only sites with freq of minor allele >= minfreq.
  - -r rho nsites     (rho here is 4Nc)
  fprintf(stderr,"\t\t -c f track_len   (f = ratio of conversion rate to rec rate. tracklen is mean length.) 
  fprintf(stderr,"\t\t\t if rho = 0.,  f = 4*N0*g, with g the gene conversion rate.
  - -G alpha  ( N(t) = N0*exp(-alpha*t) .  alpha = -log(Np/Nr)/t
  - -I npop n1 n2 ... [mig_rate] (all elements of mig matrix set to mig_rate/(npop-1) 
    - -m i j m_ij    (i,j-th element of mig matrix set to m_ij.)
    - -ma m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)
    - -n i size_i   (popi has size set to size_i*N0 
    - -g i alpha_i  (If used must appear after -M option.)
  The following options modify parameters at the time 't' specified as the first argument:
  - -eG t alpha  (Modify growth rate of all pop's.)
  - -eg t i alpha_i  (Modify growth rate of pop i.) 
  - -eM t mig_rate   (Modify the mig matrix so all elements are mig_rate/(npop-1)
  - -em t i j m_ij    (i,j-th element of mig matrix set to m_ij at time t )
  - -ema t npop  m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)
  - -eN t size  (Modify pop sizes. New sizes = size*N0 ) 
  - -en t i size_i  (Modify pop size of pop i.  New size of popi = size_i*N0 .)
  - -es t i proportion  (Split: pop i -> pop-i + pop-npop, npop increases by 1.
    - proportion is probability that each lineage stays in pop-i. (p, 1-p are admixt. proport.
    - Size of pop npop is set to N0 and alpha = 0.0 , size and alpha of pop i are unchanged.
  - -ej t i j   ( Join lineages in pop i and pop j into pop j
    - size, alpha and M are unchanged.
  - -f filename     ( Read command line arguments from file filename.)
*/
int usage();
