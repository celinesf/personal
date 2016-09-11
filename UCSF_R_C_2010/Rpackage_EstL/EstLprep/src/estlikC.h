/*! \file estlikC.h
  \brief Header specific to estlikC.c.

  - main() function
  -  cmd line recovery
  -  parameter recovery
  -  recovers statistics to consider
  - ARG and data simulations
  - Covariance matrix 
  - calculation of composite-likelihood
  - IO error and check for celine
  - calculs on matrices
*/


#define PI 3.14159265358979323846 //!< Pi constant
#define NLSTATS 15    //!< Constant for 4 branch types: L1, L2, Ls, Lf 
#define NPARAM 7 //!< # parameters in the model: theta_1, 2, A, T, M, T_c, M_c
unsigned maxsites = SITESINC ;//!< Min number of segsites (i.e., max memory assignment until now)




/* int main( */
/*          int argc,//!< # of word in cmd line */
/*          char *argv[] //!< cmd line */
/*          ); */

/******** CMD LINE RECOVERY *************/
/***/
/*! \brief Parameter recovery from cmd line 

  - Init all matrices and arrays in structure pars 
  - recover all tags from cmd line
  - Calls argcheck_estL(), argcheck2(), addtoelist() and usage_estL()
  - Called by main()
*/
void getpars(
             int argc,//!< # of word in cmd line
             char *argv[] //!< cmd line
);

/***/
/*! \brief cmd line Word Check for rho characterization 

  - -1 and -2 allowed specific to recombination parameter in cmd line
  - Called by getparv()
*/
void argcheck2( 
               char* argv//!< word in cmd line 
                ) ;

/***/
/*! \brief cmd line Word Check 

 -  Checks that arg is not "-" before one expects a flag 
 - Called by getpars()
*/
void argcheck_estL( 
              int arg, //!< arg/word #
              int argc, //!< tot # of words in cmd line
              char *argv[] //!< cmd line
               );

/***/
/*!\brief List usage of options and tags.*/
/*!
  - Called by getpars(), read_ireg() and get_stat()

  usage: ./estlikeC rho 

   - The rec. rate info must be specified by either:
    - rho (rho=4N1c) 
    - 1 (w=rho_bp in info_regions) 
    - 2 mu (w=rc in info_regions)
    - -1 lambda (c/mu~exp(1/lambda) 
    - -2 nu sigma (c/mu~normal(nu,sigma) 
  Required options: 
   - -p theta1 tehat2 thetaA Ts Tc M Mc (list of parameters to estimate the likelihood)
    - -pf filename  (file name with the sets of parameters, filename=param_info by default)
   - -r filename  (file name with information on the regions, filename=info_regions by default)

  - More options: 
    - -G howmany (=2 by default. The number of simulated data sets per sets of parameters to estimate the (composit-likelihood))
    - -R nregions (R=1 by default)
    - -l filename  (file name with information on the loci for multi-locus regions, filename=info_loci by default)
    - -s s1 s2 s3 (seeds for ramdome generator)
    - -o filename (name of file with max likelihood at checkpoint step for comparison: -2(log(L)-Log(Max))<10)
    - -m m (default m=0: use mean of regions-specific information do define hypothtical region. m=1 to take into account region-specific information. )
    - -c checkpoint (default=20000 step)
    - -h (for usage information)
 */
int usage_estL(); 

/******** MODEL/ESTIMATION PARAM RECOVERY *************/
/***/
/*! \brief Read list of 7 model parameter values to estimate 

  Recovers from a line in a file or from already specified in cmd line
  - theta1 tehat2 theatA T M Tc Mc
 and updates the ms specific parameters from those model parameter values.
 call addtoelist()
*/
void getparv(
             char* line//!< string of parameter values
             ) ;

/***/
/*! \brief Recovery informations on regions and data 

  - Reads regions info and stats from the file ppfreg
  - Keeps max or mean of info depending on 
  - Returns if everything is OK
  Called main()
*/
int read_ireg( 
              FILE* ppfreg,//!< File of info on regions
              int * pnr //!< pointer on R: the number of regions
               );

/***/
/*! \brief Update region-specific information with new set of parameters 

  Called by main()
*/
void getreg_rho(
                double *phypo, //!< array of param for hypothetical region
                int count //!< 1 simulation (init) or not
                ); 

/***/
/*! \brief Free memory allocation for the model and ms parametrers 

  Called by main()
*/
void free_pars();

/********************* STATS CONSIDERED *******************/
/***/
/*! \brief Create the list of 59 possible summary statistic 

  Called by main()
*/
char ** stat_name();

/***/
/*! \brief Recover the statistics to consider and the data 

  returns boolean if data are OK
*/
int get_stat( 
             FILE *ppfreg //!< File of info on regions and data
              );


/********************* ARG AND SAMPLE SIMULATION *******************/
/***/
/*! \brief Update ms parameters for simulation 

  Called by main()
*/
void get_hypo(
              double *phypo //!< array of param for hypothetical region
              );

/*! \brief ARG and sample genrator 

  - Calls biggerlist_estL(),
  - Calls computelstat(), segtre_mig(), prtree(), mnial2(), make_gametes2(), ttime(), poisso() and locate() in streetC.c
  - Returns segsites # S
  - Called by main()
*/
int gensam_estL(
                  char **list, //!< haplotype data
                  double* pposit //!< array of seg sites positions
                  );

/***/
/*! \brief Memory reallocation for the list of haplotypes

  - Used when S>maxsites
  - Called by gensam_estL()
*/
int biggerlist_estL(
               int nsam,  //!< number of chromosomes sampled
               char **list //!< list of haplotypes
               );

/********************* COVARIANCE MATRIX *******************/
/***/
/*! \brief Initialize memory alloc for cov marix and probabilities 

   Called in main()
*/
void init_cov();

/***/
/*! \brief Calculates covariance matrix

  - Calls function getstatS2
  - Called by main()
*/
void get_cov(
             int segr, //!< # segsites
             char ** plistr, //!< haplotype data
             int *pos//!< position os segsites in bp
             );
/***/
/*! \brief Fill up array of statistics

  called by get_cov()
 */
void fill_pstat(
                int nstat,//!< stat # considered
                double val//!< Value of the statistics in simulated data
);

/***/
/*! \brief ReInitialize values for cov marix and probabilities 

  Called by main() 
*/
void reinit_cov();

/***/
/*! \brief Frees covariance matrix and region info

  Called in main()
*/
void free_cov(
              int pnr //!< tot region # R
              );

/********************* LIK ESTIMATION *******************/
/***/
/*! \brief  Computes the probability of the statistics given THETA and an ARG. 

  Called by main() 
*/
void getprobdata( 
                 double *phypo,//!< array of param for hypothetical region
                 int nr, //!< region number considered
                 int *ok//!< boolean if prob(S|THETA) is ok for this region
);

/***/
/*! \brief Calculates p(S|g,Theta) 

  Called by main()
*/
void get_liks(
              int ok //!< Booplean if prob OK
);

/***/
/*! \brief Calculates the composite-likelihood

  - returns log(LIK)
  - Called by main()
*/
long double get_tot_lik( 
                        clock_t pstartl, //!< time I started the main function
                        int nregarg, //!< # ARG and simulated data set
                        int init  //!< 1 if end of complik, 2 if when check maxlik
);

/***/
/*! \brief Compare likelihoods from max file

  - returns log(LIK(MLE_sofar))
  - Called by main()
 */
long double comp_lik();

/****************** PRINT I/O *********************/
/***/
/*! \brief Output ms results for a simulation 

  -  Only for celine check
  - Called by main()
*/
void print_haplo(
                 int nseg,//!< # segsites
                 double * pposit,//!<  array of seg sites positions
                 char ** plist//!< haplotype data
);

/***/
/*! \brief Error output 
  - 0 param file issue
  - 1 region file exist
  - 2 stats/region file issue
  - 3 region info/recombination not ok
*/
void print_error(
                 int ne,//!< error #
                 int i ,//!< val i
                 int j//!<  val j
                 );



/******************* matrices operations *****************/

/***/
/*! \brief Matrix Inverse 

  - returns matrix inv
  - Called by main() 
*/
void Inverse(
             double **a,//!<matrix to invert
             double **inv,//!< returned inverse matrix
             int n //!<# statistics in cov matrix
);

/***/
/*! \brief Matrix cofactor

  -  returns matrix b
  - Called by function Inverse()  
*/
void CoFactor(
              double **a, //!<matrix to get cofactor
              int n,//!<# statistics in cov matrix
              double **b//!< returned matrix of cofactors
              );

/***/
/*! \brief Matrix Transpose 

 - returns t(M)
 - Called by function Inverse() 
*/
void Transpose(
               double **a,//!< matrix to transpose
               double **m,//!< returned Transpose
               int n//!<# statistics in cov matrix
               );

/***/
/*! \brief Matrix determinant

 -  returns the matrix determinant
 - Called by functions Inverse() and CoFactor()
*/
double Determinant(
                   double **a ,//!< matrix
                   int n//!< matrix size
                   );


