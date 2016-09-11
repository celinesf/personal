/*! \file ms.h
  \brief File of struct declaration for ms machinary.

  - define all the structures of ms as well as of estlik
  - Lists the functions from msfun.c
  - called by msC.c, msR.c, streecC.c, streecR.c and estlikC.c
*/

/*** object/ structure definitions ***/
#define SITESINC 10 //!< min #segsites
double *posit ;//!< array of positions of segsites

/***/
/*! \brief ms list of events */
struct devent {
  double time; //!< event time
  int popi;//!< # pop1
  int popj;//!< # pop2
  double paramv;//!< param value
  double **mat ;//!< matrix of Migration rate
  char detype ;//! type of event (see usage()).
  struct devent *nextde;//!< pointer to next event
} ;

/***/
/*! \brief ms optional param */
struct c_params {
  int npop;//!< # populations in model
  int nsam;//!< # # chromosomes
  int *config;//!< # chr in each pop
  double **mig_mat;//!< pairwise migration rates
  double r;//!< rho rec rate
  int nsites;//!< # recombining sites
  double f;//!< # gene conversion rate
  double track_len;//!< # lenght of gene conversion track
  double *size;//!< matrix of ration Ni/N1
  double *alphag;//!< matrix of growth rates
  struct devent *deventlist ;//!< pointer on event
} ;

/***/
/*! \brief ms basic param */
struct m_params {
  double theta;//!< mutation rate
  int segsitesin;//!< Fixed # seg sites
  int treeflag;//!< Flag to output ARG
  int timeflag;//!< Frag to output times in ARG
  int mfreq;//!< Min allele freq considered
} ;

/***/
/*!  \brief info for statistics and composit-lik calculation.
 
  - called by estlikC.c.
*/
struct stats_lik
{
  double *type;//!< 0-anc/maf 1-un/phased 2 3 filter LD in kb
  int ns;//!< tot # of stats
  int nsp; //!< # of other freq/stats
  int nss; //!< # S stats
  int *wstat;//!< ns values: stats number+25[# info in ireg] found in liststat
  char ** liststat;//!< string name of stats
  long double *prob;//!< 0prob_param, 1prog_gen, 2prog_reg, 3sumu, 4 sumu2
  double *infoperf;  //!< 0 # total genealogies generated, 1 # Genealogies with P(D|theta)>0, 2 # gen with P(D|theta)=0, 3 # regions with P(D|theta)>0,     4 # regions with P(D|theta)=0 
  double **lstat;//!< lenght L for 10 possible statistics
  double ** pstat;//!< 0val,1 ok?,2 nok total,3 sum,4 sum2,5+ns*2 n-s1,5+ns*2+1 sum-s1, n-s2, sim-s2....
  long double ** cov; //!< covariance matrix
  long double ** icov; //!< covariance matrix
  long double ** mstat;//!< matrix of val-mu
  int ** ncov; //!< num of comparison covariance matrix
  char * maxfile;//!< name of file with MLE so far
  int checkmax;//!< step # to check for maxfile
  int meanregion;//!< 0: mean/R, 1: variation algorithm
};

/***/
/*! \brief Info for region and param files. 

  - called by estlikC.c. */
struct t_params {

  char *fparv; //!< file of parameter values  for the model
  char *freg;  //!< file of region info + statistics
  char *floc;  //!< file of loci info + statistics
  double *paramv;//!< set of parameter consider for composite likelihood
  int pf;//!< pf==1: param from a file. pf==0: param from cmd line
  int howmany;//!< G: # of ARG/data sets to simulate
  int nregions;//!< R: # regions
  double *rho;//!< info on rho0-2, 3:rho value
  double **ireg;/*!< 0r 1x 2v 3w 4n1 5n2 6zs 7ze 8m1 9m2 10Y 11S 
                //!< +1 11Z 12xvZ 13theta 14rho 15rho/(Z-1) 16rho/omega 17xvZ/xvZ_t 18rho_n(Z-1)/rho(Z_n-1)
                //!< 19rho_n(Z-1)xvZ/rho(Z_n-1)xvZ_t */
  int **nsam;//!< list of n1 n2 combo in R regions
  int totsam;//!< # of combination n1 n2
} ;

/***/
/*! \brief Parameters declaration*/
struct params { 
  struct c_params cp;//!< optional ms param
  struct m_params mp;//!< Basic ms param
  struct t_params tp;//!< Region & model param
  struct stats_lik sp;//!< Statistics and likelihoo param
  int commandlineseedflag ;//!< seeds from cmd line
  int *tableseeds;//!< random seed values
} ;
struct params pars ;//!< parameters values for the program
	
/***/
/*! \brief Node of AGR */
struct node{
  int lpop1;  //!< # lineages from pop1 as descent
  int lpop2; //!< # lineages from pop2 as descent
  int abv;//!< ancestral node
  int ndes;//!< ? used when min frequency specified?
  float time;//!< time of node/ coalecent
};


/********************FUNCTIONS***********************/

/*** Functions from ms.c ***/

/* /\***\/ */
/* /\*! \brief cmd line Word Check  */

/*   -  Checks that arg is not "-" before one expects a flag  */
/*   - Called in getpars() in msC.c and msR.c */
/* *\/ */
/* void argcheck(  */
/*               int arg, //!< arg/word # */
/*               int argc, //!< tot # of words in cmd line */
/*               char *argv[] //!< cmd line */
/*                ); */


/*** Functions from msfunc.c ***/
/***/
/*! \brief Add an event to the list

 - Order the events by their times
 - Called by getparv() in estlikC.c and getparv() in msC.c and msR.c
 - In msfunc.c 
*/
void addtoelist( 
                struct devent *pt, //!< pointer on event to add
                struct devent *elist //!< List of events
                 );

/***/
/*! \brief Free the memory allocated in the list of events 

  - Called by main() in msC.c, msR.c and estlikC.c
*/
void free_eventlist( 
                    struct devent *pt,//!< List of of events
                    int npop//!< # populations in the model
                     );

/***/
/*! \brief Fill info for node.ndes 

  Not sure what it's for/ part of original ms
*/
void ndes_setup(
                struct node *ptree,//!< pointer of the ARG 
                int nsam//!< # of chromosomes in the sample
                );

/***/
/*! \brief Randomly choose the position of segsites along a locus */
int locate(
           int n,//!< # of chromosomes
           double beg,//!< start site of a non-recombing segment*1/(Z-1)
           double len,//!< length of the non-recombing segment*1/(Z-1)
           double *ptr//!< Array of segsite positions
           );

