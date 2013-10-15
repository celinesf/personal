/*************** mimar.h ****************************************************
 * 
 *  This file contains the structures used in the rest of the program.
 *  Some of those and attributes are similar to ms program, others are specific
 *  to MIMAR.
 *
 ************************************** Changes
 
 *--- Changed 18th March 2010 ---*
  Add objects for specifying seeds in the command line.
  -----------------------------

  *--- Changed 5nd Nov. 2008 ---*
  Added objects in mcmc to relaunch MIMAR from an interrupted run.
  -----------------------------

  *--- Changed the 18th Sept. 2007 ---*
  Added lp.name.
  -----------------------------
  *
  ***************************************************************************/

/*---------------------------------------*
 *  Structure of Past demographic event  |
 *---------------------------------------*/
struct devent
{
  double time;                 // Coalescent Time
  double timep;                // Proportion of split time given by user
  int popi;                    // # chromo in pop i
  int popj;                    // # chromo in pop j
  double paramv;               // Parameter of pop change (rate growth, increase, decrease, migration...) required by coalescent
  double paramvf;              // Fixed by user
  double **mat ;               // Matrix of migration parameters (-e -ma)
  char detype ;                // Type of event G, g, N, n, M, m, ma, s, j - Needed by coalescent
  char detypef ;               // Type of event G, g, N(b), n, M, m, ma, s, j(h) - given by user
  struct devent *nextde;       // Pointer to the next event
} ;

/*--------------------------------------*
 *  Structure Complementary parameters  |
 *--------------------------------------*/
struct c_params
{
  int npop;                    // # populations
  int nsam;                    // # chromo in the sample = sample size
  int *config;                 // Record # chromo per pop in the sample
  double **mig_mat;            // Migration matrix (-I -ma)
  double r;                    // locus specific Recombination rate
  double rfix;                 // Recombination rate
  int nsites;                  // # seg sites for recombination
  double f;                    // convserion rate
  double track_len;            // Regarding gene conversion
  double *size;                // Population size (-n)
  double *alphag;              // # alpha= growth rate
  struct devent *deventlist ;  // Pointer on a list of events
  struct devent **listevent ;  // Pointer on a list of events
  char *fout;                  // output file name
  char *fin;                   // data file name
  int ngen;                    // genealogies generated per locus for a set of parameters
  int sumout;                  // Type of output
  int interfile;               // # of steps between outputs
  /* for the following arrays,
     0=theta1, 1=rho], 2=theta2, 3=time~U,4=thetaa, 5=mig_rate i to j, 6=Mji
  */
  double **uniform;
  double *oldest;              // Recorded parameters from previous step
  double *newest;              // generate new parameters for next step
  double *newparam;            // Parametersspecific for each locus
  double *trueval;             // True value if data from simulated file
  double maxval[7];            // Maximum Value in histogram of 1000 bins
  double minval[7];            // Maximum Value in histogram of 1000 bins
};

/*-----------------------------*
 *  Structure Main parameters  |
 *-----------------------------*/
struct m_params
{
  double mu;                   // Per bp mutation rate
  double theta;                // Theta1 specific for each locus
  double thetafix;             // Theta1 fixed for all loci 
  /* int segsitesin;           // # seg sites given by user - can not be used
   */
  int treeflag;                // Gene tree
};


/*-----------------------------------*
 *  Structure for Locus information  |
 *-----------------------------------*/
struct locus
{
  /* Celine change 09/18/07 */
  char *name;                  // name of locus
  /*/////*/
  int li;                      // lenght of locus in bp
  double xi;                   // Ni=x*iNe (Theta, time and Mig too affected)
  double vi;                   // thetai=xi*vi*theta
  double wi;                   // rhoi=wi*rho (1 for autosomal, .5 for X, 0 for y and mtdna)
  int ni[3];                   // sample size: 0=total, 1=pop1, 2=pop2
  int Si[4];                   // S statistics 0=S1 1=S2 2=Ss 3=Sf;
};

/*-------------------------------*
 *  Structure for MCMC analysis  |
 *-------------------------------*/
struct mcmc_params
{
  /* Celine changed 11/05/2008*/ // Relaunch mimar from interrupted runs
  int relaunch;                // 1 if relaunching MIMAR from previous input file, 0 otherwise
  char *fpost;                 // file name for posterior distribution started from
  double *start;               // param starting point for relauncinh mimar
  /*/////*/
  int inprior;                 // 1 if all parameters are within their priors, 0 otherwise
  char type;                   // Type of criteria to stop MCMC: #minutss or #steps
  double burnin;               // time in min or # chains for burnin
  long double lenght;          // time in min or # chains before stop MCMC
  int interrec;                // # steps interval between recording 
  double *sigma;               // Variance for update distributions
  /* usually not used */
  int simul;                   // 0 if parameters updated simulatenously, 1 otherwise
  int totparam;                // Total # param to updat if update one parameter at a time
  int *counter;                // # associated for each parameter if update one parameter at a time
  int numparam;                // counter on parameters if update one parameter at a time
  double revjump;              // Prob of reversible jump for migration (see options)

};

/*------------------------*
 *  Structure Parameters  |
 *------------------------*
 | List all the parameter values, as given in argument
 | by the user.
*/
struct params
{ 
  struct mcmc_params mcmcp;    // MCMc analysis parameters
  struct locus *lp;            // Locus specific parameters
  struct c_params cp;          // Complementary parameters
  struct m_params mp;          // Main parameter
  /* Celine changed 03/18/2010 */
  int commandlineseedflag ;
  int *tableseeds;
  /*/////*/
};
