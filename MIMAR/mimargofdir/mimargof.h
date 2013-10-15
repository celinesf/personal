/****************************************************************************
 *  mimargof.h - Description    
 *  Date: 03/07/07   
 *  Author: Celine Becquet

   -------------------------
   This file contains the structures used in  the rest of the program.
   Some of those and attributes are similar to ms program, others are specific
   to MIMAR or MIMARGOF.
   *
   ************************************** Changes
 
   *--- Changed 18th March 2010 ---*
  Add objects for specifying seeds in the command line.
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
  // char *fout;               // output file name
  char *fin;                   // data file name
  int sumout;                  // Type of output
  int nsim;                    // # of simulations per set of parameters
  /* for the following arrays,
     0=theta1, 1=rho], 2=theta2, 3=time~U,4=thetaa, 5=mig_rate i to j, 6=Mji
  */
  double **uniform;
  double *oldest;              // Recorded parameters from previous step
  double *newest;              // generate new parameters for next step
  double *newparam;            // Parametersspecific for each locus
  double maxval[7];            // Maximum Value in histogram of 1000 bins
  double minval[7];            // Maximum Value in histogram of 1000 bins
  /* For the following arrays, Sum over samples of S statistics and mean Fst, Hw and TD:
     0 S1, 1 S2, 2 Ss, 3 Sf, 4Fst; 5 Hw1, 6 Hw2, 7 Td1, 8 Td2, 9p1,10p2
  */
  double  SiFst[11];           // Mean values for Goodness of fit test
  double  lSiFst[11];          // Sample specific values
  double  sSiFst[11];          // nsim Simulations specific values
};

/*-----------------------------*
 *  Structure Main parameters  |
 *-----------------------------*/
struct m_params
{
  double mu;                   // Per bp mutation rate
  // double g;                 // Number of years per generation
  double theta;                // Theta1 specific for each locus
  double thetafix;             // Theta1 fixed for all loci 
  /* int segsitesin;           // # seg sites given by user - can not be used
   */
  int treeflag;                 // Gene tree
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
  int S[4];                    // population specific S and freq 
  /* For the following arrays: locus specific values
     0 nsam1, 1 Hw1, 2 nsam2, 3 Hw2, 4 nsam total, 5 Hb, 6 fst,7 7p1,8 p2 for locus
  */
  double H[9];                 // Final value for a locus
  double tpH[9];               // Temporary array, record before check it works
};


/*------------------------*
 *  Structure Parameters  |
 *------------------------*
 | List all the parameter values, as given in argument
 | by the user.
*/
struct params
{ 
  struct locus *lp;            // Locus specific parameters
  struct c_params cp;          // Complementary parameters
  struct m_params mp;          // Main parameter
  /* Celine changed 03/18/2010 */
  int commandlineseedflag ;
  int *tableseeds;
  /*/////*/
};
