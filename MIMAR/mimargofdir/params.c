/************ param.c - mimargof *******************************************
 *  Contains function that updates the parameter values of the model.
 *
 *  - getpars is as in ms, with changes to include parameters required
 *  for MIMAR.
 *  - changeparams changes priors taking them from the prior distributions.
 *  - changeparamslocus updates the parameters for a particular locus (using the
 *  scalar x, v and w.
 *
 ************************************** Changes

 *--- Changed 18th March 2010 ---*
  Add objects/functions to specify seeds in the command line.
  -----------------------------

  *--- Changed Dec, --, 2009: ---*
  - the locus specific rec rate is now rho_y=rho*w_y*(L-1)
  instead of rho_y=rho*w_y*(L).
  - change calloc(1... by malloc(... .
  - added recombination options.
  - commented useless lines.
  -----------------------------

  *--- Changed March 3rd, 2009: ---*
  Corrected a bug on memory allocation in function getpars.
  Thanks Susan J. Miller for mentionning the bug to me.
  -----------------------------

  *--- Changed 5nd Nov. 2008 ---*
  Usage for -eh and -es uncommented and changed. Other flags commented.
  -----------------------------

  *--- Changed May 29, 2008 ---*
  In the function getpars, I removed several check fags to allow greater
  freedom to the user.
  Thanks Armando Geraldes, for mentioning the problem to me.
  -----------------------------

  *--- Changed the 18th Sept. 2007 ---*
  In function getpars, the option -lf was changed so that the locus
  name can start by any character (before the name could not start by a number).
  -----------------------------
  *
  *****************************************************************/

#include <stdio.h>  // Input output Library
#include <stdlib.h> // File gestion Library
#include <assert.h> // Verify program assertion
#include <math.h>
#include "mimargof.h"

/*------------------------------------*
 *  Declaration of external functions |
 *------------------------------------*
 | Calls for random number generators functions
*/
extern double ran1();
extern double RandomReal(double low, double high);
extern int RandomInteger(int low, int high);
extern double rnd();
extern double unif(double min, double max, int zero);
extern double randexp(double lambda);
extern double RNormal(double mu, double sd);
extern double snorm();

#define NSEGVAL 4 // Constant=number of statistics or Lenghts: L1 L2 Ls Lf
#define NPARAM 5  // 0 theta1, 1 rho or Theta2, 2 Theta2 or tcoal, 3 Tcoal or tgen, 4 t=T4N1, 5 Na. 6 m12, 7 m21
#define NPOP 2    // Fixed number of populations to 2

/*-----------------------------------*
 *  Function chamge params per locus |
 *-----------------------------------*
 |  Updates locus specific parameter values
*/
void changeparamslocus(struct params *pp, int nlocus)
{
  int i=0, j=0;
  double etime=0.0, migr=0.0;
  double tp=0;                                /* Celine changed 12/28/2009 */
  void addtoelist(struct devent *, struct devent *);
  //// From rand1.c ////
  double unif(double, double, int), randexp(double);
  struct devent *pt=NULL, *ptemp=NULL;
  void messageerr(char [1000]);
  pp->cp.deventlist=NULL;
  //--- Theta ---//
  pp->cp.nsites=pp->lp[nlocus].li;            // nsites
  pp->mp.theta=pp->mp.thetafix*pp->cp.nsites*pp->lp[nlocus].xi*pp->lp[nlocus].vi; // theta_y=theta*Z_y*x_y

  //--- Rho ---//                             // if from r prior: newlocus specific rho_y=w_y*Z*c/mu
  if(pp->cp.uniform[0][1]==1)                 // pick ln(r)=(c/mu) from uniform distribution
    pp->cp.r=pp->mp.thetafix*pp->lp[nlocus].wi*(pp->cp.nsites-1)*exp(unif(pp->cp.uniform[1][1], pp->cp.uniform[2][1], 1));
  else if(pp->cp.uniform[0][1]==2)            // pick r=(c/mu) from exp distribution of mean lambda -> rho_y=w_y*theta*c/mu
    pp->cp.r=pp->mp.thetafix*pp->lp[nlocus].wi*(pp->cp.nsites-1)*randexp(pp->cp.uniform[1][1]);
  /* Celine changed 12/28/2009 */
  else if(pp->cp.uniform[0][1]==3)            // pick r=(c/mu) from normal distribution
    {
      tp=RNormal(pp->cp.uniform[1][1], pp->cp.uniform[2][1]);
      while(tp<0)                             // rho>0
        tp=RNormal(pp->cp.uniform[1][1], pp->cp.uniform[2][1]);
      pp->cp.r=pp->mp.thetafix*pp->lp[nlocus].wi*(pp->cp.nsites-1)*tp;
    }
  else if(pp->cp.rfix==1)                     // w_y=rho_y_bp=s_y*4N_1 c_y -> rho_y=w_y*(Z_y-1)
    pp->cp.r=pp->lp[nlocus].wi*(pp->cp.nsites-1);
  else if(pp->cp.rfix==2)                     // w_y=s_y*c_y -> rho_y=w_y*theta/mu*w_y*(Z_y-1)
    pp->cp.r=(pp->mp.thetafix/pp->mp.mu)*pp->lp[nlocus].wi*(pp->cp.nsites-1);
  /*/////*/
  else                                        // if rho=4Noc given rho_y=*rho* Z-1 *w_y
    pp->cp.r=pp->cp.rfix*pp->lp[nlocus].wi*(pp->cp.nsites-1);

  //--- update # chromosomes for each pop at the locus---//
  pp->cp.nsam=pp->lp[nlocus].ni[0];           // Update total number of chromosome
  for(i=0;i<pp->cp.npop;i++)                  // Loop for each pop
    pp->cp.config[i]=pp->lp[nlocus].ni[i+1];

  //--- N2 ---//
  pp->cp.size[1]=pp->cp.newparam[1];          // update b=theta2/theta1

  //--- T ---//
  etime=0;
  if(pp->cp.listevent[0]!=NULL)
    {
      i=0;
      pt=pp->cp.listevent[i];                 // Pointer of events
      while((pt!=NULL)&&(etime==0))           // Loop to find split time
        {
          if((pt->detypef=='N')||(pt->detypef=='j'))
            etime=(double)pp->cp.newparam[2]/pp->lp[nlocus].xi; // T_y=T/x_y
          i++;
          pt=pp->cp.listevent[i];
        }
      i=0;
      pt=pp->cp.listevent[i];
      pt->nextde=NULL;
      while(pt!=NULL)                         // Loop to update events times or paramv
        {
          if((pt->detypef=='N')||(pt->detypef=='j'))
            {                                 // Updates only if event split or ancestral size
              pt->time=etime;
              //--- Na ---//
              if(pt->detype=='N')             // update a=thetaA/theta1
                pt->paramv=(double) pp->cp.newparam[4];
            }
          else                                // Other events
            {
              if(pt->detypef=='h')            // Event 'h' is a split event 'j'
                pt->detype='j';
              if(pt->detypef=='b')            // Event 'b' is a Change of size event 'N'
                pt->detype='N';
              pt->time=pt->timep*etime;       // update t=split time/ * x
            }

          if(pp->cp.deventlist==NULL)         // Point on this new event
            pp->cp.deventlist=pt;
          else if(pt->time<pp->cp.deventlist->time)
            {                                 // Order event in the list
              ptemp=pp->cp.deventlist;
              pp->cp.deventlist=pt;
              pt->nextde=ptemp;
            }
          else                                // Or add to event list
            addtoelist(pt, pp->cp.deventlist);
          i++;
          pt=pp->cp.listevent[i];
        }
      //---- Loop to check weird model events ----//
      pt=pp->cp.deventlist;
      char c='\0';
      double evt=0;
      while(pt!=NULL)                         // Loop to update events times or paramv
        {
          if(pt->time==evt)
            {
              if((pt->detype=='j')&&(c=='s')) messageerr("Option -ej needs to be before -es.\n");
              if((pt->detype=='j')&&(c=='m')) messageerr("Option -ej needs to be before -em.\n");
              if((pt->detype=='j')&&(c=='M')) messageerr("Option -ej needs to be before -eM.\n");
            }
          c=pt->detype;
          evt=pt->time;
          pt=pt->nextde;
        }
    }
  //--- Sym Migration ---//
  // unif: 1 for test & uniform distribution
  //       2 for uniform prior on M12=4N1m12, 5 ln(M)~u
  //       3 for uniform prior on m12 (IM comparison)
  migr=0;
  if((((pp->cp.uniform[0][5]==2)||(pp->cp.uniform[0][5]==1)||(pp->cp.uniform[0][5]==5))&&((pp->cp.uniform[0][6]==2)||(pp->cp.uniform[0][6]==1)||(pp->cp.uniform[0][6]==5)))&&(pp->cp.uniform[2][6]==0))
    {                                         // locus specific Migr=mig* inheritance scalar
      migr=pp->cp.newparam[5]*pp->lp[nlocus].xi;
      for(i=0;i<pp->cp.npop;i++)              // loop: migr information in matrix for all pop
        for(j=0;j<pp->cp.npop;j++) pp->cp.mig_mat[i][j]=migr/(pp->cp.npop-1);
      for(i=0;i<pp->cp.npop;i++) pp->cp.mig_mat[i][i]=migr;
    }// symetrical cases & no test
  else if(((pp->cp.uniform[0][5]==3)&&((pp->cp.uniform[0][6]==3)))&&(pp->cp.uniform[2][6]==0))
    {                                         // M=theta * m/mu
      migr=pp->cp.newparam[5]*pp->mp.thetafix*pp->lp[nlocus].xi/pp->mp.mu;
      for(i=0;i<pp->cp.npop;i++)              // loop: migr information in matrix for all pop
        for(j=0;j<pp->cp.npop;j++) pp->cp.mig_mat[i][j]=migr/(pp->cp.npop-1);
      for(i=0;i<pp->cp.npop;i++) pp->cp.mig_mat[i][i]=migr;
    }// symetrical case & test
  else
    {
      //-- M12 --//
      if((pp->cp.uniform[0][5]==1)||(pp->cp.uniform[0][5]==5)||(pp->cp.uniform[0][5]==2)||
         (pp->cp.uniform[0][5]<=0))           // M12 test or no test cases
        pp->cp.mig_mat[0][1]=pp->cp.newparam[5]*pp->lp[nlocus].xi;
      else if((pp->cp.uniform[0][5]>=3)&&
              (pp->cp.uniform[0][5]!=5))      // m12 no test case - M12=m12*theta/mu
        pp->cp.mig_mat[0][1]=pp->cp.newparam[5]*pp->mp.thetafix*pp->lp[nlocus].xi/pp->mp.mu;

      //-- M21 --//
      if((pp->cp.uniform[0][6]==1)||(pp->cp.uniform[0][6]==2)||(pp->cp.uniform[0][6]==5)||
         (pp->cp.uniform[0][6]<=0))           // M21 test or no test cases
        pp->cp.mig_mat[1][0]=pp->cp.newparam[6]*pp->lp[nlocus].xi;
      else if((pp->cp.uniform[0][6]>=3)&&
              (pp->cp.uniform[0][6]!=5))      // pick m21 no test case - M21=m21*theta/mu
        pp->cp.mig_mat[1][0]=pp->cp.newparam[6]*pp->mp.thetafix*pp->lp[nlocus].xi/pp->mp.mu;
      for(i=0;i<pp->cp.npop;i++)              // loop: migr information in matrix for all pop
        {
          pp->cp.mig_mat[i][i]=0.0;
          for(j=0;j<pp->cp.npop;j++)
            if(j!=i) pp->cp.mig_mat[i][i]+=pp->cp.mig_mat[i][j];
        }
    }// Non symetrical cases
}// Locus specific parameters


/*-------------------------------*
 *  Function Prior distributions |
 *-------------------------------*
 |  Proposes parameter values from prior distributions.
 |  Called at initial step of MCMC
*/
void changeparams(struct params *pp)
{
  double unif(double, double, int), randexp(double);
  //--- Theta ---//
  if(pp->cp.uniform[0][0]==1)                 // theta in bp ~u[]
    pp->cp.oldest[0]=pp->cp.newest[0]=unif(pp->cp.uniform[1][0], pp->cp.uniform[2][0], 0);

  //--- Theta2 ---//
  if(pp->cp.uniform[0][2]==2)                              // ln(b)~u[]
    pp->cp.oldest[1]=pp->cp.newest[1]=(double)exp(unif(pp->cp.uniform[1][2], pp->cp.uniform[2][2], 0));
  if((pp->cp.uniform[0][2]==1)||(pp->cp.uniform[0][2]==3)) // theta2~u[]
    pp->cp.oldest[1]=pp->cp.newest[1]=(double)unif(pp->cp.uniform[1][2], pp->cp.uniform[2][2], 0);
  else if(pp->cp.uniform[0][2]==4)                         // theta2 fixed
    pp->cp.oldest[1]=pp->cp.newest[1]=(double) pp->cp.uniform[1][2];
  else if((pp->cp.uniform[0][2]<=0)&&(pp->cp.sumout>=2))   // b=1, theta2 fixed=Theta1
    pp->cp.oldest[1]=pp->cp.newest[1]=pp->cp.newest[0];

  //--- T ---//
  if(pp->cp.listevent[0]!=NULL)
    {
      if((pp->cp.uniform[0][3]==1))           // Tcoal ~ u[]
        pp->cp.oldest[2]=pp->cp.newest[2]=(double) unif(pp->cp.uniform[1][3], pp->cp.uniform[2][3], 0);
      if(pp->cp.uniform[0][3]==3)             // Tgen ~ u[]
        pp->cp.oldest[3]=pp->cp.newest[3]=(double) unif(pp->cp.uniform[1][3], pp->cp.uniform[2][3], 0);
      else if(pp->cp.uniform[0][3]==4)        // Tgen fixed
        pp->cp.oldest[3]=pp->cp.newest[3]=(double) pp->cp.uniform[1][3];
      if(pp->cp.sumout>=2)                    // Tgen fixed, compute Tcoal=Tgen*mu/theta
        pp->cp.oldest[2]=pp->cp.newest[2]=pp->cp.newest[3]*pp->mp.mu/pp->cp.newest[0];
      else                                    // Tcoal fixed, compute Tgen=Tcoal*that/mu
        pp->cp.oldest[3]=pp->cp.newest[3]=pp->cp.newest[2]*pp->cp.newest[0]/pp->mp.mu;

      //--- ThetaA ---//
      if(pp->cp.uniform[0][4]==2)                              // ln(a)~U[]
        pp->cp.oldest[4]=pp->cp.newest[4]=(double)exp(unif(pp->cp.uniform[1][4], pp->cp.uniform[2][4], 0));
      if((pp->cp.uniform[0][4]==1)||(pp->cp.uniform[0][4]==3)) // ThetaA~u[]
        pp->cp.oldest[4]=pp->cp.newest[4]=(double)unif(pp->cp.uniform[1][4], pp->cp.uniform[2][4], 0);
      else if(pp->cp.uniform[0][4]==4)                         // thetaA fixed
        pp->cp.oldest[4]=pp->cp.newest[4]=(double) pp->cp.uniform[1][4];
      else if((pp->cp.uniform[0][4]<=0)&&(pp->cp.sumout>=2))   // a=1, thetaA fixed=Theta1
        pp->cp.oldest[4]=pp->cp.newest[4]=pp->cp.newest[0];
    }
  //--- Migration rates ---//
  if((pp->cp.uniform[0][NPARAM]>0)||(pp->cp.uniform[0][NPARAM+1]>0))
    {
      if(((pp->cp.uniform[0][NPARAM]==1)&&((pp->cp.uniform[0][NPARAM+1]==1)))&&
         (pp->cp.uniform[2][NPARAM+1]==0))      // Test case M~u[]
        pp->cp.newest[5]=pp->cp.oldest[5]=pp->cp.newest[6]=pp->cp.oldest[6]=unif(pp->cp.uniform[1][NPARAM], pp->cp.uniform[2][NPARAM], 2);
      else if((((pp->cp.uniform[0][NPARAM]==2)||(pp->cp.uniform[0][NPARAM]==3))&&((pp->cp.uniform[0][NPARAM+1]==2)||(pp->cp.uniform[0][NPARAM+1]==3)))&&
              (pp->cp.uniform[2][NPARAM+1]==0)) // M~u[]
        pp->cp.newest[5]=pp->cp.oldest[5]=pp->cp.newest[6]=pp->cp.oldest[6]=unif(pp->cp.uniform[1][NPARAM], pp->cp.uniform[2][NPARAM], 1);
      else if(((pp->cp.uniform[0][NPARAM]==5)&&(pp->cp.uniform[0][NPARAM+1]==5))&&
              (pp->cp.uniform[2][NPARAM+1]==0)) // ln(M)~[]
        pp->cp.newest[5]=pp->cp.oldest[5]=pp->cp.newest[6]=pp->cp.oldest[6]=exp(unif(pp->cp.uniform[1][NPARAM], pp->cp.uniform[2][NPARAM], 1));
      else
        {
          //-- M12--//
          if(pp->cp.uniform[0][NPARAM]==1)        // M12~u[] test case
            pp->cp.newest[5]=pp->cp.oldest[5]=unif(pp->cp.uniform[1][NPARAM], pp->cp.uniform[2][NPARAM], 2);
          else if((pp->cp.uniform[0][NPARAM]==3)||
                  (pp->cp.uniform[0][NPARAM]==2)) // M12~u[]
            pp->cp.newest[5]=pp->cp.oldest[5]=unif(pp->cp.uniform[1][NPARAM], pp->cp.uniform[2][NPARAM], 1);
          else if(pp->cp.uniform[0][NPARAM]==5)   // ln(M12)~[]
            pp->cp.newest[5]=pp->cp.oldest[5]=exp(unif(pp->cp.uniform[1][NPARAM], pp->cp.uniform[2][NPARAM], 0));
          else if(pp->cp.uniform[0][NPARAM]==4)   // fixed m12 update M12 when theta1 changes
            pp->cp.newest[5]=pp->cp.oldest[5]=(double) pp->cp.uniform[1][NPARAM];

          //-- M21 --//
          if(pp->cp.uniform[0][NPARAM+1]==1)        // M21~u[] test case
            pp->cp.newest[6]=pp->cp.oldest[6]=unif(pp->cp.uniform[1][NPARAM+1], pp->cp.uniform[2][NPARAM+1], 2);
          else if((pp->cp.uniform[0][NPARAM+1]==3)||
                  (pp->cp.uniform[0][NPARAM+1]==2)) // M21~u[]
            pp->cp.newest[6]=pp->cp.oldest[6]=unif(pp->cp.uniform[1][NPARAM+1], pp->cp.uniform[2][NPARAM+1], 1);
          else if(pp->cp.uniform[0][NPARAM+1]==5)   // ln(M21)~[] - ln(M21)~U[a, b]
            pp->cp.newest[6]=pp->cp.oldest[6]=exp(unif(pp->cp.uniform[1][NPARAM+1], pp->cp.uniform[2][NPARAM+1], 1));
          else if(pp->cp.uniform[0][NPARAM+1]==4)   // fixed m21 update M21 when theta1 changes
            pp->cp.newest[6]=pp->cp.oldest[6]=(double) pp->cp.uniform[1][NPARAM+1];
        }// Non symetrical case
    }// Migration rate
}// Function prior distributiom


/*--------------------------------------*
 *  Function get parameters information |
 *--------------------------------------*
 | Gets the arguments in command line and information about loci from input file.
 | Records information in structure of parameters.
*/
struct params getpars(int argc, char *argv[], int *phowmany)
{
  int arg=0, i=0, j=0, pop=0, argstart=0, npop=0, npop2=0, pop2=0, nevent=0, k=0;
  int atoi(const char*);
  double mij=0.0, migr=0.0, psize=0.0, palpha=0.0;
  /* Celine changed 12/28/2009 */
  double xw=0, nw=0; /*/////*/
  double atof(const char*);
  long double strtold(const char*, char**);
  char ch3='\0', ch4='\0', st[200]="";
  struct params p;
  struct devent *pt=NULL;
  void argcheck(int, int, char *v[]);
  void argcheck2(int, int, char *v[]);
  void format(), usage();
  void messageerr(char [1000]);
  void messageerrs(char [1000], char [1000], char [1000]);
  void messageerrd(char [1000], int, char [1000]);
  void messageerrl(char [1000], double, char [1000]);
  FILE *pf=NULL, *lf=NULL, *fopen(const char *, const char *);
  /* Celine changed 03/18/2010 */
  int commandlineseed( char **, struct params * ), nseeds=0 ;
  p.commandlineseedflag = 0 ;
  /* added celine */
  p.tableseeds = (int *) malloc( (unsigned)(3 *sizeof( int)) );
  /*/////*/

  /* Commented for Celine's use - special case */
  /* if(argc<4) messageerr("Too few command line arguments\n");
     p.cp.nsim=atoi(argv[1]);                    // Number of simulations per set of parameters
     if(p.cp.nsim<=0) messageerr("number of sim <=0. \n");
     *phowmany=atoi(argv[2]);                    // Howmany is # of loci
     *//////
  p.cp.nsim=1;
  if(argc<3) messageerr("Too few command line arguments\n");
  *phowmany=atoi(argv[1]);                    // Howmany is # of loci
  if(*phowmany<=0) messageerr("number of loci <=0. \n");

  //--- Initialization and memory allocation ---//
  //-- ms parameters --//
  /* Commented for Celine's use - special case */
  /*arg=3;
   *//////
  arg=2;
  p.cp.r=p.mp.theta=p.cp.rfix=p.mp.thetafix=p.mp.mu=p.cp.f=0.0;
  p.cp.track_len=0.;
  p.cp.npop=npop=1;
  p.cp.nsites=2;
  p.mp.treeflag=0;
  p.cp.mig_mat=(double**)malloc((unsigned)sizeof(double*));
  p.cp.mig_mat[0]=(double*)malloc((unsigned)sizeof(double)); /*--- Celine changed 03/03/2009---*/ // remove *
  p.cp.mig_mat[0][0]=0.0;
  p.cp.config=(int*)malloc((unsigned)((p.cp.npop+1)*sizeof(int)));
  (p.cp.config)[0]=p.cp.nsam=0;
  p.cp.size=(double*)malloc((unsigned)(p.cp.npop*sizeof(double)));
  (p.cp.size)[0]=1.0;
  p.cp.alphag=(double*)malloc((unsigned)((p.cp.npop)*sizeof(double)));
  (p.cp.alphag)[0]=0.0;
  p.cp.listevent=(struct devent**)malloc((unsigned)(10*sizeof(struct devent*)));
  p.cp.deventlist=NULL;
  nevent=0;
  for(i=0;i<10;i++) p.cp.listevent[i]=NULL;

  //-- MIMAR specific parameters --//
  p.cp.uniform=(double**)malloc((unsigned)(3*sizeof(double*)));
  /*0: tag
    0 theta1: 1 theta1 u[], 
    1 r=rho/theta: 1 ln(r)~u[], 2 r~exp()
    2 theta2: 1 b~u[], 2 ln(b)~u[] (option 0-1), 3 theta2~u[] (option 2-3, 4 b=1 (update theta2, option 1), 5 theta2 given and fixed (update b, option >=2)
    3 T: 1 Tcoal~u[] (option 0-1), Tgen~u[] (option 2-3), 4 Tgen fixed (update Tcoal, option >=2)
    4 theta A: 1 a~u[], 2 ln(a)~u[] (option 0-1), 3 thetaA~u[] (option 2-3), 4 a=1 (update thetaA, option 1), 5 thetaA given and fixed (update a, option >=2)
    5 and 6 M: 1=rev jump(test of migration), 2=uniform, 3 m~u[] (update M option 3), 4 m fixed (update M option 3), 5 ln(M)~u[]
    1: lower limit/mean/fixed value
    2: upper limit or empty
  */
  /**** folowing arrays of size NPARAM+POP : 0 Theta, 1 theta2, 2 Tcoal, 3 Tgen, 4 thetaA, 5 M12, 6 M21 ****/
  p.cp.newparam=(double*)malloc((unsigned)(NPARAM+NPOP)*sizeof(double));  // Coalescent parameters updated foreach locus
  p.cp.oldest=(double*)malloc((unsigned)(NPARAM+NPOP)*sizeof(double));    // estimates recorded at previous step
  p.cp.newest=(double*)malloc((unsigned)(NPARAM+NPOP)*sizeof(double));    // Proposed estimates
  p.lp=(struct locus*)malloc((unsigned)(*phowmany*sizeof(struct locus))); // locus list
  /* p.cp.fout=(char*) malloc((unsigned)((50)*sizeof(char)));             // output file name
     p.cp.fout="";*/
  p.cp.fin="";
  p.cp.sumout=2;
  p.mp.mu=p.mp.theta=p.mp.thetafix=0; // p.mp.g=
  /* sumout: kind of mcmc output
     0 for Theta1, b, tcoal, a, M12, M21 (test mig)
     1 for theta1, theta2, Tcoal, thetaA, M12, M21 (te s t mig)
     2 for Theta1, Theta2, Tgen, tcoal, ThetaA, M12, M 2 1 (priors)
     3 for Theta1, Theta2, Tgen, tcoal (prior), ThetaA, m 12, m21
  */
  for (i=0;i<11;i++)                          // Reset array of values for Goodness of fit test
    p.cp.SiFst[i]=0;

  for (i=0;i<3;i++)
    {
      p.cp.uniform[i]=(double*)calloc(3, (unsigned)(NPARAM+NPOP)*sizeof(double*));
      for(j=0;j<NPARAM+NPOP;j++)
        p.cp.uniform[i][j]=-1.;               // 0.theta1, 1.rho, 2. theta2, 3. T, 4.thetaA, 5.M12, 6.M21
    }
  for(j=0;j<NPARAM+NPOP;j++)
    {
      p.cp.maxval[j]=-1.;
      p.cp.minval[j]=-1;
      p.cp.newparam[j]=0;
      p.cp.oldest[j]=0;
      p.cp.newest[j]=0;
    }
  //--- Loop on other arguments ---//
  while(arg<argc)
    {
      if(argv[arg][0] !='-')
        messageerrs(" argument should be -", argv[arg], " ?\n");
      switch (argv[arg][1])
        {
          case 'f' :                          // read command line from a file
            arg++;
            argcheck(arg, argc, argv);
            pf=fopen(argv[arg], "r");
            if(pf==NULL) {fprintf(stderr, " no parameter file %s\n", argv[arg]); exit(0);}
            arg++;
            argc+=1;
            argv=(char**)malloc((unsigned)((argc)+1)*sizeof(char*));
            argv[arg]=(char*)malloc((unsigned)(20*sizeof(char)));

            argstart=arg;
            while(fscanf(pf, " %s", argv[arg]) !=EOF)
              {
                arg++;
                argc+=1;
                argv=(char**)realloc(argv, (unsigned)argc*sizeof(char*));
                argv[arg]=(char*)malloc((unsigned)(20*sizeof(char)));
              }
            fclose(pf);
            argc-=1;
            arg=argstart;
            break;
          case 'r' :                          // Recombination
            arg++;
            argcheck(arg, argc, argv);
            ch3=argv[arg++][0];
            if(ch3=='u')
              {
                p.cp.uniform[0][1]=1;                 // r=rho/theta: lnr~U[a, b]
                argcheck2(arg, argc, argv);
                p.cp.uniform[1][1]=atof(argv[arg++]); // a
                argcheck2(arg, argc, argv);
                p.cp.uniform[2][1]=atof(argv[arg++]); // b
                if(p.cp.uniform[1][1]>=p.cp.uniform[2][1]) messageerr("with -r u a b option b must be larger than a\n");
              }
            else if(ch3=='e')
              {
                argcheck(arg, argc, argv);
                p.cp.uniform[0][1]=2;                 // r=rho/theta: r~exp(lambda)
                argcheck(arg, argc, argv);
                p.cp.uniform[1][1]=atof(argv[arg++]); // lambda
                if(p.cp.uniform[1][1]<=0) messageerr("with -r e option must specify lambda > 0\n");
              }
            /* Celine changed 12/28/2009 */
            else if(ch3=='n')
              {
                argcheck(arg, argc, argv);
                p.cp.uniform[0][1]=3;                 // r=rho/theta: r~normal(mean, sigma)
                argcheck(arg, argc, argv);
                p.cp.uniform[1][1]=atof(argv[arg++]); // mean
                argcheck2(arg, argc, argv);
                p.cp.uniform[2][1]=atof(argv[arg++]); // sigma
                if(p.cp.uniform[1][1]<=0) messageerr("with -r n option must specify mean > 0\n");
                if(p.cp.uniform[2][1]<=0) messageerr("with -r n option must specify sigma > 0\n");
              }
            /*/////*/
            else p.cp.rfix=atof(argv[arg-1]); // User gives Rho=4N0C here
            break;
          case 'c' :                          // Conversion
            arg++;
            argcheck(arg, argc, argv);
            p.cp.f=atof(argv[arg++]);
            argcheck(arg, argc, argv);
            p.cp.track_len=atof(argv[arg++]);
            if(p.cp.track_len<1.)
              messageerr("with -c option must specify both f and track_len>0\n");
            break;
          case 't' :                          // Theta1
            arg++;
            argcheck(arg, argc, argv);
            ch3=argv[arg++][0];
            if(ch3=='u')
              {
                p.cp.uniform[0][0]=1;                 // Theta1 ~ U[a, b
                argcheck(arg, argc, argv);
                p.cp.uniform[1][0]=atof(argv[arg++]); // a
                argcheck(arg, argc, argv);
                p.cp.uniform[2][0]=atof(argv[arg++]); // b
                if(p.cp.uniform[1][0]>=p.cp.uniform[2][0]) messageerr("with -t u a b option b must be larger than a\n");
                p.cp.maxval[0]=p.cp.uniform[2][0];    // Max and min of histogram of 1000 bins
                p.cp.minval[0]=p.cp.uniform[1][0];
              }
            else                              // Theta1 fixed by user
              {
                p.mp.thetafix=atof(argv[arg-1]);
                p.cp.maxval[0]=p.cp.minval[0]=p.mp.thetafix;
              }
            break;
          case 's' :                          
            arg++;
            argcheck(arg, argc, argv);
            /* Celine changed 03/18/2010 */
            if( argv[arg-1][2] == 'e' ){  /* command line seeds */
              p.commandlineseedflag = 1 ;
              for(i=0;i<3;i++)
                argcheck( arg+i, argc, argv);
              nseeds = commandlineseed(argv+arg ,&p  );
              arg += nseeds ;
            } 
            /* else // segsites - This option can not be used
               p.mp.segsitesin=atoi(argv[arg++]);*/
            /*/////*/
            break;
          case 'T' :                          // output genetree
            p.mp.treeflag=1;
            arg++;
            break;
          case 'm' :                          // Migration
            if(argv[arg][2]=='a')             // matrix of migration rates
              {
                arg++;
                for(pop=0; pop<npop; pop++)
                  for(pop2=0; pop2<npop; pop2++)
                    {
                      argcheck(arg, argc, argv);
                      ch3=argv[arg++][0];
                      if(pop!=pop2)           // Non diagonal therms
                        {
                          if((ch3=='u')||(ch3=='t'))
                            {
                              if(ch3=='t') p.cp.uniform[0][NPARAM+pop]=1;    // tag for Mij=4Nomij~U[a, b] & test & reversible jump
                              if(ch3=='u') p.cp.uniform[0][NPARAM+pop]=2;    // tag for Mij=4Nomij~U[a, b] or mij (option3)
                              argcheck(arg, argc, argv);
                              p.cp.uniform[1][NPARAM+pop]=atof(argv[arg++]); // a
                              argcheck(arg, argc, argv);
                              p.cp.uniform[2][NPARAM+pop]=atof(argv[arg++]); // b
                              if(p.cp.uniform[1][NPARAM+pop]>p.cp.uniform[2][NPARAM+pop]) messageerr("with -ma x u a b option b must be larger than a\n");
                              if((pop==0)&&(p.cp.uniform[2][NPARAM+pop]==0)) messageerr("with -ma x u a b option must use b>0\n");
                              p.cp.maxval[NPARAM+pop]=p.cp.uniform[2][NPARAM+pop]; // Max values for mij
                              p.cp.minval[NPARAM+pop]=p.cp.uniform[1][NPARAM+pop];
                            }
                          else if(ch3=='l')
                            {
                              p.cp.uniform[0][NPARAM+pop]=5;                 // tag for ln(Mij)~U[a1, a2]
                              argcheck2(arg, argc, argv);
                              p.cp.uniform[1][NPARAM+pop]=atof(argv[arg++]); // a
                              argcheck2(arg, argc, argv);
                              p.cp.uniform[2][NPARAM+pop]=atof(argv[arg++]); // b
                              if(p.cp.uniform[1][NPARAM+pop]>p.cp.uniform[2][NPARAM+pop]) messageerr("with -ma x l a b option, b must be larger than a\n");
                              /* Celine changed 12/28/2009 */
                              /*  if((pop==0)&&(p.cp.uniform[2][NPARAM+pop]==0)) messageerr("with -ma x l a b option must use b>0\n");
                               */////
                              p.cp.maxval[NPARAM+pop]=(p.cp.uniform[2][NPARAM+pop]); // Max and min values for histogram
                              p.cp.minval[NPARAM+pop]=(p.cp.uniform[1][NPARAM+pop]);
                            }
                          else                // -ma x m12 m21 x fixed
                            {
                              p.cp.mig_mat[pop][pop2]=atof(argv[arg-1]);
                              p.cp.maxval[NPARAM+pop]=p.cp.minval[NPARAM+pop]=p.cp.mig_mat[pop][pop2];
                            }
                        }
                      else if(ch3!='x') messageerr("with -ma option, the diagonal therms need to be indicated by an 'x'.\n");
                    }
                for(pop=0; pop<npop; pop++)
                  {
                    p.cp.mig_mat[pop][pop]=0.0;
                    for(pop2=0; pop2<npop; pop2++)
                      if(pop2 !=pop) p.cp.mig_mat[pop][pop]+=p.cp.mig_mat[pop][pop2];
                  }
              }
            else                              // migration rates
              {
                arg++;
                argcheck(arg, argc, argv);
                i=atoi(argv[arg++]) -1;
                argcheck(arg, argc, argv);
                j=atoi(argv[arg++]) -1;
                if((j==i)||((i>=p.cp.npop)||(j>=p.cp.npop))||((j<0)||(i<0))) messageerr("with -m i j m_ij option i and j must be different and 1 or 2\n");
                argcheck(arg, argc, argv);
                ch3=argv[arg++][0];
                if((ch3=='u')||(ch3=='t'))    // Mij=4Nomij~U[a, b]
                  {
                    if(ch3=='t') p.cp.uniform[0][NPARAM+i]=1;    // tag for Mij=4Nomij~U[a, b]+test & rev jump
                    if(ch3=='u') p.cp.uniform[0][NPARAM+i]=2;    // tag for Mij=4Nomij~U[a, b] or mij (option3)
                    argcheck(arg, argc, argv);
                    p.cp.uniform[1][NPARAM+i]=atof(argv[arg++]); // a
                    argcheck(arg, argc, argv);
                    p.cp.uniform[2][NPARAM+i]=atof(argv[arg++]); // b
                    if(p.cp.uniform[1][NPARAM+i]>p.cp.uniform[2][NPARAM+i]) messageerr("with -m i j u a b option b must be larger than a\n");
                    if((i==0)&&(p.cp.uniform[2][NPARAM+i]==0)) messageerr("with -m 1 2 u a b option must use b>0\n");
                    p.cp.maxval[NPARAM+i]=p.cp.uniform[2][NPARAM+i]; // values for histogram
                    p.cp.minval[NPARAM+i]=p.cp.uniform[1][NPARAM+i];
                  }
                else if(ch3=='l')
                  {
                    p.cp.uniform[0][NPARAM+i]=5;                 // tag for ln(Mij)~U[a1, a2]
                    argcheck2(arg, argc, argv);
                    p.cp.uniform[1][NPARAM+i]=atof(argv[arg++]); // a
                    argcheck2(arg, argc, argv);
                    p.cp.uniform[2][NPARAM+i]=atof(argv[arg++]); // b
                    if(p.cp.uniform[1][NPARAM+i]>p.cp.uniform[2][NPARAM+i]) messageerr("with -m i j l a b option, b must be larger than a\n");
                    /* Celine changed 12/28/2009 */
                    /* if((i==0)&&(p.cp.uniform[2][NPARAM+i]==0)) messageerr("with -m 1 2 l a b option must use b>0\n");
                       if(p.cp.uniform[2][NPARAM+i]!=0)
                       {*/                                          // Max values for Mig
                    p.cp.maxval[NPARAM+i]=(p.cp.uniform[2][NPARAM+i]);
                    p.cp.minval[NPARAM+i]=(p.cp.uniform[1][NPARAM+i]);                        
                    /* }
                       else
                       p.cp.maxval[NPARAM+i]=p.cp.minval[NPARAM+i]=p.cp.uniform[2][NPARAM+i];
                    *//////
                  }
                else                          // -m i j Mij fixed
                  {
                    mij=atof(argv[arg-1]);
                    p.cp.mig_mat[i][i]+=mij - p.cp.mig_mat[i][j];
                    p.cp.mig_mat[i][j]=mij;
                    p.cp.maxval[NPARAM+i]=p.cp.minval[NPARAM+i]=p.cp.mig_mat[i][j];
                  }
              }
            break;
            /* Celine changed 05/29/2008 */
          case 'M' :
            arg++;
            argcheck(arg, argc, argv);
            ch3=argv[arg][0];
            arg++;
            i=0;
            if((ch3=='u')||(ch3=='t'))        // Mij=4Nomij~U[a, b]
              {
                if(ch3=='t') p.cp.uniform[0][NPARAM+i]=p.cp.uniform[0][NPARAM+i+1]=1; // tag for Mij=4Nomij~U[a, b]+test & rev jump
                if(ch3=='u') p.cp.uniform[0][NPARAM+i]=p.cp.uniform[0][NPARAM+i+1]=2; // tag for Mij=4Nomij~U[a, b] or mij (option3)
                argcheck(arg, argc, argv);
                p.cp.uniform[1][NPARAM+i]=atof(argv[arg++]);               // a
                argcheck(arg, argc, argv);
                p.cp.uniform[2][NPARAM+i]=atof(argv[arg++]);               // b
                p.cp.uniform[2][NPARAM+i+1]=p.cp.uniform[1][NPARAM+i+1]=0; // sym u[0, 0]
                if(p.cp.uniform[1][NPARAM+i]>p.cp.uniform[2][NPARAM+i]) messageerr("with -M u a b option b, must be larger than a\n");
                if((i==0)&&(p.cp.uniform[2][NPARAM+i]==0)) messageerr("with -M u a b option, must use b>0\n");
                p.cp.maxval[NPARAM+i]=p.cp.uniform[2][NPARAM+i];           // values for histogram
                p.cp.minval[NPARAM+i]=p.cp.uniform[1][NPARAM+i];
                p.cp.maxval[NPARAM+i+1]=p.cp.uniform[2][NPARAM+i+1];       // values for histogram
                p.cp.minval[NPARAM+i+1]=p.cp.uniform[1][NPARAM+i+1];
              }
            else if(ch3=='l')
              {
                p.cp.uniform[0][NPARAM+i]=p.cp.uniform[0][NPARAM+i+1]=5;   // tag for ln(Mij)~U[a1, a2]
                argcheck2(arg, argc, argv);
                p.cp.uniform[1][NPARAM+i]=atof(argv[arg++]);               // a
                argcheck2(arg, argc, argv);
                p.cp.uniform[2][NPARAM+i]=atof(argv[arg++]);               // b
                p.cp.uniform[2][NPARAM+i+1]=p.cp.uniform[1][NPARAM+i+1]=0; // sym u[0, 0]
                if(p.cp.uniform[1][NPARAM+i]>p.cp.uniform[2][NPARAM+i]) messageerr("with -M l a b option, b must be larger than a\n");
                /* Celine changed 12/28/2009 */               
                /* if((i==0)&&(p.cp.uniform[2][NPARAM+i]==0)) messageerr("with -M l a b option must use b>0\n");
                   if(p.cp.uniform[2][NPARAM+i]!=0)
                   {*/
                p.cp.maxval[NPARAM+i]=(p.cp.uniform[2][NPARAM+i]);     // Max values for Mig
                p.cp.minval[NPARAM+i]=(p.cp.uniform[1][NPARAM+i]);
                /* }
                   else
                   p.cp.maxval[NPARAM+i]=p.cp.minval[NPARAM+i]=p.cp.uniform[2][NPARAM+i];
                */////
                p.cp.maxval[NPARAM+i+1]=p.cp.uniform[2][NPARAM+i+1];       // values for histogram
                p.cp.minval[NPARAM+i+1]=p.cp.uniform[1][NPARAM+i+1];
              }
            else                              // -m i j Mij fixed
              {
                migr=atof(argv[arg-1]);
                for(i=0;i<p.cp.npop;i++)
                  for(j=0;j<p.cp.npop;j++) p.cp.mig_mat[i][j]=migr/(p.cp.npop-1);
                for(i=0;i<p.cp.npop;i++) p.cp.mig_mat[i][i]=migr;
                p.cp.maxval[NPARAM]=p.cp.minval[NPARAM]=p.cp.maxval[NPARAM+1]=p.cp.minval[NPARAM+1]=migr; // Sym migration rate
              }
            break;
            /*/////*/
          case 'n' :                          // population size different for pop i
            arg++;
            /* Celine changed 05/29/2008 */
            /* argcheck(arg, argc, argv);
               pop=atoi(argv[arg++]) -1;
               if((pop>=p.cp.npop)||(pop<0)) messageerr("with -n i size option i must be 1 or 2\n");
            *//////
            pop=2 -1;
            argcheck(arg, argc, argv);
            ch3=argv[arg++][0];
            if(ch3=='l')
              {
                p.cp.uniform[0][2]=2;                   // tag for b=N2/N0 lnb~U[a1, a2]
                argcheck2(arg, argc, argv);
                p.cp.uniform[1][2]=atof(argv[arg++]);   // a1
                argcheck2(arg, argc, argv);
                p.cp.uniform[2][2]=atof(argv[arg++]);   // a2
                if(p.cp.uniform[1][2]>=p.cp.uniform[2][2]) messageerr("with -n l a b option b must be larger than a\n");
                p.cp.maxval[1]=exp(p.cp.uniform[2][2]); // Max values // for N2
                p.cp.minval[1]=exp(p.cp.uniform[1][2]);
              }
            else if(ch3=='u')
              {
                p.cp.uniform[0][2]=1;                 // tag for b~U[a1, a2](option 0-1), or theta2~u[] (option 2-3)
                argcheck(arg, argc, argv);
                p.cp.uniform[1][2]=atof(argv[arg++]); // a1
                argcheck(arg, argc, argv);
                p.cp.uniform[2][2]=atof(argv[arg++]); // a2
                if(p.cp.uniform[1][2]>=p.cp.uniform[2][2]) messageerr("with -n u a b option b must be larger than a\n");
                p.cp.maxval[1]=p.cp.uniform[2][2];    // values for histogram
                p.cp.minval[1]=p.cp.uniform[1][2];
              }
            else                              // else b (option 0-1) or theta2 (option 2-3) fixed
              {
                psize=atof(argv[arg-1]);
                p.cp.size[pop]=psize;
                p.cp.maxval[1]=p.cp.minval[1]=psize;
                p.cp.uniform[0][2]=5;
              }
            break;
          case 'g' :                          // Growth
            if(npop<2) messageerr("Must use -l option first.\n");
            arg++;
            argcheck(arg, argc, argv);
            pop=atoi(argv[arg++]) -1;
            if((pop>=p.cp.npop)||(pop<0)) messageerr("with -g i alpha option i must be 1 or 2\n");
            if(arg>=argc) messageerr("Not enough arg's after -g.\n");
            palpha=atof(argv[arg++]);
            p.cp.alphag[pop]=palpha;
            break;
          case 'G' :                          // Growth rate for all pop
            arg++;
            if(arg>=argc) messageerr("Not enough arg's after -G.\n");
            palpha=atof(argv[arg++]);
            for(i=0;i<p.cp.npop;i++) p.cp.alphag[i]=palpha;
            break;
            /* Celine changed 05/29/2008 */
          case 'N' : argv[arg]="-eN";
            break;
            /*/////*/
            //--- Event ---//
          case 'e' :                          // Timing of event
            pt=(struct devent *)malloc((unsigned) sizeof(struct devent));
            pt->time=pt->timep=pt->paramv=pt->paramvf=0.0;
            pt->popi=pt->popj=0;
            pt->nextde=NULL;
            pt->detype=pt->detypef='\0';
            pt->mat=NULL;
            pt->detype=argv[arg][2];          // Type of event
            pt->detypef=pt->detype;
            ch4=argv[arg][3];
            arg++;
            if(pt->detype!='N')               /* Celine changed 05/29/2008 */
              {
                argcheck(arg, argc, argv);
                ch3=argv[arg++][0];
              }
            else argv[arg-1]="-N";            /* Celine changed 05/29/2008 */

            if(ch3=='u')                      // Time ~U[a, b]
              {
                /* if((pt->detype=='N')||(pt->detype=='j'))*/ /* Celine changed 05/29/2008 */
                if(pt->detype=='j')
                  {
                    double a, b;
                    a=b=0;
                    argcheck(arg, argc, argv);
                    a=atof(argv[arg++]);      // a
                    argcheck(arg, argc, argv);
                    b=atof(argv[arg++]);      // b
                    if(p.cp.uniform[0][3]<=0) // Set time priors
                      {
                        p.cp.uniform[0][3]=1; // tag for Time ~U[a, b]
                        p.cp.uniform[1][3]=a; // a
                        p.cp.uniform[2][3]=b; // b
                        if(p.cp.uniform[1][3]>=p.cp.uniform[2][3]) messageerrs("with -e", &pt->detype, " u a b option, b must be larger than a\n");
                      }
                    /* Celine changed 05/29/2008 */ // No use any more
                    /* else                      // eN and ej need same priors
                       {
                       if((p.cp.uniform[1][3]!=a)||(b!=p.cp.uniform[2][3])) messageerr("with otpions -ej and -N priors of times required of events must be indentical\n");
                       }
                    *//////
                    p.cp.maxval[2]=p.cp.uniform[2][3]; // Values for histogram
                    p.cp.minval[2]=p.cp.uniform[1][3];
                  }
                /* Celine changed 05/29/2008 */
                /* else if(pt->detype!='N')
                   messageerr("with otpions not -ej and -N, time of the event needs to be fixed. Can't use priors.\n");
                *//////
              }
            else
              {
                /* if((pt->detype=='N')||(pt->detype=='j'))*/ /* Celine changed 05/29/2008 */
                if(pt->detype=='j')
                  {
                    pt->time=atof(argv[arg-1]); // Set time for eN or ej
                    p.cp.maxval[2]=p.cp.minval[2]=pt->time;
                  }
                else if(pt->detype!='N')
                  /* else*/ /* Celine changed 05/29/2008 */
                  pt->timep=atof(argv[arg-1]); // Set fraction of split time for other events
              }
            p.cp.listevent[nevent]=pt;
            nevent++;
            if(nevent==11) messageerr("Can not have more than 10 otpions -e. If want to use more, increase allication of p.cp.listevent above.\n");
            //--- Type of event ---//
            switch(pt->detype)
              {
                case 'N' :                    // New ancestral population (time event=time of split of option -ej)
                  argcheck(arg, argc, argv);
                  ch3=argv[arg++][0];
                  if(ch3=='l')
                    {
                      p.cp.uniform[0][4]=2;                   // tag for a=Na/N0 from: lna~U[b1, b2]
                      argcheck2(arg, argc, argv);
                      p.cp.uniform[1][4]=atof(argv[arg++]);   // b1
                      argcheck2(arg, argc, argv);
                      p.cp.uniform[2][4]=atof(argv[arg++]);   // b2
                      if(p.cp.uniform[1][4]>=p.cp.uniform[2][4])
                        messageerrs("with -e", &pt->detype, " t l a b option b must be larger than a\n");
                      p.cp.maxval[4]=exp(p.cp.uniform[2][4]); // Value for historgram
                      p.cp.minval[4]=exp(p.cp.uniform[1][4]);
                    }
                  else if(ch3=='u')
                    {
                      p.cp.uniform[0][4]=1;                 // tag a=Na/N0 a~U[b1, b2] (option 0-1), or thetaA~U[b1, b2] (option 2-3)
                      argcheck(arg, argc, argv);
                      p.cp.uniform[1][4]=atof(argv[arg++]); // b1
                      argcheck(arg, argc, argv);
                      p.cp.uniform[2][4]=atof(argv[arg++]); // b2
                      if(p.cp.uniform[1][4]>=p.cp.uniform[2][4])
                        messageerrs("with -e", &pt->detype, " t u a b option b must be larger than a\n");
                      p.cp.maxval[4]=p.cp.uniform[2][4];    // Values for histogram
                      p.cp.minval[4]=p.cp.uniform[1][4];
                    }
                  else
                    {
                      pt->paramv=atof(argv[arg-1]);
                      p.cp.maxval[4]=p.cp.minval[4]=pt->paramv;
                      p.cp.uniform[0][4]=5;
                    }
                  break;
                case 'b' :                    // Change all populations sizes at time event
                  argcheck(arg, argc, argv);
                  pt->paramv=pt->paramv=atof(argv[arg++]);
                  break;
                case 'G' :                    // Populations growth rates change at time event
                  if(arg>=argc) messageerr("Not enough arg's after -eG.\n");
                  pt->paramv=atof(argv[arg++]);
                  break;
                case 'M' :                    // Migration matrix at time event
                  argcheck(arg, argc, argv);
                  pt->paramv=atof(argv[arg++]);
                  break;
                case 'n' :                    // pop i size changes at time of event
                  argcheck(arg, argc, argv);
                  pt->popi=atoi(argv[arg++]) -1;
                  if((pt->popi>=p.cp.npop)||(pt->popi<0)) messageerr("with -en t i size option i must be 1 or 2\n");
                  argcheck(arg, argc, argv);
                  pt->paramv=pt->paramv=atof(argv[arg++]);
                  break;
                case 'g' :                    // Rate growth of popi changes at time event
                  argcheck(arg, argc, argv);
                  pt->popi=atoi(argv[arg++]) -1;
                  if((pt->popi>=p.cp.npop)||(pt->popi<0)) messageerr("with -eg t i alpha option i must be 1 or 2\n");
                  if(arg>=argc) messageerr("Not enough arg's after -eg.\n");
                  pt->paramv=atof(argv[arg++]);
                  break;
                case 's' :                    // Split pop i into two pops (=admixture in forward simulation)
                  /* Celine changed 05/29/2008 */
                  /* argcheck(arg, argc, argv);
                     pt->popi=atoi(argv[arg++]) -1;
                  *//////
                  pt->popi=1 -1;              // split pop1 (only accepteable after 2-1 -eh case
                  if((pt->popi>=p.cp.npop)||(pt->popi<0)) messageerr("with -es t i proportion option i must be 1 or 2\n");
                  argcheck(arg, argc, argv);
                  pt->paramv=atof(argv[arg++]); // Parameter=proportion of chromosome in each new population
                  break;
                case 'm' :                    // Change 4Nomij at time event
                  if(ch4=='a')                // Migration matrix case
                    {
                      pt->detype='a';
                      argcheck(arg, argc, argv);
                      npop2=atoi(argv[arg++]);
                      pt->mat=(double**)malloc((unsigned)npop2*sizeof(double*));
                      for(pop=0; pop<npop2; pop++)
                        {
                          (pt->mat)[pop]=(double*)malloc((unsigned)npop2*sizeof(double));
                          for(i=0;i<npop2;i++)
                            {
                              if(i==pop) arg++;
                              else
                                {
                                  argcheck(arg, argc, argv);
                                  (pt->mat)[pop][i]=atof(argv[arg++]);
                                }
                            }
                        }
                      for(pop=0; pop<npop2; pop++)
                        {
                          (pt->mat)[pop][pop]=0.0;
                          for(pop2=0; pop2<npop2; pop2++){
                            if(pop2 !=pop) (pt->mat)[pop][pop]+=(pt->mat)[pop][pop2];
                          }
                        }
                    }// End of Migration matrix case
                  else                        // Information for each population
                    {
                      argcheck(arg, argc, argv);
                      pt->popi=atoi(argv[arg++]) -1;
                      argcheck(arg, argc, argv);
                      pt->popj=atoi(argv[arg++]) -1;
                      if((pt->popj==pt->popi)||((pt->popi>=p.cp.npop)||(pt->popj>=p.cp.npop))||((pt->popj<0)||(pt->popi<0)))
                        messageerr("with -em t i j m_ij option i and j must be different and 1 or 2\n");
                      argcheck(arg, argc, argv);
                      pt->paramv=atof(argv[arg++]);
                    }
                  break;
                case 'j' :                    // Collapse lineages from i to j at time event=time of split (=split of populations in forward simulation)
                  /* Celine changed 05/29/2008 */
                  /* argcheck(arg, argc, argv);
                     pt->popi=atoi(argv[arg++]) -1;
                     argcheck(arg, argc, argv);
                     pt->popj=atoi(argv[arg++]) -1;
                  *//////
                  pt->popi=2 -1;
                  pt->popj=1-1;
                  /* Celine changed 12/28/2009 */ // no use
                  /* if((pt->popj==pt->popi)||((pt->popj>=p.cp.npop)||(pt->popj>=p.cp.npop))||((pt->popj<0)||(pt->popi<0)))
                     messageerr("with -ej t i j m_ij option i and j must be different and 1 or 2\n");
                  *//////
                  break;
                case 'h' :                    // Collapse lineages from i to j at time event (!=to time of split)
                  /* Celine changed 05/29/2008 */
                  /* argcheck(arg, argc, argv);
                     pt->popi=atoi(argv[arg++]) -1;
                     argcheck(arg, argc, argv);
                     pt->popj=atoi(argv[arg++]) -1;
                  *//////
                  pt->popi=2 -1;
                  pt->popj=1-1;
                  break;
                default: messageerr("e event\n");
                  pt=pt->nextde;
                  free(pt);
              }// End switch on type of events
            break;
          case 'l' :                          // Multiple loci information
            p.cp.npop=NPOP;
            npop=p.cp.npop;
            p.cp.config=(int*) realloc(p.cp.config, (unsigned)((p.cp.npop+1)*sizeof(int)));      // Allocate memory of arrays for: -the two sample sizes
            p.cp.mig_mat=(double**)realloc(p.cp.mig_mat, (unsigned)(p.cp.npop*sizeof(double*))); // -Migration matrix
            p.cp.mig_mat[0]=(double*)realloc(p.cp.mig_mat[0], (unsigned)(p.cp.npop*sizeof(double)));
            for(i=1;i<p.cp.npop;i++)
              p.cp.mig_mat[i]=(double*)malloc((unsigned)(p.cp.npop*sizeof(double)));
            p.cp.size=(double*) realloc(p.cp.size, (unsigned)(p.cp.npop*sizeof(double)));         // -Size pop (option -n)
            p.cp.alphag=(double*) realloc(p.cp.alphag, (unsigned)(p.cp.npop*sizeof(double)));     // -Growth rate (option -g)
            for(i=1;i<p.cp.npop;i++)
              {
                (p.cp.size)[i]=(p.cp.size)[0];
                (p.cp.alphag)[i]=(p.cp.alphag)[0];
              }
            for(i=0;i<*phowmany;i++)
              {
                p.lp[i].li=0;
                for(j=0;j<4;j++)
                  {
                    if(j<3)
                      p.lp[i].ni[j]=0;
                    p.lp[i].Si[j]=0;
                  }
                p.lp[i].xi=p.lp[i].vi=p.lp[i].wi=0;
              }
            if(argv[arg][2]=='f')             // Multiple loci information from a file
              {
                arg++;
                argcheck(arg, argc, argv);
                p.cp.fin=argv[arg];
                lf=fopen(argv[arg], "r");     // open file
                if(lf==NULL) messageerrs(" No multi loci information file named ", argv[arg], "\n");
                arg++;
                fscanf(lf, "%s ", st);        // Ignore beginning of file until find "**" or "//"
                while((!feof(lf))&&(!((st[0]=='/')&&(st[1]=='/')))&&(!((st[0]=='*')&&(st[1]=='*'))))
                  fscanf(lf, "%s ", st);
                while((!feof(lf))&&(!((st[0]=='/')&&(st[1]=='/'))))
                  fscanf(lf, "%s ", st);

                if((st[0]=='/')&&(st[1]=='/')) // The locus information starts after //
                  {/* Celine changed 12/28/2009 */
                    xw=0;                     // max w_y
                    nw=10000;                 // min w_y
                    /*/////*/
                    for(i=0;i<*phowmany;i++)
                      {
                        /* Celine changed 09/18/2007 */
                        p.lp[i].name=(char*)malloc((unsigned)50*sizeof(char));
                        fscanf(lf, "%s ", p.lp[i].name); //name locus;
                        fscanf(lf, "%s\t", st); // lenght
                        p.lp[i].li=atoi(st);

                        if(p.lp[i].li<2){fprintf(stderr, "with -lf option nsites should be >1. error at loci# %d (\"%s\")\n", i+1, p.lp[i].name); format();}
                        fscanf(lf, "%s\t", st); // scalar xi
                        p.lp[i].xi=atof(st);
                        if(p.lp[i].xi<=0) {fprintf(stderr, "with -lf option xi (inheritence scalar) should be >0. error at loci# %d (\"%s\")\n", i+1, p.lp[i].name);format();}
                        fscanf(lf, "%s\t", st); // scalar vi
                        p.lp[i].vi=atof(st);
                        if(p.lp[i].vi<=0) {fprintf(stderr, "with -lf option vi (mutation rate variation scalar) should be >0. error at loci# %d (\"%s\")\n", i+1, p.lp[i].name);format();}
                        fscanf(lf, "%s\t", st); // scalar wi
                        p.lp[i].wi=atof(st);
                        /* Celine changed 12/28/2009 */
                        if(p.lp[i].wi>xw) xw=p.lp[i].wi;
                        if((p.lp[i].wi<nw)&&(p.lp[i].wi>0)) nw=p.lp[i].wi;
                        /*/////*/
                        if(p.lp[i].wi<0) {fprintf(stderr, "with -lf option wi (inheritence scalar for recombination) should be >0. error at loci# %d (\"%s\")\n", i+1, p.lp[i].name);format();}
                        for(j=0;j<npop;j++)   // sample sizes: n12 and n21
                          {
                            fscanf(lf, "%d\t", &p.lp[i].ni[j+1]);
                            if(p.lp[i].ni[j+1]<1) {fprintf(stderr, "with -lf option ni1 and ni2 should be >1. error at ni%d at loci# %d (\"%s\")\n", j+1, i+1, p.lp[i].name); format();}
                          }
                        p.lp[i].ni[0]=p.lp[i].ni[1]+p.lp[i].ni[2]; // Total sample size=ni1+ni2
                        for(j=0;j<4;j++)      // Init array of locus specific S stats
                          {
                            fscanf(lf, "%d\t", &p.lp[i].Si[j]);
                            if(p.lp[i].Si[j]<0) {fprintf(stderr, "with -lf option Si should be >0. error at Si%d at loci# %d (\"%s\")\n", j+1, i+1, p.lp[i].name); format();}
                            p.lp[i].H[j]=0;   // Init array of locus statistics for gof
                            p.lp[i].S[j]=0;   // init population specific S & freq
                          }
                        p.lp[i].S[4]=0;
                        for(j=4;j<9;j++)
                          p.lp[i].H[j]=0;
                        if((feof(lf))&&(i+1<*phowmany)) {fprintf(stderr, "Not enough loci in file \"%s\": found %d, needed %d \n", argv[arg-1], i+1, *phowmany);format();}
                        else if(!(feof(lf))&&(i+1>=*phowmany)) {fprintf(stderr, "Too many loci in file \"%s\": found >%d, needed %d \n", argv[arg-1], i+1, *phowmany);format();}
                        /*/////*/
                      }// End loop on loci
                  }// End if found loci information in file
                else {fprintf(stderr, " No multi loci information in file %s\n", argv[arg-1]); format();}
                fclose(lf);
              }// End -lf option
            /* Celine commented 12/28/2009 - No use */
            /* else
               {
               arg++;
               for(i=0;i<*phowmany;i++)
               {
               // Celine changed 09/18/2007 //
               p.lp[i].name=(char*)malloc((unsigned)50*sizeof(char));
               sprintf(p.lp[i].name, "Locus%d", i+1);
               /////
               argcheck(arg, argc, argv); // li
               p.lp[i].li=atoi(argv[arg++]);
               if(p.lp[i].li<2) messageerrd("with -l option nsites should be >1. error at Locus# ", i+1, "\n");
               argcheck(arg, argc, argv); // xi
               p.lp[i].xi=atof(argv[arg++]);
               if(p.lp[i].xi<=0) messageerrd("with -l option xi (inheritence scalar) should be>0. error at Locus# ", i+1, "\n");
               argcheck(arg, argc, argv); // vi
               p.lp[i].vi=atof(argv[arg++]);
               if(p.lp[i].vi<=0) messageerrd("with -l option vi (mutation rate variation scalar) should be>0. error at Locus# ", i+1, "\n");
               argcheck(arg, argc, argv); // wi
               p.lp[i].wi=atof(argv[arg++]);
               if(p.lp[i].wi<0) messageerrd("with -l option wi (inheritence scalar for recombination) should be>=0. error at Locus# ", i+1, "\n");
               argcheck(arg, argc, argv); // ni1
               p.lp[i].ni[1]=atoi(argv[arg++]);
               if(p.lp[i].ni[1]<1) messageerrd("with -l option ni1 and ni2 should be >1. error at ni1 at Locus# ", i+1, "\n");
               argcheck(arg, argc, argv); // ni2
               p.lp[i].ni[2]=atoi(argv[arg++]);
               if(p.lp[i].ni[1]<1) messageerrd("with -l option ni1 and ni2 should be >1. error at ni2 at Locus# ", i+1, "\n");
               p.lp[i].ni[0]=p.lp[i].ni[1]+p.lp[i].ni[2]; // nsami=ni1+ni2
               }// End loop along loci
               }// End option -l
            *//////
            break;
            /* case 'o' :                         // No use here
               arg++;
               argcheck(arg, argc, argv);
               p.cp.fout=argv[arg++];
               break;*/
          case 'u' :                          // mutation rate per bp
            arg++;
            argcheck(arg, argc, argv);
            p.mp.mu=atof(argv[arg++]);
            if(p.mp.mu>=.0001) messageerrl("with -u mu option must specify a mutation rate PER BASE PAIR (~e-8). :", p.mp.mu, " is too large\n");
            /* argcheck(arg, argc, argv);        // No use here
               p.mp.g=atof(argv[arg++]);                  if(p.mp.g<=0) messageerr("with -u mu g option must specify a # of generation per year (e.g.:1/15 for drosophila, 15-20 for Chimp):\n");
            */
            break;
            /* Celine changed 12/28/2009 - Uncomment for Celine's use */
            /* case 'q' :            // Option for output
               arg++;
               argcheck(arg, argc, argv);
               p.cp.sumout=atoi(argv[arg++]);
               break;
            *//////
          default: messageerrs(" option default: Unknown flag \"", argv[arg], "\"\n");
        }// End switch to next argument
    }// End While Look at all arguments

  if((p.cp.maxval[2]>0)&&(p.cp.maxval[0]>0))  // compute T in generations
    {
      p.cp.maxval[3]=p.cp.maxval[2]*p.cp.maxval[0]/p.mp.mu;
      p.cp.minval[3]=p.cp.minval[2]*p.cp.minval[0]/p.mp.mu;
    }
  //--- Check for main arguments ---//
  if((p.mp.thetafix==0.0)&&(p.cp.uniform[0][0]!=1)) messageerr(" -t option must be used. \n");
  else if(p.mp.mu==0.0) messageerr("-u mu option must be used with 0<mu<1E-4.\n");
  else if(p.lp[0].li<2) messageerrd("-lf option must be used with ", *phowmany, " loci\n");
  /* Celine changed 12/28/2009 - Uncomment for Celine's use */
  /* else if((p.cp.sumout<0)||(p.cp.sumout>3)) messageerr(" -q option 0, 1, 2 or 3.\n");
   *//////

  /* Celine changed 12/28/2009 */
  //-- Rho genomic fixed or 1 or 2 --//
  if(p.cp.uniform[0][1]<0)                    // Fixed rho values
    {
      if((p.cp.rfix>.1)&&(p.cp.rfix!=1)&&(p.cp.rfix!=2)) messageerrl(" -r rho -> give value for rho <.1 when rho is fixed. \n           -> Or rho=1 when, for locus y, w_y=rho_y_bp=s_y*4N_1*c_y: an estimate of the locus-specific population recombination rate per bp from linkage disequilibrium.\n           -> Or rho=2 when, for locus y, w_y=s_w*c_y=c_sex-average, where c_y is an estimate of the locus-specific cross-over rate per bp from pedigree analyses and w is the recombination scalar (i.e. c_y=c_female and s_y=1/2 for an X-linked locus). The value is now: rho=", p.cp.rfix, " \n");
      if(p.cp.rfix==1)                        // w_y=4N_y c_y >1e-4 <.1
        {
          if(xw>.1) messageerrl(" -r 1 -> the recombination scalar for locus y, w_y=s_y*4N_1*c_y <.1. The value is now: max w_y=", xw, " \n");
          if((nw>0)&&(nw<.0001)) messageerrl(" -r 1 -> the recombination scalar for locus y, w_y=s_y*4N_1*c_y > 1e-4. The value is now: min w_y=", nw, " \n");
        }
      else if(p.cp.rfix==2)                   // w_y=w c_y <1e45 <.1
        {
          if(xw>.0001) messageerrl(" -r 2 -> the recombination scalar for locus y, w_y=s_y*c_y <1e-4. The value is now: max w_y=", xw, " \n");
        }
      else
        {
          if(xw>10) messageerrl(" -r rho -> the recombination scalar for locus y, w_y<1 (may be larger if # recombining sites >> # mutating sites). The value is now: max w_y=", xw, " \n");
          if((nw>0)&&(nw<.0001)) messageerrl(" -r rho -> the recombination scalar for locus y, w_y>1e-4. The value is now: min w_y=", nw, " \n");
        }
    }
  /*/////*/

  if(&p.cp.listevent[0]==NULL) messageerr(" No Event, -ej option must be used\n");
  else
    {
      pt=p.cp.listevent[0];
      double tj=0;
      i=j=k=0;
      while((pt !=NULL)&&(j==0))
        {
          if(pt->detype=='j')
            {
              j+=2;
              tj=pt->time;
            }
          k++;
          pt=p.cp.listevent[k];
        }
      i=k=0;
      pt=p.cp.listevent[0];
      double tN=0;
      while(pt !=NULL)
        {
          if(pt->detype=='N')
            {
              j+=1;
              /* tN=pt->time;*/ /* Celine changed 05/29/2008 */
              tN=pt->time=tj;
              if((p.cp.sumout>=2)&&
                 (p.cp.uniform[0][4]==5))     // thetaA given and fixed to its value
                {
                  /* Celine changed 11/05/2008 - Uncomment for Celine's use */
                  /* if(pt->paramv>.1) messageerrd(" -N thetaA -> give value for thetaA (<.1) when -q ", p.cp.sumout, " flag is used.\n");
                   *//////
                  if(p.cp.uniform[0][0]<=0)                     // Case theta1 fixed as well
                    {
                      p.cp.newest[4]=p.cp.oldest[4]=pt->paramv; // fix estimate
                      pt->paramv/=p.mp.thetafix;                // a fixed to its value
                    }
                  else                                          // Theta1 from prior
                    {
                      p.cp.uniform[0][4]=4;                     // Need to update a if theta1 from prior
                      p.cp.uniform[1][4]=pt->paramv;
                      p.cp.newest[4]=p.cp.oldest[4]=pt->paramv; // fix a to its value
                    }
                }
            }
          else
            {
              if((pt->timep>0))               // Time fixed for other event
                {
                  if(p.cp.uniform[0][0]<=0)   // if theta1 fixed as well
                    {
                      if(p.cp.sumout>=2)
                        {
                          if(p.cp.uniform[0][3]<=0)  // split time fixed - Tcoal fixed to its value
                            pt->time=pt->timep*tj*p.mp.mu/p.mp.thetafix;
                        }
                      else if(p.cp.uniform[0][3]<=0) // split time fixed - Tcoal fixed to its value
                        pt->time=pt->timep*tj;
                    }
                }
            }
          k++;                                // New position in the list of event
          pt=p.cp.listevent[k];
        }
      if((j==0)||(j==1)) messageerr(" No Divergence Event, -ej option must be used\n");
      /* Celine changed 05/29/2008 */
      /* else if((j==3)&&(tN!=tj)) messageerr(" with -N tN size and -ej tj i j options, must use tN=tj\n");*//////
    }

  //--- Default param for a=thetaA/Theta1 and b=theta2/theta1 ---//
  if((p.cp.uniform[0][2]<=0)&&(p.cp.maxval[1]<=0)) // b=1
    p.cp.maxval[1]=p.cp.minval[1]=1;
  if((p.cp.uniform[0][4]<=0)&&(p.cp.maxval[4]<=0)) // a=1
    p.cp.maxval[4]=p.cp.minval[4]=1;
  /* Celine change 12/28/2009 */
  /* if(((p.cp.uniform[2][NPARAM+1]<=0)&&(p.cp.uniform[0][NPARAM+1]>1))&&(p.cp.uniform[0][NPARAM]<=0)) // Check for symetrical migration rate taken from priors
     messageerr("with -m 2 1 u 0 0 option must be use when the option -m 1 2 u a b to model a symatrical migration rate taken from priors\n");
     if((p.cp.uniform[0][NPARAM]>0)&&((p.cp.uniform[0][NPARAM+1]>0)&&(p.cp.uniform[0][NPARAM+1]!=p.cp.uniform[0][NPARAM])))
     messageerr("with -m 1 2 and -m 2 1 options, the priors must be similar (i.e. either u or t)\n");
  *//////
  //--- Change parameter values in option >=1 : b->theta2, a->thetaA, and tcoal or Tgen need to be calculated
  //-- Theta 1 --//
  if(p.cp.uniform[0][0]<0)                    // Fixed theta1 value
    p.cp.newest[0]=p.cp.oldest[0]=p.mp.thetafix;;

  //-- Theta 2 --//
  if(p.cp.uniform[0][2]==1)                   // Check that variance given for kernel if b taken from prior
    {
      if(p.cp.sumout>=2)
        p.cp.uniform[0][2]=3;                 // Option 2 and 3: take theta2 from prior (and not b=theta2/Theta2)
    }
  else if((p.cp.uniform[0][2]==2)&&(p.cp.sumout>=2))
    /* Uncommented for Celine's use */
    /* messageerr(" -n theta2 -> theta2 from uniform prior or b fixed when -q 2 or 3 flag is used.\n");
     *//////
    messageerr(" -n theta2 -> theta2 must be taken from uniform prior (with u flag) or b fixed.\n");
  else if(p.cp.sumout>=2)                     // Case: Theta 2 given or not specified & option >=2
    {
      if(p.cp.uniform[0][0]<=0)               // Case theta1 fixed
        {
          if((p.cp.maxval[1]==p.cp.minval[1])&&
             (p.cp.maxval[1]==1))             // b=1 so theta2=theta1
            p.cp.newest[1]=p.cp.oldest[1]=p.cp.maxval[1]=p.cp.minval[1]=p.mp.thetafix;
          else                                // Theta2 fixed, so b=theta2/theta1
            {
              p.cp.newest[1]=p.cp.oldest[1]=p.cp.size[1];
              p.cp.size[1]/=p.mp.thetafix;
            }
        }
      else                                    // Case theta1 taken from prior
        {
          if((p.cp.uniform[0][2]!=5)&&(p.cp.maxval[1]==p.cp.minval[1])&&(p.cp.maxval[1]==1))
            {                                 // b=1, update theta2=b*theta1 at each time
              p.cp.maxval[1]=p.cp.maxval[0];
              p.cp.minval[1]=p.cp.minval[0];
            }
          else                                // theta2 given, update b=theta2/theta1 every steps
            {
              p.cp.uniform[0][2]=4;
              p.cp.uniform[1][2]=p.cp.size[1];
              p.cp.newest[1]=p.cp.oldest[1]=p.cp.size[1];
            }
        }
    }
  else                                        // Case: b fixed & option <2
    {
      /* Celine changed 11/05/2008 - Uncomment for Celine's use */
      /* if(p.cp.maxval[1]<.01) messageerrd(" -n b -> give value for b (>.01) when -q", p.cp.sumout, " flag is used.\n");
       *//////
      /* Celine changed 12/28/2009 */ // added p.cp.sumout<2
      if(((p.cp.uniform[0][2]<0)||(p.cp.uniform[0][2]==5))&&(p.cp.sumout<2))
        p.cp.maxval[1]=p.cp.minval[1]=p.cp.newest[1]=p.cp.oldest[1]=p.cp.size[1];
      /*/////*/
    }
  /* Celine changed 11/05/2008 - Uncomment for Celine's use */
  /* if((p.cp.maxval[1]>.1)&&(p.cp.sumout>=2)) messageerrd(" -n theta2 -> give value for theta2 (<.1) when -q", p.cp.sumout, " flag is used.\n");
     if((p.cp.maxval[1]<.01)&&(p.cp.sumout<2)) messageerrd(" -n theta2 -> give value for b (>.01) when -q", p.cp.sumout, " flag is used.\n");*//////

  //-- Time --//
  // Depending of the option, the user gives Tgen or Tcoal - Below checks and initiation for both parameters.
  /* Uncommented for Celine's use */
  /* if((p.cp.maxval[2]<100)&&(p.cp.sumout>=2)) messageerrd(" -ej T in generations when -q ", p.cp.sumout, " flag is used.\n");
     if((p.cp.maxval[2]>100)&&(p.cp.sumout<2)) messageerrd(" -ej T in coalescent unit when -q ", p.cp.sumout, " flag is used.\n");
  *//////
  if(p.cp.uniform[0][3]==1)                   // Check that variance given for kernel if split time (tcoal) taken from prior
    {
      if(p.cp.sumout>=2)                      // Option 2 and 3: take Tgen from prior (and not tcoal)
        p.cp.uniform[0][3]=3;
    }
  if(p.cp.sumout>=2)                          // Case: Tgen given or not specified & option >=2
    {
      if(p.cp.uniform[0][3]<0)                // Case: Theta1 may vary: update tcoal every step
        {
          p.cp.newest[3]=p.cp.oldest[3]=p.cp.maxval[2];
          p.cp.uniform[0][3]=4;
          p.cp.uniform[1][3]=p.cp.maxval[2];
        }
      p.cp.maxval[3]=p.cp.maxval[2];
      p.cp.minval[3]=p.cp.minval[2];
      if(p.cp.minval[0]>0)                    // Special care for edges of range of Tcoal
        {
          p.cp.maxval[2]*=(double)(p.mp.mu/p.cp.minval[0]);
          p.cp.minval[2]*=(double)(p.mp.mu/p.cp.maxval[0]);
        }
      else
        {
          p.cp.minval[2]*=(double)(p.mp.mu/p.cp.maxval[0]);
          p.cp.maxval[2]*=(double)(p.mp.mu/8e-5);
        }
    }
  else if(p.cp.sumout<2)                      // case Tcoal & options 0, 1
    {
      p.cp.minval[3]=(double) p.cp.minval[2]*p.cp.minval[0]/p.mp.mu;
      p.cp.newest[3]=p.cp.maxval[3]=(double) p.cp.maxval[2]*p.cp.maxval[0]/p.mp.mu;
      if((p.cp.uniform[0][3]<0))              // Case Tcoal fixed
        p.cp.newest[2]=p.cp.oldest[2]=p.cp.maxval[2]=p.cp.minval[2];
    }

  //-- theta A--//
  if(p.cp.uniform[0][4]==1)                   // Check that variance given for kernel if a taken from prior
    {
      if(p.cp.sumout>=2)
        p.cp.uniform[0][4]=3;                 // option >=2: ThetaA taken from ptior (not a=ThetaA/theta1)
    }
  /* Celine changed 11/05/2008 - Uncomment for Celine's use */
  /* else if((p.cp.uniform[0][4]==2)&&(p.cp.sumout>=2))
     messageerrd(" -N thetaA -> give value for thetaA (<.1) when -q", p.cp.sumout, " flag is used.\n");
  *//////
  else if((p.cp.sumout>=2))                   // Case thetaA given or not specified & option >=2
    {
      if(p.cp.uniform[0][0]<=0)               // Case Theta1 fixed
        {                                     // a=1, ThetaA=Theta1
          if(((p.cp.maxval[4]==p.cp.minval[4])&&(p.cp.maxval[4]==1)))
            p.cp.newest[4]=p.cp.oldest[4]=p.cp.maxval[4]=p.cp.minval[4]=p.mp.thetafix;
        }
      else                                    // Case Theta1 taken from prior
        {                                     // a=1, update ThetaA=a*Theta1 each step
          if((p.cp.uniform[0][4]!=5)&&(p.cp.maxval[4]==p.cp.minval[4])&&(p.cp.maxval[4]==1))
            {
              p.cp.maxval[4]=p.cp.maxval[0];
              p.cp.minval[4]=p.cp.minval[0];
            }
          else                                // thetaA, update a=thetaA/Theta1 each step
            {
              p.cp.uniform[0][4]=4;
              p.cp.uniform[1][4]=p.cp.oldest[4]=p.cp.newest[4];
            }
        }
    }
  else if(((p.cp.uniform[0][4]<0)||(p.cp.uniform[0][4]==5))&&
          (p.cp.sumout<2))                    // Case a fixed & option <=1
    p.cp.newest[4]=p.cp.oldest[4]=p.cp.minval[4]=p.cp.maxval[4];
  /* Celine changed 11/05/2008 - Uncomment for Celine's use */
  /* if((p.cp.maxval[4]>.1)&&(p.cp.sumout>=2)) messageerrd(" -N thetaA -> give value for thetaA (<.1) when -q", p.cp.sumout, " flag is used.\n");
     if((p.cp.maxval[4]<.01)&&(p.cp.sumout<2)) messageerrd(" -N thetaA -> give value for a (>.01) when -q", p.cp.sumout, " flag is used.\n");
  *//////

  //-- M12 --//
  if((p.cp.uniform[0][5]<0)&&(p.cp.maxval[5]>0)) // fix M12
    p.cp.newest[5]=p.cp.oldest[5]=p.cp.maxval[5];
  if((p.cp.uniform[0][6]<0)&&(p.cp.maxval[6]>0)) // fix M12
    p.cp.newest[6]=p.cp.oldest[6]=p.cp.maxval[6];
  /* Celine changed 12/28/2009 - Uncomment for Celine's use */
  /* if(p.cp.sumout==3)                          // option to compare with IM: m given instead of M
     {
     if((p.cp.uniform[0][5]==1)||(p.cp.uniform[0][5]==5))
     messageerrd(" -m 1 2 m12 taken uniform prior (option u) without test when -q ", p.cp.sumout, " flag is used.\n");
     else if(p.cp.uniform[0][5]==2)          // Case m12 from prior
     p.cp.uniform[0][5]=3;
     else if(p.cp.maxval[5]>0)               // Case m12 fixed: update M12 each steap
     {
     p.cp.uniform[0][5]=4;
     p.cp.uniform[1][5]=p.cp.maxval[5];
     }
     if((p.cp.uniform[0][6]==5)||(p.cp.uniform[0][6]==1))
     messageerrd(" -m 2 1 m21 taken uniform prior (option u) without test when -q ", p.cp.sumout, " flag is used.\n");
     else if(p.cp.uniform[0][6]==2)          // Case m21 from prior
     p.cp.uniform[0][6]=3;
     else if(p.cp.maxval[6]>0)               // Case m21 fixed: update M21 each steap
     {
     p.cp.uniform[0][6]=4;
     p.cp.uniform[1][6]=p.cp.maxval[6];
     }
     if(p.cp.maxval[5]>.005)
     messageerrd(" -m 1 2 m12=M12*mu/theta1 (so very small of the order e-8) when -q ", p.cp.sumout, " flag is used.\n");
     if(p.cp.maxval[6]>.005)
     messageerrd(" -m 2 1 m21=M21*mu/theta1 (so very small of the order e-8) when -q ", p.cp.sumout, " flag is used.\n");
     }// End option 3
  *//////
    /* Celine changed 11/05/2008 - Uncomment for Celine's use */
    /* else
       {
       if((p.cp.maxval[5]>0.)&&(p.cp.maxval[5]<.005)) messageerrd(" -m 1 2 M12=4N1m12 (>.005) when -q ", p.cp.sumout, " flag is used.\n");
       if((p.cp.maxval[6]>0.)&&(p.cp.maxval[6]<.005)) messageerrd(" -m 2 1 M21=4N1m21 (>.005) when -q ", p.cp.sumout, " flag is used.\n");
       }
    *//////
  return p;
}// End Function Get Parameters


/*---------------- Functions for error messages ----------------*
  |  Allow to output different error messages with different kind of informations
*/
//--- 1 string ---//
void messageerr(char st[1000])
{
  void usage();
  fprintf(stderr, "-- %s\n", st);
  usage();
  exit(1);
}
//--- 3 strings ---//
void messageerrs(char st[1000], char st2[1000], char st3[1000])
{
  void usage();
  fprintf(stderr, "-- %s%s%s\n", st, st2, st3);
  usage();
  exit(1);
}
//--- 2 strings and integer ---//
void messageerrd(char st[1000], int i, char st3[1000])
{
  void usage();
  fprintf(stderr, "-- %s%d%s\n", st, i, st3);
  usage();
  exit(1);
}
//--- 2 strings and double ---//
void messageerrl(char st[1000], double i, char st3[1000])
{
  void usage();
  fprintf(stderr, "-- %s%lg%s\n", st, i, st3);
  usage();
  exit(1);
}

/*---------------------------*
 *  Function argument check  |
 *---------------------------*
 |  Checks that number of arguments is ok for the different options.
*/
void argcheck(int arg, int argc, char *argv[])
{
  if((arg>=argc)||(argv[arg][0]=='-'))
    {
      fprintf(stderr, "not enough arguments after %s %s and before %s\n", argv[arg-2], argv[arg-1], argv[arg]);
      fprintf(stderr, "For usage type: ms<return>\n");
      exit(1);
    }
}
/*----------------------------*
 *  Function argument check 2 |
 *----------------------------*
 |  Checks ok number of arguments for prior distribitions' definitions.
*/
void argcheck2(int arg, int argc, char *argv[])
{
  if((arg>=argc))
    {
      fprintf(stderr, "not enough arguments after %s %s and before %s\n", argv[arg-2], argv[arg-1], argv[arg]);
      fprintf(stderr, "For usage type: ms<return>\n");
      exit(1);
    }
}

/*--------------------------------*
 *  Function Format of input file |
 *--------------------------------*
 |  Called if error in the input file.
*/
void format()
{
  fprintf(stderr, " infile format for Y loci: \n");
  fprintf(stderr, " Required: \n\t\"//\n");
  fprintf(stderr, "\tName1 Z1 x1 v1 w1 n11 n21\n");
  fprintf(stderr, "\tName2 Z2 x2 v2 w2 n12 n22\n");
  fprintf(stderr, "\t...\n");
  fprintf(stderr, "\tNameY ZY xY vY wY n1Y n2Y\"\n");
  fprintf(stderr, " More Options / Information - before required part: \n");
  /* Uncommented for Celine's use */
  /* fprintf(stderr, "\tTrue values:\n\t\t\"** \n\t\tTheta N2 T Na m12 m21\"\n");
   *//////
  fprintf(stderr, "\tAny other information can preceed the required information and need to end by \"//\" \n");
  exit(1);
}

/*-----------------*
 *  Function Usage |
 *-----------------*
 | Lists the informations relative to the arguments
 | given in the command line.
*/
void usage()
{
  /* Uncommented for Celine's use - special case When nsim specified for each set of param */
  /* fprintf(stderr, "usage: mimargof nsim Y \n nsim=# of simulations to perform for each step analyzed, Y=# loci\n");
   */
  fprintf(stderr, "usage: mimargof Y \n Y=# loci\n");
  fprintf(stderr, " Required Options: -t theta -u mu -lf input ... -ej time \n");
  fprintf(stderr, " More Options / Information: \n");

  fprintf(stderr, "---- The following options must be used\n");
  fprintf(stderr, "\t -t theta    (theta=4*N1*mu, mu=per site mutation rate)\n");

  fprintf(stderr, "\t -u mu       (mu=per site mutation rate)\n");
  /* Uncommented for Celine's use */
  /* fprintf(stderr, "\t -l Z1 x1 v1 w1 n11 n12...\n");*//////
  fprintf(stderr, "\t -lf input   (File of loci information, see documentation for format)\n");

  fprintf(stderr, "\t -ej Tgen    (Split into two populations at time Tgen generation. In coalescent, join lineages in pop 2 and pop 1 into pop 1)\n");
  fprintf(stderr, "\t\t  [size, alpha and M are unchanged. for the populations]\n");

  fprintf(stderr, "---- The following options are optional\n");
  fprintf(stderr, "\t -n theta2   (Change the population mutation rate in pop 2, theta2=4N2mu. Default, theta2=theta1.\n");
  fprintf(stderr, "\t -N thetaA   (Modify ancestral pop mutation rate thetaA=4Namu. Default, thetaA=theta1. \n");
  fprintf(stderr, "\t -M M        (Symmetrical migration rate (i, j-th element of mig matrix set to M) .\n");
  fprintf(stderr, "\t -m i j M_ij (i, j-th element of mig matrix set to M_ij. i and j in [1, 2])\n");
  fprintf(stderr, "\t -r rho      (rho=4N1c, c=per site recombination rate)\n");
  /* Celine changed 12/28/2009 */
  fprintf(stderr, "\t\t -r 1         (w=rho=s*4N_1*c estimated from LD)\n");
  fprintf(stderr, "\t\t -r 2         (w=s*c; c estimated from pedigree analyses)\n");
  fprintf(stderr, "\t\t    [s=recombination scalar (i.e. c=c_female and s=1/2 for X-linked loci)]\n");
  fprintf(stderr, "\t\t -r e lambda  (r=c/mu per locus drawn from Exponential(lambda)=> E(r)=1/lambda)\n");
  /*/////*/
  /* Celine changed 12/28/2009 - Uncommented for Celine's use */
  /* fprintf(stderr, "\t\t -r u a b (r=c/mu: ln(r) per locus drawn from Uniform[a, b])\n");
     fprintf(stderr, "\t\t -r n a b (r=c/mu per locus drawn from Normal(mu=a, sigma=b))\n");
     fprintf(stderr, "\t -q sumout (type of output: default is 2. 0 for theta1, b, t, a, test Mig; 1 :theta1, theta2, T, thetaA, test mig; 2 : as 1 but theta2, thetaA and T(gen) M12 M21 taken from uniform priors or 3 as in 2 but m12 m21 taken from uniform instead)\n");
  *//////
  fprintf(stderr, "---- The following options are conserved from original ms. They are usually not used and should work, but I have not done any tests on them.\n");
  /* Celine changed 03/18/2010 */
  fprintf(stderr,"\t -se seed1 seed2 seed3     ( specify seeds from command line.)\n");/*/////*/
  fprintf(stderr, "\t -f filename    (Read command line arguments from file filename)\n");
  fprintf(stderr, "\t -c f track_len (f=ratio of conversion rate to rec rate. track_len is mean length) \n");
  fprintf(stderr, "\t -G alpha       (N(t)=N0*exp(-alpha*t). alpha=-log(Np/Nr)/t\n");
  fprintf(stderr, "\t -g i alpha_i \n");
  /* Uncommented for Celine's use */
  /* fprintf(stderr, "\t -T (Output gene tree)\n");
   *//////
  fprintf(stderr, "--- The following options modify other parameters at the time 'time*t', where time is the split time given in option -ej and t is specified as the first argument - Can not use more than 10 arguments - The user needs to be careful to specify times that are compatible with the model. Note -ej can be used only once.\n");
  fprintf(stderr, "\t -eG t alpha     (Modify growth rate of all pop's)\n");
  fprintf(stderr, "\t -eg t i alpha_i (Modify growth rate of pop i) \n");
  fprintf(stderr, "\t -eb t size      (Like -N, but size is the ratio of N_A/N_1)\n");
  fprintf(stderr, "\t -en t i size_i  (Modify pop size of pop i. New size of popi=size_i*N0 . Size_i=Ni/N0)\n");
  fprintf(stderr, "\t -eM t mig_rate  (Modify the mig matrix so all elements are mig_rate/(npop-1)\n");
  fprintf(stderr, "\t -em t i j m_ij  (i, j-th element of mig matrix set to m_ij at time t)\n");
  /* Uncommented for Celine's use */
  /* fprintf(stderr, "\t -ema t npop m_11 m_12 m_22 m_21 (Assign values to elements of migration matrix)\n");
     fprintf(stderr, "--- Be very careful with two following options since they may disrupt the model. Use eh and eb to create split and admixture of populations. But be careful in using those two options that only two populations remain at present.\\n");
     fprintf(stderr, "\t -es t p (Split: pop1 -> pop1+pop2, npop DOES NOT increase by 1.\n");
     fprintf(stderr, "\t\t p is probability that each lineage stays in pop-i. (p, 1-p are admixture proportions.\n");
     fprintf(stderr, "\t\t Size of pop npop is set to N0 and alpha=0.0, size and alpha of pop i are unchanged.\n");
     fprintf(stderr, "\t -eh t (see ej)\n");
  *//////
  fprintf(stderr, " See MIMARdoc.pdf and MIMARsimdoc.pdf for explanation of these parameters.\n");
  exit(1);
}

/*-----------------------*
 *  Function add to list |
 *-----------------------*
 | Adds a new event in the list of events (recorded in
 | the pointer of parameters).
*/
void addtoelist(struct devent *pt, struct devent *elist)
{
  struct devent *plast=NULL, *pevent=NULL, *ptemp=NULL;
  pevent=elist;
  while((pevent !=NULL)&&(pevent->time<=pt->time))
    {
      plast=pevent;
      pevent=pevent->nextde;
    }
  ptemp=plast->nextde;
  plast->nextde=pt;
  pt->nextde=ptemp;
}

/*-----------------------*
 *  Function Gen maximum |
 *-----------------------*/
double getmax(double fix, double newval)
{
  if(newval>fix)
    return fix;
  else
    return newval;
}

/*-----------------------*
 *  Function Gen minimum |
 *-----------------------*/
double getmin(double fix, double newval)
{
  if(newval<fix)
    return fix;
  else
    return newval;
}
