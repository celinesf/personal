/*******************************      mimar.c      ***********************************
 *
 *      Estimates parameters of the isolation-migration model using summary
 *  statistics.
 *
 * Usage is shown by typing mimar without arguments.
 usage: mimar nsteps bsteps L -l -u -t -ej -o [options]

 nsteps is the number of steps (or minutes) until end of MCMC.
 bsteps is the number of burnin steps.
 L is the number of loci considered.
 -l (or lf) followed by information on Y loci, with they S statistics.
 -u gives the mutation rate per base pair, mu.
 -t gives the population mutation rate per bp, 4N1*mu.
 -ej set the time of split between the two populations.
 -o provides the output file, in which a histogram of 1000 bins is printed.

 Other options: See mimardoc.pdf or after downloading and compiling, type mimar<CR>.


 **  Arguments of the options are explained here:
 *** Locus information (see documentation for details):
 nlj: the sample size from the j^th subpopulation at the l^th locus (all must be
 specified.) 
 Zl, xl, vl, wl: length and inheritance or recombination scalars for l^th locus.
 S1l S2l Ssl Sfl: S statistics for locus l.
 *** Prior distributions information
 .u, t, l, e: stand for uniform, test(reversible jump for migration rate),
 ln(param) from uniform, exponential.
 i,j: population 1 or 2.
 a, b: lower and upper limit of prior distributions.
 lambda: 1/mean of exponential prior.
 *** Other information
 mu: mutation rate per bp
 theta: 4N1*mu times the neutral mutation rate
 time: time of split.
 var: variances for kernel distributions
 ngen: # genealogies simulated per locus.
 sumout: type of output, default is 2.
 osteps:# of steps or minutes between outputs of histograms.
 size_2, size: theta2 and thetaA (by default=to theta1). For options -q<=1,
 these are the ratio theta2/theta1 and thetaA/theta1.
 int: # steps between recorded steps.
 M_ij: mutation rate: fraction in popi replaced by migrant from pop2
 rho: recombination rate per base per, 4N1c
 f: ratio of conversion rate to recombination rate. (Wiuf and Hein model.)
 track_len: mean length of conversion track in units of sites. 
 alpha: growth rate.
 mig_rate: migration rate: the fraction of each subpop made up of
 migrants times 4N1.
 p: admixture proportion.

 Note: In the above definitions, N1 is the total diploid population in the first
 population.
 A seed file called "seedmimar" will be created if it doesn't exist. The
 seed(s) in this file will be modified by the program. So subsequent runs will
 produce new output. The initial contents of seedmimar will be printed on the second line of the output.
 Output in the summary files (-o option) consists of one line with the command line arguments and one
 line with the seed(s).
 The output list the sets of parameters every recorded steps of the MCMC with
 p(Data|parameters,set of genealogies).

 To compile: cc -o mimar mimar.c params.c streec.c rand1.c -lm
 - rand1.c can be replace by rand2.c, rand1t.c or rand2t.c
 (rand1.c uses drand48(), rand1t.c uses clock for seeding and
 drand48(), rand2.c uses rand() and rand2t.c uses clock for seeding and rand()).
 - (Of course, gcc would be used instead of cc on some machines. And -O3 or
 some other optimization switches might be usefully employed with some
 compilers.)
 - I also provided a makefile that contains the compilation line with rand1.c
 use "make" to compile.

 * --- The following changed have been made in ms program and may or may not
 * still apply here.
 * Modifications made to combine ms and mss on 25 Feb 2001
 * Modifications to command line options to use switches 25 Feb 2001
 * Modifications to add // before each sample 25 Feb 2001
 * Modifications to add gene conversion 5 Mar 2001
 * Added demographic options -d 13 Mar 2001
 * Changed ran1() to use rand(). Changed seed i/o to accommodate this change. 20 April.
 * Changed cleftr() to check for zero rand(). 13 June 2001
 * Move seed stuff to subroutine seedit() 11 July 2001
 * Modified streec.c to handle zero length demographic intervals 9 Aug 2001
 * Corrected problem with negative growth rates (Thanks to D. Posada and C. Wiuf) 13 May 2002
 * Changed sample_stats.c to output thetah - pi rather than pi - thetah. March 8 2003.
 * Changed many command line options, allowing arbitrary migration matrix, and
 * subpopulation sizes. Also allows parameters to come from a file. Option to
 * output trees. Option to split and join subpopulations. March 8, 2003. (Old versions saved in msold.tar).
 *!!! Fixed bug in -en option. Earlier versions may have produced garbage when -en ... used. 9 Dec 2003
 * Fixed bug which resulted in incorrect results for the case where rho=0.0
 * and gene conversion rate > 0.0. This case was not handled correctly in early
 * versions of the program. 5 Apr 2004. (Thanks to Vincent Plagnol for pointing out this problem.)
 *---

 *--- Changed for MIMAR between 2003 and 2007 ---*
 ------------------------------------------
 Get param from user with restriction: I added as many checks as I could to stop
 the user of using wrong parameters or options.
 * The user is forced to use the isolation-migration model: There are two
 populations (not more option -I), that split Tgen ago (option -ej).
 * For multiple loci: provide information per locus (can be in a file, option
 -lf): length L in bp, xi (inheritance scalar), vi (mutation rate variation),
 wi (recombination inheritance scalar), sample sizes n1, n2 and the summary
 statistics of the polymorphism: S1 S2 Ss Sf. For a locus l, the locus specific parameters are:
 rho_l=rho*L*wi.
 theta_l=theta*L*xi*vi
 tcoal_l=tcoal/vi
 Migr_l=Migr*xi
 * The population mutations rates and recombination rates are given in bp.
 * Add option -u mu. The user needs to provide an independent estimate of mu (required to compute Tgen=tcoal*mu/theta1).
 * The parameters of interest are theta1 theta2 time (both in coalescent units and generations) thetaA and M.
 * Parameters of interest (and ln(M)) can be sampled from uniform prior
 distributions.
 The time information specified in -ej and -eN need to be the same.
 To have symmetrical migration rates m12=m21, use "-m 1 2 u a b -m 2 1 u 0 0".
 * r=rho/theta can be taken from exponential(lambda)
 * Parameters of interest updated by MCMC using normal kernel distributions
 (ln(M)~normal as well): the user needs to provide the variances for each parameter (option -v).
 * The MCMC run for bsteps steps of burnin and nsteps additional steps.
 * add option -y (c or t) to specify nsteps in minutes or # of steps.
 * Add option -x ngen. Generate ngen (default is 1) genealogies per locus.
 * Add number lineages in pop1 and 2 in structure node.
 * MCMC algorithm:
 0 Initial step: Get set of parameters from prior.
 1 If not at THETA, propose to mode to THETA' using the kernel distributions.
 2 For a locus l:
 - Generate ngen given the set of parameters THETA'. Calculate branch length L1, l2, Ls, Lf.
 - Calculate P(S|THETA', G): p(Sk|Lk)~Po(Sk, Lk*theta_l)
 - Calculate the mean over ngen genealogies of P(S|THETHA, G).
 - If mean is 0 go back to 1, else go back to 2 for locus l+1
 3 Calculate the product of likelihoods over L loci (assumed independent), Y'.
 calculate the hasting ratio h=min(1, Y'/Y) [see Becquet and Przeworski 2007 for details].
 4 Record THETA' with probability h, else record THETA.
 * Allow writing summary while the program is running. (check every 100 steps, whether mimarrun file is present and starts with "yes"
 * Add Option -L: allow to write summary output file every L steps (or minutes, option
 -y t).
 * Output is the recorded sets of parameters and they likelihood.
 * Summary output file (option -o) contains information about the CPU time and number of
 steps, as well as the histograms with 1000 of the marginal posterior
 distributions of the parameters of interest. Summaries of the posterior are also
 indicated at the end.

 *--- Information for myself on different tests in the code:
 * When simulating data sets, the true value can be given in input file after "**".
 Used to calculate performance of estimations:
 Calculate bias as theta^hat/theta of mode and mean, and Absolute difference of median
 * If reversible jump for test of migration, two additional summary files: for
 M=0 and M>0. 
 * Option -q0, for ms prior: theta1, b, tcoal, a, M12, M21
 * Option -q1, takes same input as -q0, but outputs theta1, theta2, Tgen, thetaA, M12, M21 (can use test on M)
 * Default Option -q2, takes and outputs theta1, theta2, Tgen, thetaA, M12, M21.
 * Option -q3, for IM priors. Takes and outputs theta1, theta2, Tgen, thetaA, m12, m21.
 Add summary file with comparative parameters to IM.
 * Option -d, to update parameters consequentially or simultaneously (default)
 * -es -eh for admixture and split. 
 * Option -ma works as -m, with "x" for diagonal terms.(not in -ema).
 * ****** commented ********
 * Can output prior distribution (was used in IS case).
 * Prior distribution for IM. Write info changed.
 *
 *
 ****************************** Changes in different verions
 
 *--- Changed 18th March 2010 ---*
  Add objects/functions to specify seeds in the command line.
  -----------------------------

  *--- Changed Dec. 28th 2009: ---*
  % param.c was changed: 
  - the locus specific rec rate is now rho_y=rho*w_y*(L-1)
  instead of rho_y=rho*w_y*(L).
  - changed calloc(1... by malloc(... .
  - added recombination options.
  - commented useless lines.
  -----------------------------

  *--- Changed Nov. 27th 2009: ---*
  % streec.c was changed: Corrected bugs on memory allocation.
  % mimar.c: 
  - changed calloc(1... by malloc(... .
  - correction on memort leaks.
  Thanks to the user who mentionned the bug to me.
  -----------------------------

  *--- Changed March 3rd 2009: ---*
  % params.c was changed: Corrected a bug on memory allocation in function getpars.
  % mimar.c was changed: minor changes in grid initiation for output.
  Thanks Susan J. Miller for mentionning the bug to me.
  -----------------------------

  *--- Changed May 29th 2008: ---*
  % params.c was changed: I removed several check flags to allow greater freedom
  to the user. Thanks Armando Geraldes, for mentionning the problem to me.
  -----------------------------

  *--- Changed Sept. 18th 2007 ---*
  % mimar.h was changed: Added lp.name.
  % params.c was changed: In function getpars, the option -lf was changed so that the locus
  name can start by any character (before the name could not start by a number).
  -----------------------------

**********************************************************************************/

#include <stdio.h>                            // Input output Library
#include <stdlib.h>                           // File gestion Library
#include <math.h>                             // Math Library
#include <assert.h>                           // Verify program assertion
#include <time.h>                             // Library to get CPU time
#include <string.h>
#include "mimar.h"

#define NLSTATS 4                             // Constant for 4 branch types: L1, L2, Ls, Lf
#define NPARAM 7                              // Number of parameters=theta1, theta2, tcoal, Tgen, thetaA, M12, M21
static unsigned MAX_LINE=1000;                // Number of bins in histogram

/*------------------*
 *  Structure Node  |
 *------------------*/
struct node
{
  int abv;                                    // Node above=accestor node
  int lpop1;                                  // # lineages from pop1 as descent
  int lpop2;                                  // # lineages from pop2 as descent
  float time;                                 // Time of the node
};

/*---------------------*
 *  Structure segment  |
 *---------------------*
 | Segment of chromosome issue from recombination
 | events.
*/
struct segl
{
  int beg;                                    // starting point of the segment i
  struct node *ptree;                         // Pointer on a tree the segment come from
  int next;                                   // Index number of the next segment
};

/*------------------*
 *  Structure grid  |
 *------------------*
 | Records parameter values in histogram 
*/
struct grid
{
  double paramv;                              // Value average min+(max-min)/2
  double min;                                 // Min value of bin
  double max;                                 // Max value of bin
  long double prob;                           // P(data|THETA, G)
  long double u;                              // # time recorded bin
  long double p;                              // # time value occur in prior or update
};


/*********************************************************
 *                                                       *
 * Main Function                                         *
 * -------------                                         *
 *                                                       *
 *********************************************************
 | Takes arguments and launches the MCMC analysis.       |
 *-------------------------------------------------------*/
int main(int argc, char *argv[])              // Array of char=arguments line
{
  //--- Declarations & Call Function ---//
  int i=0, howmany=0, burn=0, tag=0;          // howmany=number loci
  
  FILE *info0=NULL, *pf=NULL, *infom=NULL, *outf=NULL, *outIM=NULL, *fopen(const char*, const char*);
  /* *pfMIMAR, *pfIM (prior output) */        // Pointers on Files for outputs and input
  long double cpumin=0.0, h=0.0;
  long double *sumpu=NULL, *sump2u=NULL, *sumpu0=NULL, *sump2u0=NULL, *sumpum=NULL, *sump2um=NULL, *sumprior=NULL, *sumprior2=NULL, *infoperf=NULL;
  
  double cpusec=0.0, interval=0.0, *lstat=NULL;
  char st[200]="";
  clock_t interm=0.0;
  //// from rand1.c ////
  /* Celine changed 03/18/2010 */
  void seedit(char*, FILE*, struct params *);/*/////*/
  double ran1();
  //// From mimar.c //// 
  void writeinfo(FILE*, FILE*, struct params*, FILE*);
  void initgrid(struct grid*, struct params, int);
  void getprobdata(struct params*, double*, int, long double*, struct grid*, int, FILE*);
  long double geth(long double, long double);
  void recordmcmctest(struct params, long double*, struct grid**, struct grid**, struct grid**, struct grid*, long double*, long double*, long double*, long double*, long double*, long double*, int*, double*l, FILE*);
  void recordmcmc(struct params, long double*, struct grid**, struct grid*, long double*, long double*, int*, double*, FILE*);
  void updateparam(struct grid*, struct grid**, struct params, long double*, long double*);
  void gettime(clock_t*, long double*, double*);
  void writemimarrun(struct params, struct grid**, long double*, int, char*[], long double*, long double*, long double, double);
  int writeinterval(struct params, struct grid**, long double*, int , char*[], long double*, long double*, long double, double, int);
  void writeintervaltest(struct params, struct grid**, struct grid**, long double*, int, char*[], long double*, long double*, long double*, long double*, long double, double, int);
  void writeend(FILE*, struct params, struct grid**, long double*, int, char*[], long double*, long double*, long double, double, int);
  void writeendtest(FILE*, struct params , struct grid**, long double*, int , char*[], long double*, long double*, long double, double, int, long double);
  //// From params.c //// 
  void changeparams(struct params*);
  void proposeparamsnormal(struct params*);
  struct params getpars(int, char*[], int*);

  //--- Structure declaration---//
  struct grid **gridval=NULL, **gridval0=NULL, **gridvalm=NULL, *sumvect0=NULL, *sumvect1=NULL;
  struct params param;

  //---------- Allocation memory ------------------//
  // For the arrays bellow: 0 theta1, 1 theta2, 2 Tcoal, 3 Tgen, 4 ThetaA, 5 M12, 6 M21

  sumvect0=(struct grid*)malloc((unsigned)NPARAM*sizeof(struct grid));  // THETA: previous param
  sumvect1=(struct grid*)malloc((unsigned)NPARAM*sizeof(struct grid));  // HETA': proposed move
  gridval=(struct grid**)malloc((unsigned)NPARAM*sizeof(struct grid*)); // Histogram of accepted param.
  sumpu=(long double*)malloc((unsigned)NPARAM*sizeof(long double));     // sum u*paramv
  sump2u=(long double*)malloc((unsigned)NPARAM*sizeof(long double));    // sum u*paramv^2
  sumprior=(long double*)malloc((unsigned)NPARAM*sizeof(long double));  // prior sum u*paramv
  sumprior2=(long double*)malloc((unsigned)NPARAM*sizeof(long double)); // prior sum u*paramv^2
  /* infoperf is a table with informations about analysis:
     0 sumu , 1 sumu2, 2 # sets with P(D|theta)>0, 3 Total # sets of param,
     4 # bad sets of param:P(D|theta)=0, 5 # total genealogies generated,
     6 # Genealogies with P(D|theta)>0, 7 # gen with P(D|theta)=0, 8 # loci with P(D|theta)>0,
     9 # loci with P(D|theta)=0, 10 # steps generated,
     11 # step where param accepted with H>1, 12 ... accepted with H<1,
     13 rejected with H<1, 14 rejected with Prob=0,
     15 # stept with all param within priors, 16 outside priors,
     17 # stepts with M0, 18M>0 [test case]
  */
  infoperf=(long double*)malloc((unsigned)20*sizeof(long double));
  lstat=(double*)malloc((unsigned)NLSTATS*sizeof(double));              // 0 L1, 1 L2, 2 Ls, 3 lf
  //--- Get arguments ---//
  param=getpars(argc, argv, &howmany);                                  // Get input by user for parameters
 
  //--- Open output files and write information on analysis ---//
  if((param.cp.uniform[0][5]==1)||(param.cp.uniform[0][6]==1)) // test of migration
    {
      sprintf(st, "%sM0", param.cp.fout);                      // Output for M=0
      info0=fopen(st, "w"); 
      sprintf(st, "%sM", param.cp.fout);                       // Output for M>0
      infom=fopen(st, "w");
      gridval0=(struct grid**)malloc((unsigned)NPARAM*sizeof(struct grid*)); // Initialized extra tables for M=0 and M>0 cases
      gridvalm=(struct grid**)malloc((unsigned)NPARAM*sizeof(struct grid*));
      sumpu0=(long double*)malloc((unsigned)NPARAM*sizeof(long double));
      sump2u0=(long double*)malloc((unsigned)NPARAM*sizeof(long double));
      sumpum=(long double*)malloc((unsigned)NPARAM*sizeof(long double));
      sump2um=(long double*)malloc((unsigned)NPARAM*sizeof(long double));
    }
  /* Celine changed 11/05/2008 */ // Normal mimar
  if(param.mcmcp.relaunch==0)
    {/*/////*/
      /* sprintf(st, "%sprior", param.cp.fout); // Output for priors
         pfMIMAR=fopen(st, "w");
      */
      outf=fopen(param.cp.fout, "w");           // Results summary
      /* sprintf(st, "res%s", param.cp.fout);   // sdtoutput
         pf=fopen(st, "w");                     // redirect stdout for specific cluster
      */
      pf=stdout;
      if(param.cp.sumout==3)                    // Option for comparison with IM
        {
          sprintf(st, "%sIM", param.cp.fout);
          outIM=fopen(st, "w");
          /* sprintf(st, "%sIMprior", param.cp.fout); // Output for IM priors
             pfIM=fopen(st, "w"); 
          */
        }
      for(i=0;i<argc;i++)                       // Write argument line in summary output
        {
          fprintf(outf, "%s ", argv[i]);
          if(param.cp.sumout==3)fprintf(outIM, "%s ", argv[i]); 
        }
      /* Celine changed 03/18/2010 */
      if( !param.commandlineseedflag ) seedit("s", outf, &param);// WRITE seeds in summary output file 
      fprintf(outf,"\n%d %d %d\n",param.tableseeds[0],param.tableseeds[1],param.tableseeds[2]);
      /*/////*/
      /* writeinfo(outf, outIM, pfIM, &param); */ // Use in -q3 cases when priors output
      writeinfo(outf, outIM, &param, pf); 

      fclose(outf);
      /* Celine changed 11/05/2008 */ // relaunch mimar from interrupted posterior
    }
  else                                          // relaunch mimar from interrupted
    {
      outf=fopen(param.cp.fout,"a");            // Results summary appending
      pf=stdout;
      if(param.cp.sumout==3)                    // Option for comparison with IM
        {
          sprintf(st,"%sIM",param.cp.fout);
          outIM=fopen(st,"a");
        }
      for( i=0;i<argc;i++)                      // Write argument line in summary output
        {
          fprintf(outf,"%s ",argv[i]);
          if(param.cp.sumout==3) fprintf(outIM, "%s ", argv[i]); 
        }
      /* Celine changed 03/18/2010 */
      if( !param.commandlineseedflag ) seedit( "s",outf,&param); // WRITE seeds in summary output file
      fprintf(outf,"\n%d %d %d\n",param.tableseeds[0],param.tableseeds[1],param.tableseeds[2]);
      /*/////*/
      writeinfo(outf, outIM, &param, pf );

      fclose(outf);
    }
  /*/////*/
  //--- End Output Info ---//
  //---------- Initialisations ------------------//
  for(i=0;i<NPARAM;i++)
    {
      sumpu[i]=0.0;
      sump2u[i]=0.0;
      sumprior[i]=0.0;
      sumprior2[i]=0.0;
      gridval[i]=(struct grid *)calloc(NPARAM, (unsigned)MAX_LINE*sizeof(struct grid));
      if((param.cp.uniform[0][5]==1)||(param.cp.uniform[0][6]==1))          // Migration rate case
        {
          sumpu0[i]=0.0;
          sump2u0[i]=0.0;
          sumpum[i]=0.0;
          sump2um[i]=0.0;
          gridval0[i]=(struct grid*)calloc(NPARAM, (unsigned)MAX_LINE*sizeof(struct grid));
          gridvalm[i]=(struct grid*)calloc(NPARAM, (unsigned)MAX_LINE*sizeof(struct grid));
        }
      /* Celine Changed 03/03/2009 */ // added ||param.cp.maxval[i]!=param.cp.minval[i]
      /* Celine changed 12/28/2009 */ // removed if 
      /* if((param.cp.maxval[i]!=param.cp.minval[i])) // Case estimated parameter
         {*/
      initgrid(gridval[i], param, i);
      if((param.cp.uniform[0][5]==1)||(param.cp.uniform[0][6]==1)) // Case Mig uniform
        {
          initgrid(gridval0[i], param, i);
          initgrid(gridvalm[i], param, i);
        }
      /* }*/////    
    }// End loop on parameters
  tag=1;                  // # steps since last recorder set of parameter
  cpusec=cpumin=0;        // CPU Second and minute 
  interm=clock();         // Time interval: start of time measure
  for(i=0;i<20;i++)       // init performance table
    infoperf[i]=0.0;
  param.mcmcp.numparam=0; // num param in case -d 1 updates
  burn=0;                 // # steps of burn in, no steps recorded
  if(param.mcmcp.burnin==0)
    burn=1;
  //-- End initialisation --//
  //-------------------- Main Loop: MCMC ---------------------------//
  fprintf(pf, "Step #\ttheta1\ttheta2\tTcoal\tTgen\tthetaA\tM12\tM21\tL\n"); // Header in output
  /* Celine changed 11/05/2008 */
  if(param.mcmcp.relaunch==1)
    {
      fclose(pf);
      pf=fopen(param.mcmcp.fpost,"a");
      burn=1;
      infoperf[10]=infoperf[3]=param.mcmcp.start[0];// start number of step
      for(i=0;i<NPARAM;i++)
        {
          sumvect0[i].prob=sumvect1[i].prob=param.mcmcp.start[8];// probaility
          param.cp.newest[i]=param.cp.oldest[i]=sumvect0[i].paramv=sumvect1[i].paramv=param.mcmcp.start[i+1];
        }
    }
  /*/////*/
  while(((param.mcmcp.type=='c')&&(infoperf[10]<param.mcmcp.lenght))||
        ((param.mcmcp.type=='t')&&(cpumin<param.mcmcp.lenght)))                 // Stop when # steps=lenght, or time=lenght
    {
      infoperf[3]+=1.;               // total number of steps+1
      interval++;                    // # steps in interval. Reset after burnin and recording step in output (and histogram)
      if(infoperf[2]==0) 
        changeparams(&param);        // Sample estimated parameters from priors
      else
        proposeparamsnormal(&param); // Propose new parameters from kernel distributions
  
      if(param.mcmcp.inprior)        // Case all param values within priors: get p(D|THETA')
        {
          updateparam(sumvect1, gridval, param, sumprior, sumprior2);           // Update grid of priors
          getprobdata(&param, lstat, howmany, infoperf, sumvect1, howmany, pf); // get p(D|THETA')
          infoperf[15]+=1;                                                      // +1 set param within prior
        }
      else                           // Case a param value outside prior, record THETA
        {
          infoperf[16]+=1;
          sumvect1[0].prob=0;
        }
      if((infoperf[10]==0)&&(infoperf[2]==1)) // Initial steps: THETA0 from priors
        {
          if(sumvect1[0].prob>0)    // Accept when p(D|THETA)>0 for all loci
            {
              infoperf[12]+=1;      // 1+ accepted set of param
              if((param.cp.uniform[0][5]==1)||(param.cp.uniform[0][6]==1)) // Case test of migration
                recordmcmctest(param, infoperf, gridval, gridval0, gridvalm, sumvect1, sumpu, sump2u, sumpu0, sump2u0, sumpum, sump2um, &burn, &interval, pf);
              else                                                         // Other cases
                recordmcmc(param, infoperf, gridval, sumvect1, sumpu, sump2u, &burn, &interval, pf);
              for(i=0;i<NPARAM;i++)
                {
                  sumvect0[i]=sumvect1[i];                                  // THETA=THETA0
                  param.cp.newest[i]=param.cp.oldest[i]=sumvect1[i].paramv;
                }
            }
          interval=0;               // Reset interval 
        }
      else if(infoperf[10]>0)       // Now at THETA, propose THETA'
        {
          if(sumvect1[0].prob==0)   // Case p(D|THETA')=0, record THETA
            {
              infoperf[14]+=1;      // 1+ reject set
              if((param.cp.uniform[0][5]==1)||(param.cp.uniform[0][6]==1))  // Case test of migration
                recordmcmctest(param, infoperf, gridval, gridval0, gridvalm, sumvect0, sumpu, sump2u, sumpu0, sump2u0, sumpum, sump2um, &burn, &interval, pf);
              else                                                          // Other cases
                recordmcmc(param, infoperf, gridval, sumvect0, sumpu, sump2u, &burn, &interval, pf);
              for(i=0;i<NPARAM;i++)
                {
                  param.cp.newest[i]=param.cp.oldest[i]=sumvect0[i].paramv; // THETA=THETA
                  sumvect1[i]=sumvect0[i];
                }
            }
          else                       // Case p(D|THETA')>0 for all loci, compute h
            {
              if(sumvect0[0].prob>0) // Case p(D|THETA)>0, compute h
                h=geth(sumvect0[0].prob, sumvect1[0].prob);
              else                   // Case p(D|THETA)=0 (useless?)
                {
                  if(sumvect1[0].prob==0) h=0; // Case p(D|THETA')=0 (should never happen), record THETA
                  else h=1;                    // Case p(D|THETA')>0, record THETA'
                }
              if(h>=1)                         // Case h>=1, record THETA'
                {
                  infoperf[11]+=1;             // 1+ accepted set of param
                  if((param.cp.uniform[0][5]==1)||(param.cp.uniform[0][6]==1))  // Case test of migration 
                    recordmcmctest(param, infoperf, gridval, gridval0, gridvalm, sumvect1, sumpu, sump2u, sumpu0, sump2u0, sumpum, sump2um, &burn, &interval, pf);
                  else                                                          // Other cases
                    recordmcmc(param, infoperf, gridval, sumvect1, sumpu, sump2u, &burn, &interval, pf);
                  for(i=0;i<NPARAM;i++)
                    {
                      param.cp.newest[i]=param.cp.oldest[i]=sumvect1[i].paramv; // THETA=THETA'
                      sumvect0[i]=sumvect1[i];
                    }
                }// End case h>=1
              else                                    // Case h<1, accept THETA' with prob h
                {
                  if((long double) ran1()<=h)         // Case accept THETA'
                    {
                      infoperf[12]+=1;                // 1+ accepted set of param
                      if((param.cp.uniform[0][5]==1)||(param.cp.uniform[0][6]==1)) // Case test of migration
                        recordmcmctest(param, infoperf, gridval, gridval0, gridvalm, sumvect1, sumpu, sump2u, sumpu0, sump2u0, sumpum, sump2um, &burn, &interval, pf);
                      else                                                         // Other cases
                        recordmcmc(param, infoperf, gridval, sumvect1, sumpu, sump2u, &burn, &interval, pf);
                      for(i=0;i<NPARAM;i++)
                        {
                          param.cp.newest[i]=param.cp.oldest[i]=sumvect1[i].paramv;// THETA=THETA'
                          sumvect0[i]=sumvect1[i];
                        }
                    }
                  else                                // Case reject THETA'
                    {
                      infoperf[13]+=1;                // +1 rejected set of param
                      if((param.cp.uniform[0][5]==1)||(param.cp.uniform[0][6]==1))  // Case test of migration
                        recordmcmctest(param, infoperf, gridval, gridval0, gridvalm, sumvect0, sumpu, sump2u, sumpu0, sump2u0, sumpum, sump2um, &burn, &interval, pf);
                      else                                                          // Other cases
                        recordmcmc(param, infoperf, gridval, sumvect0, sumpu, sump2u, &burn, &interval, pf);
                      for(i=0;i<NPARAM;i++)
                        {
                          param.cp.newest[i]=param.cp.oldest[i]=sumvect0[i].paramv; // THETA=THETA
                          sumvect1[i]=sumvect0[i];
                        }
                    }// End case reject THETA'
                }// End case prob h<1
            }// End Case p(D|THETA')>0, get h
        }// End Case propose THETA'
      else interval=0;                                 // Case take THETA0 from prior again
        
      gettime(&interm, &cpumin, &cpusec);              // Compute CPU time
      if((infoperf[10]>1)&&(!((int)infoperf[10]%100))) //--- mimarrun file ---//
        writemimarrun(param, gridval, infoperf, argc, argv, sumpu, sump2u, cpumin, cpusec);
      if(param.cp.interfile>0)                         //--- Write summary in file#tag every L steps (option -L) ---//
        if(((param.mcmcp.type=='t')&&((unsigned  long int)cpumin==tag*param.cp.interfile))||((param.mcmcp.type=='c')&&(infoperf[10]==tag*param.cp.interfile)))
          {
            tag=writeinterval(param, gridval, infoperf, argc, argv, sumpu, sump2u, cpumin, cpusec, tag); 
            if((param.cp.uniform[0][5]==1)||(param.cp.uniform[0][6]==1)) // Case test of migration: write posterior for M=0 and M>0
              writeintervaltest(param, gridval0, gridvalm, infoperf, argc, argv, sumpu0, sump2u0, sumpum, sump2um, cpumin, cpusec, tag-1);
          }
    }// End Loop or MCMC analysis
  //--- output summary ---//
  outf=fopen(param.cp.fout, "a");
  gettime(&interm, &cpumin, &cpusec);  // Compute CPU time 
  writeend(outf, param, gridval, infoperf, argc, argv, sumpu, sump2u, cpumin, cpusec, 1);           // Summary output
  /* writeend(pf, param, gridval, infoperf, argc, argv, sumprior, sumprior2, cpumin, cpusec, 0); */ // Output of priors
  if(param.cp.sumout>2)                // Case -q3
    {
      writeend(outIM, param, gridval, infoperf, argc, argv, sumpu, sump2u, cpumin, cpusec, 3);
      /* writeend(pfIM, param, gridval, infoperf, argc, argv, sumprior, sumprior2, cpumin, cpusec, 2);   */ // -q3 case
    }
  if((param.cp.uniform[0][5]==1)||(param.cp.uniform[0][6]==1))// Test of migration case, summaries in file for M=0 and M>0
    {
      fprintf(info0, "No Mig\t%Lg\tMig\t%Lg\n", infoperf[17], infoperf[18]);
      fprintf(infom, "No Mig\t%Lg\tMig\t%Lg\n", infoperf[17], infoperf[18]);
      writeendtest(info0, param, gridval0, infoperf, argc, argv, sumpu0, sump2u0, cpumin, cpusec, 1, infoperf[17]);
      writeendtest(infom, param, gridvalm, infoperf, argc, argv, sumpum, sump2um, cpumin, cpusec, 1, infoperf[18]);
    }
  //--- end of program - Free memories ---// 
  /* Celine changed 03/18/2010 */
  seedit("end", outf, &param);                 // in randx.c, flag[0]!="s" so create/rewrite seed in seedmimar
  /*/////*/

  for(i=0;i<NPARAM;i++)
    free(gridval[i]);
  free(gridval);
  free(sumvect0);
  free(sumvect1);
  free(sump2u);
  free(sumpu);
  free(sumprior);
  free(sumprior2);
  free(infoperf);
  free(lstat);
  if((param.cp.uniform[0][5]==1)||(param.cp.uniform[0][6]==1))
    {
      for(i=0;i<NPARAM;i++)
        {
          free(gridval0[i]);
          free(gridvalm[i]);
        }
      free(gridval0);
      free(gridvalm);
      fclose(info0);
      fclose(infom);
      free(sump2u0);
      free(sumpu0);
      free(sump2um);
      free(sumpum);
    } 
  ///////// FREE PARAM ///////
  for(i=0;i<param.cp.npop;i++)
    free(param.cp.mig_mat[i]);
  free(param.cp.mig_mat);
  free(param.cp.config);
  /* Celine changed 11/27/2009 */
  for(i=9;i>=0;i--)
    if(param.cp.listevent[i]!=NULL && param.cp.listevent[i]->nextde!=NULL)
      free(param.cp.listevent[i]->nextde);
  if(param.cp.listevent!=NULL)
    free(param.cp.listevent);
  /*/////*/
  free(param.cp.deventlist);
  free(param.cp.size);
  free(param.cp.alphag);
  for(i=0;i<3;i++)
    free(param.cp.uniform[i]);
  free(param.cp.uniform);
  free(param.cp.oldest);
  free(param.cp.newest);
  free(param.cp.newparam);
  free(param.cp.trueval);
  free(param.mcmcp.sigma);
  free(param.mcmcp.counter);
  /* Celine changed 11/27/2009 */
  for(i=0;i<howmany;i++)
    free(param.lp[i].name);
  /*/////*/
  free(param.lp);
  fclose(outf);
  fclose(pf);
  if(param.cp.sumout>2)
    {
      // fclose(pfIM);              // IM prior output
      fclose(outIM);
    };
  /* Celine changed 03/18/2010 */
  free( param.tableseeds); /*/////*/
  return 1;
}// End of main function

/*-------------------*
 * Function initgrid |
 *-------------------*
 | Initializes the histogram for parameters distributions.
*/
void initgrid(struct grid *pgrid, struct params p, int np)
{
  int j=0;
  double min=0.0, max=0.0;
  min=p.cp.minval[np];
  max=p.cp.maxval[np];
  if(max==min)                      // Case parameter fixed. Range center of value
    {
      /*max=max+(double)max/2;      // alternative grid
        max=min+(MAX_LINE/2)*(double)min/(MAX_LINE-1);
        min=min-(MAX_LINE/2)*(double)min/(MAX_LINE-1);
      */
      max=min+(MAX_LINE/2+.05)*(double)min/(MAX_LINE-1);
      min=min-(MAX_LINE/2-.05)*(double)min/(MAX_LINE-1);
    }

  if(((np==5)||(np==6))&&(p.cp.uniform[0][np]==1))      // Case test of migration:
    {
      pgrid[0].paramv=0;            // First bin is 0 when RJ
      pgrid[0].u=0.0;
      pgrid[0].p=0.0;
      for(j=1;j<MAX_LINE;j++)
        {
          pgrid[j].max=min+(j)*(max-min)/(MAX_LINE-1);
          pgrid[j].min=min+(j-1)*(max-min)/(MAX_LINE-1);
          pgrid[j].paramv=pgrid[j].min+(double)(pgrid[j].max-pgrid[j].min)/2;
          pgrid[j].u=0.0;
          pgrid[j].p=0.0;
        }
    }
  else if(((np==5)||(np==6))&&(p.cp.uniform[0][np]==5)) // Case for ln(M)~uniform
    {
      for(j=0;j<MAX_LINE;j++)
        {
          pgrid[j].max=exp(min+(j+1)*(max-min)/MAX_LINE); // Bin exponentially distributed
          pgrid[j].min=exp(min+(j)*(max-min)/MAX_LINE);
          pgrid[j].paramv=pgrid[j].min+(double)(pgrid[j].max-pgrid[j].min)/2;
          pgrid[j].u=0.0;
          pgrid[j].p=0.0;
        }
    }
  else                                                  // All other cases
    {
      for(j=0;j<MAX_LINE;j++)
        {
          /*pgrid[j].paramv=min+(j)*(max-min)/(MAX_LINE-1); // alternative grid
            pgrid[j].max=pgrid[j].paramv+(max-min)/((MAX_LINE-1)*2);
            pgrid[j].max=pgrid[j].paramv-(max-min)/((MAX_LINE-1)*2);
          */
          pgrid[j].max=min+(j+1)*(max-min)/(MAX_LINE-1);
          pgrid[j].min=min+(j)*(max-min)/(MAX_LINE-1);
          pgrid[j].paramv=pgrid[j].min+(double)(pgrid[j].max-pgrid[j].min)/2;
          pgrid[j].u=0.0;
          pgrid[j].p=0.0;
        }
    }
}// End Initgrid


/*************************** Function for likelihood estimation ***************************/
/*-----------------*
 * Function getsam |
 *-----------------*
 | Generates a gene genealogy for a locus.
*/
void gensam(struct params *pparam, double *plstat)
{
  int nsegs=0, i=0, k=0, seg=0, ns=0, start=0, end=0, len=0 , nsites=0, nsam=0, nsam1=0, nsam2=0;
  /* nsegs is the number of segments the gametes were broken into.
   * The history of the gametes is tracked back. The histories of these segments are
   * passed back to the calling function in the array of structures seglst[]
   */
  struct segl *seglst=NULL;
  double theta=0.0, *l=NULL;
  struct segl *segtre_mig(struct c_params *p, int *nsegs);
  void computelstat(int, int, int, struct node*, double*);
  //// From streec.c //// 
  void prtree(struct node*, int);

  nsites=pparam->cp.nsites;                              // Locus lenght in bp
  seglst=segtre_mig(&(pparam->cp), &nsegs);              // Get ARGs, record # segments produced by recombination
  nsam=pparam->cp.nsam;                                  // Total # chromosomes for this locus
  nsam1=pparam->cp.config[0];                            // nsam1 and 2: # seq in pop 1 and 2
  nsam2=pparam->cp.config[1];

  l=(double *)malloc((unsigned)NLSTATS *sizeof(double)); // Array of branch lenght L1, L2, Ls, lf for a segment

  if(pparam->mp.treeflag)                                // Case output tree.
    {
      theta=pparam->mp.theta;
      ns=0;
      for(seg=0, k=0;k<nsegs;seg=seglst[seg].next, k++)
        {
          if((pparam->cp.r>0.0)||(pparam->cp.f>0.0))
            {
              end=(k<nsegs-1 ? seglst[seglst[seg].next].beg-1:nsites-1);
              start=seglst[seg].beg;
              len=end-start+1;
              fprintf(stdout, "%d\t", len);
            }
          prtree(seglst[seg].ptree, nsam);
          if((theta==0.0)) 
            free(seglst[seg].ptree);
        }
    }
  for(seg=0, k=0;k<nsegs;seg=seglst[seg].next, k++)                //--- Loop along all segments ---//
    {
      end=(k<nsegs-1 ? seglst[seglst[seg].next].beg-1:nsites-1); // End of segment
      start=seglst[seg].beg;                                       // Beginning of segment
      len=end-start+1;                                             // Lengh of segment
      computelstat(nsam, nsam1, nsam2, seglst[seg].ptree, l);      // Get L statistics for this segment
      for(i=0;i<NLSTATS;i++)                                       //--- Compute total L_k: mean weighted on segment lenghts ---//
        plstat[i]+=(double)(l[i]*len)/nsites;
      free(seglst[seg].ptree);
    }
  free(l);
}// End Gensam



/*----------------------*
 * Function getprobdata |
 *----------------------*
 | Computes the probability of the statistics given THETA and a set of genealogies.
*/
void getprobdata(struct params *pparam, double *plstat, int nloci, long double *pinfoperf, struct grid *psumvect , int phowmany, FILE *ppf)
{
  int i=0, j=0, nl=0, oki=0, oks=0; 
  long double u_g=0.0, lnprob=0.0, logl(long double), expl(long double);
  double lnps=0, lnpo(double, int);
  void changeparamslocus(struct params*, int);
  void getsam(struct params*, double*);
  void updatemainparams(struct params*, struct grid *);
  updatemainparams(pparam, psumvect); // Update the coalescent parameter
  oki=1;
  pinfoperf[8]=0.0, pinfoperf[9]=0.0, lnprob=0.0;
  for(nl=0;nl<nloci;nl++)
    {
      changeparamslocus(pparam, nl);  // Get locus specific parameters
      u_g=0.0;
      for(j=0;j<pparam->cp.ngen;j++)  //--- Loop on genealogies for a locus ---//
        {
          for(i=0;i<NLSTATS;i++) plstat[i]=0; // Initialization of Lstats
          gensam(pparam, plstat);     // Generate a new ARG
          pinfoperf[5]+=1.;           // 1+ genealogy
          lnps=0;
          oks=1;
          for(i=0;i<NLSTATS;i++)
            {
              if((plstat[i]==0)&&(pparam->lp[nl].Si[i]!=0)) // Case when ARG doesn't fit the data
                {
                  oks=0;
                  break;
                }
              else                                          // Case, ARG fit data: get P(Data|THETHA, G)
                lnps+=lnpo(pparam->mp.theta*plstat[i], pparam->lp[nl].Si[i]);
            }
          if(lnps<-10000)             // Test: make sure Likelihood large enough
            {
              oks=0;
              break;
              fprintf(ppf, "! a u too small: lnps\t%lg\n", lnps);
            }
          if(oks)                     // Case: L_ks>0, get mean over genealgies for the locus
            {
              u_g+=(long double) expl((long double) lnps)/pparam->cp.ngen;
              pinfoperf[6]+=1.;       // +1 ARG fit data
            }
          else                        // Case: one L_k=0. add 0 to likelihood
            pinfoperf[7]+=1.;         // 1+ ARG not fit
        }// End of loop on ngen
      if((!oks)&&(u_g==0))            // Case: ngen ARG not fitting
        {
          pinfoperf[9]+=1.;           // +1 locus with unfit data
          oki=0;                      // Tag: one locus not fitting, stop here
          break;
        }
      else                            // Case: ug>0
        {
          lnprob+=logl(u_g);
          pinfoperf[8]+=1.;           // +1 fitting locus
        }
      if(lnprob<-10000)               // Test likelihood large enough 
        {
          oki=0;
          break;
        }
    }// End loop on loci

  if((!oki)&&((lnprob==0)||(pinfoperf[8]!=phowmany))) // Case one locus not fitting, reject THETA
    {
      pinfoperf[4]+=1.;
      for(i=0;i<NPARAM;i++)
        psumvect[i].prob=0;
    }
  else                                // All likelihood>0 for all loci, get P(D|THETA, G)
    {
      pinfoperf[2]+=1.;
      for(i=0;i<NPARAM;i++)
        psumvect[i].prob=expl(lnprob);
    }
}// End getprobdata function


/*----------------*
 * Function lnpo  |
 *----------------*
 | Computes log of prob S from Poisson distribution with mean L*theta. log(prob(S=y|lamda=u))
*/
double lnpo(double u, int y)
{
  int i=0;
  double lni=0.0, lnp=0.0, log(double);
  for(i=1;i<=y;i++)
    lni+=log(i);
  if(u>0)
    lnp=-u+y*log(u)-lni;
  else if((u==0) &&(y==0)) lnp=0;
  return(lnp);
}// End function lnpo


/*------------------------*
 * Function computelstat  |
 *------------------------*
 | Computes the branch lenghts given rise to the 4 statistics in an ARG. Called in
 | gensam function.
*/
void computelstat(int nsam, int nsam1, int nsam2, struct node *pptree, double *pl)
{
  int i=0, abv=0, lpop1=0, lpop2=0, lpop1a=0, lpop2a=0;
  double timd=0.0;
  for(i=0;i<NLSTATS;i++) pl[i]=0;
  for(i=0;i<2*nsam-2;i++)           //--- Loop along all nodes of tree ---//
    {
      lpop1=pptree[i].lpop1;
      lpop2=pptree[i].lpop2;
      abv=pptree[i].abv;
      lpop2a=pptree[abv].lpop2;
      lpop1a=pptree[abv].lpop1;
      timd=pptree[abv].time-pptree[i].time; // Branch lenght
      if((lpop1==0)||(lpop2==0))
        {
          if(((lpop1==nsam1)&&(lpop2a<=nsam2))||((lpop2==nsam2)&&(lpop1a<=nsam1)))
            pl[3]+=timd;            // Case Lf 
          else                      // Case L1, L2
            {
              if((lpop2==0)&&(lpop1<nsam1)) pl[0]+=timd;
              else pl[1]+=timd;
            }
        }
      else
        pl[2]+=timd;                // Case Ls
    }// End loop Lcompute 
}// End computeLstat


/************************* Functions for updating coalescent parameters and estimates *************************/
/*---------------------------*
 * Function updatemainparam  |
 *---------------------------*
 | Updates the main parameters required for coalescent.
*/
void updatemainparams(struct params *pp, struct grid *pvect)
{
  int j=0;
  for(j=0;j<NPARAM;j++)
    pp->cp.newparam[j]=0;

  pp->mp.thetafix=pvect[0].paramv;             //--- Theta1 ---//
  pp->cp.newparam[0]=pvect[0].paramv;

  if(pp->cp.sumout>=1)                         //--- Theta2 and ThetaA ---//
    {
      pp->cp.newparam[1]=pvect[1].paramv/pp->mp.thetafix;         // b=theta2/theta1
      pp->cp.newparam[4]=(double)pvect[4].paramv/pp->mp.thetafix; // a=thetaA/theta1
    }
  else
    {
      pp->cp.newparam[1]=pvect[1].paramv;                         // b estimated
      pp->cp.newparam[4]=pvect[4].paramv;                         // a estimated
    }
  //--- T ---//
  pp->cp.newparam[2]=pvect[2].paramv; // Tcoalescent
  pp->cp.newparam[3]=pvect[3].paramv; // Tgen
  //-- Migration --//
  pp->cp.newparam[5]=pvect[5].paramv; // M12
  pp->cp.newparam[6]=pvect[6].paramv; // M21
}// End change Coalescent param


/*-----------------------*
 * Function updateparam  |
 *-----------------------*
 | Updates parameter estimates in histograms. Get prior probability as well.
*/
void updateparam(struct grid *psumvect, struct grid  **pgridval, struct params p, long double *psumprior, long double *psumprior2)
{
  int i=0, j=0;
  void transdoublePtoO(double *, struct params);
  double *prec=NULL;

  prec=(double*)malloc((unsigned)(NPARAM)*sizeof(double)); // Proposed estimates;

  for(i=0;i<NPARAM;i++)
    prec[i]=p.cp.newest[i];
  if(p.cp.sumout==1)                // Option -q1 prec in coal, output in pop parameters
    transdoublePtoO(prec, p);

  for(i=0;i<NPARAM;i++)             //--- Histogram of parameters ----//
    {
      psumvect[i].paramv=prec[i];   // Get param, param*u, and param^2*u for mean and var
      psumprior[i]+=psumvect[i].paramv;
      psumprior2[i]+=psumvect[i].paramv*psumvect[i].paramv;
      if(p.cp.maxval[i]>0)          // Put value into corresponding histogram bin
        {
          j=0;
          while((psumvect[i].paramv>pgridval[i][j].max))
            {
              j++;
              if(j>=MAX_LINE) break;
            }
          if(j<MAX_LINE) pgridval[i][j].p+=1.;
          else                      // Case last bin
            {
              j--;
              pgridval[i][j].p+=1.;
            }
        }
    }
  free(prec);
}// End update param in prior grid

/************************* Functions for MCMC recording *************************/
/*-----------------*
 * Function get h  |
 *-----------------*
 | Computes Hasting Ratio
*/
long double geth(long double p0, long double p1)
{
  if((long double) p1/p0<1)
    return(long double) p1/p0;
  else
    return 1.;
}// End get h function


/*---------------------*
 * Function recordmcmc |
 *---------------------*
 | Records set of param in intervals after burnin in histogram and output, and in summary file every L steps.
*/
void recordmcmc(struct params pparam, long double *pinfoperf, struct grid **pgridval, struct grid *psumvect, long double *psumpu, long double *psump2u, int *pburn, double *pinterval, FILE *ppf)
{
  void recordparam(struct params, long double*, struct grid**, struct grid*, long double*, long double*, int , FILE*);
  if((!*pburn)&&(*pinterval==pparam.mcmcp.burnin))    // Case: end burnin: 1st recording
    {
      *pburn=1;
      *pinterval=0;                                     // reset interval
      recordparam(pparam, pinfoperf, pgridval, psumvect, psumpu, psump2u, 1, ppf); 
    }
  else if((*pburn)&&((int)*pinterval==pparam.mcmcp.interrec)) // Record every i steps (-i option)
    {
      *pinterval=0;                                     // reset interval
      recordparam(pparam, pinfoperf, pgridval, psumvect, psumpu, psump2u, 1, ppf); 
    }
  else if(*pburn)                                     // Case, within Burnin: record in prior only
    recordparam(pparam, pinfoperf, pgridval, psumvect, psumpu, psump2u, 0, ppf);
  else
    pinfoperf[10]+=1.;                                 // 1+ sample analyzed
}// End record mcmc


/*----------------------*
 * Function recordparam |
 *----------------------*
 | Records set of param in every interval after burnin in histogram.
*/
void recordparam(struct params pparam, long double *pinfoperf, struct grid **pgridval, struct grid *psumvect, long double *psumpu, long double *psump2u, int rec, FILE *ppf)
{
  int i=0, j=0;
  pinfoperf[10]+=1.;                // 1+ sample analyzed
  pinfoperf[0]+=psumvect[0].prob;   // sum ui
  pinfoperf[1]+=(long double)psumvect[0].prob*psumvect[0].prob; // sum ui^2
  if(rec) fprintf(ppf, "%Lg\t", pinfoperf[10]);                 // Record in output file
  for(i=0;i<NPARAM;i++)             // Loop on param: record in histogram
    {
      psumvect[i].u=1;              // 1+ accepted param in this bin
      psumpu[i]+=psumvect[i].paramv;
      psump2u[i]+=psumvect[i].paramv*psumvect[i].paramv;
      if(rec) fprintf(ppf, "%lg\t", psumvect[i].paramv);
      /* Celine change 12/28/2009 */
      /*  if(pparam.cp.maxval[i]>0)
          {*/
      j=0;
      while((psumvect[i].paramv>pgridval[i][j].max))
        {
          j++;
          if(j>=MAX_LINE) break;
        }
      if(j<MAX_LINE) pgridval[i][j].u+=psumvect[i].u;
      else
        {
          j--;
          pgridval[i][j].u+=psumvect[i].u;
        }
      /* }*//////
    }
  if(rec) fprintf(ppf, "%Lg\n", psumvect[0].prob);
}// End record param in histograms


/*********************** Functions for writing informations on analysis ************************/
/*--------------------*
 * Function writeinfo |
 *--------------------*
 | Writes information on analysis
*/
/* void writeinfo(FILE *poutf, FILE *poutIM, FILE *ppfIM, struct params *pparam) */ // Use for -q3
void writeinfo(FILE *poutf, FILE *poutIM, struct params *pparam, FILE *ppf )
{
  void  transdoublePtoO(double [NPARAM], struct params);
  int i=0;
  char *namep="";
  fprintf(ppf, "Analysis of file %s\n", pparam->cp.fin);
  //--- Output Information for different options ---//
  if(pparam->cp.sumout==0)          // theta1, b, Tcoal, Tgen, a, M12, M21
    {
      fprintf(ppf, "\nPrior distributions ranges/parameter values - Option 0 (-1 means that no information was provided. If both values are the same, the parameter was fixed by the user)\n");
      for(i=0;i<NPARAM;i++)
        {
          if(i==0) namep="theta1";
          if(i==1) namep="b=theta2/theta1";
          if(i==2) namep="Tcoal";
          if(i==3) namep="Tgen";
          if(i==4) namep="a=thetaA/theta1";
          if(i==5) namep="M12";
          if(i==6) namep="M21";
          fprintf(ppf, "%s\t%lg\t%lg\n", namep, pparam->cp.minval[i], pparam->cp.maxval[i]); 
        }
      /* Uncommented for Celine */
      /*if(pparam->cp.trueval[0]>=0)
        {
        fprintf(ppf, "\nSimulations values\n");
        for(i=0;i<NPARAM;i++)       // theta1, b, Tcoal, Tgen, a, M12, M21
        fprintf(ppf, "%lg\t", pparam->cp.trueval[i]);
        }
      *//////
      fprintf(ppf, "\n");
    }
  if(pparam->cp.sumout==1)          // Case theta, b->theta2, Tcoal, Tgen, a->thetaA, M12, M21
    {
      fprintf(ppf, "\nPrior distributions ranges/parameter values provided by user - Option 1 (-1 means that no information was provided. If both values are the same, the parameter was fixed by the user)\n");
      for(i=0;i<NPARAM;i++)         // theta1, b, Tcoal, Tgen, a, M12, M21
        {
          if(i==0) namep="theta1";
          if(i==1) namep="b=theta2/theta1";
          if(i==2) namep="Tcoal";
          if(i==3) namep="Tgen";
          if(i==4) namep="a=thetaA/theta1";
          if(i==5) namep="M12";
          if(i==6) namep="M21";
          fprintf(ppf, "%s\t%lg\t%lg\n", namep, pparam->cp.minval[i], pparam->cp.maxval[i]);
        }
      transdoublePtoO(pparam->cp.maxval, *pparam);
      transdoublePtoO(pparam->cp.minval, *pparam);

      fprintf(ppf, "\nPrior distributions ranges/parameter values for estimates- Option 1 (-1 means that no information was provided. If both values are the same, the parameter was fixed by the user)\n");
      for(i=0;i<NPARAM;i++)         // theta1, theta2, Tcoal, Tgen, thetaA, M12, M21
        {
          if(i==0) namep="theta1";
          if(i==1) namep="theta2";
          if(i==2) namep="Tcoal";
          if(i==3) namep="Tgen";
          if(i==4) namep="thetaA";
          if(i==5) namep="M12";
          if(i==6) namep="M21";
          fprintf(ppf, "%s\t%lg\t%lg\n", namep, pparam->cp.minval[i], pparam->cp.maxval[i]);
        }
      /* Uncommented for Celine */
      /*if(pparam->cp.trueval[0]>=0)
        {
        fprintf(ppf, "\nSimulation Values\n");
        for(i=0;i<NPARAM;i++)                          // theta1, b, Tcoal, Tgen, a, M12, M21
        fprintf(ppf, "%lg\t", pparam->cp.trueval[i]);
        transdoublePtoO(pparam->cp.trueval, *pparam); // theta, b->theta2, Tcoal, Tgen, a->thetaA, M12, M21
        }
      *//////
      fprintf(ppf, "\n");
    }
  if(pparam->cp.sumout==2)          // Case theta, theta2, Tcoal, Tgen, thetaA, M12, M21
    {
      fprintf(ppf, "\nPrior distributions ranges/parameter values (-1 means that no information was provided. If both values are the same, the parameter was fixed by the user)\n");
      for(i=0;i<NPARAM;i++)
        {
          if(i==0) namep="theta1";
          if(i==1) namep="theta2";
          if(i==2) namep="Tcoal";
          if(i==3) namep="Tgen";
          if(i==4) namep="thetaA";
          if(i==5) namep="M12";
          if(i==6) namep="M21";
          fprintf(ppf, "%s\t%lg\t%lg\n", namep, pparam->cp.minval[i], pparam->cp.maxval[i]);
        }
      /* Uncommented for Celine */
      /* if(pparam->cp.trueval[0]>=0)
         {
         fprintf(ppf, "\nSimulation Values\n");
         for(i=0;i<NPARAM;i++)                            // theta1, b, Tcoal, Tgen, a, M12, M21
         fprintf(ppf, "%lg\t", pparam->cp.trueval[i]);
         // transdoublePtoO(pparam->cp.trueval, *pparam); // theta1, b->theta2, Tcoal, Tgen, a->thetaA, M12, M21
         }
      *//////
      fprintf(ppf, "\n");
    }
  if(pparam->cp.sumout>=3)          // Case (IM comp) theta, theta2, Tcoal, Tgen, thetaA, m12, m21
    {
      double mulx[NPARAM];
      for(i=0;i<NPARAM;i++)
        {
          mulx[i]=1;                                // Multiplicative factore for Im comparison
          if((i==0)||(i==1)||(i==4)) mulx[i]=pparam->lp[0].li;             // pop mutation rates per locus
          else if(i==3) mulx[i]=pparam->mp.mu*pparam->lp[0].li;            // time
          else if((i==5)||(i==6)) mulx[i]/=pparam->lp[0].li*pparam->mp.mu; // m12 and m21
        }
      /* fprintf(ppfIM, "Prior\n"); */
      /* Uncommented for Celine */
      /*if(pparam->cp.trueval[0]>=0) 
        {
        fprintf(ppf, "\nSimulation Values\n");
        for(i=0;i<NPARAM;i++)                            // theta1, b, Tcoal, Tgen, a, M12, M21
        fprintf(ppf, "%lg\t", pparam->cp.trueval[i]);
        // transdoublePtoO(pparam->cp.trueval, *pparam); // theta, b->theta2, Tcoal, Tgen, a->thetaA, M12, M21}
        *//////
      fprintf(ppf, "\nPrior distributions ranges/parameter values - Option 3 (-1 means that no information were given. If both value equal, the parameter was fixed by the user)\n");
      fprintf(poutIM, "Prior distributions ranges/parameter values in IM units - Option 3 (<=1: no information given)\n");
      for(i=0;i<NPARAM;i++)         // theta1*l, theta2*l, Tcoal, Tim, thetaA->*l, M12-?m12, M21-?m21
        {
          if(i==0) namep="theta1";
          if(i==1) namep="theta2";
          if(i==2) namep="Tcoal";
          if(i==3) namep="Tgen";
          if(i==4) namep="thetaA";
          if(i==5) namep="m12";
          if(i==6) namep="m21";
          fprintf(ppf, "%s\t%lg\t%lg\n", namep, pparam->cp.minval[i], pparam->cp.maxval[i]);

          if(i==0) namep="q1";
          if(i==1) namep="q2";
          if(i==2) namep="Tcoal";
          if(i==3) namep="Tcoal*(mu*l)";
          if(i==4) namep="qA";
          if(i==5) namep="m12/(mu*l)";
          if(i==6) namep="m21/(mu*l)";
          fprintf(poutIM, "%s\t%lg\t%lg\n", namep, pparam->cp.minval[i]*mulx[i], pparam->cp.maxval[i]*mulx[i]);
          /* fprintf(ppfIM, "%d\t%lg\t%lg\n", i, pparam->cp.minval[i]*mulx[i], pparam->cp.maxval[i]*mulx[i]);*/
        }
      /* Uncommented for Celine */
      /*if(pparam->cp.trueval[0]>=0)
        {
        fprintf(ppf, "\nSimulation Values\n");                // theta1, b, Tcoal, Tgen, a, M12, M21
        fprintf(poutIM, "\nSimulation Values in IM units\n"); // q1, q2, Tcoal, Tgen, qa, m12, m21
        for(i=0;i<NPARAM;i++)
        {
        fprintf(ppf, "%lg\t", pparam->cp.trueval[i]);
        fprintf(poutIM, "%lg\t", pparam->cp.trueval[i]*mulx[i]);
        // fprintf(ppfIM, "%lg\t", pparam->cp.trueval[i]*mulx[i]); // IM prior
        }
        }
      *//////
      fprintf(ppf, "\n");
     
    }
  /* Uncommented for Celine */
  /*if(pparam->cp.trueval[0]>=0)
    {
    fprintf(poutf, "\nmu=%lg **\n", pparam->mp.mu);
    fprintf(ppf, "\nSimulation Value - in estimation unit\n");
    for(i=0;i<NPARAM;i++) 
    {
    fprintf(ppf, "%lg\t", pparam->cp.trueval[i]);
    fprintf(poutf, "%lg\t", pparam->cp.trueval[i]);
    }
    fprintf(poutf, "\n");
    fprintf(ppf, "\n");
    }
  *//////
}// End write info file

/*********************** Functions for writing performance of run ************************/
/*-------------------*
 * Function gettime  |
 *-------------------*
 | Computes CPU time of run
*/
void gettime(clock_t *pinterm, long double *pmin, double *psec)
{
  clock_t end=0.0;
  end=clock();
  *psec=((double)(end-*pinterm))/CLOCKS_PER_SEC;
  if(*psec>=60)                     // count minutes elapsed
    {
      *pinterm=clock();
      *pmin+=(long double)*psec/60;
      *psec-=(60*(int)*psec/60);
    }
}

/*-------------------*
 * Function getESS  |
 *-------------------*
 | Computes ESS. Was used in IS method
*/
long double getess(long double *pinfoperf)
{
  long double avu=0.0, cv2=0.0;

  avu=pinfoperf[0]/pinfoperf[3];    // mean u
  cv2=pinfoperf[1];
  cv2*=pinfoperf[3];
  cv2-=pinfoperf[0]*pinfoperf[0];
  cv2/=pinfoperf[3]*(pinfoperf[3]-1);
  cv2/=avu*avu;
  return(long double) pinfoperf[3]/(1.0+cv2);
}


/*-----------------------------------------*
 * Function write performance of analysis  |
 *-----------------------------------------*/
void writeperf(FILE *poutf, long double *pinfoperf, long double pess, long double min, double sec)
{
  /* fprintf(poutf, "\nTotRecord\tAcceptH>1\tAcceptH<1\trate\tReject<1\tRejProb==0\tGoodP\tBadP\tRate\tsumu\tsumu2\tGoodProb\tTotSampled\tBadSample\tfrac\t#gen\tGoodgen\tBadGen\tfracg\tESS\t\tcpumin\tcpusec");
     fprintf(poutf, "\n%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%.15Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%lu\t%Lg\t%lg\n", pinfoperf[10], pinfoperf[11], pinfoperf[12], (pinfoperf[11]+pinfoperf[12])/pinfoperf[10], pinfoperf[13], pinfoperf[14], pinfoperf[15], pinfoperf[16], pinfoperf[16]/pinfoperf[10], pinfoperf[0], pinfoperf[1], pinfoperf[2], pinfoperf[3], pinfoperf[4], pinfoperf[4]/pinfoperf[3], pinfoperf[5], pinfoperf[6], pinfoperf[7], pinfoperf[7]/pinfoperf[5] , pess, (unsigned long int)pess, min, sec);
  */
  fprintf(poutf, "\n# generated steps.\tAccepted H>1.\tAccepted H<1.\tAcceptance rate.\tRejected H<1.\tRejected L==0.\tRejected outside prior.\tcpumin.\tcpusec.");
  fprintf(poutf, "\n%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%lg\n", 
          pinfoperf[10], pinfoperf[11], pinfoperf[12], (pinfoperf[11]+pinfoperf[12])/pinfoperf[10],
          pinfoperf[13], pinfoperf[14]-pinfoperf[16], pinfoperf[16] , min, sec);
}// End writeperf function


/*********************** Functions for writing summary histograms of posteriors ************************/
/*--------------------------*
 * Function write mimarrun  |
 *--------------------------*
 | Checks if user wants to output summary file
*/
void writemimarrun(struct params pparam, struct grid **pgridval, long double *pinfoperf, int pargc, char *pargv[], long double *psumpu, long double *psump2u, long double min, double sec)
{
  int i=0;
  char st[200]="";
  FILE *ismlf=NULL, *fopen(const char*, const char*);
  void writeend(FILE*, struct params, struct grid**, long double*, int, char*[], long double*, long double*, long double, double, int);
  ismlf=fopen("mimarrun", "r");
  if(ismlf!=NULL)                   // Check mimrrun present and ready to write in
    {
      fscanf(ismlf, "%s", st);
      if(st[0]=='y')                // write summary histograms in file
        {
          fclose(ismlf);
          ismlf=fopen("mimarrun", "w");
          fprintf(ismlf, "No\n"); 
          for(i=0;i<pargc;i++)
            fprintf(ismlf, "%s ", pargv[i]);
          writeend(ismlf, pparam, pgridval, pinfoperf, pargc, pargv, psumpu, psump2u, min, sec, 1);
        }
      fclose(ismlf);
    }
}// End write in mimarrun function

/*--------------------------*
 * Function write interval  |
 *--------------------------*
 | Writes histograms of posteriors every L steps (defines by option -L)
*/
int writeinterval(struct params pparam, struct grid **pgridval, long double *pinfoperf, int pargc, char *pargv[], long double *psumpu, long double *psump2u, long double min, double sec, int ta)
{
  int i=0;
  long double ess=0.0;
  char st[200]="";
  FILE *ismlf=NULL, *fopen(const char*, const char*);
  void writesum(FILE*, struct grid**, long double, long double, long double*, long double*, struct params, int);
  void writeperf(FILE*, long double*, long double, long double, double);
  long double getess(long double*);

  ess=getess(pinfoperf);
  sprintf(st, "%s-%d", pparam.cp.fout, ta);
  ismlf=fopen(st, "w");
  for(i=0;i<pargc;i++)
    fprintf(ismlf, "%s ", pargv[i]);
  writeperf(ismlf, pinfoperf, ess, min, sec);
  writesum(ismlf, pgridval, pinfoperf[10]-pparam.mcmcp.burnin, pinfoperf[10]-pparam.mcmcp.burnin, psumpu, psump2u, pparam, 1);
  fclose(ismlf);
  return ta+1;                      // Tag: number of time output between L steps (option -L)
}// End write interval

/*---------------------*
 * Function write end  |
 *---------------------*
 | Writes Final histograms of posteriors at the end of analysis
*/
void writeend(FILE *poutf, struct params pparam, struct grid **pgridval, long double *pinfoperf, int pargc, char *pargv[], long double *psumpu, long double *psump2u, long double min, double sec, int type)
{
  long double pess=0.0;
  void writesum(FILE*, struct grid**, long double, long double, long double*, long double*, struct params, int);
  void writeperf(FILE*, long double*, long double, long double, double);
  long double getess(long double*);
  pess=getess(pinfoperf);
  writeperf(poutf, pinfoperf, pess, min, sec);
  if((type==1)||(type==3))
    writesum(poutf, pgridval, pinfoperf[10]-pparam.mcmcp.burnin, pinfoperf[10]-pparam.mcmcp.burnin, psumpu, psump2u, pparam, type);
  else
    writesum(poutf, pgridval, pinfoperf[3], pinfoperf[3], psumpu, psump2u, pparam, type);
}// End writeend

/*---------------------*
 * Function writesum  |
 *---------------------*
 | Writes summary histograms of posteriors.
*/
void writesum(FILE *sumf, struct grid **pgridv, long double pcount, long double psumu, long double *psumpu, long double *psump2u, struct params pparam, int type)
{
  int i=0, j=0, k=0, nperc[NPARAM];
  double *mulx=NULL, *mode=NULL, **percvec=NULL, alpha[NPARAM];
  long double tempu=0.0, nomig12=0.0, nomig21=0.0;
  long double *mean=NULL, *var=NULL, *maxu=NULL, *sumui=NULL;

  /* Celine changed 12/28/2009 - calloc(1, to malloc( */
  sumui=(long double*)malloc((unsigned)NPARAM*sizeof(long double)); // sum ui for percentiles
  mean=(long double*)malloc((unsigned)NPARAM*sizeof(long double));  // mean estimates
  mode=(double*)malloc((unsigned)NPARAM*sizeof(double));            // mode estimates

  maxu=(long double*)malloc((unsigned)NPARAM*sizeof(long double));  // sum ui for mode 
  var=(long double*)malloc((unsigned)NPARAM*sizeof(long double));   // variances
  mulx=(double*)malloc((unsigned)NPARAM*sizeof(double));            // multiplicative value for tranformations (-q1 and 3)
  percvec=(double**)malloc((unsigned)NPARAM*sizeof(double*));       // .05 .5 .95 percentiles
  for(i=0;i<NPARAM;i++)            //--- Init arrays ---//
    {
      mode[i]=pgridv[i][0].paramv; // Init mode: first bin
      if((type==1)||((type==3)))   // u for mode if prior or not
        tempu=pgridv[i][0].u;             
      else
        tempu=pgridv[i][0].p;
      maxu[i]=tempu;               // init u for mode update 
      mean[i]=0;
      var[i]=0.;
      sumui[i]=0;
      percvec[i]=(double*)calloc(NPARAM, (unsigned) 3*sizeof(double));
      for(k=0;k<3;k++)
        percvec[i][k]=-1;

      mulx[i]=1;
      if(type>=2)
        {
          if((i==0)||(i==1)||(i==4)) mulx[i]=pparam.lp[0].li;            // thetas-> qs
          else if(i==3) mulx[i]=pparam.mp.mu*pparam.lp[0].li;            // t
          else if((i==5)||(i==6)) mulx[i]/=pparam.lp[0].li*pparam.mp.mu; // M into m12 and m21
        }  
      nperc[i]=0;
      alpha[i]=.025;
    }

  fprintf(sumf, "\ntheta1\t\ttheta2\t\tTcoal\t\tTgen\t\tthetaA\t\tM12\t\tM21\t\n"); // Header for histograms

  nomig12=nomig21=0;                 // sum u for M=0 and M>0
  for(j=0;j<MAX_LINE;j++)            // Loop on 1000 bins
    {
      for(i=0;i<NPARAM;i++) 
        {
          if((type==1)||((type==3))) // case posterior
            tempu=pgridv[i][j].u;
          else                       // case prior
            tempu=pgridv[i][j].p;
          if(pgridv[i][0].paramv!=pgridv[i][1].paramv)// Do not output M21 if symmatrical M
            fprintf(sumf, "%.15f\t%Lg\t", pgridv[i][j].paramv*mulx[i], (long double) tempu/psumu);
          else
            fprintf(sumf, "\t\t");
          sumui[i]+=tempu;
          if(tempu>maxu[i])          //--- get the mode ---//
            {
              mode[i]=pgridv[i][j].paramv;
              maxu[i]=tempu; 
            }
          if((i==5)&&(pgridv[i][j].paramv==0.0)) nomig12=sumui[i]/psumu; // get sum ui for M=0
          if((i==6)&&(pgridv[i][j].paramv==0.0)) nomig21=sumui[i]/psumu;
      
          if((type==1)||((type==3)))tempu=pgridv[i][j+1].u;              // Get p ui for next bin
          else tempu=pgridv[i][j+1].p;
      
          if((sumui[i]+tempu)/psumu>alpha[i])                            // Find percentile
            {
              percvec[i][nperc[i]]=pgridv[i][j].paramv;
              alpha[i]+=.475;                                            // Go to next percentile
              nperc[i]++;
              while((double)(sumui[i]+tempu)/psumu>alpha[i])
                {
                  percvec[i][nperc[i]]=pgridv[i][j].paramv;
                  alpha[i]+=.475;
                  if(alpha[i]>=1) break;
                  nperc[i]++;
                }
              if((sumui[i]/psumu<.975)&&(percvec[i][0]==percvec[i][2])&&((sumui[i]+tempu)/psumu>=.975))
                {
                  for(k=0;k<3;k++)   // Special case, all weight in 1 bin
                    percvec[i][k]=pgridv[i][j+1].paramv;
                }
            }
        }
      fprintf(sumf, "\n");
    }

  fprintf(sumf, "\n");
  /* Uncommented for Celine's uses */
  // fprintf(sumf, "\n\tmean\tmode\tvar\tper.05\tperc.5\tperc.95\tin\tme\tmod\tmed\terror_Mean\terror_Mode\terror_Med\tbias_Mean\tbias_Mode\tbias_Med\tNomig\tmig\tHo\tok?\t");
 
  fprintf(sumf, "\nparameters\tmean\tmode\t var\tper.05\tperc.5\tperc.95\t");

  for(i=0;i<NPARAM;i++)
    {
      if(i==0) 
        {
          if(type>=2) fprintf(sumf, "\nq1");
          else fprintf(sumf, "\nTheta1");
        }
      if(i==1) 
        {
          if(type>=2) fprintf(sumf, "\nq2");
          else if(pparam.cp.sumout==0) fprintf(sumf, "\nb=N2/N1");
          else fprintf(sumf, "\nTheta2");
        }
      if(i==2) fprintf(sumf, "\nTcoal=Tgen/4N1");
      
      if(i==3)
        {
          if(type>=2) fprintf(sumf, "\nTim=Tcoal*(mu*l)");
          else fprintf(sumf, "\nTgen");
        }
      if(i==4)
        {
          if(type>=2) fprintf(sumf, "\nqA");
          else if(pparam.cp.sumout==0) fprintf(sumf, "\na=Na/N1");
          else fprintf(sumf, "\nThetaA");
        }
      if(i==5)
        {
          if((type>=2)&&(pparam.cp.sumout==3)) fprintf(sumf, "\nm12/(mu*l)");
          else if(pparam.cp.sumout==3) fprintf(sumf, "\nm12");
          else fprintf(sumf, "\nM12");
        } 
      if(i==6) 
        {
          if((type>=2)&&(pparam.cp.sumout==3)) fprintf(sumf, "\nm21/(mu*l)");
          else if(pparam.cp.sumout==3) fprintf(sumf, "\nm21");
          else fprintf(sumf, "\nM21");
        }
      /* Celine changed 12/28/2009 */ // replace from "pparam.cp.maxval[i]>0"
      if((pparam.cp.maxval[i]!=pparam.cp.minval[i])||(pparam.cp.maxval[i]!=0))                            // Get point Estimates
        {
          mean[i]=(long double) psumpu[i]/psumu;
          long double temp=(long double) psump2u[i]/psumu; // Different tests on precisions of Variances
          long double tp=(long double)(mean[i]*mean[i]);
          if((double) tp==(double) temp) var[i]=(long double) pcount/(pcount-1)*((double)temp-(double)tp);
          else if((double)temp-(double)tp<1e-15) var[i]=0;
          else var[i]=((long double)pcount/(pcount-1)*(temp-tp));
          //--- Write mean mode var perc.05 median perc.95 ---//
          fprintf(sumf, "\t%Lg\t%lg\t%Lg\t%lg\t%lg\t%lg", mean[i]*mulx[i], mode[i]*mulx[i], var[i]*mulx[i], percvec[i][0]*mulx[i], percvec[i][1]*mulx[i], percvec[i][2]*mulx[i]);
          /**** Performance vs. true value in simulation ****/
          /* if(pparam.cp.trueval[i]>-1)                                                          // Case, information on true value for simulations
             {
             if((pparam.cp.trueval[i]>=percvec[i][0])&&(pparam.cp.trueval[i]<=percvec[i][2])) // True value within 90% percentile
             fprintf(sumf, "\t1");
             else fprintf(sumf, "\t0");
             //--- Estimates within [true/2, true*2] ---//
             double min=(double)pparam.cp.trueval[i]/2;
             double max=(double)pparam.cp.trueval[i]*2;
             if((max==min)||(max==0)) max+=.000000000000001;

             if((mean[i]>min)&&(mean[i]<=max))                   // Mean
             fprintf(sumf, "\t1");
             else if((mean[i]>=min)&&(mean[i]<max))
             fprintf(sumf, "\t1");
             else fprintf(sumf, "\t0");

             if((mode[i]>min)&&(mode[i]<=max))                   // Mode
             fprintf(sumf, "\t1");
             else if((mode[i]>=min)&&(mode[i]<max))
             fprintf(sumf, "\t1");
             else fprintf(sumf, "\t0");

             if((percvec[i][1]>min)&&(percvec[i][1]<=max))       // Median
             fprintf(sumf, "\t1"); 
             else if((percvec[i][1]>=min)&&(percvec[i][1]<max))
             fprintf(sumf, "\t1");
             else fprintf(sumf, "\t0");

             if((double)mean[i]==pparam.cp.trueval[i])           // Square error
             fprintf(sumf, "\t0");
             else
             fprintf(sumf, "\t%Lg", (mean[i]-(long double)pparam.cp.trueval[i]));                 // mean
             fprintf(sumf, "\t%Lg", (long double)(mode[i]-pparam.cp.trueval[i]));                 // mode
             fprintf(sumf, "\t%Lg", (long double)((double)(percvec[i][1]-pparam.cp.trueval[i]))); // median(abs difference, not squared error)


             if((double)mean[i]==pparam.cp.trueval[i])              // Biase: true/estimates
             fprintf(sumf, "\t1");
             else if((long double)pparam.cp.trueval[i]!=0)
             fprintf(sumf, "\t%Lg", ((long double)(mean[i]/(long double)pparam.cp.trueval[i])));      // mean est/true
             else fprintf(sumf, "\t%Lg", ((long double)(mean[i]-(long double)pparam.cp.trueval[i]))); // mean est-true(M=0)

             if((double)mode[i]==pparam.cp.trueval[i])             // Mode
             fprintf(sumf, "\t1");
             else if((long double)pparam.cp.trueval[i]!=0)
             fprintf(sumf, "\t%Lg", ((long double)(mode[i]/(long double)pparam.cp.trueval[i])));      // mode est/true
             else fprintf(sumf, "\t%Lg", ((long double)(mode[i]-(long double)pparam.cp.trueval[i]))); // mode est-true(M=0

             if((double)percvec[i][1]==pparam.cp.trueval[i])       // Median
             fprintf(sumf, "\t1");                                 // mean 1
             else if((long double)pparam.cp.trueval[i]!=0)
             fprintf(sumf, "\t%Lg", ((long double)(percvec[i][1]/(long double)pparam.cp.trueval[i])));      // median est/true
             else fprintf(sumf, "\t%Lg", ((long double)(percvec[i][1]-(long double)pparam.cp.trueval[i]))); // median est-true(M=0)

             //--- Test on migration ---//
             int ok=-1;
             if(i==5)
             {
             ok=(int)(((pparam.cp.trueval[i]==0) &&(nomig12>.05))||((pparam.cp.trueval[i]>0) &&(nomig12<=.05)));
             fprintf(sumf, "\t%Lg\t%Lg\t%d\t%d", nomig12 , (long double)1-nomig12, (int)(nomig12>.05) , ok); // p(M=0) and p(M>0)
             }
             else if(i==6)
             {
             ok=(int)(((pparam.cp.trueval[i]==0) &&(nomig21>.05))||((pparam.cp.trueval[i]>0) &&(nomig21<=.05)));
             fprintf(sumf, "\t%Lg\t%Lg\t%d\t%d", nomig21, (long double)1-nomig21, (int)(nomig21>.05), ok);   // p(M=0) and p(M>0)
             }
             else
             fprintf(sumf, " ");
             }
             else                   // test of migration case: pvalue for M=0 and M>0 in case when no true value available
             {
             if(i==5)
             fprintf(sumf, "\t%Lg\t%Lg\t%d", nomig12 , (long double)1-nomig12, (int)(nomig12>.05));
             if(i==6
             fprintf(sumf, "\t%Lg\t%Lg\t%d", nomig21, (long double)1-nomig21, (int)(nomig21>.05));
             }
          *//////
        }
    }
  for(i=0;i<NPARAM;i++)
    free(percvec[i]);
  free(percvec);
  free(mean);
  free(mode);
  free(maxu);
  free(mulx);
  free(var);
  free(sumui);
}// End writesum


/******************** Functions to transform parameters between IM and MIMAR types ********************/
/*--------------------------*
 * Function transdoublePtoO |
 *--------------------------*
 | Transform ratios of mutation rates into population mutation rates. For option
 | -q3, transform M into m for IM comparison. 
*/
void transdoublePtoO(double parray[NPARAM], struct params p)
{
  parray[1]*=parray[0];   // theta2=theta1*b
  parray[4]*=parray[0];   // thetaA=Theta1*a
  if(p.cp.sumout>2)       // Case option -q3
    {
      if(parray[0]>0)
        {
          parray[5]*=p.mp.mu/parray[0];  // m12=M12*mu/theta1
          parray[6]*=p.mp.mu/parray[0];  // m21=M21*mu/theta1
        }
      else parray[5]=parray[6]=0;
    }
}
/*--------------------------*
 * Function transgridPtoO   |
 *--------------------------*
 | In historams, transform ratios of mutation rates into population mutation rates. for option
 | -q3, transform M into m for IM comparison. 
*/
void transgridPtoO(struct grid parray[NPARAM], struct params p)
{
  parray[1].paramv*=parray[0].paramv;   // theta2=theta1*b
  parray[4].paramv*=parray[0].paramv;   // thetaA=Theta1*a
  if(p.cp.sumout>2)                     // Case option -q3
    {
      if(parray[0].paramv>0)
        {
          parray[5].paramv*=p.mp.mu/ parray[0].paramv; // m12=M12*mu/theta1
          parray[6].paramv*=p.mp.mu/ parray[0].paramv; // m21=M2*mu/theta1
        }
      else parray[5].paramv=parray[6].paramv=0;
    }
}


/******************** Functions for case with test of migration ********************/
/*-------------------------*
 * Function recordmcmctest |
 *-------------------------*
 | Records set of param in every interval after burnin in histogram and output, and summary file every L steps, in case of test of migration
*/
void recordmcmctest(struct params pparam, long double *pinfoperf, struct grid **pgridval, struct grid **pgridval0, struct grid **pgridvalm, struct grid *psumvect, long double *psumpu, long double *psump2u, long double *psumpu0, long double *psump2u0, long double *psumpum, long double *psump2um, int *pburn, double *pinterval, FILE *ppf)
{
  void recordparamtest(struct params, long double*, struct grid**, struct grid**, struct grid**, struct grid*, long double*, long double*, long double*, long double*, long double*, long double*, int, FILE*);
  if((!*pburn)&&(*pinterval==pparam.mcmcp.burnin)) // Case: end burnin, first recording
    {
      *pburn=1;
      *pinterval=0;                                // reset
      recordparamtest(pparam, pinfoperf, pgridval, pgridval0, pgridvalm, psumvect, psumpu, psump2u, psumpu0, psump2u0, psumpum, psump2um, 1, ppf);
    }
  else if((*pburn)&&((int)*pinterval==pparam.mcmcp.interrec)) // Record steps every i stepts (option -i)
    {
      *pinterval=0;                                // reset
      recordparamtest(pparam, pinfoperf, pgridval, pgridval0, pgridvalm, psumvect, psumpu, psump2u, psumpu0, psump2u0, psumpum, psump2um, 1, ppf);
    }
  else if(*pburn)                                  // Case in burnin, record in prior only
    recordparamtest(pparam, pinfoperf, pgridval, pgridval0, pgridvalm, psumvect, psumpu, psump2u, psumpu0, psump2u0, psumpum, psump2um, 0, ppf);
  else
    pinfoperf[10]+=1.;                             // 1+ step analyzed
}// End record mcmc with test of migration

/*--------------------------*
 * Function recordparamtest |
 *--------------------------*
 | Records set of param in grid for M=0 and M>0 every int (option -i), in case of test of migration
*/
void recordparamtest(struct params pparam, long double *pinfoperf, struct grid **pgridval, struct grid **pgridval0, struct grid **pgridvalm, struct grid *psumvect, long double *psumpu, long double *psump2u, long double *psumpu0, long double *psump2u0, long double *psumpum, long double *psump2um, int rec, FILE *ppf)
{
  int i=0, j=0;
  pinfoperf[10]+=1.;      // 1+ step analyzed
  if((psumvect[5].paramv==0)&&(psumvect[6].paramv==0))
    pinfoperf[17]+=1.;    // NoMig
  else
    pinfoperf[18]+=1.;    // Mig
  pinfoperf[0]+=psumvect[0].prob;                               // sum ui
  pinfoperf[1]+=(long double)psumvect[0].prob*psumvect[0].prob; // sum ui^2
  if(rec) fprintf(ppf, "%Lg\t", pinfoperf[10]);                 // Record every i steps
  for(i=0;i<NPARAM;i++)
    {
      psumvect[i].u=1;                                          // +1 param in this bin
      psumpu[i]+=psumvect[i].paramv;
      psump2u[i]+=psumvect[i].paramv*psumvect[i].paramv;
      if((psumvect[5].paramv==0)&&(psumvect[6].paramv==0))
        {
          psumpu0[i]+=psumvect[i].paramv;
          psump2u0[i]+=psumvect[i].paramv*psumvect[i].paramv;
        }
      else
        {
          psumpum[i]+=psumvect[i].paramv;
          psump2um[i]+=psumvect[i].paramv*psumvect[i].paramv;
        }
      if(rec) fprintf(ppf, "%lg\t", psumvect[i].paramv);
      if(pparam.cp.maxval[i]>0)
        {
          while((psumvect[i].paramv>pgridval[i][j].max))
            {
              j++;
              if(j>=MAX_LINE) break;
            }
          if(j<MAX_LINE)
            {
              pgridval[i][j].u+=psumvect[i].u;
              if((psumvect[5].paramv==0)&&(psumvect[6].paramv==0))
                pgridval0[i][j].u+=psumvect[i].u;        // p(M=0)
              else
                pgridvalm[i][j].u+=psumvect[i].u;        // p(M>0)
            }
          else                                           // Last bin
            {
              j--;
              pgridval[i][j].u+=psumvect[i].u;
              if((psumvect[5].paramv==0)&&(psumvect[6].paramv==0))
                pgridval0[i][j].u+=psumvect[i].u;
              else
                pgridvalm[i][j].u+=psumvect[i].u;
            }
        }
    }
  if(rec) fprintf(ppf, "%Lg\n", psumvect[0].prob);
}// End record parameters in histogram with test of migration

/*------------------------------*
 * Function write interval test |
 *------------------------------*
 | Writes histograms of posteriors every L steps (defines by option -L), in case
 | of test of migration.
*/
void writeintervaltest(struct params pparam, struct grid **pgridval0, struct grid **pgridvalm, long double *pinfoperf, int pargc, char *pargv[], long double *psumpu0, long double *psump2u0, long double *psumpum, long double *psump2um, long double min, double sec, int ta)
{
  int i=0;
  long double ess=0.0;
  char st[200]="";
  void writesum(FILE*, struct grid**, long double, long double, long double*, long double*, struct params, int);
  void writeperf(FILE*, long double*, long double, long double, double);
  long double getess(long double*);
  FILE *ismlf0=NULL, *ismlfm=NULL;

  ess=getess(pinfoperf);
  sprintf(st, "%s-%s-%d", pparam.cp.fout, "M0", ta); // Create file for M=0 case
  ismlf0=fopen(st, "w"); 
  sprintf(st, "%s-%s-%d", pparam.cp.fout, "M", ta);  // Crate file for M>0 case
  ismlfm=fopen(st, "w"); 

  for(i=0;i<pargc;i++)
    {
      fprintf(ismlf0, "%s ", pargv[i]);
      fprintf(ismlfm, "%s ", pargv[i]);
    }
  writeperf(ismlf0, pinfoperf, ess, min, sec);
  fprintf(ismlf0, "No Mig\t%Lg\tMig\t%Lg\n", pinfoperf[17], pinfoperf[18]);
  writesum(ismlf0, pgridval0, pinfoperf[17], pinfoperf[17], psumpu0, psump2u0, pparam, 1);
  fclose(ismlf0);
  writeperf(ismlfm, pinfoperf, ess, min, sec); 
  fprintf(ismlfm, "No Mig\t%Lg\tMig\t%Lg\n", pinfoperf[17], pinfoperf[18]);
  writesum(ismlfm, pgridvalm, pinfoperf[18], pinfoperf[18], psumpum, psump2um, pparam, 1);
  fclose(ismlfm);
}// End write interval test

/*------------------------*
 * Function writeendtest  |
 *------------------------*
 | Writes final histogram in summary file for M=0 and M>0
*/
void writeendtest(FILE *poutf, struct params pparam, struct grid **pgridval, long double *pinfoperf, int pargc, char *pargv[], long double *psumpu, long double *psump2u, long double min, double sec, int type, long double psumu)
{
  long double pess=0.0;
  long double getess(long double *);
  void writesum(FILE*, struct grid**, long double, long double, long double*, long double*, struct params, int);
  void writeperf(FILE*, long double*, long double, long double, double);

  pess=getess(pinfoperf);
  writeperf(poutf, pinfoperf, pess, min, sec); 
  if((type==1)||(type==3))
    writesum(poutf, pgridval, psumu, psumu, psumpu, psump2u, pparam, type);
  else
    writesum(poutf, pgridval, psumu, psumu, psumpu, psump2u, pparam, type);
}// End writeendtest function
