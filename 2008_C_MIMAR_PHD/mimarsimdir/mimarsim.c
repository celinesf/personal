/*******************************      mimarsim.c      ***********************************
 *
 *      Simulates gene genealogies (or ARGs) using the isolation-migration model.
 *   Computes summaries S statistics for each locus.
 *
 * Usage is shown by typing mimarsim without arguments.
 usage: mimarsim Y -lf -u -t -ej [options]

 Y is the number of loci considered.
 -lf followed by the input file name containing information on loci.
 -u gives the mutation rate per base pair, mu.
 -t gives the population mutation rate per bp, 4N1*mu.
 -ej set the time of split between the two populations.

 Other options: See mimardoc.pdf or after downloading and compiling, type mimargof<CR>.

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

 To compile: cc -o mimarsim mimarsim.c make_gametes.c params.c streec.c rand1.c -lm
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
 * More "user friendly" command line + error check that complex models work.
 * Get param from user with restriction: I added as many checks as I could to stop
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
 * Add summary statics as computed in mimargof
 *
 *
 ****************************** Changes in different verions

 *--- Changed 18th March 2010 ---*
  Add objects/functions to specify seeds in the command line.
  -----------------------------

  *--- Changed Jan 11 2010: ---*
  - Added make_gametes_noanc.c to assume that ancestral states are unknow and use allele 
    with the minor frequencies to obtain the summary staristics. 
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
  % mimarsim.c: 
  - changed calloc(1... by malloc(... .
  - correction on memort leaks.
  Thanks to the user who mentionned the bug to me.
  -----------------------------

  *--- Changed March 3rd 2009: ---*
  % params.c was changed: Corrected a bug on memory allocation in function getpars.
  Thanks Susan J. Miller for mentionning the bug to me.
  -----------------------------

  *--- Changed Nov. 5th 2008 ---*
  params.c was changed: usage for -eh and -es uncommented and changed.
  -----------------------------

  *--- Changed May 29th 2008: ---*
  % params.c was changed: I removed several check flags to allow greater freedom
  to the user. Thanks Armando Geraldes, for mentionning the problem to me.
  -----------------------------

  *--- Changed Oct. 9th  2007 ---*
  % mimarsim.c was changed: 
  - all the references to tajd.c have been commented.
  - all the parts related to the file reporting the statistics relevant to
  the goodness of fit test have been commented.
  % makefile was changed: All the references to tajd.c have been commented.
  -----------------------------

  *--- Changed Sept. 18th 2007 ---*
  % mimarsim.c was changed: Added lp.name.
  % mimarsim.h was changed: Added lp.name.
  % params.c was changed: In function getpars, the option -lf was changed so that the locus
  name can start by any character (before the name could not start by a number).
  ----------------------------- 

**********************************************************************************/

#include <stdio.h>                            // Input output Library
#include <stdlib.h>                           // File gestion Library
#include <math.h>                             // Math Library
#include <assert.h>                           // Verify program assertion
#include <string.h>
#include "mimarsim.h"

#define SITESINC 10                           // Constant min allocation for haplotypes
#define NPARAM 7                              // Number of parameters=theta1, theta2, tcoal, Tgen, thetaA, M12, M21
unsigned maxsites=SITESINC;

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

//--- Common Declarations ---//
double *posit;                                // Array of positions
int *typeseg;                                 // Type of seg sites

/*********************************************************
 *                                                       *
 * Main Function                                         *
 * -------------                                         *
 *                                                       *
 *********************************************************
 | Takes arguments and launches the simulations.      |
 *-------------------------------------------------------*/
int main(int argc, char *argv[])              // Array of char=arguments line
{
  //--- Declarations & Call Function ---//
  int i=0, j=0, k=0, count=0, howmany=0, segsites=0, okim=0;
  int **statseg=NULL, **nbvariant=NULL, S[5];
  FILE *pf=NULL, *outstat=NULL, *fopen(const char*, const char*); // Pointer on File for outputs (pf) and IM input file
  char **list=NULL, st[200]="";               // Haplotype list
  /*- IM variable -*/
  FILE *outIM=NULL; 
  int l=0;
  /* Celine changed 10/09/2007 - Uncomment for Celine's use */ /* for gof */
  /* double tajd();
     double bval=0, time=0; struct devent *pevent=NULL; *//////
  void updatemainparams(struct params*);
  int gensam(struct params*, char**, int**, int*);
  int **imatrix(int, int);
  //// From rand1.c ////
  /* Celine changed 03/18/2010 */
  void seedit(char*, FILE*, struct params *);/*/////*/
  char **cmatrix(int, int);
  //// From params.c ////
  void changeparams(struct params*);
  void changeparamslocus(struct params*, int);
  struct params getpars(int, char*[], int*);
  //--- Structure declaration---//
  struct params param;

  //--- Get arguments ---//
  param=getpars(argc, argv, &howmany);        // Get input by user for parameters
  //--- open output files ---//
  if(param.cp.fout[0]!='-')
    {
      /* IM file */
      outIM=fopen(param.cp.fout, "w");        // IM input file
      fprintf(outIM, "\nPopulation1\tPopulation2\n%d\n", howmany); // 2nd line of comment
      sprintf(st, "%s-s", param.cp.fout);
      outstat=fopen(st, "w");
    }
  pf=stdout;                                  // Mimar input
  for(i=0;i<argc;i++)                         // Information on simulation
    {
      fprintf(pf, "%s ", argv[i]);
      if(param.cp.fout[0]!='-')
        fprintf(outstat, "%s ", argv[i]);
    }
  /* Celine changed 03/18/2010 */
  if( !param.commandlineseedflag ) seedit("s", pf, &param);// WRITE seeds in summary output file 
  fprintf(pf,"\n%d %d %d\n",param.tableseeds[0],param.tableseeds[1],param.tableseeds[2]);
  /*/////*/

  //--- Initialisation & Memory allocation ---//
  posit=(double*)malloc((unsigned)(maxsites*sizeof(double)));
  nbvariant=imatrix(param.cp.npop+1, maxsites);           // array of nb of frequency spectrum
  typeseg=(int*)malloc((unsigned)(maxsites*sizeof(int))); // type of sites
  statseg=imatrix(param.cp.npop+3, howmany);              // Records locus specific S1, S2, Ss, Sf,
  changeparams(&param);                                   // Change estimates parameters from priors
  updatemainparams(&param);                               // Update parameters for the coalescent
  count=okim=0; // can write in IM if there is a seg site
  //--- Loop along the # sample in the simulation ---//
  while((howmany-count++))
    {
      if(okim==0)
        {
          for(i=0;i<maxsites;i++)                  // Init positins
            posit[i]=0.0;

          changeparamslocus(&param, count-1);      // Get locus specific parameters;
          list=cmatrix(param.cp.nsam, maxsites+1); // Allocate list of haplotypes
          if(count==1)
            {
              /* Uncomment for Celine's use */     // Begining info on parameters info
              /* if(param.cp.deventlist!=NULL)       // find time and a
                 {
                 pevent=param.cp.deventlist;
                 bval=1;
                 while(pevent!=NULL) // i==0)    // print Na
                 {
                 if(pevent->detype=='j')
                 time=pevent->time;        // T, t=T4N1

                 if(pevent->detype=='N')
                 bval=pevent->paramv;
                 pevent=pevent->nextde;
                 }
                 }
                 //- True values theta N2 T Na m12 m21 -//
                 fprintf(pf, "\nN1\tN2\tNa\ttyear\tM12\tM21\n"); // N1, N2, Na, Tyear 
                 fprintf(pf, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n", param.mp.thetafix/(4*param.mp.mu), (param.mp.thetafix*param.cp.size[1])/(4*param.mp.mu), (param.mp.thetafix*bval)/(4*param.mp.mu), time*param.mp.thetafix*param.mp.g/param.mp.mu, param.cp.mig_mat[0][1], param.cp.mig_mat[1][0]);

                 fprintf(pf, "q1\tq2\tqa\ttime\tm1\tm2\n");      // N1, N2, Na, Tyear 
                 fprintf(pf, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n", param.mp.thetafix*param.lp[0].li, param.mp.thetafix*param.lp[0].li*param.cp.size[1], param.mp.thetafix*param.lp[0].li*bval, time*param.mp.thetafix*param.lp[0].li , param.cp.mig_mat[0][1]/(param.mp.thetafix*param.lp[0].li), param.cp.mig_mat[1][0]/(param.mp.thetafix*param.lp[0].li)); // q1, q2, qa, tim, m1, m2
                 fprintf(pf, "mu\t%lg\tG\t%lg\n", param.mp.mu, param.mp.g);
                 fprintf(pf, "\nTheta\tN2\tT\tt\tNa\tMig12\tMig21");
                 // true values for this simulation in ms unit //
                 fprintf(pf, "\nms param Parameters\n%lg\t%lg\t", param.mp.thetafix, param.cp.size[1]); // mu G**, theta, N2
                 fprintf(pf, "%lg\t%lg\t", time, time*param.mp.thetafix/param.mp.mu);//  T, t=T4N1
                 fprintf(pf, "%lg\t%lg\t%lg\n", bval, param.cp.mig_mat[0][1], param.cp.mig_mat[1][0]);  // migration
                 fprintf(pf, "**Estimates\n");
              *////// //  End info on parameters info
              fprintf(pf, "Parameter values\n");   // Value for general user
              for(i=0;i<NPARAM;i++)
                fprintf(pf, "%lg\t", param.cp.newest[i]);
              fprintf(pf, "\n");
              fprintf(pf, "Name\tlenght\tx_y\tv_y\tw_y\tn_1\tn_2\tS_1\tS_2\tS_s\tS_f\t//\n");
              /* Uncomment for Celine's Use */ /* for gof */
              /* if(param.cp.fout[0]!='-')
                 {
                 fprintf(outstat, "\nParameters\n");
                 for(i=0;i<NPARAM;i++)
                 fprintf(outstat, "%lg\t", param.cp.newest[i]);
                 fprintf(outstat, "\n");
                 fprintf(outstat, "Name\tlenght\tx_y\tv_y\tw_y\tn_1\tn_2\tS_1\tS_2\tS_s\tS_f\tF_st\tPi_1\tPi_2\tD_1\tD_2\tp_1\tp_2//\n");
                 }*//////
            }
          for(i=0;i<5;i++) S[i]=0; 
          segsites=gensam(&param, list, nbvariant, S); // Allocate list of haplotypes
          statseg[0][count-1]=segsites;                // Total number of seg sites in sample
          for(i=1;i<3+param.cp.npop;i++)
            statseg[i][count-1]=0;
          if((segsites>0))                           // Case segsite < lenght of locus
            {
              /* Celine change 09/18/07 */
              fprintf(pf, "%s\t", param.lp[count-1].name); // Locus name
              /**/
              fprintf(pf, "%d\t%lg\t%lg\t%lg\t%d\t%d\t", param.lp[count-1].li, param.lp[count-1].xi, param.lp[count-1].vi, param.lp[count-1].wi, param.lp[count-1].ni[1], param.lp[count-1].ni[2]);

              if(segsites<=param.cp.nsites)                // Case segsite < lenght of locus
                {
                  if(param.cp.fout[0]!='-')                // Case IM input 
                    {
                      /*---------- Used to prepare IM into file  -----------*/
                      for(j=0;j<segsites;j++)              // Loop on seg sites to compute their positions
                        {
                          i=k=0; /* Celine Changed 03/17/2009 */ // init
                          i=(int)(posit[j]*param.cp.nsites);
                          k=(int)(posit[j+1]*param.cp.nsites);;
                        
                          if((i==k)&&(i==0))
                            posit[j+1]+=(double)1/param.cp.nsites;
                          else if(i==k)
                            {
                              l=(int)(posit[j-1]*param.cp.nsites);
                              if(i==k)
                                {
                                  if((i-1<=l)&&(k+1<param.cp.nsites))
                                    posit[j+1]+=((double)1/param.cp.nsites);
                                  else if((i-1>l))
                                    posit[j]-=((double)1/param.cp.nsites);
                                  else                     // Case Position > locus lenght
                                    {
                                      okim=1;
                                      printf("\nProblem :last segregating site conflict at locus# %d, inflile.txt aborted.\n", count);
                                      break;
                                    }
                                }
                            }
                          else if(k<i)
                            posit[j+1]+=((double)2/param.cp.nsites);
                        }// End loop segsites for position
                      
                      if(okim==0) // case position OK 
                        {
                          //--- write loci info in IM input file ---//
                          /* Celine changed 09/18/07 */
                          fprintf(outIM, "%s %d %d %d I %lg %lg\n", param.lp[count-1].name, param.lp[count-1].ni[1], param.lp[count-1].ni[2], param.lp[count-1].li, param.lp[count-1].xi, param.cp.nsites*param.mp.mu/param.mp.g);
                          /*/////*/
                          j=0;                             // position of segsite
                          for(k=0;k<param.cp.nsam;k++)     // Loop alog haplotypes
                            {
                              if(k<param.lp[count-1].ni[1]) fprintf(outIM, "pop1_%d" , k+1);
                              else fprintf(outIM, "pop2_%d" , k+1);
 
                              if(k<9) l=6; //--- put only 10 char for sequence name: IM input sequence starts at #11---//
                              else if(k<99) l=7;
                              else l=8;
                              while(l<10)                    // Loop to fill sequence name
                                {
                                  fprintf(outIM, " ");
                                  l++;
                                }
                              j=0;
                              for(i=0;i<param.cp.nsites;i++) // Loop fill sites
                                {
                                  if((int)(posit[j]*param.cp.nsites)==i)
                                    {
                                      if(list[k][j]=='0')
                                        fprintf(outIM, "a");
                                      else
                                        fprintf(outIM, "g");
                                      j++;
                                    }
                                  else
                                    fprintf(outIM, "c");
                                }
                              fprintf(outIM, "\n");
                              if(j<segsites)
                                {
                                  if(param.cp.uniform[0][5]!=5)
                                    {
                                      printf("\nProblem at locus# %d, only %d sites writen while %d was expected. inflile.txt aborted.\n", count, j, segsites);
                                      okim=1;
                                      break;
                                    }
                                }
                            }// End loop haplotypes
                          
                          if(okim==0)                      // Case Im input filled properly
                            {
                              /*fprintf(pf, "\n");
                                if(segsites>0) 
                                for(i=0;i<param.cp.nsam;i++) fprintf(pf, "%s\n", list[i]); // Write haplotypes for the sample
                              */
                              for(i=0;i<segsites;i++)                     // Calulate S statistics for the locus
                                {
                                  if(typeseg[i]<0) statseg[3][count-1]++; // shared
                                  else if(typeseg[i]<param.cp.npop+1) statseg[typeseg[i]][count-1]++; // population specific
                                  else statseg[4][count-1]++;             // fixed
                                }
                              for(i=1;i<5;i++)                            //--- Output S1 S2 Ss Sf ---//
                                fprintf(pf, "%d\t", statseg[i][count-1]);
                              fprintf(pf, "\n");
                              /* Celine changed 10/09/2007 - Uncomment for Celine's use */ /* for gof */
                              /* for(k=0;k<param.cp.nsam;k++)// loop on all sequences in the locus
                                 {
                                 for(i=k+1;i<param.cp.nsam;i++)
                                 {
                                 if((k<param.lp[count-1].ni[1])&&(i<param.lp[count-1].ni[1])) // pop1 
                                 {
                                 param.lp[count-1].H[0]++;                                // # chromosomes
                                 for(j=0;j<segsites;j++)
                                 {
                                 if(list[k][j]!=list[i][j])                                           
                                 param.lp[count-1].H[1]++;                          // # seg sites
                                 }
                                 }
                                 else if((k>=param.lp[count-1].ni[1])&&(i>=param.lp[count-1].ni[1])) // pop2
                                 {
                                 param.lp[count-1].H[4]++;                                       // # chromosomes
                                 for(j=0;j<segsites;j++)
                                 {
                                 if(list[k][j]!=list[i][j])
                                 param.lp[count-1].H[5]++;                                 // # seg sites
                                 }
                                 }
                                 else                                // total sample
                                 {
                                 param.lp[count-1].H[8]++;       // Totsal sample size
                                 for(j=0;j<segsites;j++)
                                 {
                                 if(list[k][j]!=list[i][j])
                                 param.lp[count-1].H[9]++; // total S - hb
                                 }
                                 }
                                 }// Loop on chromosome
                                 }// Loop along all sampled sequence for the locus
                                 param.lp[count-1].H[1]/=param.lp[count-1].H[0];     // Hw1
                                 param.lp[count-1].H[5]/=param.lp[count-1].H[4];     // Hw2
                                 param.lp[count-1].H[9]/=param.lp[count-1].H[8];     // Hb
                                 param.lp[count-1].H[10]=1-((param.lp[count-1].H[1]+param.lp[count-1].H[5])/2)/ param.lp[count-1].H[9]; // loci specific Fst
                             
                                 if(S[0]>0)                   // Locus popymorphic in pop1
                                 param.lp[count-1].H[2]=tajd(param.lp[count-1].ni[1], S[0], param.lp[count-1].H[1]);
                                 if(S[1]>0)                   // Locus popymorphic in pop2
                                 param.lp[count-1].H[6]=tajd(param.lp[count-1].ni[2], S[1], param.lp[count-1].H[5]);
                             
                                 if(statseg[1][count-1]>0)    // Locus popymorphic private in pop1
                                 param.lp[count-1].H[3]+=(double)  S[2]/(statseg[1][count-1]*param.lp[count-1].ni[1]*2); // p(1) for pop1
                                 if(statseg[2][count-1]>0)    // Locus popymorphic private in pop2
                                 param.lp[count-1].H[3]+=(double) S[3]/(statseg[2][count-1]*param.lp[count-1].ni[2]*2);  // p(1) for pop2
                                 if(statseg[3][count-1]>0)    // Locus popymorphic private in pop2
                                 param.lp[count-1].H[7]+=(double) S[4]/(statseg[3][count-1]*(param.lp[count-1].ni[2]+param.lp[count-1].ni[1])); // p(2)
                              *//////
                              for(i=0;i<param.cp.nsam;i++) // Free memory for this locus
                                free(list[i]);
                              free(list);
                            }// En IM input ok
                        }// End positions <lenght
                    }// End Case in IM inpout file
                  else if(okim==0)                     // Case only MIMAR input file
                    {
                      for(i=0;i<segsites;i++)                     // Calulate S statistics for the locus
                        {
                          if(typeseg[i]<0) statseg[3][count-1]++; // shared
                          else if(typeseg[i]<param.cp.npop+1) statseg[typeseg[i]][count-1]++; //population specific
                          else statseg[4][count-1]++;             // fixed
                        }
                      for(i=1;i<5;i++)                            //--- Output S1 S2 Ss Sf ---//
                        fprintf(pf, "%d\t", statseg[i][count-1]);
                      fprintf(pf, "\n");
                      
                      for(i=0;i<param.cp.nsam;i++)                // Free memory for this locus
                        free(list[i]);
                      free(list);
                    }// Case MIMAR input file
                  /* For the following arrays: locus specific values
                     0 nsam1, 1 Hw1, 2 TD1, 3 p1, 4 nsam2, 5 Hw2, 6 TD2, 7p2, 8 nsam total, 9 Hb, 10 fst, for locus
                  */
                  /* Celine changed 10/09/2007 - Uncomment for Celine's use */ /* for gof */
                  /* if(param.cp.fout[0]!='-')
                     {
                     // Celine changed 09/18/07 //
                     // fprintf(outstat, "%s\t", param.lp[count-1].name); // Locus name
                     //////
                     fprintf(outstat, "%d\t%lg\t%lg\t%lg\t%d\t%d\t", param.lp[count-1].li, param.lp[count-1].xi, param.lp[count-1].vi, param.lp[count-1].wi, param.lp[count-1].ni[1], param.lp[count-1].ni[2]);
                     
                     for(i=1;i<5;i++)                                     //--- Output S1 S2 Ss Sf ---//
                     fprintf(outstat, "%d\t", statseg[i][count-1]);
                     
                     fprintf(outstat, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", param.lp[count-1].H[10], param.lp[count-1].H[1], param.lp[count-1].H[5], param.lp[count-1].H[2], param.lp[count-1].H[6], param.lp[count-1].H[3], param.lp[count-1].H[7]);
                     
                     }*//////
                }// End case segsite< locus lenght
              else
                {
                  okim=1;
                  printf("\nCan Not Right infile.txt cause to many seg sites at locus #%d (S=%d, length=%d bp)\n", count, segsites, param.cp.nsites);
                  break;
                }
            }// End locus polymorphic
          else
            {
              okim=1;
              printf("\nLocus #%d had no seg site- Simulation avorted.\n", count);
              break;
            }
        }// End Sample ok until now
    }// End Loop on loci
  
  if(param.cp.fout[0]!='-') // free IM input file
    {
      fclose(outstat);
      fclose(outIM);
    }
  /* Celine changed 03/18/2010 */
  seedit("end", pf, &param);                 // in randx.c, flag[0]!="s" so create/rewrite seed in seedmimar
  /*/////*/
  fclose(pf);
 
  free(posit);
  free(typeseg);
  for(i=0;i<param.cp.npop+3;i++)
    {
      if(i<param.cp.npop+1)
        free(nbvariant[i]);
      free(statseg[i]);
    }
  free(nbvariant);
  free(statseg);
    
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
  /* Celine changed 11/27/2009 */
  for(i=0;i<howmany;i++)
    free(param.lp[i].name);
  /*/////*/
  free(param.lp);
  /* Celine changed 03/18/2010 */
  free( param.tableseeds); /*/////*/
  exit(0);
}// End main function

/*---------------------------*
 * Function updatemainparam  |
 *---------------------------*
 | Updates the main parameters required for coalescent.
*/
void updatemainparams(struct params *pp)
{
  int j=0;

  for(j=0;j<NPARAM;j++)
    pp->cp.newparam[j]=0;

  pp->mp.thetafix=pp->cp.newest[0];   //--- Theta 1 ---//
  pp->cp.newparam[0]=pp->cp.newest[0];

  if(pp->cp.sumout>=1)                //--- Theta2 and ThetaA ---//
    {
      pp->cp.newparam[1]=pp->cp.newest[1]/pp->mp.thetafix;          // b=theta2/theta
      pp->cp.newparam[4]=(double) pp->cp.newest[4]/pp->mp.thetafix; // a=thetaA/theta1
    }
  else
    {
      pp->cp.newparam[1]=pp->cp.newest[1]; // b estimated
      pp->cp.newparam[4]=pp->cp.newest[4]; // a estimated
    }
  //--- T ---//
  pp->cp.newparam[2]=pp->cp.newest[2];  // Tcoalescent
  pp->cp.newparam[3]=pp->cp.newest[3];  // Tgen
  //--  Migration --//
  pp->cp.newparam[5]=pp->cp.newest[5]; // M12
  pp->cp.newparam[6]=pp->cp.newest[6]; // M21
}// End change Coalescent param

/*-----------------*
 * Function getsam |
 *-----------------*
 | Generates a gene genealogy for a locus.
*/
int gensam(struct params *pparam, char **list, int **pnbvariant, int *pS) 
{
  int nsegs=0, k=0, seg=0, ns=0, start=0, end=0, len=0, segsit=0, nsites=0, nsam=0;
  /* The function returns ns, the number of sgregating sites.
   * nsegs is the gametes were broken into in tracing back the history of the gametes.
   * The histories of these segments are  passed back to the calling function in the array of structures seglst[]
   */
  struct segl *seglst=NULL;
  double nsinv=0.0, tseg=0.0, tt=0.0, theta=0.0;
  void make_gametes(int, struct node*, double, int, int, char**, int, int*, int**, int*, int*);

  void inititable(int, int, int*);
  void initimatrix(int, int, int, int**);
  void biggerimatrix(int, int**);
  void biggerlist(int, char**);
  void locate(int, double, double, double*);
  //// From make_gametes ////
  void prtree(struct node*, int);
  double ttime(struct node*, int);
  int poisso(double);
  //// From streec.c ////
  struct segl*segtre_mig(struct c_params*, int*);

  nsites=pparam->cp.nsites;                 // Locus lenght available for recombination
  nsinv=1./nsites;
  seglst=segtre_mig(&(pparam->cp), &nsegs); // Generate Xove and ARGr, record # segments
  nsam=pparam->cp.nsam;                     // # chromo in the sample
  theta=pparam->mp.theta;                   // mutation rate given by user
  ns=0;                                     // # seg sites
  if(pparam->mp.treeflag)                   // Case output tree.
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
              fprintf(stdout, "[%d]", len);
            }
          prtree(seglst[seg].ptree, nsam);
          if((theta==0.0)) 
            free(seglst[seg].ptree);
        }
    }
  //--- Loop along all segments ---//
  for(seg=0, k=0;k<nsegs; seg=seglst[seg].next, k++)
    {
      end=(k<nsegs-1 ? seglst[seglst[seg].next].beg-1:nsites-1); // End of segment
      start=seglst[seg].beg;                                     // beginning of segment
      len=end-start+1;                                           // Lengh of segment
      tseg=len*(theta/nsites);                                   // Mutation rate along the segment
      tt=ttime(seglst[seg].ptree, nsam);                         // Total time in the tree for this segment
      segsit=poisso(tseg*tt);                                    // get # segsites along genealogy for the segment

      //--- Realloc memory if # segsite bigger than max previously define ---//
      if((segsit+ns)>=maxsites)
        {
          maxsites=segsit+ns+SITESINC;                          // Bigger Max sites 
          posit=(double*)realloc(posit, maxsites*sizeof(double)); // Increase # position for the sample
          biggerlist(nsam, list);                               // Increase list of haplotype
          typeseg=(int*)realloc(typeseg, maxsites*sizeof(int)); // Increase # seg sites in typeseg
          biggerimatrix(pparam->cp.npop+1, pnbvariant);         // Increase matrix of variant
        }

      //--- Initialize table of variant and segType ---//
      initimatrix(ns, ns+segsit, pparam->cp.npop+1, pnbvariant);
      inititable(ns, segsit+ns, typeseg);
      make_gametes(nsam, seglst[seg].ptree, tt, segsit, ns, list, pparam->cp.npop, pparam->cp.config, pnbvariant, typeseg, pS); // Make gametes (Put segsites on gene genealogy)
      free(seglst[seg].ptree); 
      locate(segsit,start*nsinv, len*nsinv,posit+ns);            // Fill Position of seg sites for this segment (recal nsinv=1/nsites)                                  // Free memory
      ns+=segsit;                                                // Total # seg sites
    }
  for(k=0;k<nsam;k++) list[k][ns]='\0';                          // End of haplotype strings
  return(ns);                                                    // Total # segsites
}// End Gensam

/*---------------------------*
 * Function Initiate table   |
 *---------------------------*/
void inititable(int cstart, int ncol, int ptab[])
{
  int i=0;
  for(i=cstart;i<ncol;i++)
    ptab[i]=0;
}

/******************** Functions for matrix of int ********************/
/*-------------------*
 * Function imatrix  |
 *-------------------*
 | Initial allocation of memory to the matrix of variant (int)
*/
int **imatrix(int nline, int ncol)
{
  int i=0;
  int **m=NULL;

  if(!(m=(int**)malloc((unsigned) nline*sizeof(int*))))
    perror("alloc error in imatrix");
  for(i=0;i<nline;i++) 
    {
      if(!(m[i]=(int*)malloc((unsigned) ncol*sizeof(int))))
        perror("alloc error in imatric. 2");
    }
  return(m);
}
/*---------------------------*
 * Function Initiate matrix  |
 *---------------------------*/
void initimatrix(int cstart, int ncol, int nline, int **list)
{
  int i=0, j=0;
  for(i=0;i<nline;i++)
    {
      for(j=cstart;j<ncol;j++)
        list[i][j]=0;
    }
}
/*-------------------------*
 * Function bigger matrix  |
 *-------------------------*
 | Allocates memory to the matrix of variant (int) when there's
 | more segsites than previously allocated.
*/
void biggerimatrix(int ncol, int **list)
{
  int i=0;
  for(i=0;i<ncol;i++)
    {
      list[i]=(int*)realloc(list[i], maxsites*sizeof(int));
      if(list[i]==NULL) perror("realloc error. biggerimatrix");
    }
}

/******************** Functions for matrix of char ********************/
/*-------------------*
 * Function cmatrix  |
 *-------------------*
 | Initial memory allocation for the list of haplotype (strings)
*/
char **cmatrix(int nsam, int len)
{
  int i=0;
  char **m=NULL;

  if(!(m=(char**)malloc((unsigned)nsam*sizeof(char*))))
    perror("alloc error in cmatrix");
  for(i=0;i<nsam;i++) 
    {
      if(!(m[i]=(char*)malloc((unsigned) len*sizeof(char))))
        perror("alloc error in cmatric. 2");
    }
  return(m);
}

/*-----------------------*
 * Function bigger list  |
 *-----------------------*
 | Allocates memory to the list of haplotype (strings) when there's
 | more segsites than previously allocated.
*/
void biggerlist(int nsam, char **list)
{
  int i=0;
  for(i=0;i<nsam;i++)
    {
      list[i]=(char*)realloc(list[i], maxsites*sizeof(char));
      if(list[i]==NULL) perror("realloc error. bigger");
    }
}

/*------------------*
 * Function locate  |
 *------------------*
 | Fills the position of the segregating site of a segment.
*/
void locate(int n, double beg, double len, double *ptr) // starting point, lenght between segs sites, pointer on position
{
  int i=0;
  void ordran(int, double*);
  ordran(n, ptr);                                       // n(# seg sites) ordered random values
  for(i=0;i<n;i++)
    ptr[i]=beg+ptr[i]*len;                              // position is a random factor of the total lenght from starting point of the segment
}

