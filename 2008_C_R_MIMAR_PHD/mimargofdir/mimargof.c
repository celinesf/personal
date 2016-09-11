/*******************************      mimargof.c      ***********************************
 *
 *      Simulates gene genealogies using the isolation-migration model and estimated
 *   parameters from MIMAR. 
 *   Computes summaries of the data: sum of S statistics over loci, mean Fst, 
 *   Pi (Hw for both pops), TD (for both pops), p(1) and p(2), the average
 *   frequency of mutant allele in site private to 1 and 2 populations, accross loci.
 *
 * Usage is shown by typing mimargof without arguments.
 usage: mimargof nsim Y -lf -u -t -ej [options]

 nsim is the number of simulation per set of parameters to generate.
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

 To compile: cc -o mimargof mimargof.c make_gametes.c params.c streec.c rand1.c tajd.c -lm
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
 * Add summary statics p(1) and p(2).
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
  
  *--- Changed Dec 28th 2009: ---*
  % param.c was changed: 
  - the locus specific rec rate is now rho_y=rho*w_y*(L-1)
  instead of rho_y=rho*w_y*(L).
  - changed calloc(1... by malloc(... .
  - added recombination options.
  - commented useless lines.
  -----------------------------

  *--- Changed Nov. 27th 2009: ---*
  % streec.c was changed: Corrected bugs on memory allocation.
  % mimargof.c: 
  - changed calloc(1... by malloc(... .
  - correction on memort leaks.
  Thanks to the user who mentionned the bug to me.
  -----------------------------

  *--- Changed March 3rd 2009: ---*
  % params.c was changed: Corrected a bug on memory allocation in function getpars.
  Thanks Susan J. Miller for mentionning the bug to me.
  -----------------------------

  *--- Changed May 29th 2008: ---*
  % params.c was changed: I removed several check flags to allow greater freedom
  to the user. Thanks Armando Geraldes, for mentionning the problem to me.
  -----------------------------

  *--- Changed Sept. 18th 2007 ---*
  % mimargof.h was changed: Added lp.name.
  % params.c was changed: In function getpars, the option -lf was changed so that the locus
  name can start by any character (before the name could not start by a number).
  -----------------------------

**********************************************************************************/

#include <stdio.h>                            // Input output Library
#include <stdlib.h>                           // File gestion Library
#include <math.h>                             // Math Library
#include <assert.h>                           // Verify program assertion
#include <string.h>
#include "mimargof.h"

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
int *typeseg;                                 // Type of seg sites

/*********************************************************
 *                                                       *
 * Main Function                                         *
 * -------------                                         *
 *                                                       *
 *********************************************************
 | Takes arguments and launches the gof simulations.      |
 *-------------------------------------------------------*/
int main(int argc, char *argv[])                    // Array of char=arguments line
{
  //--- Declarations & Call Function ---//
  int i=0, j=0, k=0, count=0, howmany=0, segsites=0, okim=0, oksim=0, numsim=0, totsim=0, nokl=0, nokl0=0;
  int **statseg=NULL, **nbvariant=NULL, oks[3], okstot[3];
  FILE *pf=NULL, *fopen(const char*, const char*); // Pointer on File for outputs (pf) and IM input file
  double tajd();
  char **list=NULL;                                // Haplotype list
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
  param=getpars(argc, argv, &howmany);  // Get input by user for parameters

  pf=stdout;                            // Output
  /* Celine changed 03/18/2010 */
  if( !param.commandlineseedflag ) seedit("s", pf, &param);// WRITE seeds in summary output file 
  /*/////*/
  /* Uncommented for Celine's use */
  /* for(i=0;i<argc;i++)                // Information on simulation
     fprintf(pf, "%s ", argv[i]);
  */////\
  //---------- Initialisation & Memory allocation ------------------//
  nbvariant=imatrix(param.cp.npop+1, maxsites);            // array of nb of frequency spectrum
  typeseg=(int*)malloc((unsigned)(maxsites*sizeof(int)));  // type of sites
  statseg=imatrix(param.cp.npop+3, howmany);               // Records locus specific S1, S2, Ss, Sf, 
  changeparams(&param);                                    // Change estimates parameters from priors
  updatemainparams(&param);                                // Update parameters for the coalescent

  oksim=totsim=okstot[0]=okstot[1]=okstot[1]=okstot[2]=0; // simulation check, total #of sim, check on statistics for all simulations
  for(numsim=0; numsim<param.cp.nsim;numsim++)            // Loop along number of simulations for this set of parameters
    {
      //--- Initialization and reset of quality checks ---//
      count=nokl=nokl0=okim=oks[0]=oks[1]=oks[1]=oks[2]=0; // number of loci, # loci with ok genealogies, no set sites, #statistics ok
      for(i=0;i<11;i++)                 // Sim specific Stats
        param.cp.sSiFst[i]=0;

      //--- Loop along the loci in the simulation ---//
      while((howmany-count++))
        {
          if(okim==0)                   // Case All loci ok in the sample
            {
              for(i=0;i<11;i++)
                {
                  param.cp.lSiFst[i]=0.0;
                  if(i<9)
                    param.lp[count-1].tpH[i]=0;
                }
              changeparamslocus(&param, count-1);      // Get locus specific parameters
               
              list=cmatrix(param.cp.nsam, maxsites+1); // Allocate list of haplotypes
              segsites=gensam(&param, list, nbvariant, param.lp[count-1].S);// Generate a new gene ARG
              statseg[0][count-1]=segsites;            // Total number of seg sites in sample
              for(i=1;i<3+param.cp.npop;i++)
                statseg[i][count-1]=0;

              if((segsites>0))                         // Case segsite>0: get stats
                {
                  /*   fprintf(pf, "segsites:%d\npositions:\n",segsites); */
                  /*                      for(i=0;i<param.cp.nsam;i++) fprintf(pf, "%s\n", list[i]);  */
                  /*                      fprintf(pf, "\n"); */
                  if(segsites<=param.cp.nsites)        // Case segsite < lenght of locus
                    {
                      for(k=0;k<param.cp.nsam;k++)
                        {
                          for(i=k+1;i<param.cp.nsam;i++)
                            {
                              if((k<param.lp[count-1].ni[1])&&(i<param.lp[count-1].ni[1])) // pop1
                                {
                                  param.lp[count-1].tpH[0]++;                              // # chromosomes
                                  for(j=0;j<segsites;j++)
                                    {
                                      if(list[k][j]!=list[i][j])
                                        param.lp[count-1].tpH[1]++;                         // # seg sites
                                    }
                                }
                              else if((k>=param.lp[count-1].ni[1])&&(i>=param.lp[count-1].ni[1])) // pop2
                                {
                                  param.lp[count-1].tpH[2]++;                                     // # chromosomes
                                  for(j=0;j<segsites;j++)
                                    {
                                      if(list[k][j]!=list[i][j])                                          
                                        param.lp[count-1].tpH[3]++;                              // # seg sites
                                    }
                                }
                              else                                                               // total sample
                                {
                                  param.lp[count-1].tpH[4]++;;                                   // Totsal sample size
                                  for(j=0;j<segsites;j++)
                                    {
                                      if(list[k][j]!=list[i][j])
                                        param.lp[count-1].tpH[5]++;                              // total S
                                    }
                                }
                            }// Loop on chromosome
                        }// Loop along all sampled sequence for the locus

                      for(i=0;i<segsites;i++)                     // Calulate S statistics for the locus
                        {
                          if(typeseg[i]<0) statseg[3][count-1]++; // shared
                          else if(typeseg[i]<param.cp.npop+1) statseg[typeseg[i]][count-1]++; // population specific
                          else statseg[4][count-1]++;             // fixed
                        }
                      for(i=1;i<5;i++)                            //--- Record S1 S2 Ss Sf forthe locus ---//
                        param.cp.lSiFst[i-1]+=statseg[i][count-1];                             
                      
                      for(i=0;i<param.cp.nsam;i++)                // Free memory for this locus
                        free(list[i]);
                      free(list);
                    }// End case segsite<lenght of locus
                  else                                 // Case segsites>lenght of locus
                    {
                      okim=1;                          // Sample have a wrong locus
                      oksim=1;                         // stop this simulation
                      break;
                    }
                }// End locus polymorphic
              else                                     // Locus without seg sites
                {
                  okim=1;                              // Sample have a wrong locus (S=0)
                  oksim=2;                             // 0 for all stats
                }
            }// End Sample good until now

          if(okim==0)                                  // All loci good until now
            {
              nokl++;                                  // +1 good locus
              nokl0++;                                 // +1 polymorphic locus
              for(i=0;i<7;i++)
                {
                  if(i<4)
                    param.cp.sSiFst[i]+=param.cp.lSiFst[i];                      // sum of Sk
                  param.lp[count-1].H[i]=param.lp[count-1].tpH[i];               // locus specific stats  
                }
              param.cp.lSiFst[5]=param.lp[count-1].H[1]/=param.lp[count-1].H[0]; // Hw1
              param.cp.lSiFst[6]=param.lp[count-1].H[3]/=param.lp[count-1].H[2]; // Hw2
              param.lp[count-1].H[5]/=param.lp[count-1].H[4];                    // Hb
              param.cp.lSiFst[4]=param.lp[count-1].H[6]=1-((param.lp[count-1].H[1]+param.lp[count-1].H[3])/2)/ param.lp[count-1].H[5];// locis specific Fst
              if((param.lp[count-1].S[0]>0)&&(param.lp[count-1].ni[1]>2)) // Locus popymorphic in pop1
                {
                  param.cp.lSiFst[7]=tajd(param.lp[count-1].ni[1], param.lp[count-1].S[0], param.lp[count-1].H[1]);
                  param.cp.sSiFst[7]+=param.cp.lSiFst[7];
                  oks[0]++;             // +1 good stat for pop1
                }
              if((param.lp[count-1].S[1]>0)&&(param.lp[count-1].ni[2]>2)) // Locus popymorphic in pop2
                {
                  param.cp.lSiFst[8]=tajd(param.lp[count-1].ni[2], param.lp[count-1].S[1], param.lp[count-1].H[3]);
                  param.cp.sSiFst[8]+=param.cp.lSiFst[8];
                  oks[1]++;             // +1 good stat for pop2
                }
              if(statseg[1][count-1]>0)                                  // Locus popymorphic private in pop1
                param.lp[count-1].H[7]=param.cp.lSiFst[9]+=(double) param.lp[count-1].S[2]/(statseg[1][count-1]*param.lp[count-1].ni[1]*2); // p(1)

              if(statseg[2][count-1]>0)                                  // Locus popymorphic private in pop2
                param.lp[count-1].H[7]=param.cp.lSiFst[9]+=(double) param.lp[count-1].S[3]/(statseg[2][count-1]*param.lp[count-1].ni[2]*2); // p(1)

              param.cp.sSiFst[9]+=param.cp.lSiFst[9];                    // sum p1
              if(statseg[3][count-1]>0)                                  // Locus popymorphic private in pop2
                {
                  param.lp[count-1].H[8]=param.cp.lSiFst[10]=(double) param.lp[count-1].S[4]/(statseg[3][count-1]*(param.lp[count-1].ni[2]+param.lp[count-1].ni[1])); // p(2)
                  param.cp.sSiFst[10]+=param.cp.lSiFst[10];              // sum p2
                  oks[2]++;
                }

              param.cp.sSiFst[4]+=param.lp[count-1].H[6]; // sum Fst
              param.cp.sSiFst[5]+=param.lp[count-1].H[1]; // sum Hw1
              param.cp.sSiFst[6]+=param.lp[count-1].H[3]; // sum Hw2
            }
          else if(oksim==2)             // Case no seg site for that locus
            {
              oksim=0;                  // reset checks
              okim=0;
              nokl++;                   // 1+ locus to count in mean (all 0 values)
              for(i=0;i<9;i++)
                param.lp[count-1].H[i]=param.lp[count-1].tpH[i]; // locus specific stats
            }
        }// End loop on loci
    
      if(nokl==howmany)                 // All sample good
        {
          totsim++;                     // 1+ good simulation
          for(i=0;i<4;i++)
            param.cp.SiFst[i]+=param.cp.sSiFst[i];                    // sum of S stats along simulations
          for(i=4;i<11;i++)
            {
              if(((i<7)||(i>=9))&&(nokl0>0))
                param.cp.SiFst[i]+=(double)param.cp.sSiFst[i]/nokl0; // mean of other stats along simulations
              if((i>=7)&&(i<9)&&(oks[i-7]>0))
                {
                  param.cp.SiFst[i]+=(double)param.cp.sSiFst[i]/oks[i-7];
                  okstot[i-7]++;
                }
            }
        }
    }// End loop on simulations
  if(oksim==0)                          // All simulations worked
    {
      /* Uncommented for Celine's use */
      /* for(i=0;i<11;i++) */////
      for(i=0;i<9;i++)
        {
          if((i<7)||(i>=9))
            fprintf(pf, "%lg\t", (double) param.cp.SiFst[i]/(totsim));           // write mean of sum of S stats, Fst and Hws over simulations
          
          if((i>=7)&&(i<9))
            {
              if(oks[i-7]>0)
                fprintf(pf, "%lg\t", (double) param.cp.SiFst[i]/(okstot[i-7])); // write mean Tds if S>0 in pops 
              else fprintf(pf, "NA\t" ); 
            }
        }
      fprintf(pf, "\n");
    } 
  else                                  // Case one locus with too much seg sites
    {
      for(i=0;i<9;i++)
        fprintf(pf, "NA\t" );
      fprintf(pf, "\n");
    }
  /* Celine changed 03/18/2010 */
  seedit("end", pf, &param);                 // in randx.c, flag[0]!="s" so create/rewrite seed in seedmimar
  /*/////*/
  fclose(pf);
  
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
          biggerlist(nsam, list);                               // Increase list of haplotype
          typeseg=(int*)realloc(typeseg, maxsites*sizeof(int)); // Increase # seg sites in typeseg
          biggerimatrix(pparam->cp.npop+1, pnbvariant);         // Increase matrix of variant
        }

      //--- Initialize table of variant and segType ---//
      initimatrix(ns, ns+segsit, pparam->cp.npop+1, pnbvariant);
      inititable(ns, segsit+ns, typeseg);
      make_gametes(nsam, seglst[seg].ptree, tt, segsit, ns, list, pparam->cp.npop, pparam->cp.config, pnbvariant, typeseg, pS); // Make gametes (Put segsites on gene genealogy)
      free(seglst[seg].ptree);                                   // Free memory
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


