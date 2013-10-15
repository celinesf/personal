/************ rand1.c - mimar *******************************************
 *  Link in this file for random number generation using drand48().
 *  Uses seedms file to get seed information.
 *  Functions return non-negative, double-precision, floating-point
 *  values, uniformly distributed over the interval [0.0 , 1.0]. 
 *
 *  Addition of standart Normal distribution generator used for MCMC.
 *  Function called in param.c
 *
 ************************************** Changes
 
 *--- Changed 18th March 2010 ---*
  Add objects/functions to specify seeds in the command line.
  -----------------------------
  *
  *****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mimargof.h" /* Celine changed 03/18/2010 */

//--- Generate a random number ---//
double ran1()
{
  double drand48();                               // Function pseudo ramdom generator
  return(drand48());
}

/*-------------------*
 *  Function Seedit  |
 *-------------------*
 |  When called, the seed changes. Record new seeds so that
 |  the next simulation is different than previous one
*//* Celine changed 03/18/2010 */ // added , struct params *p
void seedit(char *flag, FILE *pf, struct params *p)
{
  unsigned short seedv[3], seedv2[3], *seed48(unsigned short *), *pseed=NULL; // The initializer function seed48() sets the value of Xi to the 48-bit value specified in the argument array
  int i=0;
  FILE *fopen(const char *, const char *), *pfseed=NULL;                      // Function FileOpen, pointer on file with seed for randomization

  if(flag[0]=='s') 
    {
      pfseed=fopen("seedmimar", "r");             // Open file seedmimar in reading
      //--- If seedms doesn't exist create it ---//
      if(pfseed==NULL) 
        {seedv[0]=3579; seedv[1]=27011; seedv[2]=59243;}
      else
        {
          seedv2[0]=3579; seedv2[1]=27011; seedv2[2]=59243;
          for(i=0;i<3;i++)
            if(fscanf(pfseed, " %hd", seedv+i)<1) // %hd for short int
              seedv[i]=seedv2[i];
          fclose(pfseed);                         // Close seedimar & free memory
        }
      seed48(seedv);
      /* Celine changed 03/18/2010 */
      //fprintf(pf, "\n%d %d %d\n", seedv[0], seedv[1], seedv[2]);
      for(i=0;i<3;i++)
        p->tableseeds[i]=seedv[i];
      /*/////*/
    }
  else                                            // Write seeds in seedmimar
    {
      pfseed=fopen("seedmimar", "w");
      pseed=seed48(seedv);
      fprintf(pfseed, "%d %d %d\n", pseed[0], pseed[1], pseed[2]);
      fclose(pfseed);
    }
}

/* Celine changed 03/18/2010 */
int
commandlineseed( char **seeds, struct params *p)
{
  unsigned short seedv[3], *seed48();
  seedv[0] = atoi( seeds[0] );
  seedv[1] = atoi( seeds[1] );
  seedv[2] = atoi( seeds[2] ); 
  int i;
  for(i=0;i<3;i++)
    p->tableseeds[i]=seedv[i];
  seed48(seedv);
  return(3);
} /*/////*/

/******************** Distributions **************/
/*------------------------*
 *  Function Random Real  |
 *------------------------*
 |  Get a double random number between low and high
*/
double RandomReal(double low, double high)
{
  double d=0.0, ran1();
  d=(double) ran1();                // ((double) RAND_MAX+1);
  return (low+d*(high-low));
}

/*----------------------------*
 *  Function Random Interger  |
 *----------------------------*
 |  Get a random interger number between low and high !INCLUSIVE
*/
int RandomInteger(int low, int high)
{
  int k=0;
  double d=0.0, ran1();
  d=(double) ran1();                // ((double) RAND_MAX+1);
  k=(int)(d*(high-low+ 1));
  return (low+k);
}

/*----------------------------*
 *  Function uniform (0, 1)    |
 *----------------------------*
 |  Get a random number between 0 and 1
*/
double rnd()
{
  double value=0.0, RandomReal(double, double);
  do
    value=RandomReal(0.0, 1.0);
  while((value==0.0)||(value==1.0));
  return value;
}

/*--------------------------*
 *  Function  Uniform [a, b] |
 *--------------------------*
 | Called when sampling from uniform priors in params.c 
*/
double unif(double min, double max, int zero)
{
  double a=0, b=0;
  double temp=0.0, ran1();
  a=min;
  b=max;
  if(zero==0)                       // zero=0 false
    while(temp==0)
      temp=ran1()*(b-a)+a;
  else if (zero==1)
    temp=ran1()*(b-a)+a;
  else
    {
      if(ran1()<.5)
        while(temp==0)
          temp=ran1()*(b-a)+a;
      else temp=0;
    }
  return temp;
}

/*------------------------*
 *  Function  exponantial |
 *------------------------*
 |  Mean is 1/lambda
*/
double randexp(double lambda)
{
  double dum=0.0, ran1();
  while(dum==0.0)
    dum=ran1();
  return (-log(dum)/lambda);
}

/*-------------------*
 *  Function Normal  |
 *--------------------*
 | Returns Normal rv with mean mu, variance sigsq.
 | Uses snorm function of Brown and Lovato. By JKP
*/
double RNormal(double mu, double sd) 
{
  double snorm();
  return (mu+ sd*snorm());
}
/*  
    (STANDARD-)  N O R M A L  DISTRIBUTION
    --------------------------------------------

    FOR DETAILS SEE:

    AHRENS, J.H. AND DIETER, U.
    EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM
    SAMPLING FROM THE NORMAL DISTRIBUTION.
    MATH. COMPUT., 27, 124 (OCT. 1973), 927 - 937.

    ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'
    (M=5) IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)

    Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
    SUNIF. The argument IR thus goes away.

    THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
    H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
*/

double snorm()
{
  static double a[32]=
    {
      0.0, 3.917609E-2, 7.841241E-2, 0.11777, 0.1573107, 0.1970991, 0.2372021, 0.2776904, 
      0.3186394, 0.36013, 0.4022501, 0.4450965, 0.4887764, 0.5334097, 0.5791322, 
      0.626099, 0.6744898, 0.7245144, 0.7764218, 0.8305109, 0.8871466, 0.9467818, 
      1.00999, 1.077516, 1.150349, 1.229859, 1.318011, 1.417797, 1.534121, 1.67594, 
      1.862732, 2.153875
    };
  static double d[31]=
    {
      0.0, 0.0, 0.0, 0.0, 0.0, 0.2636843, 0.2425085, 0.2255674, 0.2116342, 0.1999243, 
      0.1899108, 0.1812252, 0.1736014, 0.1668419, 0.1607967, 0.1553497, 0.1504094, 
      0.1459026, 0.14177, 0.1379632, 0.1344418, 0.1311722, 0.128126, 0.1252791, 
      0.1226109, 0.1201036, 0.1177417, 0.1155119, 0.1134023, 0.1114027, 0.1095039
    };
  static double t[31]=
    {
      7.673828E-4, 2.30687E-3, 3.860618E-3, 5.438454E-3, 7.0507E-3, 8.708396E-3, 
      1.042357E-2, 1.220953E-2, 1.408125E-2, 1.605579E-2, 1.81529E-2, 2.039573E-2, 
      2.281177E-2, 2.543407E-2, 2.830296E-2, 3.146822E-2, 3.499233E-2, 3.895483E-2, 
      4.345878E-2, 4.864035E-2, 5.468334E-2, 6.184222E-2, 7.047983E-2, 8.113195E-2, 
      9.462444E-2, 0.1123001, 0.136498, 0.1716886, 0.2276241, 0.330498, 0.5847031
    };
  static double h[31]=
    {
      3.920617E-2, 3.932705E-2, 3.951E-2, 3.975703E-2, 4.007093E-2, 4.045533E-2, 
      4.091481E-2, 4.145507E-2, 4.208311E-2, 4.280748E-2, 4.363863E-2, 4.458932E-2, 
      4.567523E-2, 4.691571E-2, 4.833487E-2, 4.996298E-2, 5.183859E-2, 5.401138E-2, 
      5.654656E-2, 5.95313E-2, 6.308489E-2, 6.737503E-2, 7.264544E-2, 7.926471E-2, 
      8.781922E-2, 9.930398E-2, 0.11556, 0.1404344, 0.1836142, 0.2790016, 0.7010474
    };
  static long i;
  static double snorm, u, s, ustar, aa, w, y, tt;

  u=rnd();
  s=0.0;
  if(u>0.5) s=1.0;
  u+=(u-s);
  u=32.0*u;
  i=(long) (u);
  if(i==32) i=31;
  if(i==0) goto S100;
  /* START CENTER */
  ustar=u-(double)i;
  aa=*(a+i-1);
 S40:
  if(ustar<=*(t+i-1)) goto S60;
  w=(ustar-*(t+i-1))**(h+i-1);
 S50:
  /* EXIT (BOTH CASES) */
  y=aa+w;
  snorm=y;
  if(s==1.0) snorm=-y;
  return snorm;
 S60:
  /* CENTER CONTINUED */
  u=rnd();
  w=u*(*(a+i)-aa);
  tt=(0.5*w+aa)*w;
  goto S80;
 S70:
  tt=u;
  ustar=rnd();
 S80:
  if(ustar>tt) goto S50;
  u=rnd();
  if(ustar>=u) goto S70;
  ustar=rnd();
  goto S40;
 S100:
  /* START TAIL */
  i=6;
  aa=*(a+31);
  goto S120;
 S110:
  aa+=*(d+i-1);
  i+=1;
 S120:
  u+=u;
  if(u<1.0) goto S110;
  u -=1.0;
 S140:
  w=u**(d+i-1);
  tt=(0.5*w+aa)*w;
  goto S160;
 S150:
  tt=u;
 S160:
  ustar=rnd();
  if(ustar>tt) goto S50;
  u=rnd();
  if(ustar>=u) goto S150;
  u=rnd();
  goto S140;
}
