/************ make_gametes.c - mimargof *******************************************
 * Contains functions relenvant for tree building
 * Assigns segregating sites to different branches.
 * Records what kind of site it is: S1 S2 Ss Sf.
 *
 ************************************** Changes

  *--- Changed Feb. 11st 2010: ---*
  - Changed make_gametes to assume that ancestral states are unknow and use allele 
    with the minor frequencies to obtain the summary staristics. 
  -----------------------------

 *************************************************************************/
#include <stdio.h>                  // Input output Library
#include <stdlib.h>                 // File gestion Library
#include <math.h>                   // Math Library

#define STATE1 '1'
#define STATE2 '0'

/*------------------*
 *  Structure Node  |
 *------------------*/
struct node
{
  int abv;                          // Node above=accestor
  int lpop1;                        // # lineages from pop1 as descent
  int lpop2;                        // # lineages from pop2 sample
  float time;                       // Time of the node
};

/*------------------------*
 *  Function make_Gametes |
 *------------------------*
 |  Assign states in tree for new segregating sites
 |  called by gensam
*/ 
void make_gametes(int nsam, struct node *ptree,double tt, int newsites, int ns, char **list , int npop, int *popnhap, int **freqvar, int ptyps[], int pS[])
/* newsites, ns, npop, *popnhap; // # chromo, # new sites in this segment, # seg sites until now, npop, sample size of each pop
   struct node *ptree;           // Gene genealogy
   char **list;                  // List of variant
   double tt;                    // TMRCA for the considered segment
   int **freqvar;                // Frequency spectrum at each site in all populations
   int ptyps[], pS[];            // array of Type of seg sites, of S statistics
*/
{
  int tiphap=0, ipop=0, jseg=0, kpophap=0, node=0;
  int tdesn(struct node *, int, int), pickb(int, struct node *, double);
  //--- loop along new segregating sites ---//
  for(jseg=ns;jseg<ns+newsites;jseg++)
    {
      node=pickb(nsam, ptree, tt);                       // Pick a branch at random in the tree
      tiphap=0;                                          // num of chromosome

      for(ipop=0;ipop<npop;ipop++)                       // Loop on populations
        {
          for(kpophap=0;kpophap<popnhap[ipop];kpophap++) // loop to assign variant state to descent of this node
            {
              if(tdesn(ptree, tiphap, node))
                {
                  list[tiphap][jseg]=STATE1;            // variant
                  freqvar[0][jseg]++;                   // +1 total Freq for this site
                  freqvar[ipop+1][jseg]++;              // +1 pop specific Freq for this site
                  if(ptyps[jseg]!=-1)                   // Case site not shared
                    {
                      if(freqvar[ipop+1][jseg]==0) ptyps[jseg]=0;                     // Case variant in other pop
                      else if(freqvar[0][jseg]==freqvar[ipop+1][jseg]) ptyps[jseg]=1; // Case seg site variable in this pop
                      else ptyps[jseg]=-1;                                            // Case site shared
                      if(ptyps[jseg]>0)       
                        {
                          if(freqvar[ipop+1][jseg]==popnhap[ipop]) ptyps[jseg]=npop+1; // Case: Fixed in the population
                          else
                            ptyps[jseg]=ipop+1;                                        // Case Unic in this population
                        }
                    }               
                }
              else list[tiphap][jseg]=STATE2;           // Ancestral state
              tiphap++;                                 // Check next chromosome 
            }
          if((freqvar[ipop+1][jseg]>0)&&(freqvar[ipop+1][jseg]<popnhap[ipop]))
            pS[ipop]++;                                // Add a unic site for the population
        }
      /* Celine change 02/11/2010 */// recalculate for minor allele freq
      if(freqvar[0][jseg]>(int) nsam/2 )
        {
          freqvar[0][jseg]=nsam- freqvar[0][jseg];
          for(ipop=0;ipop<npop;ipop++)
            {
              if((freqvar[ipop+1][jseg]>0)&&(freqvar[ipop+1][jseg]<popnhap[ipop]))
                pS[ipop]--; 
              freqvar[ipop+1][jseg]=popnhap[ipop]- freqvar[ipop+1][jseg];
            }

          if(freqvar[1][jseg]>0 && freqvar[2][jseg]>0) // shared
            ptyps[jseg]=-1;                                            // Case site shared
          else if(freqvar[0][jseg]==freqvar[1][jseg])
            ptyps[jseg]=1; // Case seg site variable in this pop
          else
            ptyps[jseg]=2;                                        // Case Unic in this population
       
          for(ipop=0;ipop<npop;ipop++)              
            if((freqvar[ipop+1][jseg]>0)&&(freqvar[ipop+1][jseg]<popnhap[ipop]))
              pS[ipop]++;    
        }//////
      ipop=0;
      if(ptyps[jseg]==-1)
        pS[4]+=freqvar[ipop+1][jseg]+freqvar[ipop+2][jseg]; // freq i=2 (p2) freq derived shared
      if(ptyps[jseg]==ipop+1)
        pS[ipop+2]+=freqvar[ipop+1][jseg];
      ipop=1;
      if(ptyps[jseg]==ipop+1)
        pS[ipop+2]+=freqvar[ipop+1][jseg];
    } 
}// End of make_gametes function


/*------------------*
 *  Function  ttime |
 *------------------*
 |  Returns the total time in the tree, *ptree, with nsam tips. 
 |  Called by gensam.
*/
double ttime(struct node *ptree, int nsam)
{
  double t=0.0;
  int i=0;
  t=(ptree+2*nsam-2)->time;
  for(i=nsam;i<2*nsam-1;i++)
    t+=(ptree+i)->time;
  return(t);
}

/*** The two following functions are used for genetree output ***/
/*----------*
 *  prtree  |
 *----------*
 | builts tree information for a segment
*/
void prtree(struct node *ptree, int nsam)
{
  int i=0, *descl=NULL, *descr=NULL;
  void parens(struct node *ptree, int *descl, int *descr, int noden);

  descl=(int*)malloc((unsigned)(2*nsam-1)*sizeof(int));
  descr=(int*)malloc((unsigned)(2*nsam-1)*sizeof(int));
  for(i=0;i<2*nsam-1;i++) descl[i]=descr[i]=-1;
  for(i=0;i<2*nsam-2;i++)
    {
      if(descl[(ptree+i)->abv]==-1) descl[(ptree+i)->abv]=i;
      else descr[(ptree+i)->abv]=i;
    }
  parens(ptree, descl, descr, 2*nsam-2); // Write the characteristic of the tree
  free(descl);
  free(descr);
}

/*----------*
 *  parens  |
 *----------*
 | Write the gene tree
*/
void parens(struct node *ptree, int *descl, int *descr, int noden)
{
  double time=0.0;

  if(descl[noden]==-1)
    printf("%d:%5.3lf", noden+1, (ptree+((ptree+noden)->abv))->time);
  else
    {
      printf("(");
      parens(ptree, descl, descr, descl[noden]);
      printf(", ");
      parens(ptree, descl, descr, descr[noden]);
      if((ptree+noden)->abv==0) printf(");\n");
      else
        {
          time=(ptree+(ptree+noden)->abv)->time-(ptree+noden)->time;
          printf("):%5.3lf", time);
        }
    }
}

/*---------*
 *  pickb  |
 *---------*
 | returns a random branch from the tree. The probability of picking
 | time in tree.  
*/
int pickb(int nsam, struct node * ptree, double tt)
{
  double x=0.0, y=0.0, ran1();
  int i=0;
  x=ran1()*tt;
  for(i=0, y=0;i<2*nsam-2;i++)
    {
      y+=(ptree+(ptree+i)->abv)->time-(ptree+i)->time;
      if(y >=x) return(i);
    }
  return(i);
}

/*---------*
 *  tdesn  |
 *---------*
 | returns 1 if tip is a descendant of node in *ptree, otherwise 0.
*/
int tdesn(struct node *ptree, int tip, int node)
{
  int k=0;
  for(k=tip;k<node;k=(ptree+k)->abv);
  if(k==node) return(1);
  else return(0);
}

/*---------*
 *  Pick2  |
 *---------*
 | Pick two numbers (i and j) from n
 | Called in street.c by pick2_chrom
*/
int pick2(int n, int *i, int *j)
{
  double ran1();
  *i=n*ran1();
  while((*j=n*ran1())==*i);
  return(0);
}


/*-------------------*
 *  Function ordran  |
 *-------------------*
 | Generates a table of n ordered random values
 | called by function locate
*/
void ordran(int n, double pbuf[])
{
  void ranvec(int, double*);
  void order(int, double*);
  ranvec(n, pbuf);
  order(n, pbuf);
}

/*-------------------*
 *  Function mnmial  |
 *-------------------*
 | computes # segsites for this segment 
 | called by gensam
*/
int mnmial(int n, int nclass, double p[], int rv[])  // # segsites given by user, # segments, SS: # segsites, p[] time tree for each segment
{
  double x=0.0, s=0.0, ran1();
  int i=0, j=0;

  for(i=0;i<nclass;i++) rv[i]=0;    // Initialization
  for(i=0;i<n;i++)
    {
      x=ran1();                     // time allowed for next segment
      j=0;
      s=p[0];                       // Time 1st segment
      while((x>s)&&(j<(nclass-1))) s+=p[++j]; // s=sum timing of segments until over x
      rv[j]++;                      // j it stoped to get one more seg sites
    }
  return(j);                        // return last j
}

/*------------------*
 *  Function order  |
 *------------------*
 | Orders a table of n values
*/
void order(int n, double pbuf[])    // n Size of table, pointer Buffer
{
  int gap=0, i=0, j=0;
  double temp=0.0;                  // Temporary value

  //--- Loop dicotomy to place value n ---//
  for(gap=n/2; gap>0; gap/=2)
    for(i=gap;i<n;i++)
      for(j=i-gap;j>=0 && pbuf[j]>pbuf[j+gap];j-=gap)
        {
          temp=pbuf[j];
          pbuf[j]=pbuf[j+gap];
          pbuf[j+gap]=temp;
        }
}

/*------------------*
 *  Functionranvec  |
 *------------------*
 | Generates Table of n Random values 
*/
void ranvec(int n, double pbuf[])
{
  int i=0;
  double ran1();                    // Pick Random Function as define in function ran1 (file rand1.c or rand2.c depending of compilation)

  for(i=0;i<n;i++)
    pbuf[i]=ran1();
}

/*---------------------------------*
 *  Function Poisson Distribution  |
 *---------------------------------*
 | Returns the integer value corresponding to the
 | probablity given by the poisson distribution.
*/
int poisso(double u)  //u= Mean of poisson distribution
{
  double cump=0.0, ru=0.0, p=0.0, ran1(), gasdev(double, double);
  int i=1;

  if(u>30.) return((int)(0.5+gasdev(u, u))); // Correction when mean > 30 ?? why?
  ru=ran1();
  p=exp(-u);
  if(ru<p) return(0);
  cump=p;                           // Cumulative probability function

  while(ru>(cump+=(p*=u/i)))
    i++;
  return(i);
}

/*-------------------*
 *  Function gasdev  |
 *-------------------*
 | a slight modification of crecipes version
*/
double gasdev(double m, double v)                            // m=v=u (mean poisson >30)
{
  static int iset=0;
  static float gset=0.0;
  float fac=0.0, r=0.0, v1=0.0, v2=0.0;
  double ran1();
  if(iset==0) 
    {
      do              // Compute probability r=v1^2+v2^2
        {
          v1=2.0*ran1()-1.0;
          v2=2.0*ran1()-1.0;
          r=v1*v1+v2*v2;
        } while(r>=1.0);
      fac=sqrt(-2.0*log(r)/r);
      gset=v1*fac;
      iset=1;
      return(m+sqrt(v)*v2*fac);
    } 
  else
    {
      iset=0;
      return(m+sqrt(v)*gset);
    }
}

