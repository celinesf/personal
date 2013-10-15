/*! \file statsfun.h
  \brief Header specific to statsfun.c.

  - Functions for seed generation in msR.c, msC.c and estlilC.c
  - Functions specifics to distributions.
*/


/* stats_pop only */
/********************************************************************/

/***/
/*! \brief Calculate statistics used for estimation */

/*! Called by fucntion get_cov() estlikC.c 
  Output is of the form
  0 Seg sites
  1 Freq S
  2 FST
  3 pi
  4 thetaW
  5 thetaH
  6  TajD 
  7 Hfay 
  8 D*  
  9 S_o  singletons   
  10 LD
  11 r_square
  12 nh
  13 Rm
  */
double ** getstatS2(int * popsize,int nsam,int segsites, char **list, double *typestats, int * pwstat, int* pos);

/***/
/*! Called by the main() fucntion in stat_popC.c and stat_popR.c
  Output is of the form
  0 Seg sites
  1 Freq S
  2 pi
  3 thetaW''
  4 thetaH
  5  TajD 
  6 Hfay 
  7 D*  
  8 S_o  singletons   
  9 LD
  10 r_square
  11 nh
  12 Rm
  */
double ** getstatS(int * popsize,int nsam,int segsites, char **list, double *typestats, int* pos);

/***/
/*! \brief Calculates mean pairwise differences between and within pop. and Fst.

  Called by getstatS() in stats_fun.c.
 */
double * getFST(int * popsize, int nsam, int segsites, char ** list, double *typestats);

/***/
/*! \brief Calculates the nearest neighbor statistics.

  Called by getstatS() in stats_fun.c and get_cov() in estlikC.c.
 */
double getnn2(int * popsize, int nsam, int segsites, char ** list, int typestats);
int distance(int a, int b, int segsites, char ** list);
int distance_ind(int a, int b, int segsites, char ** list);

/***/
/*! \brief Calculates Tajima's D and Fu and Li's D*.

  Called by getstatS() and getstatS2() in stats_fun.c.
 */
double* tajD_D_star(int N, int S, double k, int sing, int type);

 
double psum (int n);
int frequency( char allele,int site,int nsam,  char **list);
int myfrequency( char allele,int site,int start, int nsam,  char **list);

/***/
/*! \brief Memory allocation for haplotype list. 

  Called by getstatS() and getstatS2() in stats_fun.c and
main() functions of estlikC.c, stats_popC.c and stats_popR.c
*/
char ** cmatrix(
                int nsam ,//!< # of chromosomes
                int len//!< #number of seg sites
);


/* LD */
double freq2 (int a, int b, int start,int n, char **list);
double max2 (double a, double b);
double min2 (double a, double b);
int hapcount (char **list,int start, int nsam, int first, int last);
double * LDcalc(double p1, double p2, double p12); 
void four_gametes_test(int si,int sj, int nsam, int *pop,char ** plist, int **pgamete,int haplo);
int min_rec (int x, int size, int start, int ** pgamete, int fl);
int isgam(int c,int n1,int *pgtest, int gam);
double* get_cor(int j, int i,int psize,int ** pfreq,int ns,char **plist);

