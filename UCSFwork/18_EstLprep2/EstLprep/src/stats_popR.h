/*! \file stats_popR.h
  \brief Header file for stats_popR.c.

*/

int maxsites = 1000 ;
double **Sstat;
double *FST;
double nn2;

/***/
/*! \brief Memory reallocation for the list of haplotypes 

  - Used when S>maxsites
  - Called by main()
*/
int biggerlist_stat(
                    int nsam,  //!< number of chromosomes sampled
                    unsigned nmax , //!< max num sites until now
                    char ** list  //!< list of haplotypes
                    );
