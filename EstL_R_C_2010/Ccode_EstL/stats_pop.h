/*! \file stats_pop.h
  \brief Header file for stats_pop.c.

*/


int maxsites = 1000 ;
double **Sstat;
double *FST;
double nn2;

/***/
/*! \brief cmd line Word Check 

 -  Checks that arg is not "-" before one expects a flag 
 - Called by main()
*/
void argcheck_stat( 
              int arg, //!< arg/word #
              int argc, //!< tot # of words in cmd line
              char *argv[] //!< cmd line
               );

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

/*!\brief Lists usage of options and tags.

 Called by main()
*/	
int usage_stats();
