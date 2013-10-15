/*! \file rand1.h
  \brief Functions of random generator.

  - Functions for seed generation in msR.c, msC.c and estlilC.c
  - Functions specifics to distributions.
*/

/*** In here ***/
double ran1();      
void seedit( char *flag , struct params *p );
int  commandlineseed( char **seeds, struct params *p);
/* added for Lik stat calculation */
double lnpo(double u, int y);
double RandomReal(double low, double high);
int RandomInteger(int low, int high);
double rnd();
double unif(double min, double max, int zero);
double randexp(double lambda);
double RNormal(double mu, double sd);
double snorm();

/* external call */
double drand48();
unsigned short  *seed48();
double log(double);
