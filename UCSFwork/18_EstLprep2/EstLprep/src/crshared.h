 
/* msC and R and estlikC and R header */
#define SITESINC 10 

struct seg{
  int beg;
  int end;
  int desc;
};

struct chromo{
  int nseg;
  int pop;
  struct seg  *pseg;
};

struct segl {
  int beg;
  struct node *ptree;
  int next;
};

struct node{
  int lpop1;                                  // # lineages from pop1 as descent
  int lpop2;                                  // # lineages from pop2 as descent
  int abv;
  int ndes;
  float time;
};


double *posit ;
double segfac ;
int count, ntbs, nseeds ;
struct params pars ;	


