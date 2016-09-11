/*! \file streec.h
  \brief Struct & functions from streecC.c and streecR.c.
*/


/*! \brief Recombining segments */
struct seg{
  int beg;
  int end;
  int desc;
};
/*! \brief Recombining segments */
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


/*** called by estlikC ***/
double * computelstat(int nsam1,int nsam2 ,struct node *pptree);

int make_gametes2(int nsam, int nsam1,int nsam2, struct node *ptree, double tt, int newsites, int ns, char **list );
int mnmial2(int n,int nclass,double p[],int rv[],double tt);


/*** ARG generation ***/
struct segl * segtre_mig(struct c_params *cp, int *pnsegs ) ;

/***  recombination subroutine ***/
int xover(int nsam,int ic, int is);
int re(int nsam);
int cleftr( int nsam);
int cinr( int nsam, int nsites);
/** common ancestor subroutine */
void parens( struct node *ptree, int *descl, int *descr,  int noden); 
int ca(int nsam,int  nsites,int c1,int c2);
int isseg(int start, int c,int *psg);
int pick2_chrom(int pop ,int config[],int *pc1,int *pc2);
/*! links(c): returns the number of links between beginning and end of chrom **/
int links(int c);

/*** make_gametes.c  **/
int make_gametes(int nsam, int mfreq, struct node *ptree, double tt, int newsites, int ns, char **list );
double ttime(struct node * ptree,int nsam);
double ttimemf(  struct node * ptree,int nsam,int mfreq);
void prtree(struct node *  ptree,int nsam);
 
/***  pickb **/
int pickb(int nsam,struct node * ptree,double tt);
int pickbmf(int nsam, int mfreq,struct node * ptree,double tt );
int tdesn(  struct node *ptree,int  tip,int  node );
int pick2(int n,int *i,int *j);
int ordran(int n,double pbuf[]);
int mnmial(int n,int nclass,double p[],int rv[]);
int order(int n,double pbuf[]);
int ranvec(int n,double pbuf[]);
int poisso(double u);
/*! a slight modification of crecipes version */
double gasdev(double m,double v);
