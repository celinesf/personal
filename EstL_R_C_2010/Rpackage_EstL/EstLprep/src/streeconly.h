/*! \file streeconly.h
  \brief Header for streecC.c and streecR.c only.
*/


#define MIN(x, y) ( (x)<(y) ? (x) : (y) )
#define NL putchar('\n')
#define STATE1 '1'
#define STATE2 '0'
#define STATE3 '?'
#define size_t unsigned
#define SEGINC 80
#define NLSTATS 15  

int nchrom, begs, nsegs;
long nlinks ;
double t, cleft , pc, lnpc ;

struct node *ptree1=NULL, *ptree2=NULL;

static int *nnodes = NULL ;
static unsigned seglimit = SEGINC ;
static unsigned maxchr ;
static struct chromo *chrom = NULL ;
static struct segl *seglst = NULL ;

extern int flag;


