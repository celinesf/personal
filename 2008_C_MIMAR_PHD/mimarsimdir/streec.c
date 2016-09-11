/********** segtre_mig.c - mimarsim **********************************
 *
 *	This subroutine uses a Monte Carlo algorithm described in
 *	Hudson,R. 1983. Theor.Pop.Biol. 23: 183-201, to produce
 *	a history of a random sample of gametes under a neutral
 *	Wright-Fisher model with recombination and geographic structure. 
 *	Input parameters
 *	are the sample size (nsam), the number of sites between
 *	which recombination can occur (nsites), and the recombination
 *	rate between the ends of the gametes (r). The function returns
 *	nsegs, the number of segments the gametes were broken into
 *	in tracing back the history of the gametes. The histories of
 *	these segments are passed back to the calling function in the
 *	array of structures seglst[]. An element of this array, seglst[i], 
 * 	consists of three parts: (1) beg, the starting point of
 *	of segment i, (2) ptree, which points to the first node of the
 *	tree representing the history of the segment, (3) next, which
 *	is the index number of the next segment.
 *	   A tree is a contiguous set of 2*nsam nodes. The first nsam
 *	nodes are the tips of the tree, the sampled gametes. The other
 *	nodes are the nodes ancestral to the sampled gametes. Each node
 *	consists of an "abv" which is the number of the node ancestral to
 *	that node, an "ndes", which is not used or assigned to in this routine, 
 *	and a "time", which is the time (in units of 4N generations) of the
 *	node. For the tips, time equals zero.
 *	Returns a pointer to an array of segments, seglst.
 *------------------------------------------------------------
 |
 | The functions in these files are called in function gensam in mimar.c to
 | generate sample genealogies.
 | Most of these have not been modified from ms. The comments were added to help me
 | read through the program. 
 |
 ************************************** Changes
 
 *--- Changed Nov. 27th 2009: ---* 
 Corrected a bug on memory allocation (line 139-149).
 Thanks to the user who mentionned the bug to me.
 -----------------------------
 *
 *****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mimarsim.h"
#define NL putchar('\n')
#define size_t unsigned
#define MIN(x, y) ( (x)<(y) ? (x) : (y) )
#define ERROR(message) fprintf(stderr, message), NL, exit(1)
#define SEGINC 80 

extern int flag;

int nchrom, begs, nsegs;
long nlinks;
static int *nnodes =NULL;
double t, cleft , pc, lnpc;

static unsigned seglimit =SEGINC;
static unsigned maxchr;

struct seg
{
  int beg;
  int end;
  int desc;
}; struct chromo
{
  int nseg;          // # of segments
  int pop;           // population the chromosome belongs to
  struct seg *pseg;  // array of segments
};
static struct chromo *chrom=NULL;

/*------------------*
 *  Structure node  |
 *------------------*/
struct node
{
  int abv;            // Node above
  int lpop1;          // # descent lineages leading to population1
  int lpop2;          // # descent lineages leading to population2
  float time;         // time of the node in the tree
} *ptree1, *ptree2;

struct segl 
{
  int beg;
  struct node *ptree;
  int next;
};
static struct segl *seglst=NULL;

/*-------------------------------------------------*
 *  Function generating segments by recombination  |
 *-------------------------------------------------*
 |  all history events are added
*/
struct segl * segtre_mig(struct c_params *cp, int *pnsegs)
{
  int i=0, j=0, k=0, dec=0, pop=0, pop2=0, c1=0, c2=0, ind=0, rchrom=0;
  int migrant=0, source_pop=0, flagint=0, eflag=0, cpop=0, ic=0, nsam=0, npop=0, nsites=0;
  double sum=0.0, x=0.0, ttemp=0.0, rft=0.0, clefta=0.0, tmin=0.0, p=0.0, r=0.0, f=0.0, rf=0.0, track_len=0.0;
  double prec=0.0, cin=0.0, prect=0.0, mig=0.0, ran=0.0, coal_prob=0.0, rdum=0.0, arg=0.0;
  char event='\0'; 
  int *inconfig=NULL, *config=NULL;
  double *size=NULL, *alphag=NULL, *tlast=NULL, **migm=NULL;
  struct devent *nextevent=NULL;
  double ran1();
  int re(int), xover(int, int, int), cinr(int, int), cleftr(int), ca(int, int, int, int);
  void pick2_chrom(int, int *, int *, int *);

  nsam=cp->nsam;
  npop=cp->npop;
  nsites=cp->nsites;

  inconfig=cp->config;
  r=cp->r;
  f=cp->f;
  track_len=cp->track_len;

  migm=(double**)malloc((unsigned)npop*sizeof(double*));
  for(i=0;i<npop;i++)
    {
      migm[i]=(double*)malloc((unsigned)npop*sizeof(double));
      for(j=0;j<npop;j++) migm[i][j]=(cp->mig_mat)[i][j];
    }
  nextevent=cp->deventlist;

  /* Initialization */
  //--- Allocation for pointer on structures ---//
  if(chrom==NULL) 
    {
      maxchr=nsam+20;
      chrom=(struct chromo*)malloc((unsigned)(maxchr*sizeof(struct chromo)));
      if(chrom==NULL) perror("malloc error. segtre");
    } 
  /* Celine changed 11/27/2009 */ 
  else
    {
      if(nsam+20>maxchr)
        {
          free(chrom);
          maxchr=nsam+20;
          chrom=(struct chromo*)malloc((unsigned)(maxchr*sizeof(struct chromo)));
          if(chrom==NULL) perror("malloc error. segtre");
        }
    }/*//////*/
  if(nnodes==NULL)
    {
      nnodes=(int*) malloc((unsigned)(seglimit*sizeof(int)));
      if(nnodes==NULL) perror("malloc error. segtre_mig");
    }
  if(seglst==NULL) 
    {
      seglst=(struct segl*)malloc((unsigned)(seglimit*sizeof(struct segl)));
      if(seglst==NULL) perror("malloc error. segtre_mig.c 2");
    }
  //--- Allocate and fill internal param information ---//
  config=(int*)malloc((unsigned)((npop+1)*sizeof(int)));
  if(config==NULL) perror("malloc error. segtre.");
  size=(double*)malloc((unsigned)((npop)*sizeof(double)));
  if(size==NULL) perror("malloc error. segtre.");
  alphag=(double*)malloc((unsigned)((npop)*sizeof(double)));
  if(alphag==NULL) perror("malloc error. segtre.");
  tlast=(double*)malloc((unsigned)((npop)*sizeof(double)));
  if(alphag==NULL) perror("malloc error. segtre.");
  for(pop=0;pop<npop;pop++) 
    {
      config[pop]=inconfig[pop];
      size[pop]=(cp->size)[pop];
      alphag[pop]=(cp->alphag)[pop];
      tlast[pop]=0.0;
    }
  //--- Initialize pointer of segl ---//
  seglst[0].beg=0;
  if(!(seglst[0].ptree=(struct node*)calloc((unsigned)(2*nsam), sizeof(struct node))))
    perror("calloc error. se2");// Init ptree

  for(pop=ind=0;pop<npop;pop++)               // Loop on populations
    for(j=0;j<inconfig[pop];j++, ind++)       // loop on chromozomes per pop
      {
        chrom[ind].nseg=1;                    // Initialize at least 1 segment for this chromosome
        if(!(chrom[ind].pseg=(struct seg*)malloc((unsigned)sizeof(struct seg))))
          ERROR("calloc error. se1");
        (chrom[ind].pseg)->beg=0;             // Give size of the segment
        (chrom[ind].pseg)->end=nsites-1;
        (chrom[ind].pseg)->desc=ind;          // # of descent=# chromozome
        chrom[ind].pop=pop;                   // record pop# for this chromozome
        if(pop==0)                            // record lineage of the leaf nodes
          seglst[0].ptree[ind].lpop1+=1;
        else
          seglst[0].ptree[ind].lpop2+=1;
      }
  nnodes[0]=nsam-1;                            // # ancestor nodes
  nchrom=nsam;
  nlinks=((long)(nsam))*(nsites-1);           // nchromo* # sites
  nsegs=1;
  t=0.;
  r /=(nsites-1);                             // Recombination probability per nucleotides
  //--- values for conversion ---//
  if(f>0.0) pc=(track_len-1.0)/track_len;     // Conversion lenght
  else pc=1.0;
  lnpc=log(pc);                                // Log value
  cleft=nsam*(1.0-pow(pc,(double)(nsites-1))); // value for conv
  if(r>0.0) rf=r*f;                            // recom+ conv rate
  else rf=f/(nsites-1);                        // conv rate only
  rft=rf*track_len;                            // conv rate * lenght=lenght of expansion
  flagint=0;

  /* Main loop */
  // This loop reduces the number of chromosomes(nchom--)
  while(nchrom>1) 
    { 
      prec=nlinks*r;                          // recomb rate for the sample
      cin=nlinks*rf;                          // Conversion rate
      clefta=cleft*rft;                       // Conv prob * lenght of conv
      prect=prec+cin+clefta;                  // total # cross, conv rec event
      mig=0.0;
      for(i=0;i<npop;i++) mig+=config[i]*migm[i][i]; // # chromo * rate of migration popi to popi
      if((npop>1)&&(mig==0.0)&&(nextevent==NULL))    // Crash if no recombonbination and no other events
        {
          i=0;
          for(j=0;j<npop;j++)
            if(config[j]>0) i++;
          if(i>1)
            {
              fprintf(stderr, " Infinite coalescent time. No migration.\n");
              exit(1);
            }
        }
      eflag=0;                                 // init flag for migration / recombination / conversion events
      /* cross-over or gene conversion */
      if(prect>0.0) 
        {
          while((rdum=ran1()) ==0.0);
          ttemp=-log(rdum)/prect;             // time of recombination or conversion event
          
          if((eflag==0)||(ttemp<tmin))        // decide if cross-over occures before other timing
            {
              tmin=ttemp;
              event='r';                      // recombination event planned before any other events
              eflag=1;
            }
        }
      /* migration */
      if(mig>0.0)
        {
          while((rdum=ran1())==0.0);
          ttemp=-log(rdum)/mig;
          if((eflag==0)||(ttemp<tmin))        // Mig before recomb
            {
              tmin=ttemp;
              event='m';
              eflag=1;
            }
        }
      /* coalescent */ 
      for(pop=0; pop<npop; pop++) 
        {
          coal_prob=((double)config[pop])*(config[pop]-1.); // Probability of coalsecent event n*(n-1)
          if(coal_prob>0.0) 
            {
              while((rdum=ran1())==.0);       // Take a random #
              if(alphag[pop]==0)              // No growth case
                {
                  ttemp=-log(rdum)*size[pop] /coal_prob; // Prob coal=-log(rand*size pop)/(n*(n-1))
                  if((eflag==0) ||(ttemp<tmin))          // Coalescent event first
                    {
                      tmin=ttemp;
                      event='c';              // Coalescent event if before other events
                      eflag=1;
                      cpop=pop;               // num population on which to perform coalescent event
                    }
                }
              else                            // Growth case, increase pop size exponantially
                {
                  arg =1.-alphag[pop]*size[pop]*exp(-alphag[pop]*(t-tlast[pop]))*log(rdum)/coal_prob;
                  /*Note: arg is a ratio to decide whether coalsecent event
                    occures: if arg<=0, no coalescent event within interval */ 
                  if(arg>0.0) 
                    {                         
                      ttemp=log(arg) / alphag[pop];
                      if((eflag==0) ||(ttemp<tmin)) // + put coalsecent event
                        {
                          tmin=ttemp;
                          event='c';
                          eflag=1;
                          cpop=pop;
                        }// End if
                    }
                }// Growth case
            }// Coalescent prob>0
        }// loop on pop

      if((eflag==0) &&(nextevent==NULL))
        {
          fprintf(stderr,
                  " infinite time to next event. Negative growth rate in last time interval or non-communicating subpops.\n");
          exit(0);
        }
      //--- If(no(mig/rec/coal) and Next event exist) or(next event exist and time<tmin) DO TIMED EVENT ---//
      if(((eflag==0) &&(nextevent!=NULL))||((nextevent!=NULL) && ((t+tmin)>= nextevent->time))) // t=0 at start
        {
          t=nextevent->time;                  // Take next event timing
          switch( nextevent->detype)
            {
              case 'N' :                      // pop size change
                for(pop=0;pop<npop;pop++)
                  {
                    size[pop]=nextevent->paramv;
                    alphag[pop]=0.0;
                  }
                nextevent=nextevent->nextde;
                break;
              case 'n' :
                size[nextevent->popi]=nextevent->paramv;
                alphag[nextevent->popi]=0.0;
                nextevent=nextevent->nextde;
                break;
              case 'G' :
                for(pop=0; pop<npop; pop++)
                  {
                    size[pop]=size[pop]*exp(-alphag[pop]*(t-tlast[pop]));
                    alphag[pop]=nextevent->paramv;
                    tlast[pop]=t;
                  }
                nextevent=nextevent->nextde;
                break;
              case 'g' :
                pop=nextevent->popi;
                size[pop]=size[pop]*exp(-alphag[pop]*(t-tlast[pop]));
                alphag[pop]=nextevent->paramv;
                tlast[pop]=t;
                nextevent=nextevent->nextde;
                break;
              case 'M' :
                for(pop=0; pop<npop; pop++)
                  for(pop2=0; pop2<npop; pop2++) migm[pop][pop2]=(nextevent->paramv)/(npop-1.0);
                for(pop=0; pop<npop; pop++)
                  migm[pop][pop]=nextevent->paramv;
                nextevent=nextevent->nextde;
                break;
              case 'a' :
                for(pop=0; pop<npop; pop++)
                  for(pop2=0; pop2<npop; pop2++) migm[pop][pop2]=(nextevent->mat)[pop][pop2];
                nextevent=nextevent->nextde;
                break;
              case 'm' :
                i=nextevent->popi;
                j=nextevent->popj;
                migm[i][i]+=nextevent->paramv-migm[i][j];
                migm[i][j]=nextevent->paramv;
                nextevent=nextevent->nextde;
                break;
              case 'j' :                      /* merge pop i into pop j (join) */
                i=nextevent->popi;
                j=nextevent->popj;
                config[j]+=config[i];
                config[i]=0;
                for(ic=0;ic<nchrom;ic++) if(chrom[ic].pop==i) chrom[ic].pop=j;
                /* the following was added 19 May 2007 */
                for(k=0; k<npop; k++)
                  {
                    if(k!=i)
                      {
                        migm[k][k]-=migm[k][i];
                        migm[k][i]=0.;
                      }
                  }
                /* end addition */
                nextevent=nextevent->nextde;
                break;
              case 's' :                      /*split pop i into two;p is the proportion from pop i,and 1-p from pop n+1 */
                i=nextevent->popi;
                p=nextevent->paramv;
                /* npop++;
                   config=(int*)realloc(config,(unsigned)(npop*sizeof(int)));
                   size=(double*)realloc(size,(unsigned)(npop*sizeof(double)));
                   alphag=(double*)realloc(alphag,(unsigned)(npop*sizeof(double)));
                   tlast=(double*)realloc(tlast,(unsigned)(npop*sizeof(double)));
                   tlast[npop-1]=t;
                   size[npop-1]=1.0;
                   alphag[npop-1]=0.0;
                   migm=(double**)realloc(migm,(unsigned)(npop*sizeof(double*)));
                   for(j=0;j<npop-1;j++)
                   migm[j]=(double*)realloc(migm[j],(unsigned)(npop*sizeof(double)));
                   migm[npop-1]=(double*)malloc((unsigned)(npop*sizeof(double)));
                   for(j=0;j<npop;j++) migm[npop-1][j]=migm[j][npop-1]=0.0;
                */
                config[npop-1]=0;
                config[i]=0;
                for(ic=0;ic<nchrom;ic++)
                  {
                    if(chrom[ic].pop==i) 
                      {
                        if(ran1()<p) config[i]++;
                        else 
                          {
                            chrom[ic].pop=npop-1;
                            config[npop-1]++;
                          }
                      }
                  }
                nextevent=nextevent->nextde;
                break;
            }// End switch
        }// End if(DO TIMED EVENT)
      else
        {
          t+=tmin;                               // t is time of next event(mig, rec or coal)
          if(event=='r') 
            {
              if((ran=ran1())<(prec / prect))    /*recombination*/
                {
                  rchrom=re(nsam);               // Recombination event rchrom=tot # chromosomes
                  config[chrom[rchrom].pop]+=1;  // nsam pop+=1 ->one more chromosome
                }
              else if(ran<(prec+clefta)/(prect)) /* cleft event */
                {
                  rchrom=cleftr(nsam);           // Recombination + conversion
                  config[chrom[rchrom].pop]+=1;
                }
              else                               /* cin event */
                {     
                  rchrom=cinr(nsam, nsites);     // Conversion only
                  if(rchrom>=0) config[chrom[rchrom].pop]+=1;
                }
            }
          else if(event=='m')                    /* migration event */
            {
              x=mig*ran1();
              sum=0.0;
              for(i=0;i<nchrom;i++)
                {
                  sum+=migm[chrom[i].pop][chrom[i].pop];
                  if(x<sum) break;
                }
              migrant=i;
              x=ran1()*migm[chrom[i].pop][chrom[i].pop];
              sum=0.0;
              for(i=0;i<npop;i++)
                {
                  if(i!=chrom[migrant].pop)
                    {
                      sum+=migm[chrom[migrant].pop][i];
                      if(x<sum) break;
                    }
                }
              source_pop=i;
              config[chrom[migrant].pop]-=1;
              config[source_pop]+=1;
              chrom[migrant].pop=source_pop;
            }
          else                /* coalescent event */
            { 
              /* pick the two, c1, c2 */
              pick2_chrom(cpop, config, &c1, &c2); /* c1 and c2 are chrom's to coalesce */
              dec=ca(nsam, nsites, c1, c2);        // Coalescent
              config[cpop]-=dec;                   // remove a lineage in the pop where event occured
            }// End else coal
        }// End else if rec/coal/mig event
    } // End while nchrom>1 
  *pnsegs=nsegs;
  free(config);
  free(size);
  free(alphag);
  free(tlast);
  for(i=0;i<npop;i++) free(migm[i]);
  free(migm);
  return(seglst);
}// End segtr_migr

/*-------------------------------------------------*
 *  Recombination subroutine                       |
 *-------------------------------------------------*
 |       Picks a chromosome and splits it in two parts. If the x-over point
 |       is in a new spot, a new segment is added to seglst and a tree set up
 |       for it.
*/
int re(int nsam)
{
  struct seg *pseg=NULL;
  int  el=0, lsg=0, lsgm1=0, ic=0, is=0, spot=0;
  double ran1();
  int xover(int, int, int);
  /* First generate a random x-over spot, then locate it as to chrom and seg. */
  spot=nlinks*ran1()+1.;                  // breakpoint of recombination

  /* get chromosome #(ic) */
  for(ic=0;ic<nchrom;ic++)
    {
      lsg=chrom[ic].nseg;                 // # segments
      lsgm1=lsg-1;                        // # seg-1
      pseg=chrom[ic].pseg;                // pointer on segment
      el=((pseg+lsgm1)->end)-(pseg->beg); // size of segment
      if(spot<=el) break;                 // seg>spot: ends before next segment
      spot-=el;
    }
  is=pseg->beg+spot-1;                    // is=size seg to recombine 
  xover(nsam, ic, is);                    // ic==nchromo==nsam
  return(ic);
}

/*-------------------------------------------*
 *  Function Recombination + gene conversion |
 *-------------------------------------------*/
int cleftr(int nsam)
{
  struct seg *pseg=NULL;
  int  ic=0, is=0;
  double ran1(), x=0.0, sum=0.0, len=0.0;
  int xover(int, int, int), links(int);

  while((x=cleft*ran1())==0.0);
  sum=0.0;
  ic=-1;
  while(sum<x)
    sum+=1.0-pow(pc, links(++ic));
  pseg=chrom[ic].pseg;
  len=links(ic);
  is=pseg->beg+floor(1.0+log(1.0-(1.0-pow(pc, len))*ran1())/lnpc)-1;
  xover(nsam, ic, is);
  return(ic);
}

/*----------------------*
 *  Function conversion |
 *----------------------*/
int cinr(int nsam, int nsites)
{
  struct seg *pseg=NULL;
  int len=0, el=0, lsg=0, lsgm1=0, ic=0, is=0, spot=0, endic=0;
  double ran1();
  int xover(int, int, int), ca(int, int, int, int);

  /* First generate a random x-over spot, then locate it as to chrom and seg. */
  spot=nlinks*ran1()+1.;
  /* get chromosome #(ic) */
  for(ic=0;ic<nchrom;ic++)
    {
      lsg=chrom[ic].nseg;
      lsgm1=lsg-1;
      pseg=chrom[ic].pseg;
      el=((pseg+lsgm1)->end)-(pseg->beg);
      if(spot<=el) break;
      spot-=el;
    }
  is=pseg->beg+spot-1;
  endic=(pseg+lsgm1)->end;
  xover(nsam, ic, is);
 
  len=floor(1.0+log(ran1())/lnpc);
  if(is+len>=endic) return(ic);
  if(is+len<(chrom[nchrom-1].pseg)->beg)
    {
      ca(nsam, nsites, ic, nchrom-1);
      return(-1);
    }
  xover(nsam, nchrom-1, is+len);
  ca(nsam, nsites, ic, nchrom-1);
  return(ic);
}

/*----------------------*
 *  Function cross over |
 *----------------------*
 | ic=nchrom, is=size segment to cross over
*/
int xover(int nsam, int ic, int is)
{
  struct seg *pseg=NULL, *pseg2=NULL;
  int i=0, lsg=0, lsgm1=0, newsg=0, jseg=0, k=0, in=0;
  double ran1(), len=0.0;
  
  pseg=chrom[ic].pseg;               // pointer on last segment of chromosom
  lsg=chrom[ic].nseg;                // # seg
  len=(pseg+lsg-1)->end-pseg->beg;   // lenght total segments
  cleft-=1-pow(pc, len);
  
  /* get seg #(jseg) */
  for(jseg=0;is>=(pseg+jseg)->end;jseg++); // find seg where breakpoint is
  if(is>=(pseg+jseg)->beg) in=1;    // If breakpoint in the seg in=1
  else in=0;                        // if breakpoint at begining of seg in=0
  newsg=lsg-jseg;                   // position new segment: nseg-jseg
  
  /* copy last part of chrom to nchrom */
  nchrom++;
  if(nchrom>=maxchr) 
    {
      maxchr+=20;
      chrom=(struct chromo*)realloc(chrom,(unsigned)(maxchr*sizeof(struct chromo)));
      if(chrom==NULL) perror("malloc error. segtre2");
    }
  if(!(pseg2=chrom[nchrom-1].pseg=(struct seg*)calloc((unsigned)newsg, sizeof(struct seg))))
    ERROR(" alloc error. re1");
  chrom[nchrom-1].nseg=newsg;        // new # segments
  chrom[nchrom-1].pop=chrom[ic].pop; // put pop # for the chromosome
  pseg2->end=(pseg+jseg)->end;       // Lenght of segment
  if(in)                             // if spot in segj
    {
      pseg2->beg=is+1;               // new segment 
      (pseg+jseg)->end=is;            // cut before end of segj
    }
  else pseg2->beg=(pseg+jseg)->beg;  // else cut at beg of segj
  pseg2->desc=(pseg+jseg)->desc;     // descents new seg=descents of segj
  for(k=1;k<newsg;k++)               // loop: copy end of chromosome after new the segment
    {
      (pseg2+k)->beg=(pseg+jseg+k)->beg;
      (pseg2+k)->end=(pseg+jseg+k)->end;
      (pseg2+k)->desc=(pseg+jseg+k)->desc;
    }
  lsg=chrom[ic].nseg=lsg-newsg+in;   //=nseg-# seg cut +in
  lsgm1=lsg-1;
  nlinks-=pseg2->beg-(pseg+lsgm1)->end;
  len=(pseg+lsgm1)->end-(pseg->beg); // New lenght for next recombombination event
  cleft+=1.0-pow(pc, len);
  len=(pseg2+newsg-1)->end-pseg2->beg;
  cleft+=1.0-pow(pc, len);
  if(!(chrom[ic].pseg=(struct seg*)realloc(chrom[ic].pseg,(unsigned)(lsg*sizeof(struct seg)))))
    perror(" realloc error. re2");
  if(in)
    {
      begs=pseg2->beg;               // begs=begining new segment after cut
      for(i=0, k=0;(k<nsegs-1)&&(begs>seglst[seglst[i].next].beg-1);
          i=seglst[i].next, k++);    // Find ith segment to give tree to new segment
      if(begs!=seglst[i].beg)
        {
          /* new tree */
          if(nsegs>=seglimit)
            {
              seglimit+=SEGINC;
              nnodes=(int*)realloc(nnodes,(unsigned)(sizeof(int)*seglimit));
              if(nnodes==NULL) perror("realloc error. 1. segtre_mig.c");
              seglst=(struct segl*)realloc(seglst,(unsigned)(sizeof(struct segl)*seglimit));
              if(seglst==NULL) perror("realloc error. 2. segtre_mig.c");
              /* printf("seglimit: %d\n", seglimit); */
            }
          seglst[nsegs].next=seglst[i].next;
          seglst[i].next=nsegs;      // add a segment
          seglst[nsegs].beg=begs;    // next segment begings after cut
          if(!(seglst[nsegs].ptree=
               (struct node*) calloc((unsigned)(2*nsam), sizeof(struct node)))) perror("calloc error. re3.");
          nnodes[nsegs]=nnodes[i];    // nodes[nseg] takes pinter on node i
          ptree1=seglst[i].ptree;     // tree of segment i
          ptree2=seglst[nsegs].ptree; // tree of new seg
          nsegs++;                    // add a segment
          for(k=0; k<=nnodes[i]; k++) 
            {
              (ptree2+k)->abv=(ptree1+k)->abv; // new tree like ith tree
              (ptree2+k)->time=(ptree1+k)->time;
              (ptree2+k)->lpop1=(ptree1+k)->lpop1;
              (ptree2+k)->lpop2=(ptree1+k)->lpop2;
            }// End for k # nodes ofith segment
        }// End if begs!=beg of segment #[i]
    }// End if spot in a segment
  return(ic);
}// End xover

/*----------------------------*
 *  Function coalsecent event |
 *----------------------------*
 | Pick two chromosomes and merge them. Update trees if necessary. 
*/
int ca(int nsam, int nsites, int c1, int c2)
{
  int yes1=0, yes2=0, seg1=0, seg2=0, seg=0, tseg=0, start=0, end=0, desc=0, k=0;
  int links(int), isseg(int, int, int*);
  struct seg *pseg=NULL;
  struct node *ptree=NULL;
  seg1=0;
  seg2=0;

  if(!(pseg=(struct seg*)calloc((unsigned)nsegs, sizeof(struct seg)))) 
    perror("alloc error.ca1");

  tseg=-1;
  //--- loop along all segments ---//
  for(seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++)
    {
      start=seglst[seg].beg;        // Starting point on chromosom 
      yes1=isseg(start, c1, &seg1); // chek if segment on c1
      yes2=isseg(start, c2, &seg2); // check is segment on c2
      
      if(yes1||yes2)                // if seg on one of chromosome: do coalescent
        {
          tseg++;                           // one more segment to consider
          (pseg+tseg)->beg=seglst[seg].beg; // new begining
          end=(k<nsegs-1 ? seglst[seglst[seg].next].beg-1 : nsites-1);// define end of this segment
          (pseg+tseg)->end=end;
          
          if(yes1 && yes2)                  // if segments on both chromosomes 
            {
              nnodes[seg]++;                      // 1 more nodes in tree for this seg
              if(nnodes[seg]>=(2*nsam-2)) tseg--; // if end of tree, stop
              else                                // if ancestral node
                (pseg+tseg)->desc=nnodes[seg];    //desc of this seg=# node for seg
              ptree=seglst[seg].ptree;
            
              desc=(chrom[c1].pseg+seg1)->desc;
              (ptree+desc)->abv=nnodes[seg];      // New ancestor node for seg1 of the chromosome c1
              (ptree+nnodes[seg])->lpop1+=(ptree+desc)->lpop1;
              (ptree+nnodes[seg])->lpop2+=(ptree+desc)->lpop2;
            
              desc=(chrom[c2].pseg+seg2)->desc;
              (ptree+desc)->abv=nnodes[seg];      // Same ancestor node for seg1 of the chromosome c2
              (ptree+nnodes[seg])->lpop1+=(ptree+desc)->lpop1;
              (ptree+nnodes[seg])->lpop2+=(ptree+desc)->lpop2; 
              (ptree+nnodes[seg])->time=t;        // new time of node for this segment
            }
          else
            {
              (pseg+tseg)->desc=(yes1 ?
                                 (chrom[c1].pseg+seg1)->desc :
                                 (chrom[c2].pseg+seg2)->desc); // desc on c1 if on c1 or c2 if on c2
            }
        }
    }
  nlinks-=links(c1);                // # links between begining and end for c1
  cleft-=1.0-pow(pc,(double)links(c1));
  free(chrom[c1].pseg);
  if(tseg<0)                        // if no seg in common 
    {
      free(pseg);                   // do coalecent on whole trunck of chromo
      chrom[c1].pseg=chrom[nchrom-1].pseg;
      chrom[c1].nseg=chrom[nchrom-1].nseg;
      chrom[c1].pop=chrom[nchrom-1].pop;
      if(c2==nchrom-1) c2=c1;
      nchrom--;                     //***** nChrom-- ****//
    }
  else
    {
      if(!(pseg=(struct seg*)realloc(pseg,(unsigned)((tseg+1)*sizeof(struct seg)))))
        perror(" realloc error. ca1"); // realloc for next segment
      chrom[c1].pseg=pseg;             // coalecent on this segment
      chrom[c1].nseg=tseg+1;
      nlinks+=links(c1);
      cleft+=1.0-pow(pc,(double)links(c1));
    }
  nlinks-=links(c2);
  cleft-=1.0-pow(pc,(double)links(c2));
  free(chrom[c2].pseg);
  chrom[c2].pseg=chrom[nchrom-1].pseg; // c2 takes nchrom-1 segments
  chrom[c2].nseg=chrom[nchrom-1].nseg;
  chrom[c2].pop=chrom[nchrom-1].pop;
  nchrom--;                         //***** nChrom-- ****// 
  if(tseg<0) return(2);             /* decrease of nchrom is two */
  else return(1);
}// End ca=common ancestor

/*-----------------------------------------*
 *  Function: does is segment starts here? |
 *-----------------------------------------*
 |  Does chromosome c contain the segment on seglst which starts at
 |  start? *psg is the segment of chrom[c] at which one is to begin
 |  looking.
*/
int isseg(int start, int c, int *psg)
{
  int ns=0;
  struct seg *pseg=NULL;
  ns=chrom[c].nseg;    // nseg in chromo c
  pseg=chrom[c].pseg;  // pseg on chromo c

  /* changed order of test conditions in following line on 6 Dec 2004 */
  for(;((*psg)<ns) &&((pseg+(*psg))->beg<=start); ++(*psg))
    if((pseg+(*psg))->end>=start) return(1); // make sure end is after start
  return(0);                                 // return 0 if seg not in chromo c
}

void pick2_chrom(int pop, int config[], int *pc1, int *pc2)
{
  int c1=0, c2=0, cs=0, cb=0, i=0, count=0;
  int pick2(int, int*, int*);
  
  pick2(config[pop], &c1, &c2); // pick 2 chsomosomes for coalscent in the sample
  cs=(c1>c2) ? c2 : c1;         // cs=min c1 & c2
  cb=(c1>c2) ? c1 : c2;         // cb=max c1 & c2

  i=count=0;
  for(;;)
    {
      while(chrom[i].pop!=pop) i++; // i chromo in pop 
      if(count==cs) break;
      count++;                      // go to next chromo
      i++;                          // go to next pop
    }
  *pc1=i;                           // pinter on chromo to coalesce
  i++;
  count++;
  for(;;)                           // Find next chromosome to coalsece
    {
      while(chrom[i].pop!=pop) i++;
      if(count==cb) break;
      count++;
      i++;
    }
  *pc2=i;
}


/**** links(c): returns the number of links between beginning and end of chrom **/
int links(int c)
{
  int ns=0;
  ns=chrom[c].nseg-1;
  return((chrom[c].pseg+ns)->end-(chrom[c].pseg)->beg);
}
