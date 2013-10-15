/*! \file streetR.c
  \brief Functions to generated ARGs 

  **********  segtre_mig.c **********************************
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
  *	in tracing back the history of the gametes.  The histories of
  *	these segments are passed back to the calling function in the
  *	array of structures seglst[]. An element of this array,  seglst[i],
  * 	consists of three parts: (1) beg, the starting point of
  *	of segment i, (2) ptree, which points to the first node of the
  *	tree representing the history of the segment, (3) next, which
  *	is the index number of the next segment.
  *	     A tree is a contiguous set of 2*nsam nodes. The first nsam
  *	nodes are the tips of the tree, the sampled gametes.  The other
  *	nodes are the nodes ancestral to the sampled gametes. Each node
  *	consists of an "abv" which is the number of the node ancestral to
  *	that node, an "ndes", which is not used or assigned to in this routine,
  *	and a "time", which is the time (in units of 4N generations) of the
  *	node. For the tips, time equals zero.
  *	Returns a pointer to an array of segments, seglst.

  **************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ms.h"
#include "streec.h"
#include "streeconly.h"
#include "rand1.h"
#include <R.h>
#include <Rdefines.h>  


struct segl *
segtre_mig(struct c_params *cp, int *pnsegs ) 
{
  int i, j, k,  dec, pop, pop2, c1, c2, ind, rchrom ;
  int migrant, source_pop, *config, flagint ;
  double   sum, x,  ttemp, rft, clefta,  tmin, p  ;
  double prec, cin,  prect,  mig, ran, coal_prob, rdum , arg ;
  char  event ;
  int  eflag, cpop, ic  ;
  int nsam, npop, nsites, *inconfig=NULL ;
  double r,  f, rf,  track_len,    **migm=NULL ;
  double *size=NULL, *alphag=NULL, *tlast=NULL ;
  struct devent *nextevent=NULL ;

  cpop=0;
  tmin=0;
  event=' ';
  nsam = cp->nsam;
  npop = cp->npop;
  nsites = cp->nsites;
  inconfig = cp->config;
  r = cp->r ;
  f = cp->f ;
  track_len = cp->track_len ;
  migm = (double **)calloc( (unsigned)npop,sizeof(double *) ) ;
  for( i=0; i<npop; i++) {
    migm[i] = (double *)calloc( (unsigned)npop,sizeof( double) ) ;
    for( j=0; j<npop; j++) migm[i][j] = (cp->mig_mat)[i][j] ;
  }
  nextevent = cp->deventlist ;
	
  /* Initialization */
  if( chrom == NULL ) {
    maxchr = nsam + 20 ;
    chrom = (struct chromo *)calloc( (unsigned) maxchr,sizeof( struct chromo)) ;
    if( chrom == NULL ) perror( "malloc error. segtre");
  }
  if( nnodes == NULL ){
    nnodes = (int*) calloc((unsigned)seglimit,sizeof(int))  ;
    if( nnodes == NULL ) perror("malloc error. segtre_mig");
  }
  if( seglst == NULL ) {
    seglst = (struct segl *)calloc((unsigned)seglimit,sizeof(struct segl)) ;
    if( seglst == NULL ) perror("malloc error. segtre_mig.c 2");
  }

  config = (int *)calloc( (unsigned) (npop+1),sizeof(int) ) ;
  if( config == NULL ) perror("malloc error. segtre.");
  size = (double *)calloc( (unsigned) (npop),sizeof(double) ) ;
  if( size == NULL ) perror("malloc error. segtre.");
  alphag = (double *)calloc( (unsigned) (npop),sizeof(double) ) ;
  if( alphag == NULL ) perror("malloc error. segtre.");
  tlast = (double *)calloc( (unsigned) (npop),sizeof(double) ) ;
  if( alphag == NULL ) perror("malloc error. segtre.");
  for(pop=0;pop<npop;pop++) {
    config[pop] = inconfig[pop] ;
    size[pop] = (cp->size)[pop] ;
    alphag[pop] = (cp->alphag)[pop] ;
    tlast[pop] = 0.0 ;
  }
  seglst[0].beg = 0;
  if( !(seglst[0].ptree = (struct node *)calloc((unsigned)(2*nsam),sizeof(struct node)) ))
    perror("calloc error. se2");
  for(pop=ind=0;pop<npop;pop++)
    for(j=0; j<inconfig[pop];j++,ind++) {
			
      chrom[ind].nseg = 1;
      if( !(chrom[ind].pseg = (struct seg*)calloc((unsigned)1,sizeof(struct seg)) ))
        error("calloc error. se1");
      //ERROR("calloc error. se1");

      (chrom[ind].pseg)->beg = 0;
      (chrom[ind].pseg)->end = nsites-1;
      (chrom[ind].pseg)->desc = ind ;
      chrom[ind].pop = pop ;
      if(pop==0)                            // record lineage of the leaf nodes
        seglst[0].ptree[ind].lpop1+=1;
      else
        seglst[0].ptree[ind].lpop2+=1;
    }
 

  nnodes[0] = nsam - 1 ;
  nchrom=nsam;
  nlinks = ((long)(nsam))*(nsites-1) ;
  nsegs=1;
  t = 0.;
  r /= (nsites-1);
  if( f > 0.0 ) 	pc = (track_len -1.0)/track_len ;
  else pc = 1.0 ;
  lnpc = log( pc ) ;
  cleft = nsam* ( 1.0 - pow( pc, (double)(nsites-1) ) ) ;
  if( r > 0.0 ) rf = r*f ;
  else rf = f /(nsites-1) ;
  rft = rf*track_len ;
  flagint = 0 ;

  /* Main loop */

  while( nchrom > 1 ) {
    prec = nlinks*r;
    cin = nlinks*rf ;
    clefta = cleft*rft ;
    prect = prec + cin + clefta ;
    mig = 0.0;
    for( i=0; i<npop; i++) mig += config[i]*migm[i][i] ;
    if( (npop > 1) && ( mig == 0.0) && ( nextevent == NULL)) {
      i = 0;
      for( j=0; j<npop; j++) 
        if( config[j] > 0 ) i++;
      if( i > 1 ) {
        fprintf(stderr," Infinite coalescent time. No migration.\n");
        error("I didn't finish running ms");
        // exit(1);
      }
    }
    eflag = 0 ;

    if( prect > 0.0 ) {      /* cross-over or gene conversion */
      while( (rdum = ran1() )  == 0.0 ) ;
      ttemp = -log( rdum)/prect ;
      if( (eflag == 0) || (ttemp < tmin ) ){
        tmin = ttemp;
        event = 'r' ;
        eflag = 1;
      }
    }
    if(mig > 0.0 ) {         /* migration   */
      while( (rdum = ran1() ) == 0.0 ) ;
      ttemp = -log( rdum)/mig ;
      if( (eflag == 0) || (ttemp < tmin ) ){
        tmin = ttemp;
        event = 'm' ;
        eflag = 1 ;
      }
    }
		
    for(pop=0; pop<npop ; pop++) {     /* coalescent */
		coal_prob = ((double)config[pop])*(config[pop]-1.) ;
      if( coal_prob > 0.0 ) {
        while( ( rdum = ran1() )  == .0 ) ;
        if( alphag[pop] == 0 ){
			 ttemp = -log( rdum )*size[pop] /coal_prob ;
          if( (eflag == 0) || (ttemp < tmin ) ){
            tmin = ttemp;
            event = 'c' ;
            eflag = 1 ;
            cpop = pop;
			 }
        }
        else {
          arg  = 1. - alphag[pop]*size[pop]*exp(-alphag[pop]*(t - tlast[pop] ) )* log(rdum) / coal_prob     ;
          if( arg > 0.0 ) {                          /*if arg <= 0,  no coalescent within interval */ 
            ttemp = log( arg ) / alphag[pop]  ;
            if( (eflag == 0) || (ttemp < tmin ) ){
              tmin = ttemp;
              event = 'c' ;
              eflag = 1 ;
              cpop = pop ;
            }
          }
        }
      }		
    }

    if( (eflag == 0) && ( nextevent == NULL) ) {
      fprintf(stderr,
              " infinite time to next event. Negative growth rate in last time interval or non-communicating subpops.\n");
      error("I didn't finish running ms");
      // exit(0);
    }
    if( ( ( eflag == 0) && (nextevent != NULL))|| ( (nextevent != NULL) &&  ( (t+tmin) >=  nextevent->time)) ) {
      t = nextevent->time ;
      switch(  nextevent->detype ) {
        case 'N' :
          for(pop =0; pop <npop; pop++){
            size[pop]= nextevent->paramv ;
            alphag[pop] = 0.0 ;
          }
          nextevent = nextevent->nextde ;
          break;
        case 'n' :
          size[nextevent->popi]= nextevent->paramv ;
          alphag[nextevent->popi] = 0.0 ;
          nextevent = nextevent->nextde ;
          break;
        case 'G' :
          for(pop =0; pop <npop; pop++){
            size[pop] = size[pop]*exp( -alphag[pop]*(t - tlast[pop]) ) ;
            alphag[pop]= nextevent->paramv ;
            tlast[pop] = t ;
          }
          nextevent = nextevent->nextde ;
          break;
        case 'g' :
          pop = nextevent->popi ;
          size[pop] = size[pop]*exp( - alphag[pop]*(t-tlast[pop]) ) ;
          alphag[pop]= nextevent->paramv ;
          tlast[pop] = t ;
          nextevent = nextevent->nextde ;
          break;
        case 'M' :
          for(pop =0; pop <npop; pop++)
            for( pop2 = 0; pop2 <npop; pop2++) migm[pop][pop2] = (nextevent->paramv) /(npop-1.0) ;
          for( pop = 0; pop <npop; pop++)
            migm[pop][pop]= nextevent->paramv ;
          nextevent = nextevent->nextde ;
          break;
        case 'a' :
          for(pop =0; pop <npop; pop++)
            for( pop2 = 0; pop2 <npop; pop2++) migm[pop][pop2] = (nextevent->mat)[pop][pop2]  ;
          nextevent = nextevent->nextde ;
          break;
        case 'm' :
          i = nextevent->popi ;
          j = nextevent->popj ;
          migm[i][i] += nextevent->paramv - migm[i][j];
          migm[i][j]= nextevent->paramv ;
          nextevent = nextevent->nextde ;
          break;
        case 'j' :         /* merge pop i into pop j  (join) */
          i = nextevent->popi ;
          j = nextevent->popj ;
          config[j] += config[i] ;
          config[i] = 0 ;
          for( ic = 0; ic<nchrom; ic++) if( chrom[ic].pop == i ) chrom[ic].pop = j ;
          /*  the following was added 19 May 2007 */
          for( k=0; k < npop; k++){
            if( k != i) {
		        migm[k][k] -= migm[k][i] ;
		        migm[k][i] = 0. ;
            }
          }
          /* end addition */
          nextevent = nextevent->nextde ;
          break;
        case 's' :         /*split  pop i into two;p is the proportion from pop i, and 1-p from pop n+1  */
          i = nextevent->popi ;
          p = nextevent->paramv ;
          npop++;
          if(!(config = (int *)realloc( config, (unsigned)(npop*sizeof( int) ))))
            perror( " realloc error. conf");
          if(!(size = (double *)realloc(size, (unsigned)(npop*sizeof(double) ))))
            perror( " realloc error. size");
          if(!(alphag = (double *)realloc(alphag, (unsigned)(npop*sizeof(double) ))))
            perror( " realloc error. alpha");
          if(!(tlast = (double *)realloc(tlast,(unsigned)(npop*sizeof(double) ) ) ))
            perror( " realloc error. conf");
          tlast[npop-1] = t ;
          size[npop-1] = 1.0 ;
          alphag[npop-1] = 0.0 ;
          migm = (double **)realloc(migm, (unsigned)(npop*sizeof( double *)));
          for( j=0; j< npop-1; j++)
            migm[j] = (double *)realloc(migm[j],(unsigned)(npop*sizeof(double)));
          migm[npop-1] = (double *)calloc( (unsigned)npop,sizeof( double) ) ;
          for( j=0; j<npop; j++) migm[npop-1][j] = migm[j][npop-1] = 0.0 ;
          config[npop-1] = 0 ;
          config[i] = 0 ;
          for( ic = 0; ic<nchrom; ic++){
            if( chrom[ic].pop == i ) {
              if( ran1() < p ) config[i]++;
              else {
                chrom[ic].pop = npop-1 ;
                config[npop-1]++;
              }
            }
          }
          nextevent = nextevent->nextde ;
          break;
		}
    } 
    else {
      t += tmin ;	
      if( event == 'r' ) {   
        if( (ran = ran1()) < ( prec / prect ) ){ /*recombination*/
          rchrom = re(nsam);
          config[ chrom[rchrom].pop ] += 1 ;
        }
        else if( ran < (prec + clefta)/(prect) ){    /*  cleft event */
			 rchrom = cleftr(nsam);
			 config[ chrom[rchrom].pop ] += 1 ;
        }
        else  {         /* cin event */
			 rchrom = cinr(nsam,nsites);
			 if( rchrom >= 0 ) config[ chrom[rchrom].pop ] += 1 ;
        }
      }
      else if ( event == 'm' ) {  /* migration event */
        x = mig*ran1();
        sum = 0.0 ;
        for( i=0; i<nchrom; i++) {
          sum += migm[chrom[i].pop][chrom[i].pop] ;
          if( x <sum ) break;
        }
        migrant = i ;
        x = ran1()*migm[chrom[i].pop][chrom[i].pop];
        sum = 0.0;
        for(i=0; i<npop; i++){
          if( i != chrom[migrant].pop ){
            sum += migm[chrom[migrant].pop][i];
            if( x < sum ) break;
          }
        }
        source_pop = i;
        config[chrom[migrant].pop] -= 1;
        config[source_pop] += 1;
        chrom[migrant].pop = source_pop ;
      }
      else { 								 /* coalescent event */
        /* pick the two, c1, c2  */
        pick2_chrom( cpop, config, &c1,&c2);  /* c1 and c2 are chrom's to coalesce */
        dec = ca(nsam,nsites,c1,c2 );
        config[cpop] -= dec ;
      }
    }
  }  
  *pnsegs = nsegs ;
  free(config); 
  free( size ) ;
  free( alphag );
  free( tlast );
  for( i=0; i<npop; i++) free ( migm[i] ) ;
  free( migm ) ;
  return( seglst );
}

/******  recombination subroutine ***************************
         Picks a chromosome and splits it in two parts. If the x-over point
         is in a new spot, a new segment is added to seglst and a tree set up
         for it.   ****/


int
re( int nsam)
{
  struct seg *pseg ;
  int  el, lsg, lsgm1,  ic,  is,  spot;;

  /* First generate a random x-over spot, then locate it as to chrom and seg. */
  pseg=NULL;
  spot = nlinks*ran1() + 1.;

  /* get chromosome # (ic)  */
  for( ic=0; ic<nchrom ; ic++) {
    lsg = chrom[ic].nseg ;
    lsgm1 = lsg - 1;
    pseg = chrom[ic].pseg;
    el = ( (pseg+lsgm1)->end ) - (pseg->beg);
    if( spot <= el ) break;
    spot -= el ;
  }
  is = pseg->beg + spot -1;
  xover(nsam, ic, is);
  return(ic);
}

int
cleftr( int nsam)
{
  struct seg *pseg ;
  int    ic,  is;
  double  x, sum, len  ;
  pseg=NULL;
  while( (x = cleft*ran1() )== 0.0 )  ;
  sum = 0.0 ;
  ic = -1 ;
  while ( sum < x ) {
    sum +=  1.0 - pow( pc, links(++ic) )  ;
  }
  pseg = chrom[ic].pseg;
  len = links(ic) ;
  is = pseg->beg + floor( 1.0 + log( 1.0 - (1.0- pow( pc, len))*ran1() )/lnpc  ) -1  ;
  xover( nsam, ic, is);
  return( ic) ;
}

int
cinr( int nsam, int nsites)
{
  struct seg *pseg ;
  int len,  el, lsg, lsgm1,  ic,  is, spot, endic ;
  
  lsgm1=0;
  pseg=NULL;
  /* First generate a random x-over spot, then locate it as to chrom and seg. */

  spot = nlinks*ran1() + 1.;

  /* get chromosome # (ic)  */
  for( ic=0; ic<nchrom ; ic++) {
    lsg = chrom[ic].nseg ;
    lsgm1 = lsg - 1;
    pseg = chrom[ic].pseg;
    el = ( (pseg+lsgm1)->end ) - (pseg->beg);
    if( spot <= el ) break;
    spot -= el ;
  }
  is = pseg->beg + spot -1;
  endic = (pseg+lsgm1)->end ;
  xover(nsam, ic, is);

  len = floor( 1.0 + log( ran1() )/lnpc ) ;
  if( is+len >= endic ) return(ic) ;
  if( is+len < (chrom[nchrom-1].pseg)->beg ){
    ca( nsam, nsites, ic, nchrom-1);
    return(-1) ;
  }
  xover( nsam, nchrom-1, is+len ) ;
  ca( nsam,nsites, ic,  nchrom-1);
  return(ic);

}

int
xover(int nsam,int ic, int is)
{
  struct seg *pseg, *pseg2;
  int i,  lsg, lsgm1, newsg,  jseg, k,  in;
  double len ;

  pseg = chrom[ic].pseg ;
  lsg = chrom[ic].nseg ;
  len = (pseg + lsg -1)->end - pseg->beg ;
  cleft -= 1 - pow(pc,len) ;
  /* get seg # (jseg)  */

  for( jseg=0; is >= (pseg+jseg)->end ; jseg++) ;
  if( is >= (pseg+jseg)->beg ) in=1;
  else in=0;
  newsg = lsg - jseg ;

  /* copy last part of chrom to nchrom  */
  nchrom++;
  if( nchrom >= maxchr ) {
    maxchr += 20 ;
    chrom = (struct chromo *)realloc( chrom, (unsigned)(maxchr*sizeof(struct chromo))) ;
    if( chrom == NULL ) perror( "malloc error. segtre2");
  }
  if( !( pseg2 = chrom[nchrom-1].pseg = (struct seg *)calloc((unsigned)newsg,sizeof(struct seg)) ) )
    error(" alloc error. re1");
  //ERROR(" alloc error. re1");
  chrom[nchrom-1].nseg = newsg;
  chrom[nchrom-1].pop = chrom[ic].pop ;
  pseg2->end = (pseg+jseg)->end ;
  if( in ) {
    pseg2->beg = is + 1 ;
    (pseg+jseg)->end = is;
  }
  else pseg2->beg = (pseg+jseg)->beg ;
  pseg2->desc = (pseg+jseg)->desc ;
  for( k=1; k < newsg; k++ ) {
    (pseg2+k)->beg = (pseg+jseg+k)->beg;
    (pseg2+k)->end = (pseg+jseg+k)->end;
    (pseg2+k)->desc = (pseg+jseg+k)->desc;
  }

  lsg = chrom[ic].nseg = lsg-newsg + in ;
  lsgm1 = lsg - 1 ;
  nlinks -= pseg2->beg - (pseg+lsgm1)->end ;
  len = (pseg+lsgm1)->end - (pseg->beg) ;
  cleft += 1.0 - pow( pc, len) ;
  len = (pseg2 + newsg-1)->end - pseg2->beg ;
  cleft += 1.0 - pow(pc, len) ;
  if( !(chrom[ic].pseg = 
        (struct seg *)realloc(chrom[ic].pseg,(unsigned)(lsg*sizeof(struct seg)) )) )
    perror( " realloc error. re2");
  if( in ) {
    begs = pseg2->beg;
    for( i=0,k=0; (k<nsegs-1)&&(begs > seglst[seglst[i].next].beg-1);
		   i=seglst[i].next, k++) ;
    if( begs != seglst[i].beg ) {
      /* new tree  */

      if( nsegs >= seglimit ) {  
        seglimit += SEGINC ;
        nnodes = (int *)realloc( nnodes,(unsigned)(sizeof(int)*seglimit)) ; 
        if( nnodes == NULL) perror("realloc error. 1. segtre_mig.c");
        seglst =
          (struct segl *)realloc( seglst,(unsigned)(sizeof(struct segl)*seglimit));
        if(seglst == NULL ) perror("realloc error. 2. segtre_mig.c");
      } 
      seglst[nsegs].next = seglst[i].next;
      seglst[i].next = nsegs;
      seglst[nsegs].beg = begs ;
      if( !(seglst[nsegs].ptree = (struct node *)calloc((unsigned)(2*nsam), sizeof(struct
                                                                                   node)) )) perror("calloc error. re3.");
      nnodes[nsegs] = nnodes[i];
      ptree1 = seglst[i].ptree;
      ptree2 = seglst[nsegs].ptree;
      nsegs++ ;
      for( k=0; k<=nnodes[i]; k++) {
        (ptree2+k)->abv = (ptree1+k)->abv ;
        (ptree2+k)->time = (ptree1+k)->time;
        (ptree2+k)->lpop1=(ptree1+k)->lpop1;
        (ptree2+k)->lpop2=(ptree1+k)->lpop2;
      }
    }
  }
  return(ic) ;
}

/***** common ancestor subroutine **********************
       Pick two chromosomes and merge them. Update trees if necessary. **/

int
ca(int nsam, int nsites,int c1,int c2)
{
  int yes1, yes2, seg1, seg2, seg ;
  int tseg, start, end, desc, k;
  struct seg *pseg;
  struct node *ptree;
   
  seg1=0;
  seg2=0;

  if( !(pseg = (struct seg *)calloc((unsigned)nsegs,sizeof(struct seg) )))
    perror("alloc error.ca1");

  tseg = -1 ;

  for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
    start = seglst[seg].beg;
    yes1 = isseg(start, c1, &seg1);
    yes2 = isseg(start, c2, &seg2);
    if( yes1 || yes2 ) {
      tseg++;
      (pseg+tseg)->beg=seglst[seg].beg;
      end = ( k< nsegs-1 ? seglst[seglst[seg].next].beg-1 : nsites-1 ) ;
      (pseg+tseg)->end = end ;

      if( yes1 && yes2 ) {
        nnodes[seg]++;
        if( nnodes[seg] >= (2*nsam-2) ) tseg--;
        else
          (pseg+tseg)->desc = nnodes[seg];
        ptree=seglst[seg].ptree;
        desc = (chrom[c1].pseg + seg1) ->desc;
        (ptree+desc)->abv = nnodes[seg];
        (ptree+nnodes[seg])->lpop1+=(ptree+desc)->lpop1;
        (ptree+nnodes[seg])->lpop2+=(ptree+desc)->lpop2;

        desc = (chrom[c2].pseg + seg2) -> desc;
        (ptree+desc)->abv = nnodes[seg];
        (ptree+nnodes[seg])->lpop1+=(ptree+desc)->lpop1;
        (ptree+nnodes[seg])->lpop2+=(ptree+desc)->lpop2; 
        (ptree+nnodes[seg])->time = t;

      }
      else {
        (pseg+tseg)->desc = ( yes1 ?
                              (chrom[c1].pseg + seg1)->desc :
                              (chrom[c2].pseg + seg2)->desc);
      }
    }
  }
  nlinks -= links(c1);
  cleft -= 1.0 - pow(pc, (double)links(c1));
  free(chrom[c1].pseg) ;
  if( tseg < 0 ) {
    free(pseg) ;
    chrom[c1].pseg = chrom[nchrom-1].pseg;
    chrom[c1].nseg = chrom[nchrom-1].nseg;
    chrom[c1].pop = chrom[nchrom-1].pop ;
    if( c2 == nchrom-1 ) c2 = c1;
    nchrom--;
  }
  else {
    if( !(pseg = (struct seg *)realloc(pseg,(unsigned)((tseg+1)*sizeof(struct seg)))))
      perror(" realloc error. ca1");
    chrom[c1].pseg = pseg;
    chrom[c1].nseg = tseg + 1 ;
    nlinks += links(c1);
    cleft += 1.0 - pow(pc, (double)links(c1));
  }
  nlinks -= links(c2);
  cleft -= 1.0 - pow(pc, (double)links(c2));
  free(chrom[c2].pseg) ;
  chrom[c2].pseg = chrom[nchrom-1].pseg;
  chrom[c2].nseg = chrom[nchrom-1].nseg;
  chrom[c2].pop = chrom[nchrom-1].pop ;
  nchrom--;
  if(tseg<0) return( 2 );  /* decrease of nchrom is two */
  else return( 1 ) ;
}

/*** Isseg: Does chromosome c contain the segment on seglst which starts at
     start? *psg is the segment of chrom[c] at which one is to begin
     looking.  **/

int
isseg(int start,int c,int * psg)
{
  int ns;
  struct seg *pseg;

  ns = chrom[c].nseg;
  pseg = chrom[c].pseg;

  /*  changed order of test conditions in following line on 6 Dec 2004 */
  for(  ; ((*psg) < ns ) && ( (pseg+(*psg))->beg <= start ) ; ++(*psg) )
    if( (pseg+(*psg))->end >= start ) return(1);
  return(0);
}
	


int
pick2_chrom(int pop,int config[],int *pc1,int *pc2)
{
  int c1, c2, cs,cb,i, count;

  pick2(config[pop],&c1,&c2);
  cs = (c1>c2) ? c2 : c1;
  cb = (c1>c2) ? c1 : c2 ;
  i=count=0;
  for(;;){
    while( chrom[i].pop != pop ) i++;
    if( count == cs ) break;
    count++;
    i++;
  }
  *pc1 = i;
  i++;
  count++;
  for(;;){
    while( chrom[i].pop != pop ) i++;
    if( count == cb ) break;
    count++;
    i++;
  }
  *pc2 = i ;
  return(0);
}
	
	

/****  links(c): returns the number of links between beginning and end of chrom **/

int
links( int c)
{
  int ns;

  ns = chrom[c].nseg - 1 ;

  return( (chrom[c].pseg + ns)->end - (chrom[c].pseg)->beg);
}


	
/************ make_gametes.c  *******************************************
 *
 *
 *****************************************************************************/

#define STATE1 '1'
#define STATE2 '0'

int
make_gametes(int nsam, int mfreq, struct node *ptree, double tt, int newsites, int ns, char **list )
{
  int  tip, j,  node ;
  char tp1[ns+newsites];
  if(mfreq==0 && ns+newsites==0)// init even if no site 
    for( tip=0; tip < nsam ; tip++) 
      strcpy(list[tip],"\0");
  for(  j=ns; j< ns+newsites ;  j++ ) 
    {  
      node = pickb(  nsam, ptree, tt);
      for( tip=0; tip < nsam ; tip++) 
        { 
          if(mfreq==0 && j==0)
            strcpy(list[tip],"\0");
          if( tdesn(ptree, tip, node) )
            sprintf(tp1,"%s%c",list[tip],STATE1);
          else
            sprintf(tp1,"%s%c", list[tip], STATE2);
          strcpy(list[tip],tp1);
        }
    }
  return(0);
}

/***  ttime.c : Returns the total time in the tree, *ptree, with nsam tips. **/

double
ttime( struct node *ptree,
       int nsam)
{
  double t;
  int i;

  t = (ptree + 2*nsam-2) -> time ;
  for( i=nsam; i< 2*nsam-1 ; i++)
    t += (ptree + i)-> time ;
  return(t);
}


double
ttimemf( struct node *ptree,
         int nsam, int mfreq)
{
  double t;
  int i;

  t = 0. ;
  for( i=0;  i< 2*nsam-2  ; i++)
    if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) )
		t += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
  return(t);
}


void
prtree(  struct node *ptree,
         int nsam)
{
  int i, *descl, *descr ;

  descl = (int *)calloc( (unsigned)(2*nsam-1),sizeof( int) );
  descr = (int *)calloc( (unsigned)(2*nsam-1),sizeof( int) );
  for( i=0; i<2*nsam-1; i++) descl[i] = descr[i] = -1 ;
  for( i = 0; i< 2*nsam-2; i++){
    if( descl[ (ptree+i)->abv ] == -1 ) descl[(ptree+i)->abv] = i ;
    else descr[ (ptree+i)->abv] = i ;
  }
  parens( ptree, descl, descr, 2*nsam-2);
  free( descl ) ;
  free( descr ) ;
}

void
parens( struct node *ptree, int *descl, int *descr,  int noden)
{
  double time ;

  if( descl[noden] == -1 ) {
    Rprintf("%d:%5.3lf", noden+1, (ptree+ ((ptree+noden)->abv))->time );
  }
  else{
    Rprintf("(");
    parens( ptree, descl,descr, descl[noden] ) ;
    Rprintf(",");
    parens(ptree, descl, descr, descr[noden] ) ;
    if( (ptree+noden)->abv == 0 ) Rprintf(");\n"); 
    else {
      time = (ptree + (ptree+noden)->abv )->time - (ptree+noden)->time ;
      Rprintf("):%5.3lf", time );
    }
  }
}

/***  pickb : returns a random branch from the tree. The probability of picking
      a particular branch is proportional to its duration. tt is total
      time in tree.   ****/

int pickb(int nsam, struct node *ptree, double tt)
{
  double x, y;
  int i;

  x = ran1()*tt;
  for( i=0, y=0; i < 2*nsam-2 ; i++) {
    y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
    if( y >= x ) return( i ) ;
  }
  return( i );
}

int  pickbmf( int nsam,int mfreq,struct node *ptree,double tt)
{
  double x, y;
  int i;

  x = ran1()*tt;
  for( i=0, y=0; i < 2*nsam-2 ; i++) {
    if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) )
		y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
    if( y >= x ) return( i ) ;
  }
  return( i );
}

/* tdesn : returns 1 if tip is a descendant of node in *ptree, otherwise 0. */

int
tdesn( struct node * ptree,int tip,int node )
{
  int k;

  for( k= tip ; k < node ; k = (ptree+k)->abv ) ;
  if( k==node ) return(1);
  else return(0);
}


/* pick2()  */

int pick2( int n,int *i,int *j)
{
  *i = n * ran1() ;
  while( ( *j = n * ran1() ) == *i )
    ;
  return(0) ;
}

/**** ordran.c  ***/

int ordran( int n,
            double pbuf[])
{
  ranvec(n,pbuf);
  order(n,pbuf);
  return(0);
}


int mnmial(int n,int nclass,double p[],int rv[])
{
  double x, s;
  int i, j;
  j=0;

  for(i=0; i<nclass; i++) rv[i]=0;
  for(i=0; i<n ; i++) {
    x = ran1();
    j=0;
    s = p[0];
    while( (x > s) && ( j<(nclass-1) ) )  s += p[++j];
    rv[j]++;
  }
  return(j);
}

int
order(int n,
      double pbuf[])
{
  int gap, i, j;
  double temp;

  for( gap= n/2; gap>0; gap /= 2)
    for( i=gap; i<n; i++)
      for( j=i-gap; j>=0 && pbuf[j]>pbuf[j+gap]; j -=gap) {
        temp = pbuf[j];
        pbuf[j] = pbuf[j+gap];
        pbuf[j+gap] = temp;
      }
  return(0);
}


int
ranvec(int n,
       double pbuf[])
{
  int i;

  for(i=0; i<n; i++)
    pbuf[i] = ran1();

  return(0);
}



int
poisso(double u)
{
  double  cump, ru, ran1(), p, gasdev() ;
  int i=1;

  if( u > 30. ) return( (int)(0.5 + gasdev(u,u)) );
  ru = ran1();
  p = exp(-u);
  if( ru < p) return(0);
  cump = p;
	
  while( ru > ( cump += (p *= u/i ) ) )
    i++;
  return(i);
}


/* a slight modification of crecipes version */

double gasdev( double m,double v)
{
  static int iset=0;
  static float gset;
  float fac,r,v1,v2;

  if  (iset == 0) {
    do {
      v1=2.0*ran1()-1.0;
      v2=2.0*ran1()-1.0;
      r=v1*v1+v2*v2;
    } while (r >= 1.0);
    fac=sqrt(-2.0*log(r)/r);
    gset= v1*fac;
    iset=1;
    return( m + sqrt(v)*v2*fac);
  } else {
    iset=0;
    return( m + sqrt(v)*gset ) ;
  }
}

