
/************************************************************/

/*! \file ms C functions  */
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "ms.h"

void 
ndes_setup(struct node *ptree, int nsam )
{
  int i ;

  for( i=0; i<nsam; i++) (ptree+i)->ndes = 1 ;
  for( i = nsam; i< 2*nsam -1; i++) (ptree+i)->ndes = 0 ;
  for( i= 0; i< 2*nsam -2 ; i++)  (ptree+((ptree+i)->abv))->ndes += (ptree+i)->ndes ;
}

int
locate(n,beg,len,ptr)
     int n;
     double beg, len, *ptr;
{
  int ordran(	int,	double *);
  int i;

  ordran(n,ptr);
  for(i=0; i<n; i++)
    ptr[i] = beg + ptr[i]*len ;
  return(0);
}

void
addtoelist( struct devent *pt, struct devent *elist ) 
{
  struct devent *plast, *pevent, *ptemp  ;
  plast=NULL;
  pevent = elist ;
  while(  (pevent != NULL ) && ( pevent->time <= pt->time ) )  {
    plast = pevent ;
    pevent = pevent->nextde ;
  }

  ptemp = plast->nextde ;
  plast->nextde = pt ;
  pt->nextde = ptemp ;
}
                        
void 
free_eventlist( struct devent *pt, int npop )
{
  struct devent *next ;
  int pop ;

  while( pt != NULL){
    next = pt->nextde ;
    if( pt->detype == 'a' ) {
      for( pop = 0; pop < npop; pop++) {free( (pt->mat)[pop] ); printf("\t\t\t\tmat[i]\t\t1\n");}
      free( pt->mat );
    }
    free(pt);
    pt=NULL;
    pt = next ;
  } 
  pt=NULL;
}


