/*! \file stats_popC.c
  \brief Code for summary statistics calculations


*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rdefines.h> 
#include "statsfun.h"
#include "stats_popR.h"

/*! \brief Main function for statsPop call from R*/
SEXP statsPop(SEXP cmdline,SEXP seg,SEXP pos,SEXP haplo,SEXP type) 
{
  /*** C object and functions ***/
  int nsam, j , i,  howmany, segsites,count, length=9999   ; 
  char **list;
  int *posit;
  double *p_pos;
  /**** Protect R obecjt input *******/
  int   *p_seg, popsize[2] , npop ;
  double  *p_type,*typestats;
  char *line[1],   *dum, *tpch;
   typestats = (double *)calloc((unsigned) 4,sizeof(double));

  PROTECT(cmdline = AS_CHARACTER(cmdline));   /* get command line info*/
  line[0] = R_alloc(strlen(CHAR(STRING_ELT(cmdline, 0))), sizeof(char)); 
  strcpy(line[0], CHAR(STRING_ELT(cmdline, 0))); 
  sscanf(line[0]," %*s %d %d",  &nsam, &howmany);
  dum=(char*)calloc(100 , sizeof(char));
  tpch=(char*)calloc(100 ,sizeof(char));
  popsize[1]=i=0;                        /* obtain the popsizes */
  popsize[0]=nsam;

  tpch=NULL;
  tpch = strtok (line[0]," ");
  while(strcmp(tpch,"\n")!=0)
    {
      if(strcmp(tpch,"-I")==0)
        {
          tpch = strtok (NULL, " ");
          j=atoi(tpch);
          if(j==2)
            for(i=0;i<j;i++)
              {
                tpch = strtok (NULL, " ");
                sscanf(tpch,"%d",&popsize[i]);
              }
        }
      if(strcmp(tpch,"-r")==0)
        {
          tpch = strtok (NULL, " ");// rho
          tpch = strtok (NULL, " ");// Z-1
          length=atoi(tpch);
        }

      tpch = strtok (NULL, " ");
      if(tpch== NULL) break;
    }

  PROTECT(seg = AS_INTEGER(seg));         /* array of seg sites*/
  p_seg= INTEGER_POINTER(seg);

  PROTECT(haplo= AS_VECTOR(haplo));       /* vector of lists of haplotypes */

  PROTECT(pos= AS_VECTOR(pos));       /* vector of lists of haplotypes */

  PROTECT(type = AS_NUMERIC(type));  /* array of type*/
   p_type=NUMERIC_POINTER(type);
  /********/

  /* OUTPUT R object and their pointer */
  SEXP stat_list,v_type,v_SS,v_FS,v_FST, v_PI, v_TW, v_D, v_TH, v_H, v_SFL,v_D_star,v_D_prime, v_r2,v_nH,v_Rm, list_names ;
  double *p_FS, *p_FST, *p_PI, *p_TW, *p_D,*p_TH, *p_H ,*p_D_star, *p_D_prime, *p_r2;
  int *p_SS, *p_SFL ,*p_nH,*p_Rm ;
  PROTECT(stat_list =  allocVector(VECSXP, 15)); 
  char *names[15] = {"type","segsites","singletons","frequencies","fst","pi","thetaW", "thetaH",  "D", "H", "D_star","D_prime","r_square","nH","Rm"};// name for various object

  PROTECT(list_names = allocVector(STRSXP,15));
  for(i = 0; i < 15; i++)
    SET_STRING_ELT(list_names,i,mkChar(names[i]));

  PROTECT( v_type= allocVector(STRSXP,3));/* type of statistics*/
  for(i=0;i<4;i++)
    {
      typestats[i]=p_type[i];
      if(i==0)
        {
          strcpy(dum,"Known ancestral allele");
          if( typestats[i]==1)
            strcpy(dum,"Unknown ancestral allele");
          else typestats[i]=0;
          SET_STRING_ELT(v_type,i,mkChar(dum));
        }
      if(i==1)
        {
          strcpy(dum,"Phased data");
          if( typestats[i]==2)
            strcpy(dum,"Unphased data");
          else typestats[i]=1;
          SET_STRING_ELT(v_type,i,mkChar(dum));
        }
      if(i==2)
      {
        sprintf(dum,"Filter for pairwise LD between %lf",typestats[i] );
      }
      if(i==3)
        {
          sprintf(dum,"%s %lf KB",dum,typestats[i] );
          SET_STRING_ELT(v_type,i-1,mkChar(dum));
        }
    }
  free(dum);
  free(tpch);
  PROTECT(v_SS= allocMatrix(INTSXP,howmany,12));/* list of vectors of seg sites, Ss, S1 S2 ss Sf*/
  p_SS = INTEGER_POINTER(v_SS);
  
  PROTECT(v_SFL= allocMatrix(INTSXP,howmany,3));/* list of vectors of num of singletons */
  p_SFL =  INTEGER_POINTER(v_SFL);

  PROTECT( v_FS= allocMatrix(REALSXP,howmany,9));/* list of vectors of mean frequencies of seg sites, Ss, S1 S2 ss Sf*/
  p_FS = NUMERIC_POINTER(v_FS);

  PROTECT( v_FST= allocMatrix(REALSXP,howmany,5));/* list of vectors of FSt, mean pairwise differences Hb Hw1 Hw2, and nn2(hudson 1992) */
  p_FST = NUMERIC_POINTER(v_FST);

  PROTECT(v_PI= allocMatrix(REALSXP,howmany,3));/* list of vectors of mean nuceotide diversity pi pi1 pi2 */
  p_PI = NUMERIC_POINTER(v_PI);
  
  PROTECT(v_TW= allocMatrix(REALSXP,howmany,3));/* list of vectors of watterson theta */
  p_TW = NUMERIC_POINTER(v_TW);

  PROTECT(v_TH= allocMatrix(REALSXP,howmany,3));/* list of vectors of Th*/
  p_TH = NUMERIC_POINTER(v_TH);
 
  PROTECT( v_D= allocMatrix(REALSXP,howmany,3));/* list of vectors of Taj D*/
  p_D = NUMERIC_POINTER(v_D);

  PROTECT(v_H= allocMatrix(REALSXP,howmany,3));/* list of vectors of H*/
  p_H = NUMERIC_POINTER(v_H);
  
  PROTECT(v_D_star= allocMatrix(REALSXP,howmany,3));/* list of vectors of Fu and Li D* */
  p_D_star = NUMERIC_POINTER(v_D_star);

  PROTECT(v_D_prime= allocMatrix(REALSXP,howmany,3));/* list of vectors of watterson theta */
  p_D_prime = NUMERIC_POINTER(v_D_prime);

  PROTECT(v_r2= allocMatrix(REALSXP,howmany,3));/* list of vectors of Th*/
  p_r2 = NUMERIC_POINTER(v_r2);

  PROTECT(v_nH= allocMatrix(INTSXP,howmany,3));/* list of vectors of H*/
  p_nH = INTEGER_POINTER(v_nH);
  
  PROTECT(v_Rm= allocMatrix(INTSXP,howmany,3));/* list of vectors of Fu and Li D* */
  p_Rm = INTEGER_POINTER(v_Rm);
  /********/
  if(typestats[1]==2 && !(nsam%2==0 && popsize[0]%2==0 &&popsize[1]%2==0))
    {
      printf("Problem: I don't have even number of pseudo-haplotypes. nsam=%d n1=%d n2=%d.\n",nsam,popsize[0],popsize[1]);
  
    }
  else
    { 
      list = cmatrix(nsam,maxsites+1);
      posit = (int *)calloc((unsigned) maxsites+1,sizeof(int));
      // pposit = (double *)calloc((unsigned) maxsites+1,sizeof(double));
      count=0;
      while( howmany-count++ )
        {
          segsites=p_seg[count-1];/* add for R */
          if( segsites >= maxsites)
            {
              maxsites = segsites + 10 ;
              posit = (int *)realloc(posit, maxsites*sizeof(int) ) ;
              //   pposit = (double *)realloc(pposit, maxsites*sizeof(double) ) ;
              biggerlist_stat(nsam,maxsites, list) ;
            } 
          if(( segsites > 0)&&(popsize[0]>0))
            {
              for( i=0; i< nsam; i++)    /*** Recover haplotypes ***/
                { 
                  strcpy( list[i], CHAR(STRING_ELT(VECTOR_ELT(haplo, count-1+i*howmany), 0))) ;    
                 /*    Rprintf("%s\n",list[i]); */
                }

              p_pos= NUMERIC_POINTER(AS_NUMERIC(VECTOR_ELT(pos, count-1))) ;  /*** Recover positions ***/
              for(i=0;i<segsites;i++)
                {
                  if((int) p_pos[i]==0)
                    posit[i]=(int) (p_pos[i]*length);
                  else  posit[i]=(int) p_pos[i];  
                }
              /*** stats pop calculation ***/
              Sstat= getstatS(popsize,nsam,  segsites, list, typestats, posit);
              FST=getFST(popsize,nsam,  segsites, list, typestats);
              nn2=getnn2(popsize,nsam,  segsites, list,(int) typestats[1]);

              for (i=0; i<12; i++) /* get stats S and their freq*/
                {
                  double tp=0;
                  if(Sstat[0][i]>0)
                    tp=(double) Sstat[1][i]/(Sstat[0][i]);
                  p_SS[count-1+i*howmany]=(int) Sstat[0][i];
                  if(i<9)  p_FS[count-1+i*howmany]=tp;
                  if(i<4)  p_FST[count-1+i*howmany]=FST[i];
                  if(i==4) p_FST[count-1+i*howmany]=nn2;
                  if(i<3)
                    {
                      p_SFL[count-1+i*howmany]=Sstat[8][i];
                      p_PI[count-1+i*howmany]=Sstat[2][i];
                      p_TW[count-1+i*howmany]=Sstat[3][i];
                      p_TH[count-1+i*howmany]=Sstat[4][i];
                      p_D[count-1+i*howmany]=Sstat[5][i];
                      p_H[count-1+i*howmany]=-Sstat[6][i];
                      p_D_star[count-1+i*howmany]=Sstat[7][i];
                      p_D_prime[count-1+i*howmany]=Sstat[9][i];/* D'*/
                      p_r2[count-1+i*howmany]=Sstat[10][i]; /* r2 */
                      p_nH[count-1+i*howmany]=Sstat[11][i]; /* nH */
                      p_Rm[count-1+i*howmany]=Sstat[12][i];/* RM*/
                    }
                }// loop on 12 S statistics
              /*** Free memory ***/
              for (i=0; i<14; ++i) free (Sstat[i]);
              free (Sstat);
              free (FST);
            }// end if segsites>0
          else
            {
              for (i=0; i<12; i++) /* get stats S and their freq*/
                {
                  p_SS[count-1+i*howmany]=(int) 0;
                  if(i<9)   p_FS[count-1+i*howmany]= 0.0;
                  if(i<5) p_FST[count-1+i*howmany]=0.0;
                  if(i==0 || i==4 ) p_FST[count-1+i*howmany]=1/0.0;
                  if(i<3)
                    {
                      p_SFL[count-1+i*howmany]=0;
                      p_PI[count-1+i*howmany]=0;
                      p_TW[count-1+i*howmany]=0;
                      p_TH[count-1+i*howmany]=0;
                      p_D[count-1+i*howmany]=0;
                      p_H[count-1+i*howmany]=0;
                      p_D_star[count-1+i*howmany]=0;
                    }
                }
            }// end if no seg sites
        }// end of loop on howmany
      for(i=0;i<nsam;i++)
        free(list[i]);
      free(list);
      free(posit);
      SET_VECTOR_ELT(stat_list, 0, v_type);
      SET_VECTOR_ELT(stat_list, 1, v_SS);
      SET_VECTOR_ELT(stat_list, 2, v_SFL);
      SET_VECTOR_ELT(stat_list, 3, v_FS);
      SET_VECTOR_ELT(stat_list, 4, v_FST);
      SET_VECTOR_ELT(stat_list, 5, v_PI);
      SET_VECTOR_ELT(stat_list, 6, v_TW);
      SET_VECTOR_ELT(stat_list, 7, v_TH);
      SET_VECTOR_ELT(stat_list, 8,v_D);
      SET_VECTOR_ELT(stat_list, 9, v_H);
      SET_VECTOR_ELT(stat_list, 10, v_D_star);
      SET_VECTOR_ELT(stat_list, 11, v_D_prime);
      SET_VECTOR_ELT(stat_list, 12, v_r2);
      SET_VECTOR_ELT(stat_list, 13, v_nH);
      SET_VECTOR_ELT(stat_list, 14, v_Rm);
    }
  setAttrib(stat_list, R_NamesSymbol, list_names);
  UNPROTECT(22);

  return(stat_list);
}


int
biggerlist_stat(  int nsam ,
                  unsigned nmax ,
                  char ** list )
{
  int i;

  maxsites = nmax  ;
  for( i=0; i<nsam; i++){
    list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
    if( list[i] == NULL ) perror( "realloc error. bigger");
  }
  return(0);
}                        

