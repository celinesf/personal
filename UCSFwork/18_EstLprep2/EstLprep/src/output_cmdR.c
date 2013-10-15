/*! \file output_cmdR.c
  \brief Code to output the grid of param values for EstL from R


*/

/**** Function to write estlike command lines into files ***/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rdefines.h>

/*! \brief Main function for R



*/
SEXP outputcmd(SEXP datafile,SEXP param) 
{
  /*** C object and functions ***/
  int nline, nfile, max,j,p,q1, q2, qa, ts, tc, mp, mc,wok,ok;
  char *file[1], *outfile, *init,*cmd, *tpch, *dum;
  double *vparam, **pvalues,totline;
  FILE *fopen(), *pfout, *pfcmd; 

  /**** Protect R object input *******/
  PROTECT(datafile = AS_CHARACTER(datafile));   /* get command line info*/
  file[0] = R_alloc(strlen(CHAR(STRING_ELT(datafile, 0))), sizeof(char)); 
  strcpy(file[0], CHAR(STRING_ELT(datafile, 0))); 

  /*** recover parameter info ***/
  PROTECT(param= AS_NUMERIC(param));
  vparam= NUMERIC_POINTER(param);
  max=0;
  for(p=0;p<11;p++)
    {
      if(max<vparam[p+2*11])
        max=(int) vparam[p+2*11];
    }// end loop on 11 parameters

  /*** calculate parameter values ***/
  pvalues = (double **)calloc( (unsigned)7,sizeof( double *));
  for( p=0; p<7; p++)                      // param fix or variable
    {
      pvalues[p] = (double *)calloc( (unsigned)(max+1),sizeof( double ));
      if(vparam[p+0*11]>=0.)// should always be true?
        {
          if(vparam[p+1*11]== -1. && vparam[p+2*11]== 1.) // param fixed
           pvalues[p][0]=(double)vparam[p+0*11];
          else                     //taken from prior
            {
              if((p==0 || p==2 || p==3) && vparam[p+0*11]==0.) //theta lowe limt 0
                for( j=0; j<(int)vparam[p+2*11];j++)
                  pvalues[p][j]=(double)(vparam[p+1*11]-vparam[p+0*11])/vparam[p+2*11]+(vparam[p+1*11]-vparam[p+0*11])*j/(vparam[p+2*11]);// end theta
              else{
                if((p==5) && vparam[p+1*11]==1.) //T_c (T_c=T_s) upper lim =1
                  for( j=0; j<(int)vparam[p+2*11];j++)
                    pvalues[p][j]=(double)vparam[p+0*11]+(vparam[p+1*11]-vparam[p+0*11])*j/(vparam[p+2*11]);             
                else
                  for(j=0; j<(int)vparam[p+2*11];j++)
                    {
                    pvalues[p][j]=(double)vparam[p+0*11]+(vparam[p+1*11]-vparam[p+0*11])*j/(vparam[p+2*11]-1);
                    }
              }//end other prior
            }// end prior   Rprintf("\n");
        }// end param specified          
    }// end loop on parameters of the model
  nline=nfile=0;
  totline=(double)1/vparam[10+0*11];
  for( p=0;p<7;p++)
    if(vparam[p+2*11]>0.)
      totline*=(double) (vparam[p+2*11]);

  outfile=(char*)calloc(100 , sizeof(char));
  sprintf( outfile,"%s-c",file[0]);
  pfcmd=fopen(outfile,"w");
  if(vparam[10+0*11]>1)
    sprintf(outfile,"%s-p%d",file[0],nfile);
  else
    sprintf(outfile,"%s-p",file[0]);
  pfout=fopen(outfile,"w");

  init=(char*)calloc(1000 , sizeof(char));
  cmd=(char*)calloc(1000 , sizeof(char));
  tpch=(char*)calloc(1000 , sizeof(char));
  dum=(char*)calloc(1000 , sizeof(char));
  if(vparam[7+0*11]>0. || vparam[7+1*11]>0.)
    { 
      for(q1=0;q1<(int) vparam[0+2*11];q1++)
        {for(mp=0;mp<(int) vparam[1+2*11] ;mp++)
            {for(q2=0;q2<(int) vparam[2+2*11] ;q2++)
                {
                  if((vparam[2+0*11])==-1.) // q2 fix
                    pvalues[2][q2]=pvalues[0][q1];
                  for(qa=0;qa<(int) vparam[3+2*11] ;qa++)
                    {
                      if((vparam[3+0*11])==-1.) // qa fix
                        pvalues[3][qa]=pvalues[0][q1];
                      for(ts=0;ts<(int) vparam[4+2*11] ;ts++)
                        {
                          sprintf( init,"%lg %lg %lg %lg %lg",(double) pvalues[0][q1],(double) pvalues[2][q2],(double) pvalues[3][qa],(double) pvalues[4][ts],(double) pvalues[1][mp]); 
                          ok=1;
                          for(tc=0;tc<(int)vparam[5+2*11] ;tc++)
                            {for(mc=0;mc<(int) vparam[6+2*11] ;mc++)
                                {
                                  if(!(pvalues[1][mp]==0 && pvalues[4][ts]==0)) //Mp>0 or Ts>0
                                    {
                                      wok=1;
                                      sprintf(cmd, "%lg %lg %lg %lg %lg %lg %lg", (double) pvalues[0][q1],(double) pvalues[2][q2],(double) pvalues[3][qa],(double) pvalues[4][ts],(double) pvalues[5][tc],(double) pvalues[1][mp],(double) pvalues[6][mc]);

                                      if(pvalues[4][ts]==0. || pvalues[5][tc]==0. || (pvalues[1][mp]==pvalues[6][mc]))
                                        {
                                          sprintf(tpch,"%lg %lg %lg %lg %lg",(double) pvalues[0][q1],(double) pvalues[2][q2],(double) pvalues[3][qa],(double) pvalues[4][ts],(double) pvalues[1][mp]);
                                          if(strcmp(init,tpch)==0 && ok)
                                            {
                                              sprintf(cmd, "%lg %lg %lg %lg %lg %lg %lg", (double) pvalues[0][q1],(double) pvalues[2][q2],(double) pvalues[3][qa],(double) pvalues[4][ts],0.,(double) pvalues[1][mp],0.);
                                              ok=0;
                                            }
                                          else if(strcmp(init,tpch)==0  && !ok) wok=0 ; 
                                        }                                 
                                      if(wok)
                                        {
                                          if(nline==0)
                                            {
                                              j=0;// rho param
                                              sprintf(tpch,"./estlikC %lg", vparam[7+j*11]);
                                              j++;
                                              while( vparam[7+j*11] >0.)
                                                {  
                                                  sprintf(dum,"%s %lg",tpch, vparam[7+j*11]);
                                                  sprintf(tpch,"%s",dum);
                                                  j++;
                                                }
                                              sprintf(dum,"%s -R %d",tpch,(int) vparam[8+0*11]);
                                              sprintf(tpch,"%s",dum);
                                              sprintf(dum,"%s -G %d",tpch,(int) vparam[9+0*11]);// howmany
                                              sprintf(tpch,"%s",dum);        
                                              if(vparam[10+0*11]>1)
                                                sprintf( outfile,"%s-p%d",file[0],nfile);
                                              else
                                                sprintf(outfile,"%s-p",file[0]);
                                              sprintf(tpch,"%s -pf %s",tpch,outfile);// file of param

                                              sprintf( outfile,"%s-S",file[0]);
                                              sprintf(tpch,"%s -r %s",tpch,outfile);// file of regions
                                              if(vparam[10+0*11]>1)
                                                sprintf( outfile,"%s-lik%d",file[0],nfile);
                                              else
                                                sprintf(outfile,"%s-lik",file[0]);
                                              sprintf(tpch,"%s > %s &",tpch,outfile);// file of regions
                                              fprintf(pfcmd,"%s\n",tpch);// cmd line
                                            }
                                          fprintf(pfout,"%s\n",cmd);
                                      
                                          nline=nline+1;
                                          if(nline>totline)
                                            {
                                              nline=0;
                                              nfile=nfile+1;
                                              fclose(pfout);
                                              sprintf( outfile,"%s-p%d",file[0],nfile);
                                              pfout=fopen(outfile,"w");
                                            }
                                        }
                                    } //if mig or div
                                }// loop mc
                            }// loop tc
                        }// ts
                    }// qa
                }//q2
            }//mp
        }//q1
    }// end rec rate>0
  else
    Rprintf("Problem: No recombination rate specified. Run avorted\n");

  fclose(pfout);
  fclose(pfcmd);
  for( p=0; p<7; p++) 
    free(pvalues[p]);
  free(pvalues);
  free(outfile);
  free(cmd);
  free(init);
  free(tpch);
  free(dum);
  UNPROTECT(2);
  return(R_NilValue);
}
