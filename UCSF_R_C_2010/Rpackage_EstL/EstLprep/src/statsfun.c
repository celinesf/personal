/*! \file statsfun.c
  \brief Functions to calculate summary statistics of polymorphisms

  - Functions for seed generation in msR.c, msC.c and estlilC.c
  - Functions specifics to distributions.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "statsfun.h"

/* count stat S1 S2 ss Sf .. and average frequencies of each stats and Pi H and ThetaH for each pop and total sample - cal by estlikC*/
double ** getstatS2(int * popsize,int nsam,int segsites, char **list, double *typestats, int * pwstat, int *pos)
{
  int i,j,k,x,x1,x2, **freq=NULL , **gamete=NULL ;
  double **pstat, *tpDstats ;
  int b, piIn1, piIn2, piBet, i1, i2, j1, j2, N1, N2;
  double Hw, Hb;//,*Fst;
  char tp1[20],tp2[20];
  double *LD=NULL, *cor=NULL ;/* FOR LD */
  char **haplist=NULL; 

  /*** init ***/
  if(typestats[1]==1)
    haplist = cmatrix(nsam,segsites+1);

  if( ! ( pstat= (double **) calloc ((unsigned) 14 , sizeof (double *))))
    perror(" calloc error. stat"); 
  // 0 Seg sites,, 1 Freq S, 2 FST, 3 pi, 4 thetaW, 5 thetaH, 6  TajD , 7 Hfay , 8 D*  9 S_o  singletons   
  // 10 LD, 11 r_square, 12 nh, 13 Rm
  if( ! (freq = (int **) calloc ((unsigned) 12 , sizeof (int *))))
    perror(" calloc error. fre");// 0anc, 1derived, 2missing, 3total - same for pop1 4567 same for pop2 8 9 10 11
  for (i=0; i<14; ++i)     
    { 
      if(i<12)
        if( ! (freq[i] = (int *) calloc ((unsigned) (segsites) , sizeof (int))))
          perror(" calloc error. freqi");
      if( ! (pstat[i] = (double *) calloc ((unsigned) 12 , sizeof (double))))
        perror(" calloc error. stati");
      for (j=0; j<12; ++j)  
        pstat[i][j]=0;
    }
  
  /* Calculate freq and stats for each segsites */
  Hw = Hb = 0.;
  gamete = (int **) calloc ((unsigned)segsites , sizeof (int *));
  for( i = 0; i <segsites; i++)
    {
      gamete[i] = (int *) calloc ((unsigned)segsites , sizeof (int));
      /*** freq total sample ***/
      freq[0][i]=myfrequency('0', i,0,nsam,list);
      freq[1][i]=myfrequency('1', i,0,nsam,list);
      freq[2][i]=myfrequency('?', i,0,nsam,list);
      if(freq[0][i]+freq[1][i]==nsam)  freq[3][i]=nsam;
      else  freq[3][i]=nsam-freq[2][i];

      /*** freq pop1 ***/
      freq[4][i]=myfrequency('0', i,0,popsize[0],list);
      freq[5][i]=myfrequency('1', i,0,popsize[0],list);
      freq[6][i]=myfrequency('?', i,0,popsize[0],list);
      if(freq[4][i]+freq[5][i]==popsize[0])  freq[7][i]=popsize[0];
      else  freq[7][i]=popsize[0]-freq[6][i];

      /*** freq pop2 ***/
      freq[8][i]=myfrequency('0', i,popsize[0],nsam,list);
      freq[9][i]=myfrequency('1', i,popsize[0],nsam,list);
      freq[10][i]=myfrequency('?', i,popsize[0],nsam,list);
      if(freq[8][i]+freq[9][i]==popsize[1])  freq[11][i]=popsize[1];
      else  freq[11][i]=popsize[1]-freq[10][i];

      /* for fst */
      piIn1 = piIn2 = piBet = 0;
      i1=freq[5][i];
      i2=freq[4][i];
      j1=freq[9][i];
      j2=freq[8][i];

      if(typestats[0] ==1)/* case outgroup unknown */
        {
          if( freq[0][i]< freq[1][i])// anc>derived change
            {
              for(j=0;j<12;j+=4)
                { 
                  k=freq[j][i];
                  freq[j][i]=  freq[j+1][i];
                  freq[j+1][i]=k;
                }        
              for(j=0;j<nsam;j++)
                {
                  if(list[j][i]=='0') list[j][i]='1';
                  else if(list[j][i]=='1') list[j][i]='0';
                }
            }
          /* for fst */
          if( i1< i2)
            {
              b=i1;
              i1= i2;
              i2=b;
            }  
          if( j1< j2)
            {
              b=j1;
              j1= j2;
              j2=b;
            }/***/
        }// End if type=maf
      if(typestats[1]==2)// genotyping data
        {
          for(j=0;j<nsam;j+=2)
            {
              if(list[j][i]=='0' && list[j+1][i]=='0' ) 
                list[j][i]='0';
              else if((list[j][i]=='1'  && list[j+1][i]=='0') || (list[j+1][i]=='1'  && list[j][i]=='0')) 
                list[j][i]='1';
              else  list[j][i]='2';
                
              list[j+1][i]='0';
            }
        }

      /* for fst */
      piIn1 += i1 * i2;
      piIn2 += j1 * j2;
      piBet += i1 * j2 + j1 * i2;
      N1 = i1 + i2;
      N2 = j1 + j2;
      if (N1>1 && N2>1) 
        {
          pstat[2][2] += (double) piIn1 / (double) (N1*(N1-1.));//hw2
          pstat[2][3]+= (double) piIn2 / (double) (N2*(N2-1.));//hw1
          Hb += (double) piBet / (double) (N1*N2);
        }
      /***/

      /*************** LD CALC *******************/
      if(i>0)
        for( j = 0; j <i; j++)
          {   
            if(pos[i]-pos[j]<=typestats[3]*1000 && pos[i]-pos[j]>=typestats[2]*1000)
              {
                /*** LD CALX ***/
                /* Run 4-gamete test on all pairs of segregating sites */
                if(pwstat[56]>0||pwstat[57]>0||pwstat[58]>0)
                  four_gametes_test(i,j,nsam, popsize,list,gamete,(int) typestats[1]);

                if(typestats[1]==1)/* haplotypic LD */
                  {
                    /* LD-r2 tot sample */
                    if(pwstat[47]>0||pwstat[50]>0)
                      {
                        LD=LDcalc((double) freq[1][j]/freq[3][j],(double) freq[1][i]/freq[3][i],  freq2 (j, i,0, nsam, list))  ;
                        if(!(isnan(LD[0]) || isinf(LD[0])))
                          {
                            pstat[10][0]+=LD[0];// mean D'
                            pstat[10][3] ++;
                          }
                        if(!(isnan(LD[1]) || isinf(LD[1])))
                          {
                            pstat[11][0]+=LD[1];// mean r^2
                            pstat[11][3]++;
                          }
                        free(LD);
                      }
                    /*LD-R2 p0 */
                    if(pwstat[48]>0||pwstat[51]>0)
                      {
                        LD=LDcalc((double) freq[5][j]/freq[7][j],(double) freq[5][i]/freq[7][i],  freq2 (j, i,0, popsize[0], list))  ;
                        if(!(isnan(LD[0]) || isinf(LD[0])))
                          {
                            pstat[10][1]+=LD[0];// mean D'
                            pstat[10][4]++;
                          }
                        if(!(isnan(LD[1]) || isinf(LD[1])))
                          {
                            pstat[11][1]+=LD[1];// mean r^2
                            pstat[11][4]++;
                          }
                        free(LD);
                      }
                    /* LD-r2 p1 */
                    if(pwstat[49]>0||pwstat[52]>0)
                      {
                        LD=LDcalc((double) freq[9][j]/freq[11][j],(double) freq[9][i]/freq[11][i],  freq2 (j, i,popsize[0], nsam, list))  ;
                        if(!(isnan(LD[0]) || isinf(LD[0])))
                          {
                            pstat[10][2]+=LD[0];// mean D'
                            pstat[10][5]++;

                          }
                        if(!(isnan(LD[1]) || isinf(LD[1])))
                          {
                            pstat[11][2]+=LD[1];// mean r^2
                            pstat[11][5]++;
                          }
                        free(LD);
                      }
                  }// if haplo
                else /* if genotypes */
                  {
                    if(pwstat[50]>0||pwstat[51]>0||pwstat[52]>0)
                      {
                        cor=  get_cor(j, i,popsize[0],freq,nsam,list);
                        for(k=0;k<3;k++)
                          {
                            if(!isnan(cor[k]) && ! isinf(cor[k]) && pwstat[50+k]>0)
                              { 
                                pstat[11][k]+=cor[k]*cor[k];// mean r^2
                                pstat[11][k+3]++;
                              }
                          }
                        free(cor);
                      }
                  }
              }// end loop on j snp
            /*** end LD calculation ***/
          }// end loop on i snp
 /****** end LD ************/

      if(freq[0][i]<freq[3][i] && freq[1][i]<freq[3][i] ) /*need to be polymorphic site in sample j */
        {
          /*** various stats needed frequencies ***/
          for(j=0;j<3;j++)/* total pop1 pop2*/
            {
              if(freq[4*j+1][i]<freq[4*j+3][i] && freq[4*j+0][i]<freq[4*j+3][i] ) /*need to be polymorphic site in sample j */
                {
                  double tpf, nd ;
                  tpf=(double)freq[4*j+1][i]/freq[4*j+3][i];
                  nd = (double)(freq[4*j+3][i]/(freq[4*j+3][i]-1.0));
                  pstat[3][j] += 2.0* (tpf)*(1.0 - tpf)*nd ;   /* pi */
                  pstat[4][j] += 1. / psum (freq[4*j+3][i] ) ; /*   thetaW */
                  pstat[5][j] += (double) freq[4*j+1][i]*freq[4*j+1][i]*2/ ( freq[4*j+3][i]*(freq[4*j+3][i] -1.0));    /* thetaH */
                  pstat[7][j] += 2.0*tpf*(2.*tpf - 1.0 )*nd ;   /* H */

                  if ((int)freq[4*j+1][i]==1) // if singleton
                    pstat[9][j] += 1; /*theta Fu and Li*/
                }
            }// end loop on samples
          /* vector of seg stat (S s1 s2 ss Sss(shared by both pop), ssf1 ssf2, ssl, ssh,sf, sf1 sf2) */
          pstat[0][0]++; 
          pstat[1][0]+=(double)freq[1][i]/freq[3][i];// total seg sites
          if((freq[5][i]==0)||(freq[9][i]==0))
            {
              if(((freq[5][i]==freq[7][i])&&(freq[9][i]==0))||((freq[9][i]==freq[11][i])&&(freq[5][i]==0)))	      
                {
                  pstat[0][9]++;                // Case Lf	
                  if(freq[5][i]==freq[7][i])pstat[0][10]++;                // Case Lf1	
                  else pstat[0][11]++;                // Case Lf2	
                } 
              else                          // Case L1, L2
                {
                  if((freq[9][i]==0)&&(freq[5][i]<freq[7][i])) 
                    {
                      pstat[0][1]++;
                      pstat[1][1]+=(double)freq[5][i]/freq[7][i];
                    }	   
                  else 
                    {
                      pstat[0][2]++;
                      pstat[1][2]+=(double)freq[9][i]/freq[11][i];
                    }		     
                }// end else specific in 1 or 2
            }// If fix or specific
          else
            {             
              pstat[0][3]++;	                   // Case Ls
              pstat[1][3]+=(double)freq[1][i]/freq[3][i];         // total freq of Ss
              if((double) freq[1][i]/freq[3][i] <= .1) // S_sl
                {
                  pstat[0][7]++;
                  pstat[1][7]+=(double)(freq[9][i]+freq[5][i])/freq[3][i];
                }// fixed in 1
              else 
                {
                  pstat[0][8]++;
                  pstat[1][8]+=(double)(freq[9][i]+freq[5][i])/freq[3][i];
                }// fixed in 2
              if(freq[5][i]<freq[7][i] && freq[9][i]<freq[11][i]) // case poly in both 
                {
                  pstat[0][4]++;
                  pstat[1][4]+=(double)(freq[9][i]+freq[5][i])/freq[3][i];
                }
              else
                {
                  if(freq[5][i]==freq[7][i]) 
                    {
                      pstat[0][5]++;
                      pstat[1][5]+=(double)(freq[9][i]+freq[5][i])/freq[3][i];
                    }// fixed in 1
                  else 
                    {
                      pstat[0][6]++;
                      pstat[1][6]+=(double)(freq[9][i]+freq[5][i])/freq[3][i];
                    }// fixed in 2                 
                }
            }// shared derievd of maf alleles
        }// end if segregating        
    }// end loop on seg site to get frequencies
  /* for fst */
  if(pwstat[24]>0 ||pwstat[25]>0||pwstat[26]>0|| pwstat[27]>0)
    {
      Hw+=pstat[2][2] +pstat[2][3] ;
      pstat[2][1]=Hb;
      sprintf(tp1,"%lg",Hw);
      sprintf(tp2,"%lg",Hb);
      if(pwstat[24]>0 )
        {
          if(strcasecmp(tp1,tp2)==0)
            pstat[2][0]=0;
          else   pstat[2][0] =  (double) 1. -(Hw / Hb);
        }
    }
  /***/


 /*    for(j=0;j<nsam;j+=2) */
/*         printf("%s\n",list[j] ); */


  /*** LD statistcis ***/
  /* RM calculation */
  if(pwstat[56]>0)
  pstat[13][0]=min_rec(0, segsites, 1, gamete,1);
  if(pwstat[57]>0)
  pstat[13][1]=min_rec(0, segsites, 1, gamete,2);
  if(pwstat[58]>0)
  pstat[13][2]=min_rec(0, segsites, 1, gamete,3);
  /* mean LD and r2 */
  if(pwstat[47]>0||pwstat[48]>0||pwstat[49]>0|| pwstat[50]>0||pwstat[51]>0||pwstat[52]>0)
    {
      for(i=0;i<3;i++)
        {  
          if(pwstat[47+i]>0)// D'
            {
              if(pstat[10][i+3]>0)
                pstat[10][i]=(double) pstat[10][i]/pstat[10][i+3];
              else         
                pstat[10][i]=log(-1);
            }
          if(pwstat[50+i]>0)// R2
            {
              if(pstat[11][i+3]>0)
                pstat[11][i]=(double) pstat[11][i]/pstat[11][i+3];
              else    pstat[11][i]=log(-1);
            }
        }
    }
  /*  nhaplo */
  if(pwstat[53]>0||pwstat[54]>0||pwstat[55]>0)
    {
      if(typestats[1]==1)
        {
          for( i=0,x=0,x1=0,x2=0; i<nsam;i++) 
            { 
              strcpy(haplist[x], list[i]);
              //sprintf(haplist[x],"%s\n", list[i]);
              for (j=0, k=0; j<segsites; ++j) 
                if (haplist[x][j]=='0' || haplist[x][j]=='1' || haplist[x][j]=='2')           
                  ++k;
              if ((2*k) >= segsites)// dealing with missing data??
                {
                  if(i<popsize[0]) x1++;
                  ++x;
                }
            }
          if( pwstat[53]>0)
            {
              if (pstat[0][0]==1) // tot
                pstat[12][0]= 2;
              else if (pstat[0][0]>1 ) 
                pstat[12][0] = hapcount (haplist,0, x, 0, segsites-1);
            }
          if( pwstat[54]>0)
            {
              if (pstat[0][1]+pstat[0][4]+pstat[0][6] ==1) // pop1
                pstat[12][1]= 2;
              else if (pstat[1][1]+pstat[0][4]+pstat[0][6] >0)
                pstat[12][1] = hapcount (haplist,0, x1, 0, segsites-1);
            }
          if( pwstat[55]>0)
            {
              if (pstat[0][2]+pstat[0][4]+pstat[0][5]==1) // pop2
                pstat[12][2]= 2;
              else if (pstat[1][2]+pstat[0][4]+pstat[0][5]>0)
                pstat[12][2] = hapcount (haplist,x1, x, 0, segsites-1);
            }
        }
      else
        pstat[12][0]=pstat[12][1]=pstat[12][2]=log(-1);
    }// if nhaplo
  /*** end LD stats ***/

  /* get tajima'D and Fu and Li D-star */
  for(i=0;i<3;i++)/* total pop1 pop2*/
    {
      k=freq[4*i+3][0]; /** define max num of chromosome for tajimas' D calculation */
      for (j=1; j<segsites; ++j) 
        if (freq[4*i+3][j] > k) 
          k=freq[4*i+3][j] ;
      if(i==0 && (pwstat[29]>0 ||pwstat[32]>0||pwstat[35]>0||pwstat[38]>0 ||pwstat[41]>0 ||pwstat[44]>0 ))
        {       
          tpDstats= tajD_D_star(k,(int) pstat[0][i],pstat[3][i],pstat[9][i],(int) typestats[0] );
          pstat[6][i]=tpDstats[0];
          pstat[8][i]=tpDstats[1];
          free(tpDstats);
        }
      else if (i==1 && (pwstat[30]>0 ||pwstat[33]>0||pwstat[36]>0||pwstat[39]>0 ||pwstat[42]>0 ||pwstat[45]>0 ))       
        {         
          tpDstats= tajD_D_star(k,(int) (pstat[0][i]+pstat[0][4]+pstat[0][6]),pstat[3][i],pstat[9][i],(int) typestats[0] );
          pstat[6][i]=tpDstats[0];
          pstat[8][i]=tpDstats[1];
          free(tpDstats);
        }
      else if(pwstat[31]>0 ||pwstat[34]>0||pwstat[37]>0||pwstat[40]>0 ||pwstat[43]>0 ||pwstat[46]>0 )
        {   
          tpDstats= tajD_D_star(k,(int) (pstat[0][i]+pstat[0][4]+pstat[0][5]),pstat[3][i],pstat[9][i],(int)typestats[0] );   
          pstat[6][i]=tpDstats[0];
          pstat[8][i]=tpDstats[1];
          free(tpDstats);
        }
    }// end loop on sample to get D stats

  /* free mem */
  for (i=0; i<12; ++i) free(freq[i]);
  free (freq);
 for(i=0;i<segsites;i++)
    free(gamete[i]);
  free(gamete);
  if(typestats[1]==1)
    {
      for(i=0;i<nsam;i++)
        free(haplist[i]); 
      free(haplist); 
    }
  return pstat;
}// end getsstat2



/* count stat S1 S2 ss Sf .. and average frequencies of each stats and Pi H and ThetaH for each pop and total sample */
double ** getstatS(int * popsize,int nsam,int segsites, char **list, double* typestats, int *pos)
{
  int i,j,k,x,x1,x2, **freq=NULL, **gamete=NULL ;
  double **pstat=NULL, *tpDstats=NULL, *cor=NULL ;
  double *LD=NULL;/* FOR LD */
  char **haplist=NULL; 

  /*** init ***/
  if(typestats[1]==1)
    haplist = cmatrix(nsam,segsites+1);
  if( ! ( pstat= (double **) calloc ((unsigned) 13 , sizeof (double *))))
    perror(" calloc error. stat"); // Seg sites, Freq S, pi, thetaW,thetaH, Hfay ,  ThetaFL,  TajD ,  D*,D' r2 nh Rm'     
  if( ! (freq = (int **) calloc ((unsigned) 12 , sizeof (int *))))
    perror(" calloc error. fre");// 0anc, 1derived, 2missing, 3total - same for pop1 4567 same for pop2 8 9 10 11

  for (i=0; i<13; ++i)     
    {  
      if(i<12)
        if( ! (freq[i] = (int *) calloc ((unsigned) (segsites) , sizeof (int))))
          perror(" calloc error. freqi");
      if( ! (pstat[i] = (double *) calloc ((unsigned) 12 , sizeof (double))))
        perror(" calloc error. stati");
      for (j=0; j<12; ++j)  
        pstat[i][j]=0;
    }
  gamete = (int **) calloc ((unsigned)segsites , sizeof (int *));
  
  /* Calculate freq and stats for each segsites */
  for( i = 0; i <segsites; i++)
    {
      gamete[i] = (int *) calloc ((unsigned)segsites , sizeof (int));
      /*** freq total sample ***/
      freq[0][i]=myfrequency('0', i,0,nsam,list);
      freq[1][i]=myfrequency('1', i,0,nsam,list);
      freq[2][i]=myfrequency('?', i,0,nsam,list);
      if(freq[0][i]+freq[1][i]==nsam)  freq[3][i]=nsam;
      else  freq[3][i]=nsam-freq[2][i];
      /*** freq pop1 ***/
      freq[4][i]=myfrequency('0', i,0,popsize[0],list);
      freq[5][i]=myfrequency('1', i,0,popsize[0],list);
      freq[6][i]=myfrequency('?', i,0,popsize[0],list);
      if(freq[4][i]+freq[5][i]==popsize[0])  freq[7][i]=popsize[0];
      else  freq[7][i]=popsize[0]-freq[6][i];

      /*** freq pop2 ***/
      freq[8][i]=myfrequency('0', i,popsize[0],nsam,list);
      freq[9][i]=myfrequency('1', i,popsize[0],nsam,list);
      freq[10][i]=myfrequency('?', i,popsize[0],nsam,list);
      if(freq[8][i]+freq[9][i]==popsize[1])  freq[11][i]=popsize[1];
      else  freq[11][i]=popsize[1]-freq[10][i];
      if(typestats[0]==1)/* case outgroup unknown */
        {
          if( freq[0][i]< freq[1][i])// anc>derived change
            {
              for(j=0;j<12;j+=4)
                {                   
                  k=freq[j][i];
                  freq[j][i]=  freq[j+1][i];
                  freq[j+1][i]=k;
                }
              for(j=0;j<nsam;j++)
                {
                  if(list[j][i]=='0') list[j][i]='1';
                  else if(list[j][i]=='1') list[j][i]='0';
                }
            }
        }// End if type=maf
      if(typestats[1]==2)// genotyping data
        {
          for(j=0;j<nsam;j+=2)
            {
              if(list[j][i]=='0' && list[j+1][i]=='0' ) 
                list[j][i]='0';
              else if((list[j][i]=='1'  && list[j+1][i]=='0') || (list[j+1][i]=='1'  && list[j][i]=='0')) 
                list[j][i]='1';
              else  list[j][i]='2';
                
              list[j+1][i]='0';
            }
        }

      if(i>0)
        for( j = 0; j <i; j++)
          {   
            if(pos[i]-pos[j]<=typestats[3]*1000 && pos[i]-pos[j]>=typestats[2]*1000)
              {
                /*** LD CALX ***/
                /* Run 4-gamete test on all pairs of segregating sites */
                four_gametes_test(i,j,nsam, popsize,list,gamete,(int) typestats[1]);

                if(typestats[1]==1)/* haplotypic LD */
                  {
                    /* LD-r2 tot sample */
                    LD=LDcalc((double) freq[1][j]/freq[3][j],(double) freq[1][i]/freq[3][i],  freq2 (j, i,0, nsam, list))  ;
                    if(!(isnan(LD[0]) || isinf(LD[0])))
                      {
                        pstat[9][0]+=LD[0];// mean D'
                        pstat[9][3] ++;
                      }
                    if(!(isnan(LD[1]) || isinf(LD[1])))
                      {
                        pstat[10][0]+=LD[1];// mean r^2
                        pstat[10][3]++;
                      }
                    free(LD);
                    /*LD-R2 p0 */
                    LD=LDcalc((double) freq[5][j]/freq[7][j],(double) freq[5][i]/freq[7][i],  freq2 (j, i,0, popsize[0], list))  ;
                    if(!(isnan(LD[0]) || isinf(LD[0])))
                      {
                        pstat[9][1]+=LD[0];// mean D'
                        pstat[9][4]++;
                      }
                    if(!(isnan(LD[1]) || isinf(LD[1])))
                      {
                        pstat[10][1]+=LD[1];// mean r^2
                        pstat[10][4]++;
                      }
                    free(LD);
                    /* LD-r2 p1 */
                    LD=LDcalc((double) freq[9][j]/freq[11][j],(double) freq[9][i]/freq[11][i],  freq2 (j, i,popsize[0], nsam, list))  ;
                    if(!(isnan(LD[0]) || isinf(LD[0])))
                      {
                        pstat[9][2]+=LD[0];// mean D'
                        pstat[9][5]++;

                      }
                    if(!(isnan(LD[1]) || isinf(LD[1])))
                      {
                        pstat[10][2]+=LD[1];// mean r^2
                        pstat[10][5]++;
                      }
                    free(LD);
                  }// if haplo
                else /* if genotypes */
                  {
                    cor=  get_cor(j, i,popsize[0],freq,nsam,list);
                    for(k=0;k<3;k++)
                      {
                        if(!isnan(cor[k]) && ! isinf(cor[k]))
                          { 
                            pstat[10][k]+=cor[k]*cor[k];// mean r^2
                            pstat[10][k+3]++;
                          }
                      }
                    free(cor);
                  }
              }// end loop on j snp
            /*** end LD calculation ***/
          }// end loop on i snp

      if(freq[0][i]<freq[3][i] && freq[1][i]<freq[3][i] ) /*need to be polymorphic site in sample j */
        {
          /*** various stats needed frequencies ***/
          for(j=0;j<3;j++)/* total pop1 pop2*/
            {
              if(freq[4*j+1][i]<freq[4*j+3][i] && freq[4*j+0][i]<freq[4*j+3][i] ) /*need to be polymorphic site in sample j */
                {
                  double tpf, nd ;
                  tpf=(double)freq[4*j+1][i]/freq[4*j+3][i];
                  nd = (double)(freq[4*j+3][i]/(freq[4*j+3][i]-1.0));
                  pstat[2][j] += 2.0* (tpf)*(1.0 - tpf)*nd ;   /* pi */
                  pstat[3][j] += 1. / psum (freq[4*j+3][i] ) ; /*   thetaW */
                  pstat[4][j] += (double) freq[4*j+1][i]*freq[4*j+1][i]*2/ ( freq[4*j+3][i]*(freq[4*j+3][i] -1.0));    /* thetaH */
                  pstat[6][j] += 2.0*tpf*(2.*tpf - 1.0 )*nd ;   /* H */

                  if ((int)freq[4*j+1][i]==1) // if singleton
                    pstat[8][j] += 1; /*theta Fu and Li*/
                }
            }// end loop on samples
          /* vector of seg stat (S s1 s2 ss Sss(shared by both pop), ssf1 ssf2, ssl, ssh,sf, sf1 sf2) */
          pstat[0][0]++; 
          pstat[1][0]+=(double)freq[1][i]/freq[3][i];// total seg sites
          if((freq[5][i]==0)||(freq[9][i]==0))
            {
              if(((freq[5][i]==freq[7][i])&&(freq[9][i]==0))||((freq[9][i]==freq[11][i])&&(freq[5][i]==0)))	      
                {
                  pstat[0][9]++;                // Case Lf	
                  if(freq[5][i]==freq[7][i])pstat[0][10]++;                // Case Lf1	
                  else pstat[0][11]++;                // Case Lf2	
                } 
              else                          // Case L1, L2
                {
                  if((freq[9][i]==0)&&(freq[5][i]<freq[7][i])) 
                    {
                      pstat[0][1]++;
                      pstat[1][1]+=(double)freq[5][i]/freq[7][i];
                    }	   
                  else 
                    {
                      pstat[0][2]++;
                      pstat[1][2]+=(double)freq[9][i]/freq[11][i];
                    }		     
                }// end else specific in 1 or 2
            }// If fix or specific
          else
            {             
              pstat[0][3]++;	                   // Case Ls
              pstat[1][3]+=(double)freq[1][i]/freq[3][i];         // total freq of Ss
              if((double) freq[1][i]/freq[3][i] <= .1) // S_sl
                {
                  pstat[0][7]++;
                  pstat[1][7]+=(double)(freq[9][i]+freq[5][i])/freq[3][i];
                }// fixed in 1
              else 
                {
                  pstat[0][8]++;
                  pstat[1][8]+=(double)(freq[9][i]+freq[5][i])/freq[3][i];
                }// fixed in 2
              if(freq[5][i]<freq[7][i] && freq[9][i]<freq[11][i]) // case poly in both 
                {
                  pstat[0][4]++;
                  pstat[1][4]+=(double)(freq[9][i]+freq[5][i])/freq[3][i];
                }
              else
                {
                  if(freq[5][i]==freq[7][i]) 
                    {
                      pstat[0][5]++;
                      pstat[1][5]+=(double)(freq[9][i]+freq[5][i])/freq[3][i];
                    }// fixed in 1
                  else 
                    {
                      pstat[0][6]++;
                      pstat[1][6]+=(double)(freq[9][i]+freq[5][i])/freq[3][i];
                    }// fixed in 2                 
                }
            }// shared derievd of maf alleles
        }// end if segregating        
    }// end loop on seg site to get frequencies
        
   /*  for(j=0;j<nsam;j+=2) */
/*         printf("%s\n",list[j] ); */


  /*** LD statistcis ***/
  /* RM calculation */
  pstat[12][0]=min_rec(0, segsites, 1, gamete,1);
  pstat[12][1]=min_rec(0, segsites, 1, gamete,2);
  pstat[12][2]=min_rec(0, segsites, 1, gamete,3);
  /* mean LD and r2 */
  for(i=0;i<3;i++)
    {
      if(pstat[9][i+3]>0)
        pstat[9][i]=(double) pstat[9][i]/pstat[9][i+3];
      else         
        pstat[9][i]=log(-1);
      if(pstat[10][i+3]>0)
        pstat[10][i]=(double) pstat[10][i]/pstat[10][i+3];
      else    pstat[10][i]=log(-1);
    }
  /*  nhaplo */
  if(typestats[1]==1)
    {
      for( i=0,x=0,x1=0,x2=0; i<nsam;i++) 
        { 
          strcpy(haplist[x], list[i]);
          //sprintf(haplist[x],"%s\n", list[i]);
          for (j=0, k=0; j<segsites; ++j) 
            if (haplist[x][j]=='0' || haplist[x][j]=='1' || haplist[x][j]=='2')           
              ++k;
          if ((2*k) >= segsites)// dealing with missing data??
            {
              if(i<popsize[0]) x1++;
              ++x;
            }
        }
      if (pstat[0][0]==1) // tot
        pstat[11][0]= 2;
      else if (pstat[0][0]>1) 
        pstat[11][0] = hapcount (haplist,0, x, 0, segsites-1);
      if (pstat[0][1]+pstat[0][4]+pstat[0][6] ==1) // pop1
        pstat[11][1]= 2;
      else if (pstat[1][1]+pstat[0][4]+pstat[0][6] >0)
        pstat[11][1] = hapcount (haplist,0, x1, 0, segsites-1);
      if (pstat[0][2]+pstat[0][4]+pstat[0][5]==1) // pop2
        pstat[11][2]= 2;
      else if (pstat[1][2]+pstat[0][4]+pstat[0][5]>0)
        pstat[11][2] = hapcount (haplist,x1, x, 0, segsites-1);
    }
  else
    pstat[11][0]=pstat[11][1]=pstat[11][2]=log(-1);
  /*** end LD stats ***/

  /*** get tajima'D and Fu and Li D-star ***/
  for(i=0;i<3;i++)/* total pop1 pop2*/
    {
      k=freq[4*i+3][0]; /** define max num of chromosome for tajimas' D calculation */
      for (j=1; j<segsites; ++j) 
        if (freq[4*i+3][j] > k) 
          k=freq[4*i+3][j] ;
      if(i==0 )
        {       
          tpDstats= tajD_D_star(k,(int) pstat[0][i],pstat[2][i],pstat[8][i],(int) typestats[0] );
          pstat[5][i]=tpDstats[0];
          pstat[7][i]=tpDstats[1];
          free(tpDstats);
        }
      else if (i==1)       
        {         
          tpDstats= tajD_D_star(k,(int) (pstat[0][i]+pstat[0][4]+pstat[0][6]),pstat[2][i],pstat[8][i],(int) typestats[0] );
          pstat[5][i]=tpDstats[0];
          pstat[7][i]=tpDstats[1];
          free(tpDstats);
        }
      else
        {   
          tpDstats= tajD_D_star(k,(int) (pstat[0][i]+pstat[0][4]+pstat[0][5]),pstat[2][i],pstat[8][i],(int) typestats[0] );   
          pstat[5][i]=tpDstats[0];
          pstat[7][i]=tpDstats[1];
          free(tpDstats);
        }
    }// end loop on sample to get D stats

  /* free mem */
  for (i=0; i<12; ++i) free(freq[i]);
  free (freq);
  for(i=0;i<segsites;i++)
    free(gamete[i]);
  free(gamete);
  if(typestats[1]==1)
    {
      for(i=0;i<nsam;i++)
        free(haplist[i]); 
      free(haplist); 
    }
  return pstat;
}// end getstatS



/********************************************************************/

/**** calculate stats functions  not touched between R and C ********/
/********************************************************************/

/* Calculate Fst and mean pirwise differences*/
double * getFST(int * popsize, int nsam, int segsites, char ** list, double* typestats)
{
  int a,b, piIn1, piIn2, piBet, i1, i2, j1, j2, N1, N2;
  double Hw, Hb,*Fst;
  char tp1[20],tp2[20];

  if( ! (Fst= (double *) calloc ((unsigned) 4 , sizeof (double))))
    perror(" calloc error. fst");
  for (a=0; a<4;a++) {Fst[a]=0.0;}
  Hw = Hb = 0.;
  for (a=0; a<segsites; ++a) 
    {
      piIn1 = piIn2 = piBet = 0;
      i1=i2=j1=j2=0; 
      if(typestats[1]==1)// phased
        {
          for (b=0; b<popsize[0]; ++b) 
            {
              if (list[b][a]=='1')
                ++i1;
              else if (list[b][a]=='0')
                ++i2;
            }
          for (b=popsize[0]; b<nsam; ++b) 
            {
              if (list[b][a]=='1')
                ++j1;
              else if (list[b][a]=='0')
                ++j2;
            }
        }
      else  if(typestats[1]==2)// unphased
        {
          for (b=0; b<popsize[0]; b+=2) 
            {
              if (list[b][a]=='1')
                {
                  ++i2;
                  ++i1;
                }
              if (list[b][a]=='2')
                i1+=2;
              else if (list[b][a]=='0')
                i2+=2;
            }
          for (b=popsize[0]; b<nsam;b+=2) 
            {
              if (list[b][a]=='1')
                {
                  ++j1;
                  ++j2;;
                }
              if (list[b][a]=='2')
                j1+=2;
              else if (list[b][a]=='0')
                  j2+=2;
            }
        }
      if(typestats[0]==1)
        {
          if( i1< i2)
            {
              b=i1;
              i1= i2;
              i2=b;
            }  
          if( j1< j2)
            {
              b=j1;
              j1= j2;
              j2=b;
            }  
        }
      piIn1 += i1 * i2;
      piIn2 += j1 * j2;
      piBet += i1 * j2 + j1 * i2;
      N1 = i1 + i2;
      N2 = j1 + j2;
      
      if (N1>1 && N2>1) 
        {
          Fst[2] += (double) piIn1 / (double) (N1*(N1-1.));
          Fst[3]+= (double) piIn2 / (double) (N2*(N2-1.));
          Hb += (double) piBet / (double) (N1*N2);
        }
    }// end loop on seg sites
  Hw+=Fst[2]+Fst[3];
  Fst[1]=Hb;
  sprintf(tp1,"%lg",Hw);
  sprintf(tp2,"%lg",Hb);
  if(strcasecmp(tp1,tp2)==0)
    Fst[0]=0;
  else  Fst[0] =  (double) 1. -(Hw / Hb);
  return (Fst);
}// end getFST

/* Calculate the nearest neighboord ststistics*/
double getnn2(int * popsize, int nsam, int segsites, char ** list, int typestats)
{
  int a, b, mindist1, **distmat;
  int mindist2, d1, d2; 
  int distance (int , int , int , char **);
  
  if( ! (distmat = (int **) calloc ((unsigned)nsam , sizeof (int *))))
    perror(" calloc error. distmat");
  for (a=0; a<nsam; ++a)
    if( ! (distmat[a] = (int *)calloc ((unsigned)nsam , sizeof (int))))
      perror(" calloc error. dismat");

  d1 = d2 = 0;
  if(typestats==1)
    {
      for (a=0; a<nsam; ++a)
        for (b=0; b<nsam; ++b) 
          distmat[a][b] = distance (a, b, segsites, list);
      for (a=0; a<popsize[0]; ++a) 
        {
          mindist1 = mindist2 = 999;
          for (b=0; b<popsize[0]; ++b)
            if ((a != b) && distmat[a][b] < mindist1)
              mindist1 = distmat[a][b];
          for (b=popsize[0]; b<nsam; ++b)
            if (distmat[a][b] < mindist2)
              mindist2 = distmat[a][b];        
          d1 += mindist1; d2 += mindist2;
        }
      for (a=popsize[0]; a<nsam; ++a) 
        {
          mindist1 = mindist2 = 999;
          for (b=popsize[0]; b<nsam; ++b)
            if ((a != b) && distmat[a][b] < mindist1)
              mindist1 = distmat[a][b];
          for (b=0; b<popsize[0]; ++b)
            if (distmat[a][b] < mindist2)
              mindist2 = distmat[a][b];        
          d1 += mindist1; d2 += mindist2;
        }
    }
  else  if(typestats==2)
    {
      for (a=0; a<nsam; a+=2)
        for (b=0; b<nsam; b+=2) 
          distmat[a][b] = distance_ind(a, b, segsites, list);
      for (a=0; a<popsize[0]; a+=2) 
        {
          mindist1 = mindist2 = 999;
          for (b=0; b<popsize[0]; b+=2)
            if ((a != b) && distmat[a][b] < mindist1)
              mindist1 = distmat[a][b];
          for (b=popsize[0]; b<nsam; b+=2)
            if (distmat[a][b] < mindist2)
              mindist2 = distmat[a][b];        
          d1 += mindist1; d2 += mindist2;
        }
      for (a=popsize[0]; a<nsam; a+=2) 
        {
          mindist1 = mindist2 = 999;
          for (b=popsize[0]; b<nsam; b+=2)
            if ((a != b) && distmat[a][b] < mindist1)
              mindist1 = distmat[a][b];
          for (b=0; b<popsize[0]; b+=2)
            if (distmat[a][b] < mindist2)
              mindist2 = distmat[a][b];        
          d1 += mindist1; d2 += mindist2;
        }
    }
  for (a=0; a<nsam; ++a)
    free(distmat[a]);
  free(distmat);
  return ((double) ((double) d1 / (double)d2)); 
}// end Snn


/* Calculate the distance between two chromosomes 
 used to calc Snn
*/
int distance(int a, int b, int segsites, char ** list)
{
  int c, d;

  for (c=0, d=0; c<segsites; ++c)
    if ((list[a][c]=='0' && list[b][c]=='1') ||
        (list[a][c]=='1' && list[b][c]=='0'))
      ++d;
  return (d);
}

int distance_ind(int a, int b, int segsites, char ** list)
{
  int c, d;

  for (c=0, d=0; c<segsites; ++c)
    {
      if ((list[a][c]=='0' && list[b][c]=='1') ||
          (list[a][c]=='1' && list[b][c]=='0') || 
          (list[a][c]=='2' && list[b][c]=='1') ||
          (list[a][c]=='1' && list[b][c]=='2'))
        ++d;
      if ((list[a][c]=='0' && list[b][c]=='2') ||
          (list[a][c]=='2' && list[b][c]=='0')) 
        d+=2;
    }
  return (d);
}

/* tajD.c -- Calculates Tajima's D (1989) and Fu and Li's D* (1993).
   Freq. spectrum etc. are embedded in the source code.  */
double* tajD_D_star(int N, int S, double k, int sing, int type)
{
  int  b;
  double a1, a2, b1, b2, c1, c2, e1, e2, c, d, mu, nu,  num, den, D, sqrt();
  double num2, den2, Dstar;

  a1 = a2 = 0.0;
  for (b=1; b<N; ++b) {
    a1 += 1. / (double) b;
    a2 += 1. / (double) (b * b);
  }
  b1 = (double) (N + 1.) / (3. * (N - 1.));
  b2 = (double) (2. * N * N + 2. * N + 6.) / (9. * N * (N - 1.));
  c1 = b1 - 1. / a1;
  c2 = b2 - (double) (N + 2.) / (a1 * N) + a2 / (a1 * a1);
  e1 = c1 / a1;
  e2 = c2 / (a1 * a1 + a2);
  num = k - (double) S / a1;
  den = sqrt (e1 * S + e2 * S * (S - 1));
  D = num / den;

  c = (double) (2. * (N * a1 - 2. * (N - 1.))) / ((N - 1.) * (N - 2.));
  if(type)
    {
      d = (double) (1.5 - (2. * (a1 + 1. / N) - 3.) / (N - 2.) - 1. / N);
      d *= 2. / (N - 1.);
      d += c + (N - 2.) / (N * N - 2. * N + 1.);
      
      nu = (a2 * N * N) / ((N - 1.) * (N - 1.)) + a1 * a1 * d;
      nu -= (2. * N * a1 * (a1 +1)) / ((N - 1.) * (N - 1.));
      nu /= (a1 * a1 + a2);
      mu = N * (a1 - N / (N - 1.)) / (N - 1.) - nu;
      
      num2 = (double) N * S / (N - 1.) - a1 * sing;
    }
  else
    {      
      nu = c - (double) (N+1)/(N-1);
      nu *= (double) a1*a1 /(a2 + a1*a1);
      nu +=1; 
      num2 = (double)  S  - a1 * sing;

      mu =a1 - 1 - nu;
    }
  den2 = sqrt (mu * S + nu * S * S);
  Dstar = num2 / den2;

  double *tp;
  if( ! (tp=(double*) calloc((unsigned)2, sizeof(double))))
    perror(" calloc error. tp");
  tp[0]=D;
  tp[1]=Dstar;
  return (tp);
}  

/****************** OTHER FUNCTION ************/
/* sum 1/a - called by getstats functions*/
double psum (int n)
{
  int a;
  double sum = 0.;

  for (a=1; a<n; ++a)
    sum += 1. / (double) a;
  return (sum);
}

/*** Functions to calculate stats ***/
int
frequency( char allele,int site,int nsam,  char **list)
{
  int i, count=0;
  for( i=0; i<nsam; i++) 
    count += ( list[i][site]==allele ? 1: 0 ) ;
  
  return( count);
}        
/*
 calculate the frequency of the derived alleles at a specific site between start and nsam*/
int myfrequency( char allele,int site,int start, int nsam,  char **list)
{
  int i, count=0;

  for( i=start; i<nsam; i++) 
    {
      // printf("here %d -%s\n",i,list[i]);
      count += ( list[i][site]==allele ? 1: 0 ) ;
    }
  return( count);
} 

/* init list of haplotypes */
char **
cmatrix(nsam,len)
     int nsam, len;
{
  int i;
  char **m;

  if( ! ( m = (char **) calloc( (unsigned) nsam,sizeof( char* ) ) ) )
    perror("alloc error in cmatrix") ;
  for( i=0; i<nsam; i++) {
    if( ! ( m[i] = (char *) calloc( (unsigned) len,sizeof( char ) )))
      perror("alloc error in cmatric. 2");
  }
  return( m );
}


/***************************LDCALC*****************************/

double* get_cor(int j, int i,int psize,int ** pfreq,int ns,char ** plist)
{
  int k=0,a=0,b=0, nc[3];
  double *cor, **mean,**sigma;

  if( ! (mean=(double**) calloc((unsigned)3, sizeof(double *))))
    perror(" calloc error. mean");
  if( ! (sigma=(double**) calloc((unsigned)3, sizeof(double *))))
    perror(" calloc error. sigma");
  if( ! (cor=(double*) calloc((unsigned)3, sizeof(double ))))
    perror(" calloc error. cor");

  for(k=0; k<3;k++)
    {
      if( ! (mean[k]=(double*) calloc((unsigned)2, sizeof(double ))))
        perror(" calloc error. mean");
      if( ! (sigma[k]=(double*) calloc((unsigned)2, sizeof(double ))))
        perror(" calloc error. sigma");

      mean[k][0]=(double) pfreq[1+k*4][j]/(pfreq[3+k*4][j]/2);
      mean[k][1]=(double) pfreq[1+k*4][i]/(pfreq[3+k*4][i]/2);

      sigma[k][0]=sigma[k][1]=cor[k]=nc[k]=0;
    }
  for(k=0;k<ns;k++)
    {
      a=b=0;
      if(plist[k][j]=='1') a=1;
      else if(plist[k][j]=='2') a=2;
      if(plist[k][i]=='1') b=1;
      else if(plist[k][i]=='2') b=2;
      if( plist[k][j]!='?' &&  plist[k][i]!='?') 
        { 
          cor[0]+=(a-mean[0][0])*(b-mean[0][1]);
          sigma[0][0]+=(a-mean[0][0])*(a-mean[0][0]);
          sigma[0][1]+=(b-mean[0][1])*(b-mean[0][1]);
          nc[0]++;
          if(k<psize)
            {
              cor[1]+=(a-mean[1][0])*(b-mean[1][1]);
              sigma[1][0]+=(a-mean[1][0])*(a-mean[1][0]);
              sigma[1][1]+=(b-mean[1][1])*(b-mean[1][1]);
              nc[1]++;
            }
          else
            {
              cor[2]+=(a-mean[2][0])*(b-mean[2][1]);
              sigma[2][0]+=(a-mean[2][0])*(a-mean[2][0]);
              sigma[2][1]+=(b-mean[2][1])*(b-mean[2][1]);
              nc[2]++;
            }
        }
      k++;
    }
  for(k=0; k<3;k++)
    {
      cor[k]=((double) cor[k]/(nc[k]-1))/(sqrt((double) sigma[k][0]/(nc[k]-1)) *sqrt((double) sigma[k][1]/(nc[k]-1)));
      free(sigma[k]);
      free(mean[k]);
    }
  free(sigma);
  free(mean);
  return (cor);
}// end get_cor


int isgam(int c,int n1,int *pgtest, int gam)
{
  if(c<n1) 
    {
      ++pgtest[1];
      ++pgtest[0];
      gam++;
    }
  else  
    {
      ++pgtest[2];
      if(gam==0) ++pgtest[0];
      gam++;
      if(gam==1) gam++;// = 3 for only pop2
    }
  return gam;

}


void four_gametes_test(int si,int sj, int nsam, int *pop,char ** plist, int **pgamete,int haplo)
{
  int k=0,gtest[3],gam[4];
  gtest[0] = gtest[1]=gtest[2] = 0;
  gam[0]=gam[1]=gam[2]=gam[3]=0;// 1 for tot, 2 for p0, 3 for p3, 4 for p1&p2 
  for (k=0; k<nsam; ++k)
    {
      if(haplo==1)
        {
          if (plist[k][sj]=='0' && plist[k][si]=='0' && (gam[0]==0 || (gam[0]==1 && k>=pop[0]) ))
            gam[0]= isgam(k,pop[0],gtest, gam[0]);
          if (plist[k][sj]=='0' && plist[k][si]=='1' && (gam[1]==0 || (gam[1]==1 && k>=pop[0]) ))
            gam[1]= isgam(k,pop[0],gtest, gam[1]);

          if (plist[k][sj]=='1' && plist[k][si]=='0' && (gam[2]==0 || (gam[2]==1 && k>=pop[0]) ))
            gam[2]= isgam(k,pop[0],gtest, gam[2]);
          if (plist[k][sj]=='1' && plist[k][si]=='1' && (gam[3]==0|| (gam[3]==1 && k>=pop[0]) ))
            gam[3]= isgam(k,pop[0],gtest, gam[3]);
        }
      else// if genotyping
        {
          // 00 ->00/00 
          // 01 -> 00/01 :
          // 02 -> 01/01 
          // 10 -> 10/00:
          // 20 -> 10/10
          // 11 -> 11/00 // 10/01  ignore
          // 12 -> 01/11 
          // 21 -> 10/11 
          // 22 -> 11/11 
          if (plist[k][sj]=='0' && plist[k][si]=='0' && (gam[0]==0 || (gam[0]==1 && k>=pop[0]) ))
            gam[0]= isgam(k,pop[0],gtest, gam[0]);
          if (plist[k][sj]=='0' && plist[k][si]=='2'&& (gam[1]==0 || (gam[1]==1 && k>=pop[0]) ))
            gam[1]= isgam(k,pop[0],gtest, gam[1]);
            
         
          if (plist[k][sj]=='2' && plist[k][si]=='0' && (gam[2]==0 || (gam[2]==1 && k>=pop[0]) ))
            gam[2]= isgam(k,pop[0],gtest, gam[2]);
          if (plist[k][sj]=='2' && plist[k][si]=='2' && (gam[3]==0|| (gam[3]==1 && k>=pop[0]) ))
            gam[3]= isgam(k,pop[0],gtest, gam[3]);

          if (plist[k][sj]=='0' && plist[k][si]=='1')
            {
              if(gam[1]==0 || (gam[1]==1 && k>=pop[0]) )
                gam[1]= isgam(k,pop[0],gtest, gam[1]);
              if(gam[0]==0 || (gam[0]==1 && k>=pop[0]) )
                gam[0]= isgam(k,pop[0],gtest, gam[0]);
            }
    
          if (plist[k][sj]=='2' && plist[k][si]=='1')
            {
              if(gam[2]==0 || (gam[2]==1 && k>=pop[0]) )
                gam[2]= isgam(k,pop[0],gtest, gam[2]);
              if(gam[3]==0 || (gam[3]==1 && k>=pop[0]) )
                gam[3]= isgam(k,pop[0],gtest, gam[3]);
        
            }
   
          if (plist[k][sj]=='1' && plist[k][si]=='0')
            {
              if(gam[2]==0 || (gam[2]==1 && k>=pop[0]) )
                gam[2]= isgam(k,pop[0],gtest, gam[2]);
              if(gam[0]==0 || (gam[0]==1 && k>=pop[0]) )
                gam[0]= isgam(k,pop[0],gtest, gam[0]);
            }
          if (plist[k][sj]=='1' && plist[k][si]=='2')
            {
              if(gam[1]==0 || (gam[1]==1 && k>=pop[0]) )
                gam[1]= isgam(k,pop[0],gtest, gam[1]);
              if(gam[3]==0 || (gam[3]==1 && k>=pop[0]) )
                gam[3]= isgam(k,pop[0],gtest, gam[3]);
            }
          k+=1;
        }// end hid haplotype

      if(gtest[0]==4)// rec in total
        {
          pgamete[sj][si]=1;
          if(gtest[1]==4 && gtest[2]<4)// pop1 only 
            {
              pgamete[sj][si]=2 ;// 2 if 1 rm pop1 only
            }
          else if(gtest[2]==4 && gtest[1]<4 )
            {
              pgamete[sj][si]=3 ;
              break;
            }
          else if(gtest[2]==4 && gtest[1]==4 )
            {
              pgamete[sj][si]=4 ;
              break;
            }
        }
    }// loop on nsam/*** end of Rm calc ***/
}// end four gamnete test


double * LDcalc(double p1, double p2, double p12) 
{
  double D=0, Dmax=0;

  D = p12 - p1*p2;

  if (D<0.)
    Dmax = max2 (- p1*p2, -(1.- p1)*(1.-p2));
  else
    Dmax = min2 ( p1*(1.-p2), (1.- p1)*p2);
 
  double *tp;
  if( ! (tp=(double*) calloc((unsigned)3, sizeof(double))))
    perror(" calloc error. tp");
  tp[0]=(double) D/Dmax;// mean D'
  tp[1]=(double)(D*D)/( p1*p2*(1.- p1)*(1.-p2));// mean r^2

  return(tp);
}// end LDcalc

/*** LD statistics ***/
double freq2 (int a, int b,int start, int n, char **list )
{
  int i;
  double j=0;
  
  for (j=0,i=start; i<n; ++i)
    {
      if (list[i][a]=='1' && list[i][b]=='1')
        ++j;
    }
  return ((double) j / (double) (n-start));
}


double max2 (double a, double b)
{
  return (a>b ? a : b);
}

double min2 (double a, double b)
{
  return (a<b ? a : b);
}

int hapcount (char **list, int start, int nsam, int first, int last)
{
  int a, b, c, hlist = 0;
  for (b=start; b<nsam; ++b) 
    {
      for (c=start; c<b; ++c) 
        {
          for (a=first; a<=last && (list[b][a]==list[c][a] || (list[b][a]=='?' || list[c][a]=='?')); ++a)
            ;
          if (a==(last + 1))
            break;
        }
      if (c==b) 
        hlist++;
    }
  return (hlist);
}

int min_rec (int x, int size, int start, int ** pgamete, int fl)
{  /* Calculate min # rec. events */
  int a, b=0, c, d, flag = 0;
  if (size<2 || start >= size)
    return (0);
  for (a=start; a<size; ++a) 
    {
      for (b=0; b<(size-a); ++b)
        {
          if(fl==1)// total
            {
              if (pgamete[x+b][x+b+a] >= 1) 
                {
                  flag = 1;
                  break;
                }
            }
          else
            {
              if (pgamete[x+b][x+b+a]==fl ||pgamete[x+b][x+b+a]==4) 
                {
                  flag = 1;
                  break;
                }
            }
        }
      if (flag==1)
        break;
    }
  if (a==size)
    return (0);
  else {
    c = min_rec (x, b+1, a+1, pgamete, fl);
    d = min_rec (x+b+a, (size-b-a), a,pgamete,fl );
    return (1+c+d);
  }
}// end min_rec


