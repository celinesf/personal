/*! \file stats_popC.c
  \brief Code for summary statistics calculations


*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "statsfun.h"
#include "stats_pop.h"

/*! \brief Main function for stats_pop */
int main(int argc,char *argv[])
{
  /*** C object and functions ***/        
  int nsam, j , i,howmany, segsites, count ,arg=0 ,length=9999, type=0,nreg=1; 
  char **list=NULL ;
  int *posit=NULL;
  double tp;
  /**** add for C *******/
  FILE *fopen(), *pfin=NULL ;
  int    popsize[2] ;
  double typestats[4] ;
  char line[100001],*tpch, tpl[100001],msfile[100];
  //char *names[14] = {"segsites","singletons","frequencies","fst","pi","thetaW", "thetaH",  "D", "H", "D_star","D_prime","r_square","nH","Rm"};// name for various object  

  /*** recover parameters: -f low up [filters] -h 1/2 (genotype) -a anc/maf ***/
  /* Init default values */
  typestats[2]=5;
  typestats[3]=10;
  typestats[1]=1;
  typestats[0]=0 ;
  /***/

  if( argc < 2 ){ 
    fprintf(stderr,"Too few command line arguments\n");
    usage_stats() ;
  }
  else 
    {
      arg=1;
      while( arg < argc ){
        if( argv[arg][0] != '-' ) { fprintf(stderr," argument should be -%s ?\n", argv[arg]); usage_stats();}
        switch ( argv[arg][1] ){
          case 'm' : 
            arg++;
            argcheck_stat( arg, argc, argv);
            strcpy(msfile, argv[arg++]) ;
            break;
          case 'n' : 
            arg++;
            argcheck_stat( arg, argc, argv);
            nreg= atoi (argv[arg++]) ;
            break;
          case 'f' : 
            arg++;
            argcheck_stat( arg, argc, argv);
            typestats[2] = atoi(  argv[arg++] );
            argcheck_stat( arg, argc, argv);
            typestats[3] = atoi(  argv[arg++] );
            break;
          case 'a' : 
            arg++;
            argcheck_stat( arg, argc, argv);
            strcpy(tpl, argv[arg++]) ;
            if(strcmp(tpl,"maf")==0)
              { typestats[0]=1; }
       
            break;
          case 'p' : 
            arg++;
            argcheck_stat( arg, argc, argv);
            typestats[1] = atoi(  argv[arg++] );
            break;
          case 'h' : 
            fprintf(stderr," option default\n");  usage_stats() ;
            break;
          default: fprintf(stderr," option default\n");  usage_stats() ;
        }
      }
      // pfparv=fopen(pars.tp.fparv,"r"); 
      pfin=fopen(msfile,"r"); 
      fgets( line, 100000, pfin); /* cmdline */
      if( typestats[0]==1)
        printf("UnKnown ancestral allele\n");
      else   printf("Known ancestral allele\n");
      if( typestats[1]==2)
        printf("UnPhased data\n");
      else   printf("Phased data\n");
      printf("Filter for pairwise LD: %lf %lf kb\n", typestats[2],typestats[3] );
      /*** end recover parameters ***/
      printf("r\tn_1\tn_2\tL\tS\tS_1\tS_2\tS_s\tS_ss\tS_sf1\tS_sf2\tS_sl\tS_sh\tS_f\tS_f1\tS_f2\tS_o\tS_o1\tS_o2\tF(S)\tF(S_1)\tF(S_2)\tF(S_s)\tF(S_ss)\tF(S_sf1)\tF(S_sf2)\tF(S_sl)\tF(S_sh)\tF_st\tH_b\tH_w1\tH_w2\tSnn\tpi\tpi_1\tpi_2\ttheta_W\ttheta_W1\ttheta_W2\ttheta_H\ttheta_H1\ttheta_H2\tD\tD_1\tD_2\tH\tH_1\tH_2\tD_star\tD_star1\tD_star2\tD_prime\tD_prime1\tD_prime2\tr_square\tr_square1\tr_square2\tnH\tnH_1\tnH_2\tRm\tRm_1\tRm_2");
   

      sscanf(line," %*s  %d %d",  &nsam, &howmany);
     
      i= 0;popsize[1]=0;  /* obtain the popsizes */
      popsize[0]=nsam;

      tpch=NULL;
      tpch = strtok (line,"////");
      if(strcmp(tpch,".")==0)
        {
          tpch = strtok (NULL," ");
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
        }// from msC
      else// from msR with multiple regions
        {
          // howmany=nsam;
          type++;
        }

      if(typestats[1]==2 && !(nsam%2==0 && popsize[0]%2==0 &&popsize[1]%2==0) && type==0)
        {
          printf("Problem: I don't have even number of pseudo-haplotypes. nsam=%d n1=%d n2=%d.\n",nsam,popsize[0],popsize[1]);
  
        }
      else
        {
          fgets( line, 1000, pfin);/* seeds or empty*/
          count=1;
          while( !feof(pfin))//howmany-count++) 
            {  
              if(type==1)
                {
                  while ( (line[0] != '/')  )      
                    if( fgets( line, 100000, pfin) == NULL ) exit(0);  
                  sscanf(line," %*s  %d %d",  &popsize[0], &popsize[1]);
                  nsam=popsize[0]+popsize[1];
                  tpch=NULL;
                  tpch = strtok (line," ");
                  
                  while(strcmp(tpch,"\n")!=0)
                    {
                      if(strcmp(tpch,"-r")==0)
                        {
                          tpch = strtok (NULL, " ");// rho
                          tpch = strtok (NULL, " ");// Z-1
                          length=atoi(tpch);
                        }
                      
                      tpch = strtok (NULL, " ");
                      if(tpch== NULL) break;
                    }
                 
                  if(typestats[1]==2 && !(nsam%2==0 && popsize[0]%2==0 &&popsize[1]%2==0) && type==0)
                    {
                      printf("Problem: I don't have even number of pseudo-haplotypes. nsam=%d n1=%d n2=%d.\n",nsam,popsize[0],popsize[1]);
  
                    }
                }
             
              if(count==1)
                {
                  list = cmatrix(nsam,maxsites+1);
                  posit = (int*)calloc((unsigned) maxsites+1,sizeof( int));
                /*   if(type==0) */
/*                     { */
/*                     /\*   printf("nsam=%d\thowmany=%d\n", nsam, howmany); *\/ */
/* /\*                       printf("n1: %d \t n2: %d\tlength: %d\n", popsize[0],popsize[1],length+1); *\/ */
/*                       printf("%d\t%d\t%d\t", popsize[0],popsize[1],length+1); */
/*                     } */
                }
              while ( (line[0] != 's')  )      
                if( fgets( line, 100000, pfin) == NULL )
                  exit(0);   
              segsites=0;
              if ( (line[0] == 's')  )
                sscanf( line, "  segsites: %d", &segsites );
              if( segsites >= maxsites)
                {
                  maxsites = segsites + 10 ;
                  posit = (int *)realloc(posit, maxsites*sizeof(int) ) ;
                  biggerlist_stat(nsam,maxsites, list) ;
                }
              //printf("segsites: %d\n",segsites);
              if( segsites > 0)
                {
                  /*** recover positions and list ***/
                  fgets( line, 1000000, pfin) ;  
                  tpch=NULL;
                  tpch = strtok (line,"\n");// divide stats line
                  strcpy(tpl,tpch);
                  tpch = strtok (tpl,"positions: ");
                  for(i=0;i<segsites;i++)
                    {
                      tp=atof(tpch);
                      if((int) tp==0)
                        posit[i]=(int) (tp*length);
                      else  posit[i]=(int) tp;
                      // printf("%d\t %d\n",i,posit[i]);
                      tpch = strtok (NULL, " ");
                      if(tpch== NULL) break;
                    }
                  for( i=0; i<nsam;i++)
                    { 
                      fscanf(pfin," %s", list[i] );
                      // printf("%s\n",list[i]); 
                    } /* haplo*/   
                  /***/
                  /*** stats pop calculation ***/
                  Sstat= getstatS(popsize,nsam,  segsites, list, typestats,posit);
        
                  FST=getFST(popsize,nsam,  segsites, list, typestats);
                  nn2=getnn2(popsize,nsam,  segsites, list,(int) typestats[1]);

                  printf("\n%d\t",count);
                  //if(type==1){
                    //printf("nsam=%d\thowmany=%d\n", nsam, howmany);
                  printf("%d\t%d\t%d\t", popsize[0],popsize[1],length+1);
                    //}
                  // printf("\n%s:\t",names[0]);
                  for (i=0; i<12; i++) /* get stats S*/
                    printf("%d\t",(int)Sstat[0][i] );
           
                  //printf("\n%s:\t",names[1]);
                  for (i=0; i<3; i++) /* get stats singletonsq*/
                    printf("%d\t",(int)Sstat[8][i]);
          
                  // printf("\n%s:\t",names[2] );
                  for (i=0; i<9; i++) /*  freq*/
                    {
                      double tp=0;
                      if(Sstat[0][i]>0)
                        tp=(double) Sstat[1][i]/(Sstat[0][i]);
                      printf("%lf\t",tp );
                    }

                  //printf("\n%s:\t",names[3]);/* fst */
                  for (i=0; i<4; i++) 
                    printf("%lf\t",FST[i]);
                  printf("%lg\t",nn2);
             
                  //printf("\n%s:\t",names[4]); /* pi */
                  for (i=0; i<3; i++) 
                    printf("%lf\t",Sstat[2][i]);
             
                  // printf("\n%s:\t",names[5]); /* thetaW */
                  for (i=0; i<3; i++) 
                    printf("%lf\t",Sstat[3][i]);

                  //printf("\n%s:\t",names[6]); /* thetaH */
                  for (i=0; i<3; i++) 
                    printf("%lf\t",Sstat[4][i]);

                  //printf("\n%s:\t",names[7]); /* D */
                  for (i=0; i<3; i++) 
                    printf("%lf\t",Sstat[5][i]);
          
                  // printf("\n%s:\t",names[8]); /* H */
                  for (i=0; i<3; i++) 
                    printf("%lf\t",-Sstat[6][i]);

                  //  printf("\n%s:\t",names[9]); /* D* */
                  for (i=0; i<3; i++) 
                    printf("%lf\t",Sstat[7][i]);
               
                  //  printf("\n%s:\t",names[10]); /* D'*/
                  for (i=0; i<3; i++) 
                    printf("%lf\t",Sstat[9][i]);

                  //  printf("\n%s:\t",names[11]); /* r2 */
                  for (i=0; i<3; i++) 
                    printf("%lf\t",Sstat[10][i]);
                  //  printf("\n%s:\t",names[12]); /* nH */
                  for (i=0; i<3; i++) 
                    printf("%lf\t",Sstat[11][i]);
          
                  // printf("\n%s:\t",names[13]); /* Rm */
                  for (i=0; i<3; i++) 
                    printf("%lf\t",Sstat[12][i]);

                  /*** Free memory ***/
                  for (i=0; i<13; ++i) free (Sstat[i]);
                  free (Sstat);
                  free (FST);
                }// end if segsites>0
              else
                {
                  fgets( line, 100000, pfin);/* /n */    
                  printf("\n%d\t%d\t",count,segsites);
                }// end if no seg sites
              //printf("\n");
              count++;
            }// end of loop on howmany
       
          for(i=0;i<nsam;i++)
            {
              free(list[i]);
              /*    free(haplist[i]); */
            }
          free(list);
          /*   free(haplist); */
          free(posit);
        }// end if h=2 but note even ch
     
    }
  /* else */
/*     { */
/*       printf("nofile\n"); */
/*       usage_stats(); */
/*     } */
  fclose(pfin); 
  return(0);
}// end main


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


int
usage_stats()
{
  fprintf(stderr,"usage: ./stats_pop -m msfile\n");
  fprintf(stderr,"  Options: \n");
  fprintf(stderr,"\t -m filename\n");
  fprintf(stderr,"\t -a anc/maf (un/known ancestral allele)\n");
  fprintf(stderr,"\t -f l u (lower and upper kb lim for filtering pairwise SNP to calculate LD stats. Dafault is 5 to 10kb)\n");
  fprintf(stderr,"\t -p 1/2 (1/2 for un/phased data)\n");
  fprintf(stderr,"\t -h for help [show this]\n");
  exit(1);
}                      

void
argcheck_stat( int arg, int argc, char *argv[] )
{
  if( (arg >= argc ) || ( argv[arg][0] == '-') ) {
    fprintf(stderr,"not enough arguments after %s\n", argv[arg-1] ) ;
    fprintf(stderr,"For usage type: ms(\"./ms\")<return>\n");
    exit(0);
  }
}
