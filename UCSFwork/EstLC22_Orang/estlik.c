/*! \file estlikC.c
  \brief Code for estimation program


*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "ms.h"
#include "streec.h"
#include "estlikC.h"
#include "rand1.h"
#include "statsfun.h"

/***/
/*! \brief main function

  - call getpars()
  - open files of region/stats info and file of parameters if not in cmdline
  - call get_stat()
  - call init_cov() & read_ireg()
  - call cmatrix()
  - loop on param values
   - call getparv()
   - call reinit_cov() 
   - loop on ARG
    - call getreg_rho() (once or every steps)
    - call get_hypo() & gensam_estL()
    - if(S>0)
     - Loop on R call getprobdata()
     - Call get_cov() and get_lik()
    - If at checkmax call get_tot_lik() and comp_lik() 
   - get_tot_lik() for this set of param
  - call free_pars(), close files and free other memories
 */
int main(int argc,char *argv[])
{
  clock_t startl;
  int i, howmany,count,segsites=0, *pos=NULL;
  char **list=NULL;
  long double LIK=0,MAX=0;
  /*** addition for  thinning ***/
  int nr, oks=1, okp=1,okmax=1;
  char *line=NULL;
  double *hypo=NULL; // 0 rho_n, 1 rho_x=omega, 2Z_(maxrhoR), 3 xvZ_t, 4 max_rho_r
  FILE  *pfparv=NULL, *pfreg=NULL, *fmax=NULL;
  /*** only for multi locus ***/
  /* End declaration */

 /*    for(i=0;i<argc;i++) // print cmd line // comment */
  /*       printf("%s ",argv[i]); /\*\////\*\/ */ /*     printf("\n"); */

  /*** recover cmd line ***/
  if(!(line = (char *) calloc((unsigned)1000,sizeof(char))))
    perror(" calloc error. line");
  getpars( argc, argv); // Init pars 

  if( !pars.commandlineseedflag ) seedit( "s", &pars); /* change celine */// for write seeds in file // 
  //   printf("\n%d %d %d\n",pars.tableseeds[0],pars.tableseeds[1],pars.tableseeds[2]);// print seeds // comment

  if( pars.tp.pf)
    pfparv=fopen(pars.tp.fparv,"r"); // File of parameters
  if( pfparv!=NULL ||  pars.tp.pf==0)  /*** Parameter values file exists ***/
    {
      pfreg=fopen(pars.tp.freg,"r"); // file of regions
      if( pfreg!=NULL )  /*** regions info file exists ***/
        { 
          oks=get_stat(pfreg);
          if(pars.sp.ns>0 && oks)
            {
              init_cov();  /* Init pstats and cov */
              oks=read_ireg(pfreg,&nr);/* recover info regions */
              /*** Info regions ok ***/
              if(oks)
                { 
                  /* Init matrices and vector for regions and simulated ARG & data */
                  pars.cp.nsam = pars.cp.config[0]+ pars.cp.config[1];
                  if(pars.tp.pf)
                    if(!(pars.tp.paramv= (double *)calloc( (unsigned)NPARAM,sizeof(double ))))
                      perror(" calloc error. paramn");
                  if(!( hypo = (double *) calloc((unsigned) 5 ,sizeof(double))))
                    perror(" calloc error. hypo");
                  if(pars.tp.pf)
                    fgets( line, 1000, pfparv);
                  else
                    sprintf(line,"%s","0");

                  howmany=pars.tp.howmany; 
                  segsites=pars.mp.segsitesin;
                  //  printf("%d\n",segsites);
                  if(segsites==0)
                    {
                      list = cmatrix(pars.cp.nsam,maxsites+1);
                      posit = (double *)calloc( (unsigned) maxsites,sizeof( double) );
                    }
                  else
                    {
                      list = cmatrix(pars.cp.nsam, segsites+1 ) ;// init list for the region nr
                      if(!(posit = (double *)calloc((unsigned) segsites,sizeof( double))))
                        perror(" calloc error. posit");  /*** only for multi locus ***/
                    }
                  while(line[0]!='\0' && line[0]!='\n') /*** Loop on sets of parameters ***/
                    {  
                      if(pars.tp.pf)
                        for(i=0; i<NPARAM;i++)/* Init */
                          pars.tp.paramv[i]=0;
          
                      getparv(line);  /*** recover param values ***/

                      if(pars.tp.rho[0]==1)
                        pars.tp.rho[3]=pars.tp.rho[0]; // (4Nec)_r given 
                      if(pars.tp.rho[0]==2)
                        pars.tp.rho[3]=(double)pars.tp.paramv[0]/pars.tp.rho[1];// c_r given-> theta/mu=4N_1[*c_r]       

                      /*** end recover region info */
                      startl = clock();
                      count=0; 
                      reinit_cov();  /* REInit pstats and cov */
                      for(i=0; i<5;i++)
                        {
                          pars.sp.prob[i]=0;
                          pars.sp.infoperf[i]=0;
                        } 
                      okmax=1;
                      while( howmany-count++ ) // loop on sets of parameters : ARG-stats-lik
                        {
                          /*** Calculate reg values if -r <0 Each estimation step ***/
                          if( (pars.tp.rho[0]<0 ) ||(count==1 && pars.tp.rho[3]>0))
                            {
                              getreg_rho( hypo,count);/*/////*/
                            }
                          /*** Generate Hypothetical ARG ***////
                          get_hypo(hypo); 
                          if(  pars.mp.theta>0)
                            segsites = gensam_estL( list, posit);
                     
                          /////////////////////      print_haplo(segsites,posit,list);
                          if(segsites>0)
                            {
                              okp=1;
                              pars.sp.prob[1]=pars.sp.infoperf[3]= pars.sp.infoperf[4]=0.0; // Init num region ok
                              for(nr=0;nr<pars.tp.nregions;nr++) /*** P(S_i|ARG,THETA) + make sample w/ theta_r*/
                                { 
                                  pars.tp.ireg[nr][17]=(double) pars.tp.ireg[nr][16]/hypo[1];// rho/omega
                                  pars.tp.ireg[nr][18]=(double) pars.tp.ireg[nr][13]/hypo[3];// xvZ/xvZ_t
                                  pars.tp.ireg[nr][19]=(double) hypo[0]/pars.tp.ireg[nr][16];// rho_bp_n/rho_bp
                                  pars.tp.ireg[nr][20]=(double) pars.tp.ireg[nr][18] * ((double) hypo[0]/pars.tp.ireg[nr][16]) ;// rho_bp_n/rho_bp * xvZ/xvZ_t

                                  /*printf("\n%d \t",nr); */
                                  /*                                           int i; */
                                  /*                                           for(i=0;i<25+pars.sp.ns;i++) */
                                  /*                                             printf("%d-%lg\t",i, pars.tp.ireg[nr][i]); */
                                  /*                                           printf("\n"); */

                                  /*** rescale  ARG and Caclulate Lstat ***/ 
                                  pars.mp.theta =pars.tp.ireg[nr][14];
                                  printf("nr %d theta %lf\n",nr, pars.mp.theta);
                                  if( okp )// no need calculate for other regions if one didn't fit'
                                    getprobdata(hypo, nr,&okp);
                                  /***/////***/
                                  if(!okp) 
                                    break;
                                }// end loop on regions
                              if(okp)
                                {
                                  if( pars.sp.nsp>0)/*** Dp/cov calculation ***/
                                    {
                                      if(!(pos= (int *)calloc( (unsigned)segsites+1,sizeof(int ))))
                                        perror(" calloc error. pos");
                                      for(i=0;i<segsites;i++)
                                        pos[i]=(int) (posit[i]*hypo[2]);
                                      get_cov( segsites, list, pos);
                                      free(pos);
                                    }
                                  get_liks(okp); /* Lik of S|G and THETA) */
                                }

                              if((count==pars.sp.checkmax) )
                                {
                                  okmax=1;
                                  if( (fmax=fopen(pars.sp.maxfile,"r"))!=NULL) // File of parameters
                                    {
                                      LIK=get_tot_lik( startl,  count,2);
                                      fclose(fmax);
                                    }
                                  if(!isnan(LIK))
                                    {
                                      MAX=comp_lik();
                                      if(-2*(LIK-MAX)>10) okmax=0;
                                    }
                                  else okmax=0;
                                }// end if time to check max lik
                            }// if segsite>0
                          if(!okmax) break;
                        }// end of loop on samples per one set of parameters
                      LIK=get_tot_lik( startl, count,1);
                      if(count-1 == howmany && (LIK>MAX || MAX==0 || isnan(MAX)))
                        {
                          if( (fmax=fopen(pars.sp.maxfile,"w"))!=NULL) // File of parameters
                            {
                              for(i=0; i<NPARAM;i++)
                                fprintf(fmax,"%lg\t",pars.tp.paramv[i]);
                              fprintf(fmax,"%Lg %Lg\t%Lg\t%Lg\t%lf\t%d\n", LIK-logl(pars.sp.prob[0]),logl(pars.sp.prob[0]),LIK,(long double) expl(LIK), (double)((double)clock() - (double)startl) / (double)CLOCKS_PER_SEC, count);
                              fclose(fmax);
                            }
                        }
                      fflush(stdout);
                      free_eventlist(pars.cp.deventlist, pars.cp.npop );  
                      pars.cp.deventlist = NULL ;
                      sprintf(line,"%s","\0");// Init
                      if(pars.tp.pf)
                        {
                          fgets( line, 1000, pfparv); // Read new set of parameter
                          if(feof(pfparv))break;
                        }   
            
                    }// end loop on sets of parameters
                  /* free scaled data sets */
                  /*   for(i=0;i<pars.cp.nsam;i++)  */
                  /*                     free(listr[i]); */
                  /*   free(listr); */
                  /*                   free(positr);  */
                 
                  free(posit); 
                  for(i=0;i<pars.cp.nsam;i++)
                    free(list[i]);
                  free(list);
                  /* free */ 
                  free(hypo);
                  free(pars.tp.paramv); 
                }// end if all region recombnining and ok
              else print_error(3,nr,0); 
              //////////////////////////////////      // if(nr==pars.tp.nregions) nr--;
              free_cov(pars.tp.nregions-1); /* Free */
            }// end problem in stat and region info
          else print_error(2,0,0); 
          fclose(pfreg);
        }// end if region file opens
      else print_error(1,0,0);
      if( pars.tp.pf && pfparv!=NULL) 
        fclose(pfparv);
    }// end if paramv file opens 
  else print_error(0,0,0);
  /**** Free added by celine ***/
  free(line);
  // free_eventlist(pars.cp.deventlist, pars.cp.npop );  
  free_pars(); /* FREE PARS */ 
  if( !pars.commandlineseedflag ) seedit( "end",&pars );/* change celine */
  return(0);
}// end main C


/***********************************************************/
/*************** CMD PARAMETERS FUNCTIONS ******************/
void getpars(int argc, char *argv[])
{
  int arg, i ,npop, nseeds=0;

  if( argc < 2 ){ 
    fprintf(stderr,"Too few command line arguments\n");
    usage_estL();
  }
  pars.tp.rho=NULL;
  if(!(  pars.tp.rho= (double *)calloc( (unsigned)4,sizeof(double ))))
    perror(" calloc error. rho");
  arg=1;
  pars.tp.rho[0]=atof(argv[arg]); // Rho 
  if( pars.tp.rho[0]== 0 ) { fprintf(stderr,"First argument error. rho = 0. \n"); usage_estL();}
  else if(pars.tp.rho[0]==-1.)  // c/mu~Exp lamba
    {
      arg++;
      argcheck_estL( arg, argc, argv);// lambda */
      pars.tp.rho[1]=atof( argv[arg]);
    }
  else if(pars.tp.rho[0]==-2.)// c/mu~Normal nu sigma
    {
      arg++;
      argcheck_estL( arg, argc, argv);// nu
      pars.tp.rho[1]=atof( argv[arg]);
      arg++;
      argcheck_estL( arg, argc, argv);// sigma 
      pars.tp.rho[2]=atof( argv[arg]);
    }
  else if(pars.tp.rho[0]==2.) // c_r fixed
    {
      arg++;
      argcheck_estL( arg, argc, argv);
      pars.tp.rho[1]=atof( argv[arg]);// mu
    }
  else if(pars.tp.rho[0]!=1.) // rho_r_bp fixed
    pars.tp.rho[3]=pars.tp.rho[0];

  printf("rho %lg %lg %lg %lg \n", pars.tp.rho[0], pars.tp.rho[1], pars.tp.rho[2], pars.tp.rho[3]);

  pars.cp.nsam =1;
  pars.tp.howmany=2;
  pars.commandlineseedflag = 0 ;
  if(!(pars.tableseeds = (int *) calloc( (unsigned)3 ,sizeof( int) )))
    perror(" calloc error. seed");
  /****/
  pars.cp.r = pars.mp.theta = pars.cp.f = 0.0 ;
  pars.cp.track_len = 0. ;
  pars.cp.npop = npop = 2 ;
  if(!(pars.cp.mig_mat = (double **)calloc( (unsigned) pars.cp.npop ,sizeof( double *) )))
    perror(" calloc error. mingmat");
  if(!(pars.cp.mig_mat[0] = (double *)calloc( (unsigned)pars.cp.npop,sizeof(double ))))
    perror(" calloc error. mmi");
  pars.cp.mig_mat[0][0] = 0.0 ;
  pars.mp.segsitesin = 0 ;
  pars.mp.treeflag = 0 ;
  pars.mp.timeflag = 0 ;
  pars.mp.mfreq = 1 ;
  if(!(pars.cp.config = (int *) calloc( (unsigned)  pars.cp.npop +1  ,sizeof( int))))
    perror(" calloc error. config0");
  (pars.cp.config)[0] = pars.cp.nsam ;
  if(!(pars.cp.size= (double *) calloc( (unsigned) pars.cp.npop ,sizeof( double ) )))
    perror(" calloc error. size0");
  (pars.cp.size)[0]= (pars.cp.size)[1] = 1.0 ; 
  if(!(pars.cp.alphag = (double *) calloc( (unsigned) pars.cp.npop  ,sizeof( double ))))
    perror(" calloc error. alpha0");
  (pars.cp.alphag)[0] = 0.0 ;
  pars.cp.nsites = 2 ;

  for(i=1; i<pars.cp.npop; i++) 
    { if(!(pars.cp.mig_mat[i] =(double *)calloc( (unsigned) pars.cp.npop,sizeof(double))))
        perror(" calloc error. mmati");
      pars.cp.mig_mat[0][i] = pars.cp.mig_mat[0][0];
      (pars.cp.size)[i] = (pars.cp.size)[0] ;
      (pars.cp.alphag)[i] = (pars.cp.alphag)[0] ;
      (pars.cp.config)[i] = (pars.cp.config)[0] ;
    }
  pars.cp.deventlist = NULL ;

  /* added celine */
  if(!(pars.sp.wstat=( int *) calloc((unsigned)59,sizeof(int))))
    perror(" calloc error. wstat");// which stats to calculate
  if(!(pars.sp.type=( double *) calloc((unsigned)4,sizeof(double))))
    perror(" calloc error. type");// which stats to calculate
  pars.sp.type[0]=0;// anc
  pars.sp.type[1]=1;// phase
  pars.sp.type[2]=5;// 5klb
  pars.sp.type[3]=10;// 10kb filter
  if(!( pars.sp.infoperf=(double *)calloc((unsigned)5 ,sizeof(double))))
    perror(" calloc error. info"); // info on perf
  if(!( pars.sp.prob=(long double *)calloc((unsigned)5 ,sizeof(long double))))
    perror(" calloc error. prob"); // prob_param, prog_gen, prog_reg
  if(!( pars.tp.fparv = (char *) calloc((unsigned)100 ,sizeof(char))))
    perror(" calloc error. fparv");
  if(!( pars.tp.freg = (char *) calloc((unsigned)100 ,sizeof(char))))
    perror(" calloc error. freg");
  if(!( pars.tp.floc = (char *) calloc((unsigned)100 ,sizeof(char))))
    perror(" calloc error. nfloc");
  if(!( pars.sp.maxfile= (char *) calloc((unsigned)100 ,sizeof(char))))
    perror(" calloc error. max");
  pars.sp.checkmax=20000;
  pars.sp.meanregion=0;
  pars.tp.pf=1;
  sprintf(pars.tp.fparv,"%s", "param_file");
  sprintf(pars.tp.freg,"%s", "info_regions");
  sprintf(pars.tp.floc,"%s", "info_loci");
  pars.tp.nregions = 1;
  sprintf(pars.sp.maxfile,"%s","maxparam");
  /*** recover list of possible statistics ***/
  pars.sp.liststat=stat_name();// init list of statistics

  arg ++ ;
  while( arg < argc ){
    if( argv[arg][0] != '-' ) { fprintf(stderr," argument should be -%s ?\n", argv[arg]); usage_estL();}
    switch ( argv[arg][1] ){
      case 'G':/* howmany ARG per sets of parameters */
        arg++;
        argcheck_estL( arg, argc, argv);
        pars.tp.howmany=atoi( argv[arg++]);
        break;
      case 'p' : 
        arg++;
        argcheck_estL( arg, argc, argv);
        if( argv[arg-1][2] == 'f' ){
          sprintf( pars.tp.fparv,"%s", argv[arg++] ); 
          arg--;
        }
        else {
          if(!(pars.tp.paramv= (double *)calloc( (unsigned)NPARAM,sizeof(double ))))
            perror(" calloc error. paramn");
          pars.tp.pf=0;
          for(i=0;i<NPARAM;i++)
            {
              argcheck_estL( arg+i, argc, argv);
              pars.tp.paramv[i]=atof(argv[arg+i]);         
            }
          arg += NPARAM-1;
        }

        sprintf( pars.tp.fparv,"%s", argv[arg++] ); 
        break;	
      case 'r' : 
        arg++;
        argcheck_estL( arg, argc, argv);
        sprintf( pars.tp.freg,"%s", argv[arg++] ); 
        break;	
      case 'l' : 
        arg++;
        argcheck_estL( arg, argc, argv);
        sprintf( pars.tp.floc,"%s", argv[arg++] ); 
        break;	
      case 'R' : 
        arg++;
        argcheck_estL( arg, argc, argv);
        pars.tp.nregions = atoi( argv[arg++] ); 
        break;	
       
      case 's' :
        arg++;
        argcheck_estL( arg, argc, argv);
        pars.commandlineseedflag = 1 ;
        for(i=0;i<3;i++)
          argcheck_estL( arg+i, argc, argv);
        nseeds = commandlineseed(argv+arg ,&pars );
        arg += nseeds ;
        break;
      case 'h' :
        usage_estL(); 
        break;
      case 'T' : 
        pars.mp.treeflag = 1 ;
        arg++;
        break;
      case 'o' : 
        arg++;
        argcheck_estL( arg, argc, argv);
        pars.sp.maxfile = argv[arg++]   ;
        break;
      case 'c' : 
        arg++;
        argcheck_estL( arg, argc, argv);
        pars.sp.checkmax = atoi(argv[arg++])   ;
        break;
      case 'm' : 
        arg++;
        argcheck_estL( arg, argc, argv);
        pars.sp.meanregion = atoi(argv[arg++])   ;
        if(pars.sp.meanregion!=0 && pars.sp.meanregion!=1 )
          {
            fprintf(stderr,"-m %d option default\n",pars.sp.meanregion); usage_estL() ;
          }
        break;
      default: fprintf(stderr," option default\n"); usage_estL() ;
    }
  }
}// end getpars

/*** Function to check for rho characterization (allowed -1 -2) ***/
void argcheck2(char * argv)
{
  if( ( argv[0] == '-')&& (argv[1]!= '2') &&( argv[1]!= '1')) {
    fprintf(stderr,"not enough arguments before %s\n", argv ) ;
    fprintf(stderr,"For usage type: ms(\"./ms\")<return>\n");
    exit(0);
  }
}

int
usage_estL()
{
  fprintf(stderr,"usage: ./estlikeC rho \n");
  fprintf(stderr, "---- The rec. rate info must be specified by either:\n");
  fprintf(stderr,"\t\t rho (rho=4N1c)>e-6 \n");
  fprintf(stderr,"\t\t 1 (w=rho_bp in info_regions) \n");
  fprintf(stderr,"\t\t 2 mu (w=rc in info_regions)\n");
  fprintf(stderr,"\t\t -1 lambda (c/mu~exp(1/lambda) \n");
  fprintf(stderr,"\t\t -2 nu sigma (c/mu~normal(nu,sigma) \n");
  fprintf(stderr," Required options: \n");
  fprintf(stderr,"\t -p theta1 tehat2 thetaA Ts Tc M Mc (list of parameters to estimate the likelihood)\n");
  fprintf(stderr,"\t\t -pf filename  (file name with the sets of parameters, filename=param_info by default)\n");
  fprintf(stderr,"\t -r filename  (file name with information on the regions, filename=info_regions by default)\n");

  fprintf(stderr," More options: \n");
  fprintf(stderr,"\t -G howmany (=2 by default. The number of simulated data sets per sets of parameters to estimate the (composit-likelihood))\n");
  fprintf(stderr,"\t -R nregions (R=1 by default)\n");
  fprintf(stderr,"\t -l filename  (file name with information on the loci for multi-locus regions, filename=info_loci by default)\n");
  fprintf(stderr,"\t -s s1 s2 s3 (seeds for ramdome generator)\n");
  fprintf(stderr,"\t -o filename (name of file with max likelihood at checkpoint step for comparison: -2(log(L)-Log(Max))<10)\n");
  fprintf(stderr,"\t -m m (default m=0: use mean of regions-specific information do define hypothtical region. m=1 to take into account region-specific information. )\n");
  fprintf(stderr,"\t -c checkpoint (default=20000 step)\n");
  fprintf(stderr,"\t -h (for usage information)\n");
  fprintf(stderr," See Rmspack-vignette.pdf for more detailed on the method, input files and estimated models.\n");
  exit(1);
}


/***********************************************************/
/*************** MODEL PARAMETERS FUNCTIONS ******************/
/*** Function to reconver the model parameters ***/
void getparv(char *line )
{ 
  /* Init */
  int i,j, ok=1;
  char *tp=NULL, tpl[1000];
  struct devent *ptemp=NULL , *pt=NULL ;
  /****/

  if(pars.tp.pf)
    {
      tp = strtok (line,"\n");// divide stats line
      strcpy(tpl,tp);
      tp = strtok (tpl," ");
      ok=0;
      while (tp != NULL)
        {
          for(i=0;i<NPARAM;i++)
            { //printf("-%s-\n",tp);
              argcheck2(tp);
              pars.tp.paramv[i]=atof(tp);   
              tp = strtok (NULL, " ");
              if(tp== NULL) break;
            }
          if(i<6) {fprintf(stderr,"PROBLEM: I recorded only %d parameter values instead of seven.. Run aborted\n",i ); exit(0);}
          else {ok=1; break;}
          tp = strtok (NULL, " ");
        }// end loop on word in line
    }
  /*** Record model parameters in ms.pars structure ***/
  if(!ok){fprintf(stderr,"PROBLEM: I did not record any parameter values. Run aborted\n"); exit(0);}
  else
    {
      if(pars.tp.paramv[1] !=pars.tp.paramv[0])//  -n N2/N1
        pars.cp.size[1] = (double) pars.tp.paramv[1]/pars.tp.paramv[0];
      //  printf("-n %lf\t", pars.cp.size[1] ); 
      for( i=0; i<pars.cp.npop; i++) // -M M_p
        for( j=0; j<pars.cp.npop; j++) pars.cp.mig_mat[i][j] = pars.tp.paramv[5]/(pars.cp.npop-1) ;

      for( i=0; i< pars.cp.npop; i++) pars.cp.mig_mat[i][i] = pars.tp.paramv[5] ;
      //   printf("-M %lf ", pars.cp.mig_mat[0][0] );

      if( pars.tp.paramv[3]>0) // -ej Ts 1 2
        {
          if(!( pt = (struct devent *)calloc((unsigned)1, sizeof( struct devent) ) ))
            perror(" calloc error. ptj");
          pt->detype = 'j' ;
          pt->time = pars.tp.paramv[3] ;
          pt->popi =  pars.cp.npop-1 ;
          pt->popj = 0 ;// T_s
          pt->nextde = NULL ;
          // printf("-e%c %lf %lf\t",pt->detype, pt->time, pt->paramv );
          if( pars.cp.deventlist == NULL )
            pars.cp.deventlist = pt ;
          else if ( pt->time < pars.cp.deventlist->time ) {
            ptemp = pars.cp.deventlist ;
            pars.cp.deventlist = pt ;
            pt->nextde = ptemp ;
          }
          else
            addtoelist( pt, pars.cp.deventlist ) ;  
          if(pars.tp.paramv[2] !=pars.tp.paramv[0])// -eN T_s NA/N1
            {
              if(!(pt = (struct devent *)calloc((unsigned)1, sizeof( struct devent) ) ))
                perror(" calloc error. ptN");
              pt->detype = 'N' ;
              pt->time = pars.tp.paramv[3] ; // T_s
              pt->paramv =(double) pars.tp.paramv[2]/pars.tp.paramv[0]; // NA/N1
              // printf("-e%c %lf %lf\t",pt->detype, pt->time, pt->paramv );
              pt->nextde = NULL ;
              if( pars.cp.deventlist == NULL )
                pars.cp.deventlist = pt;
              else if ( pt->time < pars.cp.deventlist->time ) {
                ptemp = pars.cp.deventlist ;
                pars.cp.deventlist = pt ;
                pt->nextde = ptemp ;
              }
              else addtoelist( pt, pars.cp.deventlist ) ;
            }
        }// -ej event
      else if(pars.tp.paramv[5] ==0)// check coalescent possible
        {
          fprintf(stderr,"PROBLEM: M_p=%lf and the time of split=%lf. This would lead to Infinite coalescent time. Run aborted\n",pars.tp.paramv[5],pars.tp.paramv[3] ); exit(0);
        }

      if( pars.tp.paramv[4]>0 && pars.tp.paramv[3]>0 && pars.tp.paramv[5]!=pars.tp.paramv[6]) // -eM Tc M_p
        {
          if(!(pt = (struct devent *)calloc((unsigned)1, sizeof( struct devent) )))
            perror(" calloc error. ptM");
          pt->detype = 'M' ;
          pt->time = pars.tp.paramv[4]*pars.tp.paramv[3] ; // T_c= tau*T_s
          pt->paramv =pars.tp.paramv[6]; // M_c
          pt->nextde = NULL ;// printf("-e%c %lf %lf ",pt->detype, pt->time, pt->paramv );
          if( pars.cp.deventlist == NULL )
            pars.cp.deventlist = pt ;
          else if ( pt->time < pars.cp.deventlist->time ) {
            ptemp = pars.cp.deventlist ;
            pars.cp.deventlist = pt ;
            pt->nextde = ptemp ;
          }
          else
            addtoelist( pt, pars.cp.deventlist ) ;
        }
    }// Record model parameters
}// end function getpvalues

int read_ireg(FILE * ppfreg, int * pnr)
{
  int nr,  sstat=0;
  char *dum=NULL, tp[10], line[100000], tpl[100000];
  int i, out=0, n1=0, n2=0;
  
  if(!( pars.tp.ireg = (double **) calloc((unsigned)pars.tp.nregions ,sizeof(double *))))
    perror(" calloc error. ireg");
  if(!( pars.tp.nsam = (int **) calloc((unsigned)pars.tp.nregions ,sizeof(int *))))
    perror(" calloc error. ñsam");
  if(!(pars.sp.lstat=(double **) calloc((unsigned)pars.tp.nregions,sizeof(double*))))
    perror(" calloc error. lstat"); // Array of branch lenght L1, L2, Ls, lf for a segment
  
  pars.tp.totsam=0;// # combination of sample size
  for(nr=0;nr<pars.tp.nregions && !feof(ppfreg);nr++)
    { 
      if(!( pars.tp.ireg[nr] = (double *) calloc( (unsigned) (25+pars.sp.ns) ,sizeof(double ) )))
        perror(" calloc error. ireg i");
      if(!( pars.tp.nsam[nr] = (int *) calloc( (unsigned) (2) ,sizeof(int ) )))
        perror(" calloc error. nsam i");
      if(!(pars.sp.lstat[nr]=(double *) calloc((unsigned) NLSTATS+1,sizeof(double))))
        perror(" calloc error. lstat"); // Array of branch lenght L1, L2, Ls, lf for a segment
      sprintf(line,"%s","\0");//Init
      fgets( line, 1000, ppfreg);
      dum = strtok (line,"\n");// divide stats line
      if(dum !=NULL)// recover info
        {
          strcpy(tpl,dum);
          dum = strtok (tpl," \t");// divide stats line
          sscanf(dum,"%lg",&pars.tp.ireg[nr][0]);
          dum = strtok (NULL, " \t");// (name of region)
          for(i=1; i<12; i++)
            {
              dum = strtok (NULL, " \t");
              if( dum == NULL) break;
              sscanf(dum,"%lg",&pars.tp.ireg[nr][i]);
              if( (i==5 ||i==4) &&!((int) (pars.tp.ireg[nr][i])%2 ==0))
                {fprintf(stderr,"PROBLEM: I expect even numbers of chromosomes but found n_pop%d=%d at region #%d\n",i-4,(int) pars.tp.ireg[nr][i],nr+1); usage_estL() ;}
              if(i==7 && pars.tp.ireg[nr][i]-pars.tp.ireg[nr][i-1]<=1)
                {fprintf(stderr,"PROBLEM: I expect that every region is >2bp long but found Z=%d at region #%d\n",(int) (pars.tp.ireg[nr][i]-pars.tp.ireg[nr][i-1]),nr+1); usage_estL() ;}
            }
          sstat=0;
          for(i=25; i<25+pars.sp.ns; i++)
            { 
              dum = strtok (NULL, " \t");
              if( dum == NULL) break;
              sscanf(dum,"%s\n",tp);
              if(strcasecmp(tp,"NA")==0)
                pars.tp.ireg[nr][i]=log(-1);
              else
                sscanf(dum,"%lg\n",&pars.tp.ireg[nr][i]);
              
              if(i<25+pars.sp.nss) sstat+=pars.tp.ireg[nr][i];
              if(i>=25+pars.sp.nss && !isnan(pars.tp.ireg[nr][i]))// mean Dp
                {
                  
                  pars.sp.mstat[i-25-pars.sp.nss][0]+=pars.tp.ireg[nr][i];
                  pars.sp.mstat[i-25-pars.sp.nss][1]++;
                }
            }
          if(!(dum!= NULL && pars.tp.ireg[nr][3]>0 && pars.tp.ireg[nr][1]*pars.tp.ireg[nr][2]*(pars.tp.ireg[nr][7]-pars.tp.ireg[nr][6])>0))// check rho & theta>0 & all info
            out++;
          if( pars.sp.meanregion)/* use variation in n1 n2*/
            {
              if((int) pars.tp.ireg[nr][4]> pars.cp.config[0])
                pars.cp.config[0]=(int) pars.tp.ireg[nr][4];
              if((int) pars.tp.ireg[nr][5]> pars.cp.config[pars.cp.npop-1])
                pars.cp.config[pars.cp.npop-1]=(int) pars.tp.ireg[nr][5];
              
              for(i=0;i<pars.tp.totsam;i++)/* combination nsam */
                {
                  if(pars.tp.nsam[i][0]==pars.tp.ireg[nr][4] &&  pars.tp.nsam[i][1]==pars.tp.ireg[nr][5])
                    {
                      pars.tp.ireg[nr][21]=i;
                      break;
                    }
                }
              if(i== pars.tp.totsam)
                {
                  pars.tp.nsam[i][0]=pars.tp.ireg[nr][4] ;
                  pars.tp.nsam[i][1]=pars.tp.ireg[nr][5];
                  pars.tp.ireg[nr][21]=i;
                  pars.tp.totsam++;
                }/***/
            }
          else /* take mean n1 n2 */
            {
              n1+=(int) pars.tp.ireg[nr][4];
              n2+=(int) pars.tp.ireg[nr][5];
            }

          /*/*  max segsites *\/ */
          /*              if(pars.tp.ireg[nr][11]>pars.mp.segsitesin) pars.mp.segsitesin=pars.tp.ireg[nr][11]; */
          pars.tp.ireg[nr][22]=pars.tp.ireg[nr][11]-sstat;
          
        }// end n>0
      else
        out++;
      if( pars.tp.ireg[nr][4]+ pars.tp.ireg[nr][5]<2) out++;
      if(out) break;
    }// end loop recover info reg
  if(! pars.sp.meanregion)
    {
      if((int)   round((double)n1/pars.tp.nregions) % 2!=0)
        {
          if( abs(round((double)n1/pars.tp.nregions)+1 - (double) n1/pars.tp.nregions)> abs(round((double)n1/pars.tp.nregions)-1 - (double) n1/pars.tp.nregions) )
            n1=(int)  (round((double)n1/pars.tp.nregions)-1) *pars.tp.nregions;
          else
            n1=(int)  (round((double)n1/pars.tp.nregions)+1) *pars.tp.nregions;
        }
      
      if((int)   round((double)n2/pars.tp.nregions) % 2!=0)
        {
          if( abs(round((double)n2/pars.tp.nregions)+1 - (double) n2/pars.tp.nregions)> abs(round((double)n2/pars.tp.nregions)-1 - (double) n2/pars.tp.nregions) )
            n2=(int)  (round((double)n2/pars.tp.nregions)-1) *pars.tp.nregions;
          else
            n2=(int)  (round((double)n2/pars.tp.nregions)+1) *pars.tp.nregions;
        } 
      pars.cp.config[0]=(int) round((double)n1/pars.tp.nregions) ;
      pars.cp.config[pars.cp.npop-1]=(int) round((double)n2/pars.tp.nregions) ;
      
      pars.tp.nsam[pars.tp.totsam][0]=     pars.cp.config[0] ;
      pars.tp.nsam[pars.tp.totsam][1]=pars.cp.config[pars.cp.npop-1];
      pars.tp.totsam++;
    }
  *pnr=nr;
  if(nr!=pars.tp.nregions)
    return(0) ; // Exit if problem with file of region 
  else
    return(1);
}// end function read)ireg



/*** Function to calculate reg specific info and find max/min rho, theta, Z ***/
void getreg_rho( double *phypo , int count) 
{
  /* Init */
  int nr;
  if(pars.sp.meanregion || (!pars.sp.meanregion && count==1)) /* init */
    {
      phypo[0]=10000;
      phypo[2]=phypo[1]=phypo[3]=phypo[4]=0.;
      /*** Get rho for each region ***/
      for(nr=0;nr<pars.tp.nregions;nr++) 
        {    
          /*** check on rho param and scalar compatibilities ***/////////// change for small R
          if((((pars.tp.rho[0]>0 && pars.tp.rho[0]<1) || pars.tp.rho[0]==-1 || pars.tp.rho[0]==-2) && pars.tp.ireg[nr][3]>0 && (pars.tp.ireg[nr][3]<.1|| pars.tp.ireg[nr][3]>10)) || (pars.tp.rho[0]==1 && pars.tp.ireg[nr][3]>0 && (pars.tp.ireg[nr][3]>.1||pars.tp.ireg[nr][3]<1e-6)) || (pars.tp.rho[0]==2 && pars.tp.ireg[nr][3]>0 && pars.tp.ireg[nr][3]>1e-6))// rho=4Nc
            { 
              fprintf(stderr,"PROBLEM: With the recombination rate scalar: w_%d=%lg while rho is %lg %lg %lg.  Run aborted\n",nr, pars.tp.ireg[nr][3], pars.tp.rho[0],pars.tp.rho[1],pars.tp.rho[2]); 
              exit(0);
            }
          /////////////////////////////////////// CHANGE??? ASK JEFF ??? ///////////////////////////
          /*** sample rho value ***/
          if(pars.tp.rho[0]==-1 && pars.sp.meanregion) // c/mu~Exp lamba
            pars.tp.rho[3]=randexp(pars.tp.rho[1])*pars.tp.paramv[0];
          if(pars.tp.rho[0]==-2 && pars.sp.meanregion) // c/mu~Normal nu sigma
            {
              pars.tp.rho[3]=RNormal(pars.tp.rho[1],pars.tp.rho[2])*pars.tp.paramv[0];
              while(pars.tp.rho[3]<=0)
                pars.tp.rho[3]=RNormal(pars.tp.rho[1],pars.tp.rho[2])*pars.tp.paramv[0];
            }
          /********/
        
          pars.tp.ireg[nr][12]=(double)pars.tp.ireg[nr][7]-pars.tp.ireg[nr][6]; // Z
          pars.tp.ireg[nr][13]=(double)pars.tp.ireg[nr][1]*pars.tp.ireg[nr][2]*pars.tp.ireg[nr][12]; // xvZ
          pars.tp.ireg[nr][14]=(double)pars.tp.paramv[0]*pars.tp.ireg[nr][13]; // theta= theta_1 * xvZ 
          pars.tp.ireg[nr][15]=(double)pars.tp.rho[3]*pars.tp.ireg[nr][3]*(pars.tp.ireg[nr][12]-1); // rho_r=rho* w *Z-1
          pars.tp.ireg[nr][16]=(double)pars.tp.ireg[nr][15]/(pars.tp.ireg[nr][12]-1); // rho_r/Z-1
      
          if(pars.sp.meanregion )/* if consider region-specific characteristics */
            {
              if(pars.tp.ireg[nr][16]<phypo[0]) // min rho_bp-rho_n
                phypo[0]=pars.tp.ireg[nr][16];
              if(pars.tp.ireg[nr][16]>phypo[1]) // max rho_bp=omega
                phypo[1]=pars.tp.ireg[nr][16];
              if(pars.tp.ireg[nr][15]>phypo[4]) // max rho_r
                { phypo[4]=pars.tp.ireg[nr][15];// max rho* w *Z-1
                  phypo[2]=pars.tp.ireg[nr][12];// Z of max rho* w *Z-1
                }
              if(pars.tp.ireg[nr][13]>phypo[3]) // max Zvx
                phypo[3]=pars.tp.ireg[nr][13]; 
            }
          else /* hypo is mean of characteristics */
            {
              phypo[1]+=pars.tp.ireg[nr][16];  printf("%lf\n",phypo[1]);
              phypo[4]+=pars.tp.ireg[nr][15];// max rho* w *Z-1
              phypo[2]+=pars.tp.ireg[nr][12];// Z of max rho* w *Z-1
              phypo[3]+=pars.tp.ireg[nr][13];
            }
        }// end loop obtain info reg
      /*** check on c/mu ***///////// change ??
      printf("%lg before mean %lf %lf %lf %lf %lf\n",  pars.tp.rho[3],phypo[0],phypo[1],phypo[2],phypo[3],phypo[4]);
      if(!pars.sp.meanregion )/* mean region characteristics */
        {
          phypo[0]=phypo[1]=(double) phypo[1]/pars.tp.nregions;
          phypo[4]=(double)phypo[4]/pars.tp.nregions;// max rho* w *Z-1
          phypo[2]=round((double)phypo[2]/pars.tp.nregions);// Z of max rho* w *Z-1
          phypo[3]=round((double)phypo[3]/pars.tp.nregions);
        }
      printf("after mean %lf %lf %lf %lf %lf\n",phypo[0],phypo[1],phypo[2],phypo[3],phypo[4]);
    }
  else /* reinit & no need to go in loop */
    {
      phypo[0]=phypo[1]=(double) phypo[1]/pars.tp.rho[3];
      phypo[4]=(double)phypo[4]/pars.tp.rho[3];// max rho* w *Z-1
    
      if(pars.tp.rho[0]==-1 && !pars.sp.meanregion) // c/mu~Exp lamba
        pars.tp.rho[3]=randexp(pars.tp.rho[1])*pars.tp.paramv[0];
      if(pars.tp.rho[0]==-2 && !pars.sp.meanregion) // c/mu~Normal nu sigma
        {
          pars.tp.rho[3]=RNormal(pars.tp.rho[1],pars.tp.rho[2])*pars.tp.paramv[0];
          while(pars.tp.rho[3]<=0)
            pars.tp.rho[3]=RNormal(pars.tp.rho[1],pars.tp.rho[2])*pars.tp.paramv[0];
        }
      phypo[0]=phypo[1]=(double) phypo[1]*pars.tp.rho[3];
      phypo[4]=(double)phypo[4]*pars.tp.rho[3];// max rho* w *Z-1
    }
  if(phypo[1]/pars.tp.paramv[0]>200)
    { 
      fprintf(stderr,"PROBLEM: With the recombination rate scalar: the ratio c/mu is %lg. Run aborted\n", (double) phypo[1]/pars.tp.paramv[0]); 
      exit(0);
    }/////
 
}// End function getreg_rho



/***********************************************************/
/****************** STATS CONSIDERED *********************/
char ** stat_name()
{
  int i;
  char **liststat=NULL, *dum=NULL;
  char lstat[1000];

  strcpy(lstat,"S S_1 S_2 S_s S_ss S_sf1 S_sf2 S_sl S_sh S_f S_f1 S_f2 S_o S_o1 S_o2 F(S) F(S_1) F(S_2) F(S_s) F(S_ss) F(S_sf1) F(S_sf2) F(S_sl) F(S_sh) F_st H_b H_w1 H_w2 Snn pi pi_1 pi_2 theta_W theta_W1 theta_W2 theta_H theta_H1 theta_H2 D D_1 D_2 H H_1 H_2 D_star D_star1 D_star2 D_prime D_prime1 D_prime2 r_square r_square1 r_square2 nH nH_1 nH_2 Rm Rm_1 Rm_2");
  if(!(liststat=( char **) calloc((unsigned)59,sizeof(char*))))
    perror(" calloc error. liststat");
  dum = strtok (lstat," ");
  for(i=0;i<59;i++)
    {
      if(!(liststat[i]=( char *) calloc((unsigned)10,sizeof(char))))
        perror(" calloc error. lstati");
      strcpy(liststat[i],dum);
      dum = strtok (NULL," ");
    } 
  return(liststat);
}

int get_stat( FILE * ppfreg)
{
  int i,j, ok=1;
  char  *dum, *line, tpl[1000];
  dum=line=NULL;
  if(!(line = (char *) calloc((unsigned)1000 ,sizeof(char))))
    perror(" calloc error. line2");
  /*** recover statistics to consider ***///////////////////////////////////////////
  fgets( line, 1000, ppfreg); // $type
  if(strncasecmp(line,"$type",4)==0)
    { 
      for(i=0;i<3;i++)
        {
          fgets( line, 1000, ppfreg); // type
           strcpy(tpl,line);
          dum = strtok (tpl,"\t");
          if( strncasecmp(dum,"r\n",1)!=0)
            {
              if(strncasecmp(line,"unknown\n",3)==0) // type defaul =0 anc
                pars.sp.type[0]=1;
              // fgets( line, 1000, ppfreg);
              if(strncasecmp(line,"unphased\n",3)==0) // type defaul =0 anc
                pars.sp.type[1]=2;
              //fgets( line, 1000, ppfreg);
              if(strncasecmp(line,"Filter\n",3)==0) // type defaul =0 anc
                {
                  dum = strtok (line,"\n");
                  strcpy(tpl,dum);
                  dum = strtok (tpl," ");
                  dum = strtok (NULL, " ");
                  i=1;
                  while( dum != NULL && i!=5)
                    {
                      if(i==1 && strcasecmp(dum,"for")!=0) ok=0;
                      if(i==2 && strcasecmp(dum,"pairwise")!=0) ok=0;
                      if(i==3 && strcasecmp(dum,"LD")!=0) ok=0;
                      if(i==4 && strcasecmp(dum,"between")!=0) ok=0;
                      dum = strtok (NULL, " ");
                      if(ok==0) break;
                      i++;
                    }
                  if(i==5 && ok==1)
                    {
                      pars.sp.type[2]=atof(dum);
                      dum = strtok (NULL, " ");
                      pars.sp.type[3]=atof(dum);
                    }
                  if(ok==0) ok=1;
                }
            }
          else break;
          // fgets( line, 1000, ppfreg);
        }
     
    }// type specified 
  if(strncasecmp(dum,"r\n",1)!=0)
    fgets( line, 1000, ppfreg);
  dum = strtok (line,"\n");// divide stats line
  strcpy(tpl,dum);
  dum = strtok (tpl,"//");// divide stats line
  dum = strtok (NULL,"//");// take the statistics
  if(dum != NULL)
    {
      strcpy(line,dum);
      pars.sp.nsp= pars.sp.nss= pars.sp.ns=0;
      dum = strtok (line," \t");
      j=0;
      while( dum != NULL )
        {
          for(i=j; i<59;i++)
            {
              if(strcasecmp(dum,pars.sp.liststat[i])==0)
                {
                  if(! (i>=53 && i<=55 && pars.sp.type[1]==2))// nh
                    {
                      pars.sp.wstat[i]=25+pars.sp.ns;// position in ireg
                      pars.sp.ns++;
                      if(i<15) pars.sp.nss++;// Sstats
                      else 
                        pars.sp.nsp++;// other stats
                    }
                  else  {fprintf(stderr,"PROBLEM: I have unphased data so can't calculate the statistics %s.\n",pars.sp.liststat[i]); usage_estL() ;}
                  j=i;
                  break;
                }
            }
          /*** Check independence S ***///////////////  comment for user 
          if((i>0 && i<NLSTATS &&  pars.sp.wstat[i]>0 &&  pars.sp.wstat[0]>0) /*S- S*/
             || (i>=12 && i<14 && pars.sp.wstat[i]>0 &&  pars.sp.wstat[1]>0) /*S1- So*/
             ||((i==12 || i==14) && pars.sp.wstat[i]>0 &&  pars.sp.wstat[2]>0) /*S2- So*/
             ||(i>3 && i<=8 && pars.sp.wstat[i]>0 &&  pars.sp.wstat[3]>0) /*Ss- Ss_*/
             ||((i==7 || i==8) && pars.sp.wstat[i]>0 &&  (pars.sp.wstat[4]>0||pars.sp.wstat[5]>0||pars.sp.wstat[6]>0)) /*Ss Ssf- Shl*/
             ||(i>9 && i<=11 && pars.sp.wstat[i]>0 &&  pars.sp.wstat[9]>0) /*Sf- Sf_*/
             ||(i>12 && i<=14 && pars.sp.wstat[i]>0 &&  pars.sp.wstat[12]>0)) /*So- So_*/
            {ok=0;  break;  }
          /*/////*/
          if(i==59)
            { 
              printf("The statistics %s doesn't exist or is out of order. Run aborted\n", dum) ;
              ok=0;
              break;
            }
          dum = strtok (NULL, " \t");
          if( dum == NULL) break;
        }
      /***/
    }// end if S stat specified
  else ok=0;
  free(line);
  return(ok);
}// end get_stat

/***********************************************************/
/****************** ARG AND SAMPLE SIMULATION *********************/
void get_hypo(double *phypo)
{
  /*  pars.mp.theta =0;  */

  pars.mp.theta = (double) ((double)phypo[1]/phypo[0])*phypo[3]*pars.tp.paramv[0]; // theta_1*max(xvZ)*omega/rho_n
  pars.cp.r=phypo[1]*(phypo[2]-1); // rho_hypo=Omega*maxZ-1 
  pars.cp.nsites = (int) phypo[2]; // Z= Z for max(rho_r)
  printf("q %lf r %lf L %d n %d /p0 %lf p1 %lf p2 %lf p3 %lf p4 %lf\n",pars.mp.theta,pars.cp.r,pars.cp.nsites, pars.cp.nsam, phypo[0], phypo[1], phypo[2], phypo[3], phypo[4]);

}

int gensam_estL( char **list, double * pposit) 
{
  int nsegs=1, i, j,k, seg, ns, start, end, len , segsit;
  struct segl *seglst=NULL ; 
  double nsinv, tseg, tt, theta;
  double *pk=NULL;
  int *ss=NULL;
  int segsitesin,nsites;
  int nsam ;
  ns=0;
  nsites = pars.cp.nsites ;
  nsinv = 1./nsites;
  seglst =segtre_mig(&(pars.cp), &nsegs ) ;
  nsam = pars.cp.nsam;
  /* added for estimation */
  int nsam1=0, nsam2=0;
  double *l=NULL;
  nsam1=pars.cp.config[0];                            // nsam1 and 2: # seq in pop 1 and 2
  nsam2=pars.cp.config[1];

  segsitesin = pars.mp.segsitesin ;
  theta = pars.mp.theta ;
   printf("%d %lf\n",segsitesin,theta);
  if( pars.mp.treeflag ) {
    ns = 0; 
    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
      if( (pars.cp.r > 0.0 ) || (pars.cp.f > 0.0) ){
        end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
        start = seglst[seg].beg ;
        len = end - start + 1 ;
        fprintf(stdout,"[%d]", len);
      }
      prtree( seglst[seg].ptree, nsam ) ;
    }
  }
  printf("am here -1 %d %lf\n",segsitesin,theta);
  if( (segsitesin == 0) && ( theta > 0.0) ) {
    ns = 0 ;
    printf("am here 0\n");
    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++)
      {
        end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
        start = seglst[seg].beg ;
        len = end - start + 1 ;
        /*** lstats ***/
        if(pars.sp.nss>0)
          {
            for(j=0; j<pars.tp.totsam;j++)
              {
                l=computelstat(pars.tp.nsam[j][0],pars.tp.nsam[j][1], seglst[seg].ptree);      // Get L statistics for this segment
                for(i=0;i<NLSTATS;i++)  //--- Compute total L_k: mean weighted on segment lenghts ---//
                  {
                    if(k==0)// reinit at first segment
                      pars.sp.lstat[j][i]=0;
                    pars.sp.lstat[j][i]+=(double)(l[i]*len)/nsites;
                    //printf("%d %d %lg\n", j, i,l[i] );
                  }
                free(l);
              }
          }    printf("am here 1\n");
        tseg = len*(theta/nsites) ;
        tt = ttime(seglst[seg].ptree, nsam);
        segsit = poisso( tseg*tt );
        if(segsit+ ns<=pars.cp.nsites)
          { 
            if( (segsit + ns) >= maxsites ) {
              maxsites =pars.cp.nsites + SITESINC ;
              posit = (double *)realloc(posit, maxsites*sizeof(double) ) ;
              biggerlist_estL(nsam, list) ;
            }
            make_gametes(nsam,seg,seglst[seg].ptree,tt, segsit, ns, list );
            free(seglst[seg].ptree); printf("am here free\n");
            locate(segsit,start*nsinv, len*nsinv,posit+ns); // find positions for each site
          }
        else{
          ns += segsit;
          break;
        }
        ns += segsit;
      }
  }
  else if( segsitesin > 0 )
    {
      ss = (int *)calloc((unsigned)nsegs,sizeof(int));
      pk = (double *)calloc((unsigned)nsegs,sizeof(double));
      if( (pk==NULL) || (ss==NULL) ) perror("malloc error. gensam.2");
      tt=0;
      for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++)
        {
          end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
          start = seglst[seg].beg ;
          len = end - start + 1 ;
          tseg = len/(double)nsites ;
          /*** lstats ***/
          if(pars.sp.nss>0)
            {
              for(j=0; j<pars.tp.totsam;j++)
                {
                  printf("%d %d %d\n",j,pars.tp.nsam[j][0],pars.tp.nsam[j][1]);
                  l=computelstat(pars.tp.nsam[j][0],pars.tp.nsam[j][1], seglst[seg].ptree);      // Get L statistics for this segment
                  for(i=0;i<NLSTATS;i++)  //--- Compute total L_k: mean weighted on segment lenghts ---//
                    {
                      if(k==0)// reinit at first segment
                        pars.sp.lstat[j][i]=0;
                      pars.sp.lstat[j][i]+=(double)(l[i]*len)/nsites;
                      //printf("%d %d %lg\n", j, i,l[i] );
                    }
                  free(l);
                }
            }
          l=computelstat(nsam1,nsam2, seglst[seg].ptree);      // Get L statistics for total ARG
          /*/////*/
          pk[k]=(double) l[0]*tseg ;
          free(l);
          tt += pk[k] ;
        }
      if( tt > 0.0 ) mnmial2(segsitesin,nsegs,pk,ss,tt);
      else
        for( k=0; k<nsegs; k++) ss[k] = 0 ;
      ns = 0 ;
    
      for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++)
        {
          end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
          start = seglst[seg].beg ;
          len = end - start + 1 ;
          tseg = len/(double)nsites;
          make_gametes(nsam,1,seglst[seg].ptree,(double) pk[k]/tseg, ss[k], ns, list);
          free(seglst[seg].ptree) ;
          locate(ss[k],start*nsinv, len*nsinv,pposit+ns);
          ns += ss[k] ;
        }
      free(pk);
      free(ss);
    }
  return( ns ) ;
}

/***********************************************************/
/****************** COVARIANCE MATRIX *********************/
void init_cov()
{
  int i;
  if(pars.sp.nsp>0)
    {
      if(!(pars.sp.pstat= (double **)calloc( (unsigned) pars.sp.nsp,sizeof( double *))))
        perror(" calloc error. pstat");
      if(!( pars.sp.cov= (long double **)calloc( (unsigned) pars.sp.nsp,sizeof(long double *))))
        perror(" calloc error. cov"); // covariance matrix
      if(!( pars.sp.icov= (long double **)calloc( (unsigned) pars.sp.nsp,sizeof(long  double *))))
        perror(" calloc error. icov"); // covariance matrix
      if(!( pars.sp.ncov= (int **)calloc( (unsigned) pars.sp.nsp,sizeof( int*))))
        perror(" calloc error. ncov");// n ok cov
      if(!(pars.sp.mstat= (long double **)calloc( (unsigned) pars.sp.nsp,sizeof(long double *))))
        perror(" calloc error. mstat");// mean obserc=ved Dp
      for(i=0;i<pars.sp.nsp;i++)
        {
          if(!( pars.sp.pstat[i]= (double *)calloc( (unsigned) 5+pars.sp.nsp*2,sizeof( double) )))
            perror(" calloc error. pstati");
          if(!(pars.sp.cov[i]= (long double *)calloc( (unsigned) pars.sp.nsp,sizeof(long double))))
            perror(" calloc error. covi");
          if(!(pars.sp.icov[i]= (long double *)calloc( (unsigned) pars.sp.nsp,sizeof(long double))))
            perror(" calloc error. icovi");
          if(!(pars.sp.ncov[i]= (int *)calloc( (unsigned) pars.sp.nsp,sizeof( int))))
            perror(" calloc error. ncov");
          if(!(pars.sp.mstat[i]= (long double *)calloc( (unsigned) 2,sizeof(long double) )))
            perror(" calloc error. mstats");
        }
    }
}

void get_cov(int segr, char ** plistr, int *pos)
{
  int i=0,j=0,k=0,ns=0;
  double **Sstat=NULL,  nn2=0;/* vector of seg stat (S s1 s2 ss Sss(shared by both pop), ssf1 ssf2,sf, sf1 sf2) + freq +pi +thetaW +thetaH +H + tajD + thetaFK + D_star*/
  /*** stats pop calculation ***/
  Sstat= getstatS2(pars.cp.config,pars.cp.nsam,  segr, plistr, pars.sp.type, pars.sp.wstat, pos);
  ns=j=0; 

  /**************** OUTPUT TO MAKE SURE SAME CALCUL AS Sc.C*********/
 /*  k=0; */
/*   for(i=0; i<12; i++) */
/*     { printf("%d %d %s %lg\n",k,i,pars.sp.liststat[k], Sstat[0][i] ); */
/*       k++; */
/*     } */
/*   for (i=0; i<3; i++) /\* get stats singletonsq*\/ */
/*     { */
/*       printf("%d %d %s %lg\n",k,i,pars.sp.liststat[k], Sstat[9][i] ); */
/*       k++; */
/*     } */

/*   for (i=0; i<9; i++) /\*  freq*\/ */
/*     { */
/*       double toto=(double) Sstat[1][i]/(Sstat[0][i]); */
/*         printf("%d %d %s %lg %lg %lg\n",k,i,pars.sp.liststat[k], Sstat[0][i],Sstat[1][i] , toto); */
/*       k++; */
/*     } */


  k=NLSTATS;/*  freq*/ 
  if( pars.sp.wstat[k]>0 ||pars.sp.wstat[k+1]>0||pars.sp.wstat[k+2]>0||pars.sp.wstat[k+3]>0||pars.sp.wstat[k+4]>0||pars.sp.wstat[k+5]>0||pars.sp.wstat[k+6]>0||pars.sp.wstat[k+7]>0||pars.sp.wstat[k+8]>0 )
    for(i=0; i<9; i++) 
      {
        if( pars.sp.wstat[k]>0)
          {
            //  printf("%d %s %lg %lg\n",i,pars.sp.liststat[k-NLSTATS], Sstat[0][i],Sstat[1][i] );
            if(Sstat[0][i]>0)
              {  
                Sstat[j][i]=(double) Sstat[1][i]/(Sstat[0][i]);
                pars.sp.pstat[ns][1]=  1;
                pars.sp.pstat[ns][2]++;
              }
            else 
              {
                Sstat[j][i]=0; 
                pars.sp.pstat[ns][1]=  0;
              } 
            pars.sp.pstat[ns][3]+=  Sstat[j][i];
            pars.sp.pstat[ns][4]+=  Sstat[j][i]* Sstat[j][i];
            pars.sp.pstat[ns][0]=  Sstat[j][i];
            ns++;
          }
        k++;
      } // end freq
  j=2;/* fst */
  k=NLSTATS+9; 
  if( pars.sp.wstat[k]>0 ||pars.sp.wstat[k+1]>0||pars.sp.wstat[k+2]>0||pars.sp.wstat[k+3]>0)
    for(i=0; i<4; i++) 
      {
        if( pars.sp.wstat[k]>0)
          {//printf("\n%d %d %d %d %s %lg\n",k,ns,i,j,pars.sp.liststat[k],Sstat[j][i] );
            fill_pstat(ns,Sstat[j][i]);
            ns++;
          }
        k++;
      }// end fst
  k=NLSTATS+9+4;/* nn2 */
  if( pars.sp.wstat[k]>0)
    { 
      nn2=getnn2(pars.cp.config,pars.cp.nsam,  segr, plistr, (int)  pars.sp.type[1] );
      //printf("\nhere %d %d %d %d %s %lg \n",k,ns,i,j,pars.sp.liststat[k],nn2 );
      Sstat[j][i]=nn2;
      fill_pstat(ns,Sstat[j][i]);
      ns++;
    }

  k=NLSTATS+9+4+1+3;/* ThetaW */
  j=4;
  if( pars.sp.wstat[k]>0 ||pars.sp.wstat[k+1]>0||pars.sp.wstat[k+2]>0)
    for (i=0; i<3; i++)
      {
        if( pars.sp.wstat[k]>0)
          {//printf("H3 %d %d %s %lg\n",k,ns,pars.sp.liststat[k],Sstat[j][i] );
            fill_pstat(ns,Sstat[j][i]);
            ns++;
          }
        k++;
      }
  k=NLSTATS+9+4+1+3+3;/* ThetaH */
  j=5;
  if( pars.sp.wstat[k]>0 ||pars.sp.wstat[k+1]>0||pars.sp.wstat[k+2]>0)
    for (i=0; i<3; i++)
      {
        if( pars.sp.wstat[k]>0)
          {//printf("H4 %d %d %s %lg\n",k,ns,pars.sp.liststat[k],Sstat[j][i] );
            fill_pstat(ns,Sstat[j][i]);
            ns++;
          }
        k++;
      }
  k=NLSTATS+9+4+1+3+3+3;/* D */
  j=6;
  if( pars.sp.wstat[k]>0 ||pars.sp.wstat[k+1]>0||pars.sp.wstat[k+2]>0)
    for (i=0; i<3; i++)
      {
        if( pars.sp.wstat[k]>0)
          {//printf("%d %d %s %lg\n",k,ns,pars.sp.liststat[k],Sstat[j][i] );
            fill_pstat(ns,Sstat[j][i]);
            ns++;
          }
        k++;
      }
  k=NLSTATS+9+4+1+3+3+3+3;/* H */
  j=7;
  if( pars.sp.wstat[k]>0 ||pars.sp.wstat[k+1]>0||pars.sp.wstat[k+2]>0)
    for (i=0; i<3; i++)
      {
        Sstat[j][i]=- Sstat[j][i];
        if( pars.sp.wstat[k]>0)
          {//printf("%d %d %s %lg\n",k,ns,pars.sp.liststat[k],Sstat[j][i] );
            fill_pstat(ns,Sstat[j][i]);
            ns++;
          }
        k++;
      }
  k=NLSTATS+9+4+1+3+3+3+3+3;/* Dstar */
  j=8;
  if( pars.sp.wstat[k]>0 ||pars.sp.wstat[k+1]>0||pars.sp.wstat[k+2]>0)
    for (i=0; i<3; i++)
      {
        if( pars.sp.wstat[k]>0)
          {//printf("%d %d %s %lg\n",k,ns,pars.sp.liststat[k],Sstat[j][i] );
            fill_pstat(ns,Sstat[j][i]);
            ns++;
          }
        k++;
      }
  k=NLSTATS+9+4+1+3+3+3+3+3+3;/* LD*/
  j=10;
  if( pars.sp.wstat[k]>0 ||pars.sp.wstat[k+1]>0||pars.sp.wstat[k+2]>0)
    for (i=0; i<3; i++)
      {
        if( pars.sp.wstat[k]>0)
          {//printf("%d %d %s %lg\n",k,ns,pars.sp.liststat[k],Sstat[j][i] );
            fill_pstat(ns,Sstat[j][i]);
            ns++;
          }
        k++;
      }
  k=NLSTATS+9+4+1+3+3+3+3+3+3+3;/* rsquare*/
  j=11;
  if( pars.sp.wstat[k]>0 ||pars.sp.wstat[k+1]>0||pars.sp.wstat[k+2]>0)
    for (i=0; i<3; i++)
      {
        if( pars.sp.wstat[k]>0)
          {//printf("%d %d %s %lg\n",k,ns,pars.sp.liststat[k],Sstat[j][i] );
            fill_pstat(ns,Sstat[j][i]);
            ns++;
          }
        k++;
      }
  k=NLSTATS+9+4+1+3+3+3+3+3+3+3+3;/* nhaplo*/
  j=12;
  if( pars.sp.wstat[k]>0 ||pars.sp.wstat[k+1]>0||pars.sp.wstat[k+2]>0)
    for (i=0; i<3; i++)
      {
        if( pars.sp.wstat[k]>0)
          {//printf("%d %d %s %lg\n",k,ns,pars.sp.liststat[k],Sstat[j][i] );
            fill_pstat(ns,Sstat[j][i]);
            ns++;
          }
        k++;
      }
  k=NLSTATS+9+4+1+3+3+3+3+3+3+3+3+3;/* rM*/
  j=13;
  if( pars.sp.wstat[k]>0 ||pars.sp.wstat[k+1]>0||pars.sp.wstat[k+2]>0)
    for (i=0; i<3; i++)
      {
        if( pars.sp.wstat[k]>0)
          {//printf("%d %d %s %lg\n",k,ns,pars.sp.liststat[k],Sstat[j][i] );
            fill_pstat(ns,Sstat[j][i]);
            ns++;
          }
        k++;
      }
  /*** cov matrix ***/
  if(ns==pars.sp.nsp)
    for (i=0; i<pars.sp.nsp; i++)
      for (j=0; j<=i /*<pars.sp.nsp*/; j++)
        {
          if(pars.sp.pstat[i][1]==1 && pars.sp.pstat[j][1]==1)
            {
              //   printf("C %d %d %lg %lg %lg %lg %lg %lg %lg\n",i,j,pars.sp.cov[i][j], pars.sp.pstat[i][0], pars.sp.pstat[j][0], pars.sp.pstat[i][5+j*2], pars.sp.pstat[i][5+i*2], pars.sp.pstat[i][5+j*2+1], pars.sp.pstat[i][5+i*2+1]);
              pars.sp.cov[i][j]+=(long double) pars.sp.pstat[i][0]* pars.sp.pstat[j][0];
              pars.sp.ncov[i][j]++;
              pars.sp.ncov[j][i]= pars.sp.ncov[i][j];
              pars.sp.cov[j][i]= (long double)pars.sp.cov[i][j];
              /* calculate sum/ n for this pair of stats */
              pars.sp.pstat[i][5+j*2]++;
              pars.sp.pstat[j][5+i*2]=pars.sp.pstat[i][5+j*2];
              pars.sp.pstat[i][5+j*2+1]+=pars.sp.pstat[i][0];
              if(j!=i)
                pars.sp.pstat[j][5+i*2+1]+=pars.sp.pstat[j][0];
            } 
        }
  /*** Free memory ***/
  for (i=0; i<14; ++i) free (Sstat[i]);
  free (Sstat);
}//end get cov

void fill_pstat(int nstat,double val)
{
  if(!isnan( val))
    {
      pars.sp.pstat[nstat][3]+=  val; 
      pars.sp.pstat[nstat][4]+=  val* val; 
      pars.sp.pstat[nstat][0]=  val;
      pars.sp.pstat[nstat][1]=  1;
      pars.sp.pstat[nstat][2]++;
    }
  else  pars.sp.pstat[nstat][1]= 0;
}

void reinit_cov()
{
  int i,j;
  if(pars.sp.nsp>0)
    {
      for(i=0;i<pars.sp.nsp;i++)
        {
          for(j=0;j<pars.sp.nsp;j++)
            {
              pars.sp.pstat[i][j]= 0;
              pars.sp.cov[i][j]=0;
              pars.sp.icov[i][j]= 0;
              pars.sp.ncov[i][j]= 0;
            }
        }
    }
}

void free_cov(int pnr)
{
  int i;
  for(i=0;i<pnr+1;i++){
    free(pars.tp.ireg[i]);
    free(pars.tp.nsam[i]);
    free( pars.sp.lstat[i]);
  }
  free( pars.sp.lstat);
  free(pars.tp.ireg);
  free(pars.tp.nsam);
  if(pars.sp.nsp>0)
    {
      for(i=0;i<pars.sp.nsp;i++)
        {
          free(pars.sp.pstat[i]);
          free(pars.sp.cov[i]);
          free(pars.sp.ncov[i]);
          free(pars.sp.icov[i]);
          free(pars.sp.mstat[i]);
        }
      free(pars.sp.pstat);
      free(pars.sp.cov);
      free(pars.sp.icov);
      free(pars.sp.ncov);
      free(pars.sp.mstat);
    }
}


void free_pars()
{
  int i;
  for(i=0;i<pars.cp.npop;i++)
    free(pars.cp.mig_mat[i]);
  free(pars.cp.mig_mat);
  free(pars.cp.config); 
  free(pars.cp.size);
  free(pars.cp.alphag);
  free( pars.tableseeds);

  free(pars.tp.fparv);
  free(pars.tp.freg);
  free(pars.tp.floc);
  free( pars.tp.rho);

  free( pars.sp.type);
  free(pars.sp.wstat); 
  for(i=0;i<59;i++)
    free(pars.sp.liststat[i]);
  free(pars.sp.liststat);
  free( pars.sp.prob);
  free( pars.sp.infoperf);
  free( pars.sp.maxfile);
}

/***********************************************************/
/****************** LIK ESTIMATION *********************/
/*----------------------*
 * Function getprobdata |
 *----------------------*
 | Computes the probability of the statistics given THETA and a set of genealogies.
*/
void getprobdata( double *phypo, int nr, int *ok)
{
  int i=0;
  double tt=0;// total branch of tree used 
  double *rlstat=NULL; 
  char tp1[20],tp2[20];
  if(!(rlstat = (double *)calloc((unsigned) NLSTATS+1,sizeof( double))))
    perror(" calloc error. rlstat");
  pars.sp.prob[2]=0;
  *ok=1;
  for(i=0;i<NLSTATS;i++)
    {
      if(pars.sp.wstat[i]>0)
        {
          tt+=pars.sp.lstat[(int)pars.tp.ireg[nr][21]][i];
          rlstat[i]=(double) pars.sp.lstat[(int)pars.tp.ireg[nr][21]][i]*pars.tp.ireg[nr][16]/phypo[1];// L*rho_r/omega (rescaling)
          if((rlstat[i]==0)&&(pars.tp.ireg[nr][pars.sp.wstat[i]]!=0)) // Case when ARG doesn't fit the data
            {
              //printf("NOTOK %d %lg %Lg %lg %lg %lg \n", i,rlstat[i],pars.sp.prob[2],pars.tp.ireg[nr][14],rlstat[i],pars.tp.ireg[nr][pars.sp.wstat[i]]);
              *ok=0;
              break;
            }
          else{                                          // Case, ARG fit data: get P(Data|THETHA, G)
            pars.sp.prob[2]+=(long double) lnpo((double) pars.tp.ireg[nr][14]*rlstat[i],(int) (pars.tp.ireg[nr][pars.sp.wstat[i]]));
            //printf("OK %d %Lg %lg %lg %lg \n",i, pars.sp.prob[2],pars.tp.ireg[nr][14],rlstat[i],pars.tp.ireg[nr][pars.sp.wstat[i]]);
          }   
        }// end if stats for estimation
    }// end loop on STATS
 
  sprintf(tp1,"%lg",tt);
  sprintf(tp2,"%lg",pars.sp.lstat[(int)pars.tp.ireg[nr][21]][0]);
  if(strcasecmp(tp1,tp2)==0)
    tt=0;
  else
    tt= (double) ( (double) pars.sp.lstat[(int)pars.tp.ireg[nr][21]][0]-tt)*pars.tp.ireg[nr][16]/phypo[1];// L*rho_r/omega (rescaling)
  if((tt==0)&&(pars.tp.ireg[nr][22]!=0)) // Case when ARG doesn't fit the data
    {
      //printf("NOTOK2 %d %lg %Lg %lg %lg %lg \n", i,tt,pars.sp.prob[2],pars.tp.ireg[nr][14],tt,pars.tp.ireg[nr][22]);
      *ok=0;
    }
  else{ 
    if(pars.tp.ireg[nr][22]!=0)                                         // Case, ARG fit data: get P(Data|THETHA, G)
      {
        pars.sp.prob[2]+=(long double) lnpo((double) pars.tp.ireg[nr][14]*tt,(int) (pars.tp.ireg[nr][22]));
        // printf("OK2 %d %lg %Lg %lg %lg %lg \n",i, tt, pars.sp.prob[2],pars.tp.ireg[nr][14],tt,pars.tp.ireg[nr][22]);
      }
  }   
  if(pars.sp.prob[2]< -10000)             // Test: make sure Likelihood large enough
    {
      *ok=0;
      printf("! a u too small: lnps\t%Lg\n",  pars.sp.prob[2]);
    }
  
  if(*ok) // Ok region
    {
      pars.sp.infoperf[3]++;
      pars.sp.prob[1]+=pars.sp.prob[2];// prod of P(S|g,THETA) for genealogy
      //printf("adding %lg %Lg %Lg\n", pars.sp.infoperf[3] ,pars.sp.prob[1],pars.sp.prob[2]);
    }
  else{
    pars.sp.infoperf[4]++; // rejected region
    // printf("rej %lg %Lg %Lg\n", pars.sp.infoperf[4] ,pars.sp.prob[1],pars.sp.prob[2]);
  }
  free(rlstat);
}// End getprobdata function

void get_liks(int ok)
{
  /////////////////// if(pars.sp.prob[1]<-10000) ok=0; // L=0 
  pars.sp.infoperf[0]++;// genealogy generated
  if(ok) // all regions ok for this genealogy
    {
      // printf("ps  %Lg %Lg\n",pars.sp.prob[3],pars.sp.prob[1]);
      pars.sp.prob[3]+=(long double) expl((long double) pars.sp.prob[1]);// sumu for this set of param
      pars.sp.prob[4]+=(long double) expl((long double) pars.sp.prob[1])*expl((long double) pars.sp.prob[1]);// sumu2
      pars.sp.infoperf[1]++;
      //  printf("aps  %Lg %Lg\n",pars.sp.prob[3],pars.sp.prob[1]);
    }
  else  pars.sp.infoperf[2]++;// one genealogy not good
}

long double get_tot_lik(  clock_t pstartl, int nregarg, int init)
{
  int i, j,k, okp=1;
  long double LDp=0, lik=0;
  char tp1[20],tp2[20];
  long double **cov=NULL;
  FILE  *pf=NULL;
  if(init==1) pf=stdout;
  /*  else if(init==0)/\* useless? *\/ */
  /*     { */
  /*       sprintf(tp1,"%s-%s-%d", pars.tp.fparv,pars.tp.freg, nf); */
  /*       pf=fopen(tp1, "w");  */
  /*     } */

  if(!( cov= (long double **)calloc( (unsigned) pars.sp.nsp,sizeof(long double *))))
    perror(" calloc error. cov"); // covariance matrix
  for(i=0;i<pars.sp.nsp;i++)
    {
      if(!(cov[i]= (long double *)calloc( (unsigned) pars.sp.nsp,sizeof(long double))))
        perror(" calloc error. covi");
    }

  if(pars.sp.prob[3]>0 && pars.sp.infoperf[0]>0)////////////////////////////
    {
      okp=1;
      /*** Prob from cov ***/
      for (i=0; i<pars.sp.nsp; i++)
        {
          for (j=0; j <=i ; j++)
            {
              if(pars.sp.ncov[i][j]>1)
                {
                  // printf("BEF i%d j%d c%Lg \t %Lg %Lg\t%Lg\n", i, j,pars.sp.cov[i][j],(long double) pars.sp.cov[i][j]/pars.sp.ncov[i][j] ,((long double) pars.sp.pstat[i][5+j*2+1]/ pars.sp.pstat[i][5+j*2]) *((long double) pars.sp.pstat[j][5+i*2+1]/pars.sp.pstat[j][5+i*2]),(long double) pars.sp.ncov[i][j]/(pars.sp.ncov[i][j]-1)*((long double) pars.sp.cov[i][j]/pars.sp.ncov[i][j] - ((long double) pars.sp.pstat[i][5+j*2+1]/ pars.sp.pstat[i][5+j*2]) *((long double) pars.sp.pstat[j][5+i*2+1]/pars.sp.pstat[j][5+i*2])));
                  sprintf(tp1,"%Lg",(long double) pars.sp.cov[i][j]/pars.sp.ncov[i][j]);
                  sprintf(tp2,"%Lg",((long double) pars.sp.pstat[i][5+j*2+1]/ pars.sp.pstat[i][5+j*2]) *((long double) pars.sp.pstat[j][5+i*2+1]/pars.sp.pstat[j][5+i*2]));
                  if( strcasecmp(tp1,tp2)==0  ) cov[i][j]=0;
                  else
                    cov[i][j]= (long double) pars.sp.ncov[i][j]/(pars.sp.ncov[i][j]-1)*((long double) pars.sp.cov[i][j]/pars.sp.ncov[i][j] - ((long double) pars.sp.pstat[i][5+j*2+1]/ pars.sp.pstat[i][5+j*2]) *((long double) pars.sp.pstat[j][5+i*2+1]/pars.sp.pstat[j][5+i*2]));
                                 
                  if(j!=i)
                    cov[j][i]=  cov[i][j];
                  //printf("AFT %d %d n%d c%Lg %lg %lg %lg %lg %d \t %Lg\n", i, j, pars.sp.ncov[i][j],pars.sp.cov[i][j],pars.sp.pstat[i][5+j*2+1],pars.sp.pstat[i][5+j*2],pars.sp.pstat[j][5+i*2+1],pars.sp.pstat[j][5+i*2], (int) (strcasecmp(tp1,tp2) ), cov[i][j] );
                }
              else
                { 
                  print_error(5,i , j);
                  okp=0;
                  break;
                }
              if(!okp)  break;
            }
         
          if(!okp)  break;
        } /*/////*/
   /*    for (i=0; i<pars.sp.nsp; i++) */
/*         { */
/*           // printf("%d - ",i); */
/*           for (j=0; j <pars.sp.nsp ; j++) */
/*             { */
/*               printf("%01.20Lg\t",cov[i][j]); */
/*             } */
/*           printf("\n"); */
/*         } */
      /*** Like Dp|THETA ***/           
      ////////////////////////////
      /*   if(pars.sp.prob[3]>0 && pars.sp.infoperf[0]>0)//////////////////////////// */
      /*                         { */
      ////////////////////////////
      if(pars.sp.nsp>0 && okp)
        {
          Inverse(cov,pars.sp.icov,pars.sp.nsp); /* inverse cov */
          for (i=0; i<pars.sp.nsp; i++)
            {
              for (j=0; j <=i; j++)
                {
                  if(okp)
                    {
                      sprintf(tp1,"%Lg",pars.sp.icov[i][j]);
                      sprintf(tp2,"nan");
                      if(strcasecmp(tp1,tp2)!=0)
                        {
                         
                          if(j!=i)
                            LDp+= (long double) 2*((long double)pars.sp.mstat[i][0]/pars.sp.mstat[i][1] - ((long double) pars.sp.pstat[i][3]/ pars.sp.pstat[i][2]) )* pars.sp.icov[i][j]* ((long double)pars.sp.mstat[j][0]/pars.sp.mstat[j][1] - ((long double) pars.sp.pstat[j][3]/ pars.sp.pstat[j][2]) );
                          else
                            LDp+= (long double) ((long double)pars.sp.mstat[i][0]/pars.sp.mstat[i][1] - ((long double) pars.sp.pstat[i][3]/ pars.sp.pstat[i][2]) )* pars.sp.icov[i][j]* ((long double)pars.sp.mstat[j][0]/pars.sp.mstat[j][1] - ((long double) pars.sp.pstat[j][3]/ pars.sp.pstat[j][2]) );
                        }
                      else
                        {
                          LDp=pars.sp.icov[i][j];
                          okp=0;
                        }
                      //  printf("%d i %d j %d m0 %lg m1 %lg p3 %lg p2 %lg cij %lg m0 %lg m1 %lg pj3 %lg pj2 %lg\n",okp,i,j, pars.sp.mstat[i][0],pars.sp.mstat[i][1] , pars.sp.pstat[i][3], pars.sp.pstat[i][2], pars.sp.icov[i][j],pars.sp.mstat[j][0],pars.sp.mstat[j][1] , pars.sp.pstat[j][3], pars.sp.pstat[j][2]);
                    }
                 
                  if(init==1)  /* reinit */
                    { 
                      pars.sp.ncov[i][j]=pars.sp.ncov[j][i]=0;
                      if(i==pars.sp.nsp-1)
                        for(k=0; k<5;k++)
                          pars.sp.pstat[j][k]=0;
                      pars.sp.pstat[i][5+j*2]= pars.sp.pstat[j][5+i*2]= pars.sp.pstat[i][5+j*2+1]=pars.sp.pstat[j][5+i*2+1]=0;
                    }
                }
            }
          LDp=(long double) -.5*(LDp + 3*log(2*PI));/* transform into Log(L) */
        }
      /*/////*/
                
      if(okp)
        {
          /*** Lik of S|THETA) ***/// i.e. av over G
          pars.sp.prob[0]= (long double) pars.sp.prob[3]/pars.sp.infoperf[0];
          if(init!=2)
            {                 
              for(i=0; i<NPARAM;i++)
                fprintf(pf,"%lg\t",pars.tp.paramv[i]);
              fprintf(pf,"%Lg\t%Lg\t%Lg\t%Lg\t%lf\t%d\n", LDp,logl(pars.sp.prob[0]), (long double) LDp+logl(pars.sp.prob[0]),(long double) expl(LDp+logl(pars.sp.prob[0])), (double)((double)clock() - (double)pstartl) / (double)CLOCKS_PER_SEC, nregarg);
            }/***/////***/ end L(S|THETA)
          lik= (long double) LDp+logl(pars.sp.prob[0]);
        }
      else
        { 
          if(init!=2)
            { 
              for(i=0; i<NPARAM;i++)
                fprintf(pf,"%lg\t",pars.tp.paramv[i]);
              fprintf(pf,"ldp %Lg\t logl %Lg\tNA\t0\tmin %lf\t%d\n", LDp,logl( (long double) pars.sp.prob[3]/pars.sp.infoperf[0]),(double)((double)clock() - (double)pstartl) / (double)CLOCKS_PER_SEC, nregarg);
            }
          lik= (long double) log(-1);
        }
    }// if L(S)=0
  else
    {
      if(init!=2)
        { 
          for(i=0; i<NPARAM;i++)
            fprintf(pf,"%lg\t",pars.tp.paramv[i]);
          fprintf(pf,"%Lg\tNA\tNA\tNA\t%lf\t%d\n", LDp, (double)((double)clock() - (double)pstartl) / (double)CLOCKS_PER_SEC, nregarg);
        }
      lik=(long double) log(-1);
    }
  for(i=0;i<pars.sp.nsp;i++)
    free(cov[i]);
  free(cov);
  if(init==0 )
    fclose(pf);
  return (lik);
}// end get_tot_lik

long double comp_lik()
{
  int i=0,j=0;;
  FILE *fmax=NULL;
  char *tp=NULL, tpl[1000], line[1000];
  long double max=0,lik=0;
  if( (fmax=fopen(pars.sp.maxfile,"r"))!=NULL)
    {
      fgets( line, 1000, fmax); // Read new set of parameter
      j=0;
      while(line[0]!='\0' && line[0]!='\n')
        { 
          tp = strtok (line,"\n");// divide stats line
          strcpy(tpl,tp);
          tp = strtok (tpl," \t");
          for(i=0;i<9;i++)
            {
              tp = strtok (NULL, " \t");
              if(tp== NULL) break;
            }
          lik=atof(tp);
          if(max<lik ||j==0) max=lik;
          sprintf(line,"%s","\0");// Init
          fgets( line, 1000, fmax); // Read new set of parameter
          j++;
        }
      fclose(fmax);
    }
  return(max);
}// end compare

/***********************************************************/
/****************** CALLED FROM MS.H *********************/
void
argcheck_estL( int arg, int argc, char *argv[] )
{
  if( (arg >= argc ) || ( argv[arg][0] == '-') ) {
    fprintf(stderr,"not enough arguments after %s\n", argv[arg-1] ) ;
    fprintf(stderr,"For usage type: estlikC -h (\"./estlikC -h\")<return>\n");
    exit(0);
  }
}

int
biggerlist_estL(int nsam ,
                char ** list )
{
  int i;
  for( i=0; i<nsam; i++){
    list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
    if( list[i] == NULL ) perror( "realloc error. bigger");
  }
  return(0);
}


/***********************************************************/
/****************** PRINT I/O *********************/
void print_haplo(int nseg, double * pposit, char ** plist)
{
  int i;
  printf("segsites: %d\n",nseg); // Print Simul. hypo. data
  if( nseg > 0 )	printf("positions: ");
  for( i=0; i<nseg; i++)  printf("%6.4lf ",pposit[i] );// only for multi locus //
  printf("\n");
  if( nseg > 0 )
    for(i=0;i<pars.cp.nsam; i++) printf("%s\n", plist[i] ); ////
}

void print_error(int ne, int i , int j)
{
  if(ne==0)
    fprintf(stderr,"The file %s is missing or not avalaible. Run aborted\n", pars.tp.fparv ) ;
  if(ne==1)
    fprintf(stderr,"The file %s is missing or not avalaible. Run aborted\n", pars.tp.freg ) ;
  if(ne==2)
    fprintf(stderr,"There is a problem with the list of statistics in the file %s.\n Check that you used only independent S-statistics and the rigth \"words\". Run aborted\n", pars.tp.freg ) ;
  if(ne==3)
    fprintf(stderr,"Problem with the information on the regions in the file %s at region #%d.\n Check that all the information is provided for this region, that n1+n2>2 (I found n1=%d n2=%d). Check the recombination info (I found w=%lg) or mutation rate info (I found xvZ=%lg).\n Run aborted.\n", pars.tp.freg,i+1,(int) pars.tp.ireg[i][4] ,(int) pars.tp.ireg[i][5] , pars.tp.ireg[i][3], pars.tp.ireg[i][1]*pars.tp.ireg[i][2]*(pars.tp.ireg[i][7]-pars.tp.ireg[i][6])) ;
  if(ne==5)
    fprintf(stderr,"Problem with the covariance matrix the statistics #%i and #%i did not have enought simulated data. Add more simulated data sets with tag '-a'.\n Run aborted. \n", i+pars.sp.nss,j+pars.sp.nss) ;

}

/***********************************************************/
/**************************MATRIX**************************/

/****** Matrices Operation Functions ******/
void Inverse(long double **a,long double **inv,int n )
{
  int i, j;
  long double  d=0;
  CoFactor(a,n,inv); 
  Transpose(inv,inv,n); 
  d=Determinant(a,n);
  for ( i = 0;i < n;i++ )
    for ( j = 0;j < n;j++ )
      {
        inv[ i ][ j ] =(long double) inv[ i ][ j ] / d;
        // if(d<1e-15 && d>-1e-15) inv[ i ][ j ]=log(-1) ; 
        // printf("i %d %d %010.8Lf %Lg\n",i,j,inv[ i ][ j ] ,d);
      }
}

/*
  Find the cofactor matrix of a square matrix
*/
void CoFactor(long double **a,int n, long double **b)
{
  int i,j,ii,jj,i1,j1;
  long double det=0;
  long double **c;

  if( ! (c =(long double**)calloc((unsigned)(n),sizeof(double *))))
    perror(" calloc error. c");
  for (i=0;i<n;i++)
    if( ! (c[i] = calloc((unsigned)(n),sizeof(long double))))
      perror(" calloc error. ci");
  for (j=0;j<n;j++) 
    {
      for (i=0;i<n;i++) 
        { 
          /* Form the adjoint a_ij */
          i1 = 0;
          for (ii=0;ii<n;ii++) 
            {
              if (ii==i)
                continue;
              j1 = 0;
           
              for (jj=0;jj<n;jj++) 
                {   
                  if (jj==j)
                    continue;
                  c[i1][j1] = a[ii][jj];
                  // printf("c %d %d %lf\t %d %d %lf\n",ii,jj,a[ii][jj], i1, j1,c[i1][j1]  );
                  j1++;
                } 
              i1++;
            }
          /* Calculate the determinate */
          det = Determinant(c,n-1);
         
          /* Fill in the elements of the cofactor */
          b[i][j] =(long double) pow(-1.0,i+j+2.0) * det;
        }
    }
  for (i=0;i<n;i++)
    free(c[i]);
  free(c);
}

/*
  Transpose of a square matrix, do it in place
*/
void Transpose(long double **a,long double **m,int n)
{
  int i,j;
 long double tp=0;
  for (i=0;i<n;i++)
    {
      for (j=0;j<i+1;j++) 
        { 
          tp=a[i][j];
          m[i][j] =  a[j][i];
          m[j][i] =  tp;
        }
    }
}


/*
  Recursive definition of determinate using expansion by minors.
*/
long double Determinant(long double **a,int n)
{
  int i,j,j1,j2;
  long double det = 0;
  long double **m = NULL;
  if (n < 1) {}  /* Error */
  else if (n==1) /* Shouldn't get used */
    det =(long double) a[0][0];
  else if (n==2)  
    det =(long double) a[0][0] * a[1][1] - a[1][0] * a[0][1];
  else 
    {
      det = 0;
      for (j1=0;j1<n;j1++) 
        { 
          if( ! ( m = (long double **)calloc((unsigned)(n),sizeof(long double *))))
            perror(" calloc error. m");
          for (i=0;i<n;i++)
            if( ! (m[i] = (long double *)calloc((unsigned)(n),sizeof(long double))))
              perror(" calloc error. mi");
          for (i=1;i<n;i++) 
            { 
              j2 = 0;
              for (j=0;j<n;j++) 
                {
                  if (j==j1)  continue;
                  m[i-1][j2] = a[i][j];
                  j2++;
                }
            }
          det +=(long double) pow(-1.0,1.0+j1+1.0) * a[0][j1] * Determinant(m,n-1);
          for (i=0;i<n;i++)
            free(m[i]);
          free(m);
        }
    }
  return(det);
}


