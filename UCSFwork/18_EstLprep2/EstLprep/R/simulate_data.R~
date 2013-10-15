### Function to simulate data under IMc models ### 

`simulate_data` <-
  function(paramfile="simulation.par",datafile="reg", regfile="info.reg", locifile="info.loc", append=TRUE) 
{
### define stats used in estimation
  liststats=list()
  liststats[[1]]=c("S","S_1","S_2","S_s","S_ss","S_sf1","S_sf2","S_sl","S_sh","S_f","S_f1","S_f2","S_o","S_o1","S_o2","F(S)","F(S_1)","F(S_2)","F(S_s)","F(S_ss)","F(S_sf1)","F(S_sf2)","F(S_sl)","F(S_sh)","F_st","H_b","H_w1","H_w2","Snn","pi","pi_1","pi_2","theta_W","theta_W1","theta_W2","theta_H","theta_H1","theta_H2","D","D_1","D_2","H","H_1","H_2","D_star","D_star1","D_star2","D_prime","D_prime1","D_prime2","r_square","r_square1","r_square2","nH","nH_1","nH_2","Rm","Rm_1","Rm_2")
                                        #liststats=strsplit(liststats, split="\t")  
  nstats=get_liststats(liststats=liststats,summaries=liststats)
  
### get parameters from file
  listparam=c("theta_1","M_present","theta_2","theta_A","T_split","T_change","M_change","rho", "nregions","type")
  param=scan(paramfile, comment.char="#",what=c("c","n","n","n"),allowEscapes=FALSE,fill=TRUE,sep="\n",strip.white=TRUE,quiet = TRUE) #

  ret=order_param(param=param, listparam=listparam)
  check=check_param(param=ret$vparam, paramfile=paramfile, listparam=listparam,type=ret$type)
  if(check$ok) ## All the  parameters are ok
    {
      vparam=check$param
      type=check$type
      reg=read.table(regfile,skip=1,fill=TRUE)
      loci=0
      if(sum(reg[(1:vparam[9,1]),12])>vparam[9,1])
        loci=read.table(locifile,skip=1,fill=TRUE)
      check=check_reg(nregions=vparam[9,1],info_region=reg, info_loci=loci, locifile=locifile)
      reg=check$reg
      sumstats=list()
      sumstatsloci=list()
######## Get param values
      mig=time=struct=0 
      par=vector("numeric",l=7)
      for(p in 1:7)                     # param fix or variable
        {
          if(!is.na(vparam[p,1]))
            {
              if(is.na(vparam[p,2]))
                {
                  par[p]=vparam[p,1]
                }else
              {
                if((p==1 || p==3 || p==4) && vparam[p,1]==0) #thetas
                  {
                    par[p]=sample((vparam[p,2]-vparam[p,1])/vparam[p,3]+(vparam[p,2]-vparam[p,1])/vparam[p,3]*(0:(vparam[p,3]-1)), size=1)
                  }else
                {
                  if((p==6) && vparam[p,2]==1.) #T_c (T_c=T_s) upper lim =1
                    {
                      par[p]=sample(vparam[p,1]+(vparam[p,2]-vparam[p,1])/vparam[p,3]*(0:(vparam[p,3]-1)), size=1)                            
                    }else{
                      par[p]=sample(vparam[p,1]+(vparam[p,2]-vparam[p,1])/(vparam[p,3]-1)*(0:(vparam[p,3]-1)), size=1)
                      if((p==2 && par[p]==0 && is.na(vparam[5,2])) ||(p==5 && par[p]==0 && par[2]==0) ) # T_s=0 M=0
                        {
                          while(par[p]==0)
                            {
                              par[p]=sample(vparam[p,1]+(vparam[p,2]-vparam[p,1])/(vparam[p,3]-1)*(0:(vparam[p,3]-1)), size=1)
                            }
                        }
                    }
                }
              }
            }
          else
            { if( (p==3 || p==4) && par[p]==0)
                {
                  par[p]=par[1]
                }
            }
        }
######## Get param values
      nloci=1
      for(nr in 1:vparam[9,1])
        {
          if(check$ok)
            { 
              nsam=reg[nr,6]+reg[nr,7]
              Z=reg[nr,9]-reg[nr,8]     #lentgh in bp
              zxv=Z*reg[nr,3]*reg[nr,4] # inheritance/ mutation scalar
              w=reg[nr,5]  
              rho=0
              if(w>0 && !is.na(vparam[8,1]))
                {
                  tp=get_rho(param=vparam[8,],nregion=nr,Z=Z,W=w,theta=par[1])
                  rho=c(tp$rho,Z)
                  check$ok=tp$ok
                }else
                {
                  tp=c("PROBLEM: Recombination needs to be specified (for effectively non recombining regions, specify a small number for w). rho=",vparam[8,],", and the genomic region reg", nr, " has w=",reg[nr,5],"."); stp=c("'"," "," ","'","","","","","","","","","",""); error_message(tp,stp);check$ok=0;
                }
              if(check$ok && reg[nr,6]>0 && reg[nr,7]>0) # case with 2 pops ni>0
                {
                  struct=c(2,reg[nr,6],reg[nr,7])                       
                  if(par[2]>0)          # Migration constant>0
                    {mig=par[2]
                   }else{
                     if(par[5]==0)      #time not defined
                       { 
                         tp=c("PROBLEM: When some regions have data sampled from two populations M_p>0 is required to be specified with the keyword ",listparam[2],". The value is now ",listparam[2],"_c=",vparam[2,1],", while the genomic region reg", nr, " has n_1=",reg[nr,6], " and n_2=",reg[nr,6], "."); stp=c("'","'","'","","","'","","","","","","","",""); error_message(tp,stp);check$ok=0;}
                   }
                  time=c(par[5],par[5]*par[6])
                  mig=c(par[2],par[7])
                }
              if(check$ok)
                {
                  cmd=cmd_line(nsam=nsam, theta=c(par[1]*zxv,par[3]*zxv,par[4]*zxv), rho=rho,structure=struct,migration=mig, time=c(par[5],par[5]*par[6]))
            
                  msframe=ms(cmd)
                
                  msloci= output_ms(datafile=datafile,nregion=nr, msframe=msframe, info_region=reg[nr,], info_loci=loci[loci$V1==reg[nr,1],],append=append,nregions=vparam[9,1])
                  msloci$cmdline_seeds=msframe$cmdline_seeds
                  st=stats_pop(msframe=msloci, type=type)
                  sumstats=get_stats(nloci=vparam[9,1], nlocus=nr, sumstats=sumstats,statslocus=st,info=reg[nr,],summaries=nstats$stats)

                  if(reg[nr,12]>1)      #case with multiple loci
                    {
                      for(l in 1:reg[nr,12])
                        {
                          if(msloci$segsites>0)
                            { 
                              mslocus=get_loci(msframe=msloci, info_region=reg[nr,], info_loci=loci[loci$V1==reg[nr,1],], nlocus=l)
                              mslocus$cmdline_seeds=msframe$cmdline_seeds
                            }else 
                          mslocus=msframe
                          stl=stats_pop(msframe=mslocus, type=type)
                          sumstatsloci=get_stats(nloci=sum(reg[,12]), nlocus=nloci, sumstats=sumstatsloci,statslocus=stl,info=loci[loci$V1==reg[nr,1],][loci[loci$V1==reg[nr,1],]$V2==l,],summaries=nstats$stats)
                          nloci=nloci+1
                        }                       
                    } else
                  {
                    sumstatsloci=get_stats(nloci=sum(reg[,12]), nlocus=nloci, sumstats=sumstatsloci,statslocus=st,info=reg[nr,],summaries=nstats$stats)
                    nloci=nloci+1
                  }
                }
            }                           # end ok last region
        }                               ### end loop on regions
      if(check$ok){
        output_stats(datafile=datafile, sumstats=sumstats,liststats=liststats[[1]][nstats$list], regfile= regfile,locifile=NA, reg=reg)
        if(sum(reg[(1:vparam[9,1]),12])>vparam[9,1])
          output_stats(datafile=paste(datafile,"loci",sep="-"), sumstats=sumstatsloci,liststats=liststats[[1]][nstats$list], regfile= regfile,locifile=locifile, reg=reg, loci=loci)
      }
    } ### end of if param ok
}     ### end of function simulate_data


### Function to write ms outputs into a regr file. ###
`output_ms` <- function(datafile,nregion,msframe, info_region, info_loci=0,append=FALSE,nregions)
 {
   if(append==TRUE)
     {
       name=paste(datafile,"",sep="")
      
       if(nregion==1)
         {
           x="//"
           x= paste(x,nregions,sep=" ");#nregions to appends
           write(file=name,x=x)
       }
        x="\n//"
       x= paste(x,info_region[1,6],sep=" ");#n1
       x= paste(x,info_region[1,7],sep=" ");#n2
       x= paste(x,msframe$cmdline_seed[1],sep=" ");#cmd
       x= paste(x,"-seeds",sep=" ");#prep for seeds
       x= paste(x,msframe$cmdline_seed[2],sep=" ");#seeds
       write(file=name,x=x, append=TRUE)  #seed
     }
   else
     {
       name=paste(datafile,nregion,sep="")
       write(file=name,x=msframe$cmdline_seed[1])               # cmdline       
       write(file=name,x=msframe$cmdline_seed[2], append=TRUE)  #seed
       write(file=name,x="\n//", append=TRUE)
     }      
   msframe$positions[[1]]=msframe$positions[[1]]*(info_region[1,9]-info_region[1,6])
   if(info_region[1,12]>1 && msframe$segsites>0)  #case with multiple loci
     msframe=out_ms_loci(msframe, info_region, info_loci)
   if(info_region[1,10]>0)# missing in pop1
     {
       nm=sample(0:(msframe$segsites*info_region[1,6]*info_region[1,10]), size=1,replace = FALSE)
       miss1=sample(1:(msframe$segsites*info_region[1,6]), size=nm,replace = FALSE)
        if(nm>0){for(i in 1:nm)
         {
           n=round(miss1[i]/msframe$s+.49)
           p=miss1[i]-((n-1)*msframe$s)
           str=substring(msframe$haplotypes[n],c(1,p,p+1),c(p-1,p,msframe$s))
           str[2]="?"
           str=paste(paste(str[1],str[2],sep=""),str[3],sep="")
           msframe$haplotypes[n]=str
         } }      
     }
   if(info_region[1,11]>0)# missing in pop2
     {
       nm=sample(0:(msframe$segsites*info_region[1,7]*info_region[1,11]), size=1,replace = FALSE)
       miss2=sample(1:(msframe$segsites*info_region[1,7]), size=nm,replace = FALSE)
# print(c(nm))
       if(nm>0){ for(i in 1:nm)
         {
           n=round(miss2[i]/msframe$s+.49)
           p=miss2[i]-((n-1)*msframe$s)
           str=substring(msframe$haplotypes[info_region[1,6]+n],c(1,p,p+1),c(p-1,p,msframe$s))
          str[2]="?"
           str=paste(paste(str[1],str[2],sep=""),str[3],sep="")
           msframe$haplotypes[info_region[1,6]+n]=str
         }}  
     }
     
   tp=paste("segsites:",msframe$segsites,sep=" ") #seg
   write(file=name,x=tp, append=TRUE)
   if(msframe$segsites>0)
   {
     tp="positions:"
     i=1
     while(i<= msframe$segsites)
       {
         tp=paste(tp,as.integer(msframe$positions[[1]][i]),sep=" ")
         i=i+1
       }
     write(file=name,x=tp, append=TRUE)
     write(file=name,x=msframe$haplotypes, append=TRUE)
   }
   return(msframe)
 }# end ms_output


### Function accounting for missing data and multiple loci and gaps in the ms output for a particular region. ###
`out_ms_loci`<- function (msframe, info_region, info_loci)
 {
   newms=list()
   newms$segsites=0;
   newms$positions[[1]]=vector();
   newms$haplotypes=matrix("",nr=1,ncol=info_region$V6+info_region$V7); 
   k=0
   j=1

   nmiss=list()
   nmiss$a=list()
   nmiss$b=list()
   for(l in 1: info_region$V12)
     {
       nmiss$a[[l]]=nmiss$b[[l]]=0
       if(info_loci[l,3]<info_region$V6) # missing data
         {
           nmiss$a[[l]]=vector("integer",l=info_region$V6-info_loci[l,3])
           nmiss$a[[l]]=sample(1:info_region$V6, size=info_region$V6-info_loci[l,3],replace = FALSE)
         }
       if(info_loci[l,4]<info_region$V7) # missing data
         {
           nmiss$b[[l]]=vector("integer",l=info_region$V7-info_loci[l,4])
           nmiss$b[[l]]=sample(info_region$V6+1:info_region$V7, size=info_region$V7-info_loci[l,4],replace = FALSE)
         }
     }
   for(i in 1:msframe$segsites)
     {j=1
      while((as.integer(msframe$positions[[1]][i])<info_loci[j,5] || as.integer(msframe$positions[[1]][i])>=info_loci[j,6]) && j< dim(info_loci)[1])
        { j=j+1}
      if(as.integer(msframe$positions[[1]][i])>=info_loci[j,5] && as.integer(msframe$positions[[1]][i])<info_loci[j,6] )
        {
          k=k+1
          newms$segsites=  k
          newms$positions[[1]][k]=as.integer(msframe$positions[[1]][i])
          newms$ns[k]=i           #seg sites number to keep           
          for(n1 in 1:info_region$V6)# loop on sample for pop1
            {
              miss=0
              if(length(nmiss$a[[j]])>0)
                for(i in 1:length(nmiss$a[[j]]))
                  if(n1==nmiss$a[[j]][i])
                    miss=1
               
              if(miss)
                newms$haplotypes[1,n1]=paste(newms$haplotypes[1,n1],"?", "", sep="")
              else                 
                newms$haplotypes[1,n1]=paste(newms$haplotypes[1,n1], unlist(strsplit(msframe$haplotypes[1,n1], ""))[[1]], sep="")
            }      
          for(n2 in info_region$V6+1:info_region$V7)
            {
              miss=0
              if(length(nmiss$b[[j]])>0)
                for(i in 1:length(nmiss$b[[j]]))
                  if(n2==nmiss$b[[j]][i])
                    miss=1              
              if(miss)
                newms$haplotypes[1,n2]=paste(newms$haplotypes[1,n2],"?", "", sep="")
              else
                newms$haplotypes[1,n2]=paste(newms$haplotypes[1,n2],unlist(strsplit(msframe$haplotypes[1,n2], ""))[[1]], sep="")
            }  
        } ## found locus          
    }     ### end loop on segsites
   return(newms)
 } ### end out_ms_loci
