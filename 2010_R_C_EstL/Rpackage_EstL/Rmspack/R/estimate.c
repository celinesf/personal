### Function to estimate parameters of IMc models ### 

`input_estlik` <-
function(paramfile="estimation.par",datafile="datafiles/reg", regfile="datafiles/info.reg", locifile="datafiles/info.loc", outparam=TRUE,outdata=TRUE)
{
  if( file.exists(paramfile))
    {
      parfile=get_fname(parf=paramfile, regf=datafile)
        param=scan(paramfile, comment.char="#",what=c("c","n","n","n"),allowEscapes=FALSE,fill=TRUE,sep="\n",strip.white=TRUE,quiet = TRUE)  
        listparam=c("theta_1","M_present","theta_2","theta_A","T_split","T_change","M_change","rho", "nregions","type","stats","howmany", "parallel" )
        ret=order_param(param=param, listparam=listparam,outparam=outparam,outdata=outdata)
        check=check_param(param=ret$vparam, paramfile=paramfile, listparam=listparam,type=ret$type,summaries=ret$summaries)
        vparam=check$param
        vparam[is.na(vparam[])]=-1 
        vparam[which(vparam[1:7,3]==-1),3]=1
        if(outparam && check$ok) 
          {    	
            tp=c("I will output the parameter values in files starting with ",parfile,"."); stp=c("'","'","'","'",""); error_message(tp,stp);
            output_cmd(datafile=parfile,param=vparam)
          }else
          {
            for(j in 1:vparam[11,1])
              {
                cmd="./estlikC"
                  cmd=paste(cmd,vparam[8,1],sep=" ")
                  i=2
                  while(vparam[8,i]>=0) 
                    {
                      cmd=paste(cmd,vparam[8,i],sep=" ")
                      i=i+1
                    }
                cmd=paste(cmd,"-R",sep=" ")
                  cmd=paste(cmd, vparam[9,1],sep=" ");
                cmd=paste(cmd,"-G",sep=" ");
                cmd=paste(cmd, vparam[10,1], sep=" ");
                cmd=paste(cmd,"-r",   sep=" ");   
                cmd=paste(cmd,paste(parfile,"S",sep="-"),sep=" ");        
                cfile=paste(paramfile,"c",sep="-")
           
                  if(vparam[11,1]==1)
                    {
                      cmd=paste(cmd,"-pf",sep=" ");
                      pfile=paste(paramfile,"p",sep="-")
                      cmd=paste(cmd,pfile, sep=" "); 
                      cmd=paste(cmd,">",sep=" ");
                      pfile=paste(parfile,"lik",sep="-")
                      cmd=paste(cmd,pfile,sep=" ");
                    }
                  else
                    {
                      cmd=paste(cmd,"-pf",sep=" ");
                      pfile=paste(paramfile,"p",sep="-")
                        pfile=paste(pfile,j-1,sep="")
                        cmd=paste(cmd,pfile, sep=" "); 
                      cmd=paste(cmd,">",sep=" ");
                      pfile=paste(parfile,"lik",sep="-")
                        pfile=paste(pfile,j-1,sep="")
                        cmd=paste(cmd,pfile,sep=" ");
                    }
                if(j==1)
                  {write(file=cfile,x=cmd, append=FALSE)
                      }else
                  {write(file=cfile,x=cmd, append=TRUE)}
              }
            tp=c("I will output a cmd line in the file ",cfile,"."); stp=c("'","'","'","'",""); error_message(tp,stp);
          } 
      nr=vparam[9,1]
        if(outdata)
          { 
            if(file.exists(regfile))
              {
                type=check$type
                summaries=check$summaries
                tp=c("I will output the data in the file ",parfile,"-S' ."); stp=c("'","","'","'",""); error_message(tp,stp);
                data4estlik(paramfile=paramfile,datafile=datafile, regfile=regfile, locifile= locifile, summaries=summaries, type=type,nr=nr) 
              }
            else
              {tp=c("PROBLEM: the file of information on the regions ",regfile, " doesn't exist. I could not output the data."); stp=c("'","'","'","'"); error_message(tp,stp);}
          }
    }
  else
    {tp=c("PROBLEM: the parameter file ",paramfile, " doesn't exist. I could not output the grid of parameters and the cmd line."); stp=c("'","'","'","'"); error_message(tp,stp);}
}# end input_estlik

'data4estlik' <- function(paramfile="estimation.par",datafile="datafiles/reg", regfile="datafiles/info.reg", locifile="datafiles/info.loc", summaries="S_1 S_2 S_sl S_sh F_st D_1 D_2 D_star1 D_star2",type=c(0,1,5,10),nr=1)
{
  ok=1
    liststats=list()
    liststats[[1]]=c("S","S_1","S_2","S_s","S_ss","S_sf1","S_sf2","S_sl","S_sh","S_f","S_f1","S_f2","S_o","S_o1","S_o2","F(S)","F(S_1)","F(S_2)","F(S_s)","F(S_ss)","F(S_sf1)","F(S_sf2)","F(S_sl)","F(S_sh)","F_st","H_b","H_w1","H_w2","Snn","pi","pi_1","pi_2","theta_W","theta_W1","theta_W2","theta_H","theta_H1","theta_H2","D","D_1","D_2","H","H_1","H_2","D_star","D_star1","D_star2","D_prime","D_prime1","D_prime2","r_square","r_square1","r_square2","nH","nH_1","nH_2","Rm","Rm_1","Rm_2")
    summaries=clean_str(str=summaries, sep=" ")
    nstats=get_liststats(liststats=liststats,summaries=strsplit(summaries, split=" "))
  
    reg=read.table(regfile,skip=1,fill=TRUE)
    nregions=dim(reg)[1]
    if(nr!=nregions)
      {tp=c("PROBLEM: I found ",nregions, " regions/lines in the file ",regfile, " but in the file of parameters I found ",nr," regions. I didn't output the datafiles."); stp=c("'","'","'","'","'","'"); error_message(tp,stp)
      }else
      {
        loci=0
          if(sum(reg[(1:nregions),12])>nregions)
            { 
              if(file.exists(locifile))
                {loci=read.table(locifile,skip=1,fill=TRUE)}
              else  {tp=c("PROBLEM: the file of information on the loci ",locifile, " doesn't exist. I could not output the grid of parameters and the cmd line."); stp=c("'","'","'","'"); error_message(tp,stp);ok=0}
            }
        if(ok)
          {
            check=check_reg(nregions=nregions,info_region=reg, info_loci=loci, locifile=locifile)
              if(check$ok) ## if region/loci are ok
                  {
                    append=0
                    if(file.exists(datafile))
                      {
                        tp=c("I found all the data in the file ",datafile,"."); stp=c("'","'","'","'",""); error_message(tp,stp);
                        msout=get_ms_output(msoutput=datafile)
                        append=1
                        if(dim(msout$segsites)[1]!= nr)
                          {tp=c("PROBLEM: I found ",nregions, " regions in the file ",regfile, " but in the data file ",datafile," I found ",dim(msout$segsites)[1]," regions. I didn't output the datafiles."); stp=c("'","'","'","'","'","'","'","'"); error_message(tp,stp);check$ok=0}
                      }
     
                    data=list()
                    sumstats=list()
                    sumstatsloci=list()
                    nloci=1
                    nl=1
                    nsam=0
                    for(nr in 1:nregions)
                      {
                        if(check$ok)
                          {
                            data[[nr]]=list()
                            if(append==1)
                              {
                                data[[nr]]$cmdline_seeds[1]=msout$cmdline_seeds[1,nr]
                                data[[nr]]$segsites=msout$segsites[nr]
                                data[[nr]]$positions[[1]]=msout$positions[[nr]]
                                data[[nr]]$haplotypes=msout$haplotypes[[nr]]
                                tp=strsplit(data[[nr]]$cmdline_seeds[1], split=" ")
                                nsam=as.numeric(tp[[1]][2])+as.numeric(tp[[1]][3])
                              } else{
                              fdata=paste(datafile,nr, sep="")
                              if(file.exists(fdata))	
                                {
                                  tp=c("I am recovering the data for region #",nr, " in file ", fdata,"."); stp=c("","","'","'","'","'","'","'"); error_message(tp,stp);
                                  data[[nr]]=get_ms_output(msoutput=fdata)
                                  tp=strsplit(data[[nr]]$cmdline_seeds[1], split=" ")
                                  nsam=as.numeric(tp[[1]][2])
                                }
                              else{tp=c("PROBLEM: I couldn't find the data file for the region #",nr, " named ",fdata, ". I didn't output the datafiles."); stp=c("","","'","'","'","'","'","'"); error_message(tp,stp);check$ok=0}
                            }
                            if( nsam== reg[nr,6]+reg[nr,7] && check$ok)
                              {
                                nsam= reg[nr,6]+reg[nr,7] 
                                w=reg[nr,5]  
                                Z=  reg[nr,9]-reg[nr,8]
                                rho=c(w,Z-1) #lentgh in bp
                                zxv=Z*reg[nr,3]*reg[nr,4] # inheritance/ mutation scalar
                  
                                mig=struct=time=0                                          
                                if(check$ok && reg[nr,6]>0 && reg[nr,7]>0) # case with 2 pops
                                                                             { struct=c(2,reg[nr,6],reg[nr,7])}
                                if(check$ok)
                                  {
                                    cmd=cmd_line(nsam=nsam, theta=zxv, rho=rho,structure=struct,migration=1, time=time)             
                                    data[[nr]]$cmdline_seeds[1]=cmd[1]
                                    st=stats_pop(msframe= data[[nr]], type=type)
                                    reg[nr,13]=st$segsites[1,1]
                                    sumstats=get_stats(nloci=nregions, nlocus=nr, sumstats=sumstats,statslocus=st,info=reg[nr,],summaries=nstats$stats)
                      
                                    if(reg[nr,12]>1) #case with multiple loci
                                                       {
                                                         for(l in 1:reg[nr,12])
                                                           {
                                                             mslocus=get_loci(msframe=data[[nr]], info_region=reg[nr,], info_loci=loci[loci$V1==reg[nr,1],], nlocus=l)
                                                             mslocus$cmdline_seeds=data[[nr]]$cmdline_seeds
                                                             stl=stats_pop(msframe=mslocus, type=type)
                                                             sumstatsloci=get_stats(nloci=sum(reg[,12]),nlocus=nloci, sumstats=sumstatsloci,statslocus=stl,info=loci[loci$V1==reg[nr,1],][loci[loci$V1==reg[nr,1],]$V2==l,],summaries=nstats$stats)
                                                             loci[nl,7]=stl$segsites[1,1]
                                                             nl=nl+1
                                                             nloci=nloci+1
                                                           }                       
                                                       }else
                                      {
                                        sumstatsloci=get_stats(nloci=sum(reg[,12]), nlocus=nloci, sumstats=sumstatsloci,statslocus=st,info=reg[nr,],summaries=nstats$stats)
                                        nloci=nloci+1
                                      }
                                  } ## end ok
                              }     ### end if nsam match
                          }         # end ok last region
                      }             ### end loop on regions
                    if(check$ok)
                      { 
                        parfile=get_fname(parf=paramfile, regf=datafile)
                        output_stats(datafile=parfile, sumstats=sumstats, liststats=liststats[[1]][nstats$list], reg=reg,locifile=NA ,regfile=regfile)
                        if(sum(reg[(1:nregions),12])>nregions)
                          output_stats(datafile=paste(parfile,"loci",sep="-"), sumstats=sumstatsloci,liststats=liststats[[1]][nstats$list], reg=reg,loci=loci,locifile=locifile,regfile=regfile)
                            ### files with cmd lines for estlik
                      } ### end of if reg ok
                  }# end check ok 
          }# end if locifile ok
      }# end if ok # regions in paramfile and region info file
    
}# end data4estlik


