### Function to call the program ms ###
#' Call ms program
#'
#' Generate independent genetic samples using ms
#' @param cmdline the cmd line to run ms
#' @return list of ms samples
`ms` <-
  function (cmdline) 
{
  return(.Call("msMain",cmdline=cmdline, PACKAGE="msR"))  
}


### Function to write the output of the R function ms into a text file in an ms-like manner ###

`write_ms_output`<- function (msframe, msoutput="msout")
  { 
    name=msoutput
    write(file=name,x=msframe$cmdline_seed[1])               # cmdline
    write(file=name,x=msframe$cmdline_seed[2], append=TRUE)  #seed
    
    tp= strsplit(msframe$cmdline_seed[1]," ")
    howmany=as.integer(tp[[1]][3])
    
    for(count in 1:howmany)
      {
        write(file=name,x="\n//", append=TRUE)
        
        tp=paste("segsites:",msframe$segsites[count],sep=" ") #seg
        write(file=name,x=tp, append=TRUE)
        if(msframe$segsites[count]>0)
          {
            tp="positions:"
            i=1
            while(i<= msframe$segsites[count])
              {
                tp=paste(tp,msframe$positions[[count]][i],sep=" ")
                i=i+1
              }
            write(file=name,x=tp, append=TRUE)
            write(file=name,x=msframe$haplotypes[count,], append=TRUE)
          }
      }
  } ### end out_ms_loci


### Function to create a data frame in the R function ms format from data in ms output format. ###

`get_ms_output` <-
  function(msoutput="msout")
{
### get arguments from file
  msout=list()
  print("bforms read");
  msfile=read.table(msoutput,fill=TRUE, header=FALSE)
  print("done")
  msout$cmdline_seeds[1]=""  
  msout$cmdline_seeds[2]=""  
  cmd=msfile[1,1] 
  if(cmd == "./msR"|| cmd == "./ms"||cmd == "./msC" ) ## cmd line
    {
      for(i in 2:length(msfile[1,]))
        {cmd=paste( cmd,msfile[1,i],sep=" ") }
      nsample=msfile[1,2]
      msout$cmdline_seeds[1]=cmd
      tp= strsplit(msout$cmdline_seed[1]," ")
      nsam=as.integer(tp[[1]][2])
      howmany=as.integer(tp[[1]][3])
      
      cmd=msfile[2,1] ## seeds
      for(i in 2:3)
        {cmd=paste( cmd,msfile[2,i],sep=" ") }
      msout$cmdline_seeds[2]=cmd
      
      ## ms frame
      msout$segsites=matrix(0,nr=howmany,nc=1)
      msout$positions=list()
      msout$haplotypes=matrix("",nrow=howmany,ncol=nsam)
      
      a=1
      for(count in 1:howmany)
        {print(count)
          while(msfile[a,1]!="segsites:")
            {a=a+1 }
          
          if(msfile[a,1]=="segsites:")
            { msout$segsites[count,1]=msfile[a,2] }
          else{tp=c("WARNING: The ms output is not in correct format: I am looking for 'segsites:' but I found: ",as.character(msfile$V1[a])," instead."); stp=c("'","'","'","'" );  error_message(tp,stp)}
          
          if( msout$segsites[count,1]>0)
            {
              ## positions
              msout$positions[[count]]=vector(mode="numeric",l=msout$segsites[count,1])
              
              a=a+1
              if(msfile[a,1]=="positions:")
                {
                  for(i in 1:msout$segsites[count,1])
                    {
                      if(is.factor(msfile[a,i+1]))
                        { msout$positions[[count]][i]=suppressWarnings(as.numeric(levels(msfile[a,i+1]))[msfile[a,i+1]]) }
                      else
                        { msout$positions[[count]][i]=as.numeric(msfile[a,i+1]) }
                    }
                }
              else{tp=c("WARNING: The ms output is not in correct format: I am looking for 'positions:' but I found: ",as.character(msfile$V1[a])," instead."); stp=c("'","'","'","'" );  error_message(tp,stp)}
              
              
              for(i in 1:nsam)
                {
                  msout$haplotypes[count,i]=as.character(msfile[i+a,1])
                }
            }                           # end if seg>0
          else
            {
              msout$positions[[count]]=vector("integer",l=1)
              for(i in 1:nsam)
                {
                  msout$haplotypes[count,i]=as.character("")
                }
            }
          count=count+1
        } # end loop on samples
    }                                   # end its an ms output
  else
    {
      msfile2=read.csv(msoutput,fill=TRUE, header=FALSE) 
      a=1
      tp=strsplit(as.character(msfile2[a,1]),split=" ")[[1]]
      if(tp[1]=="//")
        {
          howmany=as.integer(tp[2])
          msout$cmdline_seeds=matrix("",nrow=1,ncol=howmany)
          msout$segsites=matrix(0,nr=howmany,nc=1)
          msout$positions=list()
          msout$haplotypes=list()
          maxn=minn=0
          a=2
          for(count in 1:howmany)
            {
              while(strsplit(as.character(msfile2[a,1]),split=" ")[[1]][1]!="//")# find nexr sample
                {a=a+1 }
              if(strsplit(as.character(msfile2[a,1]),split=" ")[[1]][1]=="//")# new sample
                {
                  tp=strsplit(as.character(msfile2[a,1]),split=" ")[[1]]
                  cmd=tp[1] 
                  for(i in 2:length(tp))
                    {cmd=paste( cmd,tp[i],sep=" ") }
                  nsam=as.integer(tp[2])+as.integer(tp[3])
                  msout$cmdline_seeds[count]=cmd
                  
                  while(strsplit(as.character(msfile2[a,1]),split=" ")[[1]][1]!="segsites:")
                    {a=a+1 }
                  
                  if(strsplit(as.character(msfile2[a,1]),split=" ")[[1]][1]=="segsites:")
                    { 
                      tp=strsplit(as.character(msfile2[a,1]),split=" ")[[1]]
                      msout$segsites[count,1]=as.integer(tp[2]) 
                    }else{tp=c("WARNING: The ms output is not in correct format: I am looking for 'segsites:' but I found: ",as.character(msfile$V1[a])," instead."); stp=c("'","'","'","'" );  error_message(tp,stp)}
                  
                  if( msout$segsites[count,1]>0)
                    {
                      ## positions
                      msout$positions[[count]]=vector(mode="numeric",l=msout$segsites[count,1])
                      msout$haplotypes[[count]]=matrix("",nr=1,nc=nsam)
                      a=a+1
                      if(strsplit(as.character(msfile2[a,1]),split=" ")[[1]][1]=="positions:")
                        {
                          pos=strsplit(as.character(msfile2[a,1]),split=" ")[[1]]
                          for(i in 1:msout$segsites[count,1])
                            {msout$positions[[count]][i]=as.numeric(pos[i+1])}
                          for(i in 1:nsam)
                            {msout$haplotypes[[count]][1,i]=as.character(msfile2[i+a,1])}
                        }
                      else{tp=c("WARNING: The ms output is not in correct format: I am looking for 'positions:' but I found: ",as.character(msfile$V1[a])," instead."); stp=c("'","'","'","'" );  error_message(tp,stp)}                    
                    } else{
                      msout$positions[[count]]=vector("integer",l=1)
                      msout$haplotypes[[count]]=matrix("",nr=1,nc=nsam)
                      for(i in 1:nsam)
                        {
                          msout$haplotypes[[count]][1,i]=as.character("")
                        }
                    }# end else
                }# if next sample
            }# loop regions
        }# appedned regions
    }#not ms
  return(msout)
}
