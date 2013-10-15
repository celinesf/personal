### Function to generate command line arguments for the R function ms. ###


`output_cmd` <-
  function(datafile, param)
  {
    return(.Call("outputcmd",datafile,param, PACKAGE="output_cmdR"))
  
  }


`cmd_line` <-
function(nsam,theta,howmany=1,rho=0,structure=1,time=0,migration=0,extra="NA",seeds="NA")
{
  msline="./msR"
  msline=paste(msline,nsam)
  msline=paste(msline,howmany)
  msline=paste(msline,"-t")
  msline=paste(msline,theta[1])
  ok=1 ## ok=0 if error in command line
### rec.information
  if(rho[1] >0)
    { 
    	if(length(rho)==2) ## make sure specify lenght of locus
        {
          msline=paste(msline,"-r")
          for(s in 1 :length(rho))
      		{msline=paste(msline,rho[s])}
        }
		else
        {
          tp=c(" PROBLEM: The '-r' tag requires 2 arguments. Only one value was given: ",rho, "." ); stp=c("'","'","'","'","'","'","'","'","'"); error_message(tp,stp); return(0) ;
        }
    }
### seeds
  if(seeds[1] != "NA" && ok)
    {
      msline=paste(msline,"-seeds")
      for(s in 1 :length(seeds))
        {msline=paste(msline,seeds[s]) }
    }
### Structure
  if(structure[1]==2 &&ok)
    {
		if(structure[2]+structure[3]==nsam)
        {
          if(length(migration)==1)
            {migration[2]=migration[1]}
          msline=paste(msline,"-I")     # - I tag
          for(s in 1 :length(structure))
      		{msline=paste(msline,structure[s])}
          msline=paste(msline,migration[1])
          ## thetas info
          if(length(theta)==1){theta[3]=theta[2]=theta[1]}
          else
            { 
              if(theta[1]<=0)
                {
                  tp=c(" PROBLEM: The region-specific population mutation rate should be >0: ",theta[1], "." ); stp=c("'","'","'","'","'","'","'","'","'"); error_message(tp,stp); return(0) ;
                }
              
              else
                {
                  if(theta[2]<=0)
                    {
                      tp=c(" PROBLEM: The region-specific population mutation rate for population 2 should be >0: ",theta[2], "." ); stp=c("'","'","'","'","'","'","'","'","'"); error_message(tp,stp); return(0);
                    }
                  else
                    {  
                      if(theta[2]!=theta[1])
                        {
                          msline=paste(msline,"-n 2")
                          msline=paste(msline,theta[2]/theta[1])
                        }
                    }
                }
            }           
        
          ## time info
          if(time[1]>0)                 # Tdiv
            {
              msline=paste(msline,"-ej")
              msline=paste(msline,time[1])
              msline=paste(msline,"1 2")             
              if(length(theta)==3)
                {
                  if(theta[3]<=0){
                    tp=c(" PROBLEM: The region-specific population mutation rate for the ancestral population should be >0: ",theta[3], "." ); stp=c("'","'","'","'","'","'","'","'","'"); error_message(tp,stp); return(0)
                  }
                  else
                    {
                      if(theta[3]!=theta[1])
                        {
                          msline=paste(msline,"-eN")
                          msline=paste(msline,time[1])
                          msline=paste(msline,theta[3]/theta[1])
                        }
                    }
                }
            }
          else
            {
              if(migration[1]<=0)
                {
                  tp=c(" PROBLEM: In an Island model, M_p>0: as is, I have the tag '-I",structure, migration[1], "' and no split time specified." ); stp=c(" "," "," "," ","","'","'","'","'"); error_message(tp,stp); return(0)
                }
            }
          if(length(time)==1) {time[2]=0} # no change of gene flow rate
          if(time[2]>0)                   # Tmigration specified
            {
              if(time[2]<time[1] && migration[1]!=migration[2])
                {
                  msline=paste(msline,"-eM")
                  msline=paste(msline,time[2])
                  msline=paste(msline,migration[2])
                }
              else
                {
                  if(time[1]<time[2])   # || time[1]==time[2])
                    {
                     
                      tp=c(" PROBLEM: The time of split should be > to the time of change of gene flow rates, unlike the values given: ",time[1],"<=", time[2], "." ); stp=c("'","'","'","'","'","'","'","'","'"); error_message(tp,stp); return(0)
                    }
                }
            }
        } ## end nsam=ns1 and 2
      else
        {
          tp=c(" PROBLEM:  nsam is ",nsam,", which does not add up with the numbers of chromosomes reported for the structured populations: " ,structure[2]+structure[3], "." ); stp=c("'","'","'","'","'","'","'","'","'"); error_message(tp,stp); return(0)
        }		
    }                     ## end if 2 pops
  if(extra != "NA" && ok) # if want extra complex model add tag here
    {
      msline=paste(msline,extra)
    }
  return(msline)
} ## end function

