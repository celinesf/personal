### Function to calculate rho from prior or fix non-variable part of rho, and check that w correct fromat in regfile ###


'get_rho' <-
  function(param,nregion,Z,W, theta=1)
{
  rok=list();
  rok$ok=1
  rok$rho=param[1]*(Z-1)*W
  if((param[1]>0 && param[1]<1 || param[1]==-1 || param[1]==-2) && W>0 && (W<.1||W>2)) # rho=4Nc
    { tp=c("0? PROBLEM: With the recombination rate: When I calculate the region-specific intra-region recombination rate I obtain rho_",nregion,"=rho(Z-1)*w=",  rok$rho,". The parameters on the recombination are ",param, " and w=", W, " in the genomic region region ", nregion,  "."); stp=c("","","'","'",", ",", ","'","'","'","'","'","'","'","'"); error_message(tp,stp);rok$ok=0 
    }
  if(param[1]==1)                     # 1 , 4Nc fixed
    {
     if(W>0 && (W>.1||W<1e-6))
       {
         tp=c("1PROBLEM: With the recombination rate: When I calculate the region-specific intra-region recombination rate I obtain rho_",nregion,"=(Z-1)*w=",(Z-1)*W,". The parameters on the recombination are ",param[1], " and w=", W, " in the genomic region region ", nregion,  "."); stp=c("","","'","'","'","'","'","'","'","'","'","'"); error_message(tp,stp);rok$ok=0;
       }else
     {rok$rho=(Z-1)*W}
   }
  if(param[1]==2)                     # 2 mu , c fixed
    {
     if(W>0 && W>1e-6)
       {
         tp=c("2PROBLEM: The recombination rate is too high. When I calculate the region-specific intra-region recombination rate I obtain rho_",nregion,"=theta_1*(Z-1)*w/mu=",theta*(Z-1)*W/param[2],". The parameters on the recombination are ",param[1],param[2] , " and w=", W, " in the genomic region region ", nregion,  "."); stp=c("","","'","'","'"," ","'","'","'","'","'","'"); error_message(tp,stp);rok$ok=0;
       }else
     {rok$rho=theta*(Z-1)*W/param[2]}
   }                         
  if(param[1]==-1 && !is.na(param[2]) ) ## take c/mu from exponentiel lambda
    { rok$rho=theta*rexp(1,rate=param[2])*(Z-1)*W }
  if(param[1]==-2 && !is.na(param[2])&& !is.na(param[3]) ) ## take c/mu from exponentiel lambda
    {
     tp=rnorm(1,mean=param[2],sd=param[3])
     while(tp<0)
       {tp=rnorm(1,mean=param[2],sd=param[3])}
     rok$rho=theta*tp*(Z-1)*W
   }

  if(rok$ok && rok$rho/(theta*(Z-1))>100)
    {
      tp=c("3PROBLEM: With the recombination rate: When I calculate the region-specific intra-region recombination rate I obtain rho_",nregion,"=rho(Z-1)*w=",  rok$rho,". The parameters on the recombination are ",param, " and w=", W, " in the genomic region region ", nregion,  ".", theta); stp=c("","","'","'","'",", ",", ","'","'","'","'","'",""); error_message(tp,stp);rok$ok=0 
    }
  return(rok)
}
