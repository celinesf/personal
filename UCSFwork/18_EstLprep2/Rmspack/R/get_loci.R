
### Function accounting for missing data and multiple loci and gaps in the ms output for a particular region. ###
`get_loci`<- function (msframe, info_region, info_loci, nlocus)
 {
   newms=list()
   newms$segsites=0;
   newms$positions[[1]]=vector();
   newms$haplotypes=matrix("",nr=1,ncol=info_region$V6+info_region$V7); 
   k=0
   j=nlocus
   if(msframe$segsites>0)
   {
   	for(i in 1:msframe$segsites)
     {
     
       if(as.integer(msframe$positions[[1]][i])>=info_loci[j,5] && as.integer(msframe$positions[[1]][i])<info_loci[j,6] )
         {
           k=k+1
           newms$segsites=  k
           newms$positions[[1]][k]=as.integer(msframe$positions[[1]][i])
           newms$ns[k]=i          #seg sites number to keep           
           for(n1 in 1:info_region$V6)
             newms$haplotypes[1,n1]=paste(newms$haplotypes[1,n1], unlist(strsplit(msframe$haplotypes[1,n1], ""))[i], sep="")
           
           for(n2 in info_region$V6+1:info_region$V7)
             newms$haplotypes[1,n2]=paste(newms$haplotypes[1,n2],unlist(strsplit(msframe$haplotypes[1,n2], ""))[i], sep="")
               
         } ## found locus          
     	}     ### end loop on segsites 
     }  
   return(newms)
 } ### get_loci
