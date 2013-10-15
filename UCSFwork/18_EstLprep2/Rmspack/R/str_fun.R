### Function to output error messages ###

'error_message' <- function(error_message, separator)
  {
    tp=paste(error_message[1],error_message[2], sep=separator[1])
    for(i in 3 : length(error_message))
      tp=paste(tp,error_message[i], sep=separator[i-1])
    print(tp)
  }


'clean_str' <-function(str, sep=" ")
{
  tp=strsplit(str, split=",")
  tp2=str
  if(tp[[1]][1]!="" & length(tp[[1]])>1){
  tp2=tp[[1]][1]
  for(i in 2:length(tp[[1]]))
  {
  	if(!(tp[[1]][i]=="" || tp[[1]][i]==" " || tp[[1]][i]=="\t" || tp[[1]][i]=="," || tp[[1]][i]=="."))
  	tp2=paste(tp2,tp[[1]][i],sep=sep)
  }}
  tp=strsplit(tp2, split=".")
   if(tp[[1]][1]!="" & length(tp[[1]])>1){
  tp2=tp[[1]][1]
  for(i in 2:length(tp[[1]]))
  {
  	if(!(tp[[1]][i]=="" || tp[[1]][i]==" " || tp[[1]][i]=="\t" || tp[[1]][i]=="," || tp[[1]][i]=="."))
  	tp2=paste(tp2,tp[[1]][i],sep=sep)
  }}
  tp=strsplit(tp2, split="\t")
  if(tp[[1]][1]!="" & length(tp[[1]])>1){
  tp2=tp[[1]][1]
  for(i in 2:length(tp[[1]]))
  {
  	if(!(tp[[1]][i]=="" || tp[[1]][i]==" " || tp[[1]][i]=="\t" || tp[[1]][i]=="," || tp[[1]][i]=="."))
  	tp2=paste(tp2,tp[[1]][i],sep=sep)
  }}

  tp=strsplit(tp2, split=" ")
  if(tp[[1]][1]!="" & length(tp[[1]])>1){
  tp2=tp[[1]][1]
  for(i in 2:length(tp[[1]]))
  {
  	if(!(tp[[1]][i]=="" || tp[[1]][i]==" " || tp[[1]][i]=="\t" || tp[[1]][i]=="," || tp[[1]][i]=="."))
  	tp2=paste(tp2,tp[[1]][i],sep=sep)
  }}
  return(tp2)
}


'get_fname' <-
  function(parf="estpar",regf="regions")
{
 if(is.na(strsplit(regf,"/")[[1]][2])){
	parfile= paste(parf,strsplit(regf,"/")[[1]][1],sep="-")
 }else
	{parfile= paste(parf,strsplit(regf,"/")[[1]][2],sep="-")}
return(parfile)
}