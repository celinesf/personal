# INSTRUCTIONS
# At any points during the analysis, copy the output file from 
# and rename the copy by "outputmimar1" . 
# If the analyses were not completed, open this files and delete 
# the last incomplete rows, and save.
# Then execute the following lines in R.
# The five posterior autocorrelation functions, then the 5 scalar 
# functions will be plotted for the 5 parameters of interest.
# - Note that if you recorded a lot of steps, R might take a long
# time downloading the results.
# Recording only every "int" steps (switch "-i int" in MIMAR) does
# not affect the results drasticaly but allows to read the file much faster.

rm(list=ls())                          # delete previous objects in R

########### Setting parameteres ##########
g=read.table("outputmimar1", skip=12)  # change the name of MIMAR output file if needed
##### End of parameters and option setting ######
                        
par(mfrow=c(3,2))
############### Theta1 
acf(g$V2,col="red",main=expression(theta[1]) )
############### Theta2 
acf(g$V3,col="red",main=expression(theta[2]) )
############### ThetaA 
acf(g$V6,col="red",main=expression(theta[A]) )
############### Mig 
acf(g$V7,col="red",main=expression(M) )
############### T
acf(g$V5,col="red",main=expression(T) )

par(mfrow=c(3,2))
############### Theta1 
plot(g$V2,col="red",ylab=expression(theta[1]),xlab="" )
############### Theta2 
plot(g$V3,col="red",ylab=expression(theta[2]),xlab=""  )
############### ThetaA 
plot(g$V6,col="red",ylab=expression(theta[A]),xlab=""  )
############### Mig 
plot(g$V7,col="red",ylab=expression(M),xlab=""  )
############### T
plot(g$V5,col="red",ylab=expression(T),xlab=""  )
