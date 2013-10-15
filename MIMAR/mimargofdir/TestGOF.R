# INSTRUCTIONS	
# Run perlgof on the output of MIMAR for two different seeds.
# The output files are assumed to be outputgof1 and outputgof2. 
# Some variables need to be changes:	
# - In the array named "true" change the values to the observed value: 
# sum over loci of S statistics and mean across loci of Fst, pi1, pi1, D1 and D2
# - Change the value of nbin to smooth the distributions
# Execute in R
# Output:
# - Four graphs with the distributions of the statistics and the observed values in vertical lines
# - The p-value for each statistics: if p-value<.05, the data does not fit the estimated model well.
################
rm(list=ls())

########### Setting parameteres ##########
# change the file names with the posterior distributions    
gof_file1="outputgof1"
gof_file2="outputgof2"
# change the values below to the observed value in the input file
true=c(
27,
28,
27,
0,
0.08308975,
4.788889,
4.91111125,
-0.28935375,
0.072975,
0.17240075,
0.26331175
) 
# change this number to smooth the distributions
nbin=10
                           
##### End of parameters and option setting ######

g=read.table(gof_file1,skip=1, na.strings = "NA")
gg=read.table(gof_file2,skip=1, na.strings = "NA") 
g=na.omit(g)
gg=na.omit(gg)
par(mfrow=c(2,2))
############### S1 
h1=(max(g$V3,gg$V3)-min(g$V3,gg$V3))/nbin
h2=c(min(g$V3,gg$V3)-h1,max(g$V3,gg$V3)+h1)
b=c()
for(i in 1:(nbin+2))
{
b[i]=h2[1]+(i-1)*h1
}
f=hist(g$V3,breaks=b,plot=F)
f$counts=f$counts/sum(f$counts)
ff=hist(gg$V3,breaks=b, plot=F)
ff$counts=ff$counts/sum(ff$counts)
S1=f
SS1=ff
############### S2 
h1=(max(g$V4,gg$V4)-min(g$V4,gg$V4))/nbin
h2=c(min(g$V4,gg$V4)-h1,max(g$V4,gg$V4)+h1)
b=c()
for(i in 1:(nbin+2))
{
b[i]=h2[1]+(i-1)*h1
}
f=hist(g$V4,breaks=b, plot=F)
f$counts=f$counts/sum(f$counts)
ff=hist(gg$V4,breaks=b, plot=F)
ff$counts=ff$counts/sum(ff$counts)
S2=f
SS2=ff
############### Ss 
h1=(max(g$V5,gg$V5)-min(g$V5,gg$V5))/nbin
h2=c(min(g$V5,gg$V5)-h1,max(g$V5,gg$V5)+h1)
b=c()
for(i in 1:(nbin+2))
{
b[i]=h2[1]+(i-1)*h1
}
f=hist(g$V5,breaks=b, plot=F)
f$counts=f$counts/sum(f$counts)
Ss=f
ff=hist(gg$V5,breaks=b, plot=F)
ff$counts=ff$counts/sum(ff$counts)
SSs=ff
############### Sf 
h1=(max(g$V6,gg$V6)-min(g$V6,gg$V6))/nbin
h2=c(min(g$V6,gg$V6)-h1,max(g$V6,gg$V6)+h1)
b=c()
for(i in 1:(nbin+2))
{
b[i]=h2[1]+(i-1)*h1
}
f=hist(g$V6,breaks=b, plot=F)
f$counts=f$counts/sum(f$counts)
Sf=f
ff=hist(gg$V6,breaks=b, plot=F)
ff$counts=ff$counts/sum(ff$counts)
SSf=ff
limx=c(min(true[1],S1$breaks,true[2],S2$breaks,true[3],Ss$breaks,true[4],Sf$breaks),max(true[1],S1$breaks,true[2],S2$breaks,true[3],Ss$breaks,true[4],Sf$breaks))
m=max((S1$counts+SS1$count)/2,(S2$counts+SS2$count)/2,(Ss$counts+SSs$count)/2,(Sf$counts+SSf$count)/2)
plot(S1$mids,(S1$counts+SS1$count)/2,'l',cex.lab=1.4,xlab=expression("S statistics"),ylab="",xlim=limx,ylim=c(0,m), col="blue",lwd=2)
lines(c(true[1],true[1]),c(0,m),col="blue",lwd=1)
lines(S2$mids,(S2$counts+SS2$count)/2,col="red",lwd=2)
lines(c(true[2],true[2]),c(0,m),col="red",lwd=1)
lines(Ss$mids,(Ss$counts+SSs$count)/2,col="dark grey",lwd=2)
lines(c(true[3],true[3]),c(0,m),col="dark grey",lwd=1)
lines(Sf$mids,(Sf$counts+SSf$count)/2,col="black",lwd=2)
lines(c(true[4],true[4]),c(0,m),col="black",lwd=1)
legend("topright", inset=.02,box.lty=0, lty=1,lwd=c(2,2,2,2),col=c("blue","red","dark grey","black"),  
legend=c(expression(S[1]),expression(S[2]), expression(S[s]), expression(S[f])))
############### Fst 
h1=(max(g$V7,gg$V7)-min(g$V7,gg$V7))/nbin
h2=c(min(g$V7,gg$V7)-h1,max(g$V7,gg$V7)+h1)
b=c()
for(i in 1:(nbin+2))
{
b[i]=h2[1]+(i-1)*h1
}
f=hist(g$V7,breaks=b, plot=F)
f$counts=f$counts/sum(f$counts)
Fst=f
ff=hist(gg$V7,breaks=b, plot=F)
ff$counts=ff$counts/sum(ff$counts)
FFst=ff
limx=c(min(true[5],Fst$breaks,FFst$breaks),max(true[5],Fst$breaks,FFst$breaks))
plot(Fst$mids,smooth((Fst$counts+FFst$counts)/2),cex.lab=1.4,'l',xlab=expression(F[st]),ylab="", col="black",lwd=2,xlim=limx)
lines(c(true[5],true[5]),c(0,max(FFst$counts,Fst$counts)),col="black")
############### HW 
#H1
h1=(max(g$V8,gg$V8)-min(g$V8,gg$V8))/nbin
h2=c(min(g$V8,gg$V8)-h1,max(g$V8,gg$V8)+h1)
b=c()
for(i in 1:(nbin+2))
{
b[i]=h2[1]+(i-1)*h1
}
f=hist(g$V8,breaks=b, plot=F)
f$counts=f$counts/sum(f$counts)
Hw1=f
ff=hist(gg$V8,breaks=b, plot=F)
ff$counts=ff$counts/sum(ff$counts)
HHw1=ff
#H2
h1=(max(g$V9,gg$V9)-min(g$V9,gg$V9))/nbin
h2=c(min(g$V9,gg$V9)-h1,max(g$V9,gg$V9)+h1)
b=c()
for(i in 1:(nbin+2))
{
b[i]=h2[1]+(i-1)*h1
}
f=hist(g$V9,breaks=b, plot=F)
f$counts=f$counts/sum(f$counts)
Hw2=f
ff=hist(gg$V9,breaks=b, plot=F)
ff$counts=ff$counts/sum(ff$counts)
HHw2=ff
limx=c(min(true[6],Hw1$breaks,true[7],Hw2$breaks),max(true[6],Hw1$breaks,true[7],Hw2$breaks))
m=max((Hw1$counts+HHw1$counts)/2,(Hw2$counts+HHw2$counts)/2)
plot(Hw1$mids,(Hw1$counts+HHw1$counts)/2,'l',cex.lab=1.4,xlab=expression(pi),ylab="",xlim=limx,ylim=c(0,m), col="blue",lwd=2)
lines(c(true[6],true[6]),c(0,m),col="blue")
lines(Hw2$mids,(Hw2$counts+HHw2$counts)/2,col="red",lwd=2)
lines(c(true[7],true[7]),c(0,m),col="red",lwd=1)
legend("topright", inset=.02,box.lty=0, lty=1,lwd=c(2,2),col=c("blue","red"),  
legend=c(expression(pi[1]),expression(pi[2])))
############### TD 
#TD1
h1=(max(g$V10,gg$V10)-min(g$V10,gg$V10))/nbin
h2=c(min(g$V10,gg$V10)-h1,max(g$V10,gg$V10)+h1)
b=c()
for(i in 1:(nbin+2))
{
b[i]=h2[1]+(i-1)*h1
}
f=hist(g$V10,breaks=b, plot=F)
f$counts=f$counts/sum(f$counts)
TD1=f
ff=hist(gg$V10,breaks=b, plot=F)
ff$counts=ff$counts/sum(ff$counts)
TTD1=ff
#TD2
h1=(max(g$V11,gg$V11)-min(g$V11,gg$V11))/nbin
h2=c(min(g$V11,gg$V11)-h1,max(g$V11,gg$V11)+h1)
b=c()
for(i in 1:(nbin+2))
{
b[i]=h2[1]+(i-1)*h1
}
f=hist(g$V11,breaks=b, plot=F)
f$counts=f$counts/sum(f$counts)
TD2=f
ff=hist(gg$V11,breaks=b, plot=F)
ff$counts=ff$counts/sum(ff$counts)
TTD2=ff
limx=c(min(true[8],TD1$breaks,true[9],TD2$breaks),max(true[8],TD1$breaks,true[9],TD2$breaks))
m=max((TD1$counts+TTD1$counts)/2,(TD2$counts+TTD2$counts)/2)
plot(TD1$mids,(TD1$counts+TTD1$counts)/2,'l',cex.lab=1.4,xlab=expression(D),ylab="",xlim=limx,ylim=c(0,m), col="blue",lwd=2)
lines(c(true[8],true[8]),c(0,m),col="blue")
lines(TD2$mids,(TD2$counts+TTD2$counts)/2,col="red",lwd=2)
lines(c(true[9],true[9]),c(0,m),col="red",lwd=1)
legend("topright", inset=.02,box.lty=0, lty=1,lwd=c(2,2),col=c("blue","red"),  
legend=c(expression(D[1]),expression(D[2])))
########## test
i=1;
s1val=0;
while((true[1]>S1$mids[i])&&(S1$breaks[i]<max(S1$breaks)))
{
s1val=s1val+(S1$counts[i]+SS1$counts[i])/2
i=i+1;
}
if(s1val>.5)
{
s1val=1-s1val
}
i=1;
s2val=0;
while((true[2]>S2$mids[i])&&(S2$breaks[i]<max(S2$breaks)))
{
s2val=s2val+(S2$counts[i]+SS2$counts[i])/2
i=i+1;
}
if(s2val>.5)
{
s2val=1-s2val
}

i=1;
ssval=0;
while((true[3]>Ss$mids[i])&&(Ss$breaks[i]<max(Ss$breaks)))
{
ssval=ssval+(Ss$counts[i]+SSs$counts[i])/2
i=i+1;
}
if(ssval>.5)
{
ssval=1-ssval
}

i=1;
sfval=0;
while((true[4]>Sf$mids[i])&&(Sf$breaks[i]<max(Sf$breaks)))
{
sfval=sfval+(Sf$counts[i]+SSf$counts[i])/2
i=i+1;
}
if(sfval>.5)
{
sfval=1-sfval
}
i
i=1;
fstval=0;
while((true[5]>Fst$mids[i])&&(Fst$breaks[i]<max(Fst$breaks)))
{
fstval=fstval+(Fst$counts[i]+FFst$counts[i])/2
i=i+1;
}
if(fstval>.5)
{
fstval=1-fstval
}
i=1;
pi1val=0;
while((true[6]>Hw1$mids[i])&&(Hw1$breaks[i]<max(Hw1$breaks)))
{
pi1val=pi1val+(Hw1$counts[i]+HHw1$counts[i])/2
i=i+1;
}
if(pi1val>.5)
{
pi1val=1-pi1val
}
i=1;
pi2val=0;
while((true[7]>Hw2$mids[i])&&(Hw2$breaks[i]<max(Hw2$breaks)))
{
pi2val=pi2val+(Hw2$counts[i]+HHw2$counts[i])/2
i=i+1;
}
if(pi2val>.5)
{
pi2val=1-pi2val
}
i=1;
d1val=0;
while((true[8]>TD1$mids[i])&&(TD1$breaks[i]<max(TD1$breaks)))
{
d1val=d1val+(TD1$counts[i]+TTD1$counts[i])/2
i=i+1;
}
if(d1val>.5)
{
d1val=1-d1val
}
i=1;
d2val=0;
while((true[9]>TD2$mids[i])&&(TD2$breaks[i]<max(TD2$breaks)))
{
d2val=d2val+(TD2$counts[i]+TTD2$counts[i])/2
i=i+1;
}
if(d2val>.5)
{
d2val=1-d2val
}
s1val
s2val
ssval
sfval
fstval
pi1val
pi2val
d1val
d2val