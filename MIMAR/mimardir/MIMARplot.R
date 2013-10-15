# INSTRUCTIONS
# At any points during the analysis, copy the output files from 
# two different seeds for the same input file and rename the copy by
# "outputmimar1" and "outputmimar2". 
# If the analyses were not completed, open these two files and delete 
# the last incomplete rows, and save.
# Then execute the following lines in R.
# The five posterior distributions will be printed (theta_1, theta_2, 
# T, theta_A, M) for both seeds (in red and blue).
# - You can change the number of bin in each histogram by 
# changing the value of the object "b" below.
# - Note that if you recorded a lot of steps, R might take a long time downloading the results.
# Recording only every "int" steps does not affect the results drasticaly,
# and allow to read the file much faster.
# - Note (4/24/2008 - Thanks Stephen Wright and Tanja Slotte for mentionning the problem)
# The bin for the migration rate now increases exponentially in size. 
# Thus now the histogram from the summary file and the output file looks the same for M
# ! Comment and uncoment the lines mentionned below if you use a uniform prior on M.
# This will make the bins for M uniform in size.

rm(list=ls())                          # delete previous objects in R

########### Setting parameteres ##########
# change the name of MIMAR output file if needed
g=read.table("outputmimar1", skip=12)  # seed1 results
gg=read.table("outputmimar2", skip=12) # seed2 results
b=1000                                 # number of bin in histogram
##### End of parameters and option setting ######

par(mfrow=c(3,2))
############### Theta1 
f=hist(g$V2,breaks=b, plot=F)
f$counts=f$counts/sum(f$counts)
ff=hist(gg$V2,breaks=b, plot=F)
ff$counts=ff$counts/sum(ff$counts)
q1=f
qq1=ff
m=max(smooth(q1$counts),smooth(qq1$counts))
lim=c(min(q1$mids,qq1$mids),max(q1$mids,qq1$mids))
plot(q1$mids,smooth(q1$counts),cex.lab=1.5,'l',xlab=expression(theta[1]),ylab="", col="red",ylim=c(0,m),xlim=lim)
lines(qq1$mids,smooth(qq1$counts), col="blue")
############### Theta2 
f=hist(g$V3,breaks=b, plot=F)
f$counts=f$counts/sum(f$counts)
ff=hist(gg$V3,breaks=b, plot=F)
ff$counts=ff$counts/sum(ff$counts)
q2=f
qq2=ff
m=max(smooth(q2$counts),smooth(qq2$counts))
lim=c(min(q2$mids,qq2$mids),max(q2$mids,qq1$mids))
plot(q2$mids,smooth(q2$counts),cex.lab=1.5,'l',xlab=expression(theta[2]),ylab="", col="red",ylim=c(0,m),xlim=lim)
lines(qq2$mids,smooth(qq2$counts), col="blue")
############### ThetaA 
f=hist(g$V6,breaks=b, plot=F)
f$counts=f$counts/sum(f$counts)
qa=f
ff=hist(gg$V6,breaks=b, plot=F)
ff$counts=ff$counts/sum(ff$counts)
qqa=ff
m=max(smooth(qa$counts),smooth(qqa$counts))
lim=c(min(qa$mids,qqa$mids),max(qa$mids,qqa$mids))
plot(qa$mids,smooth(qa$counts),cex.lab=1.5,'l',xlab=expression(theta[A]),ylab="", col="red",ylim=c(0,m),xlim=lim)
lines(qqa$mids,smooth(qqa$counts), col="blue")
############### Migration rate
# bm=b	                                        #uncomment if the prior on M is uniform
bm=c()	                                        # comment if the prior on M is uniform
h2=c(log(max(g$V7,gg$V7)),log(min(g$V7,gg$V7)))	# comment if the prior on M is uniform
h1=(h2[2]-h2[1])/b					            # comment if the prior on M is uniform
for(i in 1:(b+1)){bm[i]=exp(h2[2]-(i-1)*h1)}	# comment if the prior on M is uniform

f=hist(g$V7,breaks=bm, plot=F)
f$counts=f$counts/sum(f$counts)
ff=hist(gg$V7,breaks=bm, plot=F)
ff$counts=ff$counts/sum(ff$counts)
m=f
mm=ff
mt=max(smooth(m$counts),smooth(mm$counts))
lim=c(min(m$mids,mm$mids),max(m$mids,mm$mids))
plot(m$mids,smooth(m$counts),cex.lab=1.5,'l',xlab=expression(M),ylab="", col="red",ylim=c(0,mt),xlim=lim)
lines(mm$mids,smooth(mm$counts), col="blue")
############### Time in generation
f=hist(g$V5,breaks=b, plot=F)
f$counts=f$counts/sum(f$counts)
t=f
ff=hist(gg$V5,breaks=b, plot=F)
ff$counts=ff$counts/sum(ff$counts)
tt=ff
mt=max(smooth(t$counts),smooth(tt$counts))
lim=c(min((t$mids),(tt$mids)),max((t$mids),(tt$mids)))
plot(t$mids,smooth(t$counts),cex.lab=1.5,'l',xlab=expression(T[gen]),ylab="", col="red",ylim=c(0,mt),xlim=lim)
lines(tt$mids,smooth(tt$counts), col="blue")



####### mode and 90th precentiles (rows) for
# theta_1 (col 1), theta_2 (col 2), theta_A (col 3), M (col 4), and T_gen(col 5)
res=matrix(nrow=3,ncol=5)
res[1,1]=mean(c(q1$mids[order(q1$counts,decreasing=T)[1]],qq1$mids[order(qq1$counts,decreasing=T)[1]]))
res[2,1]=mean(c(quantile(g$V2,probs=.05),quantile(gg$V2,probs=.05)))
res[3,1]=mean(c(quantile(g$V2,probs=.95),quantile(gg$V2,probs=.95)))
res[1,2]=mean(c(q2$mids[order(q2$counts,decreasing=T)[1]],qq2$mids[order(qq2$counts,decreasing=T)[1]]))
res[2,2]=mean(c(quantile(g$V3,probs=.05),quantile(gg$V3,probs=.05)))
res[3,2]=mean(c(quantile(g$V3,probs=.95),quantile(gg$V3,probs=.95)))
res[1,3]=mean(c(qa$mids[order(qa$counts,decreasing=T)[1]],qqa$mids[order(qqa$counts,decreasing=T)[1]]))
res[2,3]=mean(c(quantile(g$V6,probs=.05),quantile(gg$V6,probs=.05)))
res[3,3]=mean(c(quantile(g$V6,probs=.95),quantile(gg$V6,probs=.95)))
res[1,4]=mean(c(m$mids[order(m$counts,decreasing=T)[1]],mm$mids[order(mm$counts,decreasing=T)[1]]))
res[2,4]=mean(c(quantile(g$V7,probs=.05),quantile(gg$V7,probs=.05)))
res[3,4]=mean(c(quantile(g$V7,probs=.95),quantile(gg$V7,probs=.95)))
res[1,5]=mean(c(t$mids[order(t$counts,decreasing=T)[1]],tt$mids[order(tt$counts,decreasing=T)[1]]))
res[2,5]=mean(c(quantile(g$V5,probs=.05),quantile(gg$V5,probs=.05)))
res[3,5]=mean(c(quantile(g$V5,probs=.95),quantile(gg$V5,probs=.95)))
restot=res
res=matrix(nrow=3,ncol=5)
res[1,1]=q1$mids[order(q1$counts,decreasing=T)[1]]
res[2,1]=quantile(g$V2,probs=.05)
res[3,1]=quantile(g$V2,probs=.95)
res[1,2]=q2$mids[order(q2$counts,decreasing=T)[1]]
res[2,2]=quantile(g$V3,probs=.05)
res[3,2]=quantile(g$V3,probs=.95)
res[1,3]=qa$mids[order(qa$counts,decreasing=T)[1]]
res[2,3]=quantile(g$V6,probs=.05)
res[3,3]=quantile(g$V6,probs=.95)
res[1,4]=m$mids[order(m$counts,decreasing=T)[1]]
res[2,4]=quantile(g$V7,probs=.05)
res[3,4]=quantile(g$V7,probs=.95)
res[1,5]=t$mids[order(t$counts,decreasing=T)[1]]
res[2,5]=quantile(g$V5,probs=.05)
res[3,5]=quantile(g$V5,probs=.95)
res1=res
res=matrix(nrow=3,ncol=5)
res[1,1]=qq1$mids[order(qq1$counts,decreasing=T)[1]]
res[2,1]=quantile(gg$V2,probs=.05)
res[3,1]=quantile(gg$V2,probs=.95)
res[1,2]=qq2$mids[order(qq2$counts,decreasing=T)[1]]
res[2,2]=quantile(gg$V3,probs=.05)
res[3,2]=quantile(gg$V3,probs=.95)
res[1,3]=qqa$mids[order(qqa$counts,decreasing=T)[1]]
res[2,3]=quantile(gg$V6,probs=.05)
res[3,3]=quantile(gg$V6,probs=.95)
res[1,4]=mm$mids[order(mm$counts,decreasing=T)[1]]
res[2,4]=quantile(gg$V7,probs=.05)
res[3,4]=quantile(gg$V7,probs=.95)
res[1,5]=tt$mids[order(tt$counts,decreasing=T)[1]]
res[2,5]=quantile(gg$V5,probs=.05)
res[3,5]=quantile(gg$V5,probs=.95)
res2=res


### estimates from seed1
res1
### estimates from seed2
res2
### average estimates over the two seeds
restot
