rm(list = ls())


#columns are snps (0,1,2)
#rows are individuals
 
eigenstratcel<-function(chr){
name=paste("X",chr,sep="-")	#get X matrix	
X=read.table(name,stringsAsFactors = F)
E<-eigen(X)
name=paste("SNPIND",chr,sep="-")#get info for sign
info=read.table(name,stringsAsFactors = F)
snp=info[1]
ind=info[2]		
mu<-(sqrt(snp-1)+sqrt(ind))^2/snp
sigma<-(sqrt(snp-1)+sqrt(ind))/snp*(1/sqrt(snp-1)+1/sqrt(ind))^(1/3)
for(i in 1:ind[[1]])
{
E$TW[i]<-(E$values[i]*ind/sum(E$values)-mu)/sigma
}
E$mu<-mu
E$sigma<-sigma
E$snp=snp
E$ind=ind		
class(E)<-"eigenstratcel"
E
}

plot.eigenstratcel<-function(x,...){plot(x$vectors[,1:2],...)}
print.eigenstratcel<-function(x){print(x$TW)}






chr="sum-auto"
file1=paste("Fig-PCA",chr,sep="-")


E<-eigenstratcel(chr)
eig=paste("eigen",chr,sep="-")
write(E$values,eig,n=1)
eig=paste("eigenvector",chr,sep="-")
write(t(E$vectors[,1:6]),eig,n=6,sep="\t")
P=read.table("IndvPheno-Latlon",header=T)

colo=c()
ch=c()
for(i in 1:40)
{
	if(P$V[i]==1)
	{colo=c(colo,"cyan4")
	ch=c(ch,3)	}
	else if(P$V[i]==2)
	{colo=c(colo,"slateblue")
		ch=c(ch,4)}
	else if(P$V[i]==3)
	{colo=c(colo,"darkblue")
		ch=c(ch,5)}
	#else if(i==30 |i==34)
	#{colo=c(colo,grey(.4))
	#ch=c(ch,2)}
	else
	{colo=c(colo,"black")
		ch=c(ch,1)}
}#japan

leg=c("Domesticated-V=1","Domesticated-V=2","Domesticated-V>2","Wild")
cleg=c("cyan4","slateblue","darkblue","black")
pleg=c(3,4,5,1)


par(mfrow=c(1,1))
####### 1 2

postscript("Fig1.eps",paper="special",width = 8.0, height = 8.0,)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(E$vectors[,1:2], col=colo, xlab="Eigenvector 1", ylab="Eigenvector 2",pch=ch,cex=1.5,cex.axis=1.5,cex.lab=1.5)
legend("topleft",leg,cex=1.5,col=cleg,pch=pleg, bty="n")
dev.off()


colo=c()
ch=c()
for(i in 1:40)
{
	if(P$V[i]==1)
	{colo=c(colo,"cyan4")
	ch=c(ch,3)	}
	else if(P$V[i]==2)
	{colo=c(colo,"slateblue")
		ch=c(ch,4)
	}
	else if(P$V[i]==3)
	{colo=c(colo,"darkblue")
		ch=c(ch,5)}
	else if(i==30 |i==34)
	{colo=c(colo,"black")
	ch=c(ch,6)}
	else
	{colo=c(colo,"black")
		ch=c(ch,1)}
}#japan

leg=c("Domesticated-V=1","Domesticated-V=2","Domesticated-V>3","Wild","B.m. Ziyang & Pengshan ")
cleg=c("cyan4","slateblue","darkblue","black","black")
pleg=c(3,4,5,1,6)
######### 1/4
postscript("FigS2.eps",paper="special",width = 8.0, height = 8.0,)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(E$vectors[,1],E$vectors[,4], col=colo,xlab="Eigenvector 1", ylab="Eigenvector 4",pch=ch,cex=1.5,cex.axis=1.5,cex.lab=1.5)
legend("bottomright",leg,cex=1.5,col=cleg,pch=pleg, bty="n")
dev.off()


########## 2/4
postscript("FigS4.eps",paper="special",width = 8.0, height = 8.0,)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(E$vectors[,2], E$vectors[,4],col=colo,xlab="Eigenvector 2", ylab="Eigenvector 4",pch=ch,cex=1.5,cex.axis=1.5,cex.lab=1.5)
legend("right",leg,cex=1.5,col=cleg,pch=pleg, bty="p",box.lty=2)
dev.off()







colo=c()
ch=c()
for(i in 1:40)
{
	if(P$V[i]==1)
	{colo=c(colo,"cyan4")
	ch=c(ch,3)	}
	else if(P$V[i]==2)
	{colo=c(colo,"slateblue")
		if(i==1 |i==3)
	#{colo=c(colo,grey(.4))
	{ch=c(ch,2)}
	else{ch=c(ch,4)}
	}
	else if(P$V[i]==3)
	{colo=c(colo,"darkblue")
		ch=c(ch,5)}
	#else if(i==1 |i==34)
	#{colo=c(colo,grey(.4))
	#ch=c(ch,2)}
	else
	{colo=c(colo,"black")
		ch=c(ch,1)}
}#japan

leg=c("Domesticated-V=1","Domesticated-V=2","Domesticated-V>3","Wild","J7532 & J872 ")
cleg=c("cyan4","slateblue","darkblue","black","slateblue")
pleg=c(3,4,5,1,2)

######### 1/3

postscript("FigS1.eps",paper="special",width = 8.0, height = 8.0,)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(E$vectors[,1],E$vectors[,3], col=colo,xlab="Eigenvector 1", ylab="Eigenvector 3",pch=ch,cex=1.5,cex.axis=1.5,cex.lab=1.5)
legend("topleft",leg,cex=1.5,col=cleg,pch=pleg, bty="n")
dev.off()

############ 2/3
postscript("FigS3.eps",paper="special",width = 8.0, height = 8.0,)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(E$vectors[,2:3], col=colo,xlab="Eigenvector 2", ylab="Eigenvector 3",pch=ch,cex=1.5,cex.axis=1.5,cex.lab=1.5)
legend("topright",leg,cex=1.5,col=cleg,pch=pleg, bty="n")
dev.off()


colo=c()
ch=c()
for(i in 1:40)
{
	if(P$V[i]==1)
	{colo=c(colo,"cyan4")
	ch=c(ch,3)	}
	else if(P$V[i]==2)
	{colo=c(colo,"slateblue")
		if(i==1 |i==3)
	#{colo=c(colo,grey(.4))
	{ch=c(ch,2)}
	else{ch=c(ch,4)}
	}
	else if(P$V[i]==3)
	{colo=c(colo,"darkblue")
		ch=c(ch,5)}
	else if(i==30 |i==34)
	{colo=c(colo,"black")
	ch=c(ch,6)}
	else
	{colo=c(colo,"black")
		ch=c(ch,1)}
}#japan

leg=c("Domesticated-V=1","Domesticated-V=2","Domesticated-V>3","Wild","J7532 & J872 ","B.m. Ziyang & Pengshan ")
cleg=c("cyan4","slateblue","darkblue","black","slateblue","black")
pleg=c(3,4,5,1,2,6)
######
postscript("Fig2.eps",paper="special",width = 8.0, height = 8.0,)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(E$vectors[,3:4], col=colo,xlab="Eigenvector 3", ylab="Eigenvector 4",pch=ch,cex=1.5,cex.axis=1.5,cex.lab=1.5)
legend("bottomright",leg,col=cleg,cex=1.5,pch=pleg,bty="n",box.lty=2)
dev.off()

dev.off()

