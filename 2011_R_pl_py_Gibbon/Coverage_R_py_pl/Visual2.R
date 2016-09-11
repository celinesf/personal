#just one figure per sample?
#Please tabulate: mean coverage, median coverage, proportion of sites
#with no coverage, proportion of sites with >=5x coverage and proportion
#of sites with >=10x coverage.  Thanks!


file="gibbon_coverage_vo"
OUT="VOK"
fpdf=paste("SUM_cov",OUT,sep="_")
fpdf=paste(fpdf,"pdf",sep=".")
ftext=paste("SUM_cov",OUT,sep="_")

for(i in 0:97)
{
a=scan(file,what="c",nlines=1,skip=i)
b=as.numeric(a[2:length(a)])
d=c(d,c)
}

pdf(fpdf)
br=c(0,5,10,max(d))
h=hist(d,breaks=br,main=file,xlab="coverage",plot=T)
dev.off()

x=matrix(0,nc=5,nr=1)
x[1,1]=mean(d)
x[1,2]=quantile(d,p=.5)
x[1,3]=h$counts[1]/sum(h$counts)
x[1,4]=(h$counts[2]+h$counts[3])/sum(h$counts)
x[1,5]=h$count[3]/sum(h$counts)

write(x,f=ftext)



############
br=c(0,5,10,max(b))
h=hist(b,breaks=br,main=file,xlab="coverage",plot=T)

if(i==0)	{fh=h
	}else
	{
		fh$breaks[4]=max(fh$breaks[4],h$breaks[4])
		fh$mids[3]=max(fh$mids[3],h$mids[3])
		fh$counts=fh$counts+h$counts
		fh$intensities=fh$intensities+h$intensities
		fh$density=fh$density+h$density
	}