rm(list = ls())

pdf("FSFS-auto-50.pdf")

col = c("blue", "pink")
Xw=read.table("28u-FSFS18-wild")
Xd=read.table("28u-FSFS48-dom")
mc=max(dim(Xw),dim(Xd))
M=matrix(0,ncol=,mc,nrow=2)
for(i in 1:dim(Xw)[1])
{M[1,i]=Xw$V2[i]	}

for(i in 1:dim(Xd)[1])
{M[2,i]=Xd$V2[i]	}

barplot2(M,beside = TRUE, col=col, name= Xd$V1,legend.text=c("11 Wild", "29 Domesticated"),ylab="count", xlab="Folded Site Frequencies",  main="Freq. on Autosomes, with <50% missing data")
##################### foled
col = c("blue", "pink")
Xw=read.table("28f-FSFS18-wild")
Xd=read.table("28f-FSFS18-dom")
mc=max(dim(Xw),dim(Xd))
M=matrix(0,ncol=,mc,nrow=2)
for(i in 1:dim(Xw)[1])
{M[1,i]=Xw$V2[i]	}
for(i in 1:dim(Xd)[1])
{M[2,i]=Xd$V2[i]	}
barplot2(M,beside = TRUE, col=col, name= Xd$V1,legend.text=c("11 Wild", "29 Domesticated"),ylab="count", xlab="Folded Site Frequencies",  main="Freq. on Autosomes (<50% MD), projected on 18 chr.")

############## total

Xu=read.table("28u-FSFS66-total")
Xf=read.table("28f-FSFS18-total")
mc=max(dim(Xu),dim(Xf))
M=matrix(0,ncol=,mc,nrow=2)
for(i in 1:mc)
{M[1,i]=Xu$V2[i]	}
for(i in 1:mc)
{M[2,i]=Xf$V2[i]	}

barplot2(M,beside = TRUE,col=col, name= Xu$V1,legend.text=c("40 indv.", "40 ind. projected on 18 chr."),ylab="count", xlab="Folded Site Frequencies",  main="Freq. on Autosomes, with <50% missing data")

dev.off()
##########################
#########################

rm(list = ls())
col = c("blue", "pink")
pdf("FSFS-AutoUN-50.pdf")
Xw=read.table("29u-FSFS18-wild")
Xd=read.table("29u-FSFS48-dom")
mc=max(dim(Xw),dim(Xd))
M=matrix(0,ncol=,mc,nrow=2)
for(i in 1:dim(Xw)[1])
{M[1,i]=Xw$V2[i]	}
for(i in 1:dim(Xd)[1])
{M[2,i]=Xd$V2[i]	}
barplot2(M,beside = TRUE, col=col, name= Xd$V1,legend.text=c("11 Wild", "29 Domesticated"),ylab="count", xlab="Folded Site Frequencies",  main="Freq. on Autosomes + UN, with <50% missing data")

################ folded
Xw=read.table("29f-FSFS18-wild")
Xd=read.table("29f-FSFS18-dom")
mc=max(dim(Xw),dim(Xd))
M=matrix(0,ncol=,mc,nrow=2)
for(i in 1:dim(Xw)[1])
{M[1,i]=Xw$V2[i]	}
for(i in 1:dim(Xd)[1])
{M[2,i]=Xd$V2[i]	}

barplot2(M,beside = TRUE, col=col, name= Xd$V1,legend.text=c("11 Wild", "29 Domesticated"),ylab="count", xlab="Folded Site Frequencies",  main="Freq. on Autosomes + UN (<50% MD), projected on 18 chr.")

############## total
Xu=read.table("29u-FSFS66-total")
Xf=read.table("29f-FSFS18-total")
mc=max(dim(Xu),dim(Xf))
M=matrix(0,ncol=,mc,nrow=2)
for(i in 1:mc)
{M[1,i]=Xu$V2[i]	}
for(i in 1:mc)
{M[2,i]=Xf$V2[i]	}

barplot2(M,beside = TRUE, col=col,name= Xu$V1,legend.text=c("40 indv.", "40 ind. projected on 18 chr."),ylab="count", xlab="Folded Site Frequencies",  main="Freq. on Autosomes + UN, with <50% missing data")

dev.off()