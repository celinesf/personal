rm(list = ls())

pdf("FSFS-auto-95.pdf")

col = c("blue", "pink")
Xw=read.table("28u-FSFS8-wild")
Xd=read.table("28u-FSFS24-dom")
mc=max(dim(Xw),dim(Xd))
M=matrix(0,ncol=,mc,nrow=2)
for(i in 1:dim(Xw)[1])
{M[1,i]=Xw$V2[i]	}

for(i in 1:dim(Xd)[1])
{M[2,i]=Xd$V2[i]	}

barplot2(M,beside = TRUE, col=col, name= Xd$V1,legend.text=c("11 Wild", "29 Domesticated"),ylab="count", xlab="Folded Site Frequencies",  main="Freq. on Autosomes, with <95% missing data")
##################### foled
col = c("blue", "pink")
Xw=read.table("28f-FSFS8-wild")
Xd=read.table("28f-FSFS8-dom")
mc=max(dim(Xw),dim(Xd))
M=matrix(0,ncol=,mc,nrow=2)
for(i in 1:dim(Xw)[1])
{M[1,i]=Xw$V2[i]	}
for(i in 1:dim(Xd)[1])
{M[2,i]=Xd$V2[i]	}
barplot2(M,beside = TRUE, col=col, name= Xd$V1,legend.text=c("11 Wild", "29 Domesticated"),ylab="count", xlab="Folded Site Frequencies",  main="Freq. on Autosomes (<95% MD), projected on 8 chr.")

############## total

Xu=read.table("28u-FSFS36-total")
Xf=read.table("28f-FSFS8-total")
mc=max(dim(Xu),dim(Xf))
M=matrix(0,ncol=,mc,nrow=2)
for(i in 1:mc)
{M[1,i]=Xu$V2[i]	}
for(i in 1:mc)
{M[2,i]=Xf$V2[i]	}

barplot2(M,beside = TRUE,col=col, name= Xu$V1,legend.text=c("40 indv.", "40 ind. projected on 8 chr."),ylab="count", xlab="Folded Site Frequencies",  main="Freq. on Autosomes, with <95% missing data")

dev.off()
##########################
#########################

rm(list = ls())
col = c("blue", "pink")
pdf("FSFS-AutoUN-95.pdf")
Xw=read.table("29u-FSFS8-wild")
Xd=read.table("29u-FSFS24-dom")
mc=max(dim(Xw),dim(Xd))
M=matrix(0,ncol=,mc,nrow=2)
for(i in 1:dim(Xw)[1])
{M[1,i]=Xw$V2[i]	}
for(i in 1:dim(Xd)[1])
{M[2,i]=Xd$V2[i]	}
barplot2(M,beside = TRUE, col=col, name= Xd$V1,legend.text=c("11 Wild", "29 Domesticated"),ylab="count", xlab="Folded Site Frequencies",  main="Freq. on Autosomes + UN, with <95% missing data")

################ folded
Xw=read.table("29f-FSFS8-wild")
Xd=read.table("29f-FSFS8-dom")
mc=max(dim(Xw),dim(Xd))
M=matrix(0,ncol=,mc,nrow=2)
for(i in 1:dim(Xw)[1])
{M[1,i]=Xw$V2[i]	}
for(i in 1:dim(Xd)[1])
{M[2,i]=Xd$V2[i]	}

barplot2(M,beside = TRUE, col=col, name= Xd$V1,legend.text=c("11 Wild", "29 Domesticated"),ylab="count", xlab="Folded Site Frequencies",  main="Freq. on Autosomes + UN (<95% MD), projected on 8 chr.")

############## total
Xu=read.table("29u-FSFS36-total")
Xf=read.table("29f-FSFS8-total")
mc=max(dim(Xu),dim(Xf))
M=matrix(0,ncol=,mc,nrow=2)
for(i in 1:mc)
{M[1,i]=Xu$V2[i]	}
for(i in 1:mc)
{M[2,i]=Xf$V2[i]	}

barplot2(M,beside = TRUE, col=col,name= Xu$V1,legend.text=c("40 indv.", "40 ind. projected on 8 chr."),ylab="count", xlab="Folded Site Frequencies",  main="Freq. on Autosomes + UN, with <95% missing data")

dev.off()