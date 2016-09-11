


lim=!is.na(P$La)

M=matrix(0,ncol=6,nrow=4)
M[1,1]=1
M[1,2]=kruskal.test(E$vectors[,1], P$D.W)$p.
M[1,3]=kruskal.test(E$vectors[,1][P$V!=0], P$V[P$V!=0])$p.
M[1,4]=kruskal.test(E$vectors[,1][P$M!=0], P$M[P$M!=0])$p.
M[1,5]=kruskal.test(E$vectors[,1][lim], P$Lo[lim])$p.
M[1,6]=kruskal.test(E$vectors[,1][lim], P$La[lim])$p.

M[2,1]=2
M[2,2]=kruskal.test(E$vectors[,2], P$D.W)$p.
M[2,3]=kruskal.test(E$vectors[,2][P$V!=0], P$V[P$V!=0])$p.
M[2,4]=kruskal.test(E$vectors[,2][P$M!=0], P$M[P$M!=0])$p.
M[2,5]=kruskal.test(E$vectors[,2][lim], P$Lo[lim])$p.
M[2,6]=kruskal.test(E$vectors[,2][lim], P$La[lim])$p.

M[3,1]=3
M[3,2]=kruskal.test(E$vectors[,3], P$D.W)$p.
M[3,3]=kruskal.test(E$vectors[,3][P$V!=0], P$V[P$V!=0])$p.
M[3,4]=kruskal.test(E$vectors[,3][P$M!=0], P$M[P$M!=0])$p.
M[3,5]=kruskal.test(E$vectors[,3][lim], P$Lo[lim])$p.
M[3,6]=kruskal.test(E$vectors[,3][lim], P$La[lim])$p.

M[4,1]=4
M[4,2]=kruskal.test(E$vectors[,4], P$D.W)$p.
M[4,3]=kruskal.test(E$vectors[,4][P$V!=0], P$V[P$V!=0])$p.
M[4,4]=kruskal.test(E$vectors[,4][P$M!=0], P$M[P$M!=0])$p.
M[4,5]=kruskal.test(E$vectors[,4][lim], P$Lo[lim])$p.
M[4,6]=kruskal.test(E$vectors[,4][lim], P$La[lim])$p.
write.table(M,"KW_test.xls",sep="\t", row.names=F,col.names=F)

M=matrix(0,ncol=6,nrow=4)
M[1,1]=1
M[1,2]=cor.test(E$vectors[,1], P$D.W, method="sp")$p.
M[1,3]=cor.test(E$vectors[,1][P$V!=0], P$V[P$V!=0], method="sp")$p.
M[1,4]=cor.test(E$vectors[,1][P$M!=0], P$M[P$M!=0], method="sp")$p.
M[1,5]=cor.test(E$vectors[,1][lim], P$Lo[lim], method="sp")$p.
M[1,6]=cor.test(E$vectors[,1][lim], P$La[lim], method="sp")$p.

M[2,1]=2
M[2,2]=cor.test(E$vectors[,2], P$D.W, method="sp")$p.
M[2,3]=cor.test(E$vectors[,2][P$V!=0], P$V[P$V!=0], method="sp")$p.
M[2,4]=cor.test(E$vectors[,2][P$M!=0], P$M[P$M!=0], method="sp")$p.
M[2,5]=cor.test(E$vectors[,2][lim], P$Lo[lim], method="sp")$p.
M[2,6]=cor.test(E$vectors[,2][lim], P$La[lim], method="sp")$p.

M[3,1]=3
M[3,2]=cor.test(E$vectors[,3], P$D.W, method="sp")$p.
M[3,3]=cor.test(E$vectors[,3][P$V!=0], P$V[P$V!=0], method="sp")$p.
M[3,4]=cor.test(E$vectors[,3][P$M!=0], P$M[P$M!=0], method="sp")$p.
M[3,5]=cor.test(E$vectors[,3][lim], P$Lo[lim], method="sp")$p.
M[3,6]=cor.test(E$vectors[,3][lim], P$La[lim], method="sp")$p.

M[4,1]=4
M[4,2]=cor.test(E$vectors[,4], P$D.W, method="sp")$p.
M[4,3]=cor.test(E$vectors[,4][P$V!=0], P$V[P$V!=0], method="sp")$p.
M[4,4]=cor.test(E$vectors[,4][P$M!=0], P$M[P$M!=0], method="sp")$p.
M[4,5]=cor.test(E$vectors[,4][lim], P$Lo[lim], method="sp")$p.
M[4,6]=cor.test(E$vectors[,4][lim], P$La[lim], method="sp")$p.
write.table(M,"Sperman_test.xls",sep="\t", row.names=F,col.names=F)

M=matrix(0,ncol=12,nrow=4)
M[1,1]=1
M[1,2]=cor.test(E$vectors[,1], P$D.W, method="k")$e
M[1,3]=cor.test(E$vectors[,1], P$D.W, method="k")$p.
M[1,4]=cor.test(E$vectors[,1][P$V!=0], P$V[P$V!=0], method="k")$e
M[1,5]=cor.test(E$vectors[,1][P$V!=0], P$V[P$V!=0], method="k")$p.
M[1,6]=cor.test(E$vectors[,1][P$M!=0], P$M[P$M!=0], method="k")$e
M[1,7]=cor.test(E$vectors[,1][P$M!=0], P$M[P$M!=0], method="k")$p.
M[1,8]=cor.test(E$vectors[,1][lim], P$La[lim], method="k")$e
M[1,9]=cor.test(E$vectors[,1][lim], P$La[lim], method="k")$p.
M[1,10]=cor.test(E$vectors[,1][lim], P$Lo[lim], method="k")$e
M[1,11]=cor.test(E$vectors[,1][lim], P$Lo[lim], method="k")$p.

M[2,1]=2
M[2,2]=cor.test(E$vectors[,2], P$D.W, method="k")$e
M[2,3]=cor.test(E$vectors[,2], P$D.W, method="k")$p.
M[2,4]=cor.test(E$vectors[,2][P$V!=0], P$V[P$V!=0], method="k")$e
M[2,5]=cor.test(E$vectors[,2][P$V!=0], P$V[P$V!=0], method="k")$p.
M[2,6]=cor.test(E$vectors[,2][P$M!=0], P$M[P$M!=0], method="k")$e
M[2,7]=cor.test(E$vectors[,2][P$M!=0], P$M[P$M!=0], method="k")$p.
M[2,8]=cor.test(E$vectors[,2][lim], P$La[lim], method="k")$e
M[2,9]=cor.test(E$vectors[,2][lim], P$La[lim], method="k")$p.
M[2,10]=cor.test(E$vectors[,2][lim], P$Lo[lim], method="k")$e
M[2,11]=cor.test(E$vectors[,2][lim], P$Lo[lim], method="k")$p.


M[3,1]=3
M[3,2]=cor.test(E$vectors[,3], P$D.W, method="k")$e
M[3,3]=cor.test(E$vectors[,3], P$D.W, method="k")$p.
M[3,4]=cor.test(E$vectors[,3][P$V!=0], P$V[P$V!=0], method="k")$e
M[3,5]=cor.test(E$vectors[,3][P$V!=0], P$V[P$V!=0], method="k")$p.
M[3,6]=cor.test(E$vectors[,3][P$M!=0], P$M[P$M!=0], method="k")$e
M[3,7]=cor.test(E$vectors[,3][P$M!=0], P$M[P$M!=0], method="k")$p.
M[3,8]=cor.test(E$vectors[,3][lim], P$La[lim], method="k")$e
M[3,9]=cor.test(E$vectors[,3][lim], P$La[lim], method="k")$p.
M[3,10]=cor.test(E$vectors[,3][lim], P$Lo[lim], method="k")$e
M[3,11]=cor.test(E$vectors[,3][lim], P$Lo[lim], method="k")$p.

M[4,1]=4
M[4,2]=cor.test(E$vectors[,4], P$D.W, method="k")$e
M[4,3]=cor.test(E$vectors[,4], P$D.W, method="k")$p.
M[4,4]=cor.test(E$vectors[,4][P$V!=0], P$V[P$V!=0], method="k")$e
M[4,5]=cor.test(E$vectors[,4][P$V!=0], P$V[P$V!=0], method="k")$p.
M[4,6]=cor.test(E$vectors[,4][P$M!=0], P$M[P$M!=0], method="k")$e
M[4,7]=cor.test(E$vectors[,4][P$M!=0], P$M[P$M!=0], method="k")$p.
M[4,8]=cor.test(E$vectors[,4][lim], P$La[lim], method="k")$e
M[4,9]=cor.test(E$vectors[,4][lim], P$La[lim], method="k")$p.
M[4,10]=cor.test(E$vectors[,4][lim], P$Lo[lim], method="k")$e
M[4,11]=cor.test(E$vectors[,4][lim], P$Lo[lim], method="k")$p.

write.table(M,"KE_test.xls",sep="\t", row.names=F,col.names=F)

#postscript("Cor.eps", paper="letter") 
par(mfrow=c(3,2))

############# WD
plot(E$vectors[,1], P$D.W, col=colo, main="Correlation with Wild & Domesticated",xlab="Eigenvector 1", ylab="Wild/Domesticated",yaxt="n",pch=ch)
axis(2, c(1,2),c("Domesticated","Wild"))
c=cor.test(E$vectors[,1], P$D.W, method="kendall")
t=paste("***-R=",signif(c$e,3),sep="")
t=paste(t,signif(c$p.,3),sep=" p=")
k=kruskal.test(E$vectors[,1], P$D.W)
t=paste(t,signif(k$p.,3),sep=" p=")
reg=lm(P$D.W~E$vectors[,1])
abline(reg,lty=2,col=grey(.4))
anova(reg)

legend("bottomleft",c(leg,t),col=c(cleg,"white"),pch=pleg, bty="n")

###########V/M
plot(E$vectors[,1][P$VM!=0], P$VM[P$VM!=0], col=colo, main="Correlation with Voltinism & Moulting in dom.",xlab="Eigenvector 1", ylab="V/M",yaxt="n",pch=ch)
#axis(2, 1:8,c("V1M3","V1M4","V1M5","V2M2","V2M3","V2M4","V3M3","V3M4"))
axis(2, seq(1,8,by=2),c("V1M3","V1M5","V2M3","V3M3"))
axis(4, seq(2,8,by=2),c("V1M4","V2M2","V2M4","V3M4"))
c=cor.test(E$vectors[,1][P$VM!=0], P$VM[P$VM!=0], method="kendall")
t=paste("***-R=",signif(c$e,3),sep="")
t=paste(t,signif(c$p.,3),sep=" p=")
k=kruskal.test(E$vectors[,1][P$VM!=0], P$VM[P$VM!=0])
t=paste(t,signif(k$p.,3),sep=" p=")
reg=lm(P$VM[P$VM!=0]~E$vectors[,1][P$VM!=0])
abline(reg,lty=2,col=grey(.4))
anova(reg)
legend(.05,8,c(leg[1:3],t),col=c(cleg[1:3],"white"),pch=pleg[1:3], bty="n")

######### Infered sex
plot(E$vectors[,1], P$iS, col=colo, main="Correlation with infered sex in dom.",xlab="Eigenvector 1", ylab="Sex",yaxt="n",pch=ch)
axis(2, c(1,2),c("Female","Male"))
c=cor.test(E$vectors[,1], P$iS, method="kendall")
t=paste("*-R=",signif(c$e,3),sep="")
t=paste(t,signif(c$p.,3),sep=" p=")
k=kruskal.test(E$vectors[,1], P$iS)
t=paste(t,signif(k$p.,3),sep=" p=")
reg=lm(P$iS~E$vectors[,1])
abline(reg,lty=2,col=grey(.4))
anova(reg)
legend("center",c(leg[1:3],t),col=c(cleg[1:3],"white"),pch=pleg[1:3], bty="n")

############# V
plot(E$vectors[,1][P$V!=0], P$V[P$V!=0], col=colo, main="Correlation with Voltinism in dom.",xlab="Eigenvector 1", ylab="# generations per year",yaxt="n",pch=ch)
axis(2, 1:3,c("1","2",">2"))
c=cor.test(E$vectors[,1][P$V!=0], P$V[P$V!=0], method="kendall")
t=paste("**-R=",signif(c$e,3),sep="")
t=paste(t,signif(c$p.,3),sep=" p=")
k=kruskal.test(E$vectors[,1][P$V!=0], P$V[P$V!=0])
t=paste(t,signif(k$p.,3),sep=" p=")
reg=lm(P$V[P$V!=0]~E$vectors[,1][P$V!=0])
abline(reg,lty=2,col=grey(.4))
anova(reg)
legend(.05,2.8,c(leg[1:3],t),col=c(cleg[1:3],"white"),pch=pleg[1:3], bty="n")

############# M
plot(E$vectors[,1][P$M!=0], P$M[P$M!=0], col=colo, main="Correlation with Moulting in dom.",xlab="Eigenvector 1", ylab="Number of molts",yaxt="n",pch=ch)
axis(2, 2:5,c("2","3","4","5"))
c=#cor.test(E$vectors[,1][P$M!=0], P$M[P$M!=0], method="kendall")
t=paste("*R=",signif(c$e,3),sep="")
t=paste(t,signif(c$p.,3),sep=" p=")
k=kruskal.test(E$vectors[,1][P$M!=0], P$M[P$M!=0])
t=paste(t,signif(k$p.,3),sep=" p=")
reg=lm(P$M[P$M!=0]~E$vectors[,1][P$M!=0])
abline(reg,lty=2,col=grey(.4))
anova(reg)
legend(.05,5.2,c(leg[1:3],t),col=c(cleg[1:3],"white"),pch=pleg[1:3], bty="n")

######### sex
plot(E$vectors[,1][P$sex!=0], P$sex[P$sex!=0], col=colo, main="Correlation with sex in dom.",xlab="Eigenvector 1", ylab="Sex",yaxt="n",pch=ch)
axis(2, c(1,2),c("Female","Male"))
c=cor.test(E$vectors[,1][P$sex!=0], P$sex[P$sex!=0], method="kendall")
t=paste("*-R=",signif(c$e,3),sep="")
t=paste(t,signif(c$p.,3),sep=" p=")
k=kruskal.test(E$vectors[,1][P$sex!=0], P$sex[P$sex!=0])
t=paste(t,signif(k$p.,3),sep=" p=")
reg=lm(P$sex[P$sex!=0]~E$vectors[,1][P$sex!=0])
abline(reg,lty=2,col=grey(.4))
anova(reg)
legend(.055,1.49,c(leg[1:3],t),col=c(cleg[1:3],"white"),pch=pleg[1:3], bty="n")




par(mfrow=c(2,1))
############## Long
plot(E$vectors[,1][lim], P$Lo[lim], col=colo[lim], main="Correlation with longitude",xlab="Eigenvector 1", ylab="Longitude",pch=ch[lim])
c=cor.test(E$vectors[,1][lim], P$Lo[lim], method="kendall")
t=paste("R=",signif(c$e,3),sep="")
t=paste(t,signif(c$p.,3),sep=" p=")
k=kruskal.test(E$vectors[,1][lim], P$Lo[lim])
t=paste(t,signif(k$p.,3),sep=" p=")
reg=lm(P$Lo[lim]~E$vectors[,1][lim])
#abline(reg,lty=2,col=grey(.4))
anova(reg)
legend("bottomleft",c(leg,t),col=c(cleg,"white"),pch=pleg, bty="n")

######### lat
plot(E$vectors[,1][lim], P$La[lim], col=colo[lim], main="Correlation with latitude",xlab="Eigenvector 1", ylab="Latitude",pch=ch[lim])
c=cor.test(E$vectors[,1][lim], P$La[lim], method="kendall")
t=paste("R=",signif(c$e,3),sep="")
t=paste(t,signif(c$p.,3),sep=" p=")
k=kruskal.test(E$vectors[,1][lim], P$La[lim])
t=paste(t,signif(k$p.,3),sep=" p=")
reg=lm(P$La[lim]~E$vectors[,1][lim])
#abline(reg,lty=2,col=grey(.4))
anova(reg)
legend("bottom",c(leg,t),col=c(cleg,"white"),pch=pleg, bty="n")





#dev.off()