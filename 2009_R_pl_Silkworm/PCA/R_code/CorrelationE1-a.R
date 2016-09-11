

postscript("Cor-a-E1.eps", paper="letter") 
lim=!is.na(P$La)

par(mfrow=c(3,2))

############# WD
plot(E$vectors[,1], P$D.W, col=colo, main="Correlation with Wild & Domesticated",xlab="Eigenvector 1", ylab="Wild/Domesticated",yaxt="n",pch=ch)
axis(2, c(1,2),c("Domesticated","Wild"))
c=cor.test(E$vectors[,1], P$D.W, method="kendall")
t=paste("***-R=",signif(c$e,3),sep="")
t=paste(t,signif(c$p.,3),sep=" p=")
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
reg=lm(P$VM[P$VM!=0]~E$vectors[,1][P$VM!=0])
abline(reg,lty=2,col=grey(.4))
anova(reg)
legend(.05,8,c(leg[1:3],t),col=c(cleg[1:3],"white"),pch=pleg[1:3], bty="n")


############# V
plot(E$vectors[,1][P$V!=0], P$V[P$V!=0], col=colo, main="Correlation with Voltinism in dom.",xlab="Eigenvector 1", ylab="# generations per year",yaxt="n",pch=ch)
axis(2, 1:3,c("1","2",">2"))
c=cor.test(E$vectors[,1][P$V!=0], P$V[P$V!=0], method="kendall")
t=paste("**-R=",signif(c$e,3),sep="")
t=paste(t,signif(c$p.,3),sep=" p=")
reg=lm(P$V[P$V!=0]~E$vectors[,1][P$V!=0])
abline(reg,lty=2,col=grey(.4))
anova(reg)
legend(.05,2.8,c(leg[1:3],t),col=c(cleg[1:3],"white"),pch=pleg[1:3], bty="n")


######### Infered sex
plot(E$vectors[,1], P$iS, col=colo, main="Correlation with infered sex in dom.",xlab="Eigenvector 1", ylab="Sex",yaxt="n",pch=ch)
axis(2, c(1,2),c("Female","Male"))
c=cor.test(E$vectors[,1], P$iS, method="kendall")
t=paste("*-R=",signif(c$e,3),sep="")
t=paste(t,signif(c$p.,3),sep=" p=")
reg=lm(P$iS~E$vectors[,1])
abline(reg,lty=2,col=grey(.4))
anova(reg)
legend("center",c(leg[1:3],t),col=c(cleg[1:3],"white"),pch=pleg[1:3], bty="n")
######### sex
plot(E$vectors[,1][P$sex!=0], P$sex[P$sex!=0], col=colo, main="Correlation with sex in dom.",xlab="Eigenvector 1", ylab="Sex",yaxt="n",pch=ch)
axis(2, c(1,2),c("Female","Male"))
c=cor.test(E$vectors[,1][P$sex!=0], P$sex[P$sex!=0], method="kendall")
t=paste("*-R=",signif(c$e,3),sep="")
t=paste(t,signif(c$p.,3),sep=" p=")
reg=lm(P$sex[P$sex!=0]~E$vectors[,1][P$sex!=0])
abline(reg,lty=2,col=grey(.4))
anova(reg)
legend(.055,1.49,c(leg[1:3],t),col=c(cleg[1:3],"white"),pch=pleg[1:3], bty="n")

############# M
plot(E$vectors[,1][P$M!=0], P$M[P$M!=0], col=colo, main="Correlation with Moulting in dom.",xlab="Eigenvector 1", ylab="Number of molts",yaxt="n",pch=ch)
axis(2, 2:5,c("2","3","4","5"))
c=cor.test(E$vectors[,1][P$M!=0], P$M[P$M!=0], method="kendall")
t=paste("*R=",signif(c$e,3),sep="")
t=paste(t,signif(c$p.,3),sep=" p=")
reg=lm(P$M[P$M!=0]~E$vectors[,1][P$M!=0])
abline(reg,lty=2,col=grey(.4))
anova(reg)
legend(.05,5.2,c(leg[1:3],t),col=c(cleg[1:3],"white"),pch=pleg[1:3], bty="n")





par(mfrow=c(2,1))
############## Long
plot(E$vectors[,1][lim], P$Lo[lim], col=colo[lim], main="Correlation with longitude",xlab="Eigenvector 1", ylab="Longitude",pch=ch[lim])
c=cor.test(E$vectors[,1][lim], P$Lo[lim], method="kendall")
t=paste("R=",signif(c$e,3),sep="")
t=paste(t,signif(c$p.,3),sep=" p=")
reg=lm(P$Lo[lim]~E$vectors[,1][lim])
#abline(reg,lty=2,col=grey(.4))
anova(reg)
legend("bottomleft",c(leg,t),col=c(cleg,"white"),pch=pleg, bty="n")

######### lat
plot(E$vectors[,1][lim], P$La[lim], col=colo[lim], main="Correlation with latitude",xlab="Eigenvector 1", ylab="Latitude",pch=ch[lim])
c=cor.test(E$vectors[,1][lim], P$La[lim], method="kendall")
t=paste("R=",signif(c$e,3),sep="")
t=paste(t,signif(c$p.,3),sep=" p=")
reg=lm(P$La[lim]~E$vectors[,1][lim])
#abline(reg,lty=2,col=grey(.4))
anova(reg)
legend("bottom",c(leg,t),col=c(cleg,"white"),pch=pleg, bty="n")





dev.off()