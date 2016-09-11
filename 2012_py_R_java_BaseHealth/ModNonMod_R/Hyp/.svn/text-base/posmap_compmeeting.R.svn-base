hyp <- UImap[UImap$DiseaseID==1,]
hyp2 <- read.csv("/home/pouria3/Dropbox/EngExl/T2D/T2D_SNP_UIMapping.csv")
jpeg("/home/pouria3/Dropbox/EngExl/T2D/pos.jpg", width = 1080, height = 480)
par(lty=2, lwd=2, mfrow=c(1, 2))
plot(hyp$Chromosome, hyp$Position, col="red")
plot(as.character(hyp2$Chromosome), hyp2$Position, col='blue', type="p")
axis(2, c(9:0))
dev.off()