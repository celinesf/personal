br=c(0,1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50)


pdf("Read_sizes.pdf")

a=read.table("BPjo")$V1
h=hist(a,main="Joanna",xlab="bp",plot=T)


a=read.table("BPvo")$V1
h=hist(a,main="Vok",xlab="bp",plot=T)

a=read.table("BPas")$V1
h=hist(a,main="Asia",xlab="bp",plot=T)


a=read.table("BPch")$V1
h=hist(a,main="China",xlab="bp",plot=T)


dev.off()




