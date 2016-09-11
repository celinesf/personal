#################################
#BP
#################################
bpq<-read.csv ("bpq_a.csv", header=T);
bp=cbind(as.numeric(as.character(bpq$SEQN)),as.numeric(as.character(bpq$BPQ020)),as.numeric(as.character(bpq$BPQ040A)));
bpq<-read.csv ("bpq_b.csv", header=T);
bp=rbind(bp,cbind(as.numeric(as.character(bpq$SEQN)),as.numeric(as.character(bpq$BPQ020)),as.numeric(as.character(bpq$BPQ040A))));
bpq<-read.csv ("bpq_c.csv", header=T);
bp=rbind(bp,cbind(as.numeric(as.character(bpq$SEQN)),as.numeric(as.character(bpq$BPQ020)),as.numeric(as.character(bpq$BPQ040A))));
bpq<-read.csv ("bpq_d.csv", header=T);
bp=rbind(bp,cbind(as.numeric(as.character(bpq$SEQN)),as.numeric(as.character(bpq$BPQ020)),as.numeric(as.character(bpq$BPQ040A))));
bp=data.frame(bp);
names(bp)[1:3]=c('SEQN','HBP','RX');
bpx<-read.csv ("bpx_a.csv", header=T);
bx=cbind(as.numeric(as.character(bpx$SEQN)),as.numeric(as.character(bpx$BPXDI1)),as.numeric(as.character(bpx$BPXSY1)),
as.numeric(as.character(bpx$BPXDI2)),as.numeric(as.character(bpx$BPXSY2)),as.numeric(as.character(bpx$BPXDI3)),
as.numeric(as.character(bpx$BPXSY3)),as.numeric(as.character(bpx$BPXDI4)),as.numeric(as.character(bpx$BPXSY4)));
bpx<-read.csv ("bpx_b.csv", header=T);
bx=rbind(bx,cbind(as.numeric(as.character(bpx$SEQN)),as.numeric(as.character(bpx$BPXDI1)),as.numeric(as.character(bpx$BPXSY1)),
as.numeric(as.character(bpx$BPXDI2)),as.numeric(as.character(bpx$BPXSY2)),as.numeric(as.character(bpx$BPXDI3)),
as.numeric(as.character(bpx$BPXSY3)),as.numeric(as.character(bpx$BPXDI4)),as.numeric(as.character(bpx$BPXSY4))));
bpx<-read.csv ("bpx_c.csv", header=T);
bx=rbind(bx,cbind(as.numeric(as.character(bpx$SEQN)),as.numeric(as.character(bpx$BPXDI1)),as.numeric(as.character(bpx$BPXSY1)),
as.numeric(as.character(bpx$BPXDI2)),as.numeric(as.character(bpx$BPXSY2)),as.numeric(as.character(bpx$BPXDI3)),
as.numeric(as.character(bpx$BPXSY3)),as.numeric(as.character(bpx$BPXDI4)),as.numeric(as.character(bpx$BPXSY4))));
bpx<-read.csv ("bpx_d.csv", header=T);
bx=rbind(bx,cbind(as.numeric(as.character(bpx$SEQN)),as.numeric(as.character(bpx$BPXDI1)),as.numeric(as.character(bpx$BPXSY1)),
as.numeric(as.character(bpx$BPXDI2)),as.numeric(as.character(bpx$BPXSY2)),as.numeric(as.character(bpx$BPXDI3)),
as.numeric(as.character(bpx$BPXSY3)),as.numeric(as.character(bpx$BPXDI4)),as.numeric(as.character(bpx$BPXSY4))));
bx=data.frame(bx);
names(bx)[1:9]=c('SEQN','BP1L','BP1H','BP2L','BP2H','BP3L','BP3H','BP4L','BP4H');
idx=which(is.na(bx[,2])==TRUE);
#bx[idx,2]=10000;
bx$bpmin_max=bx[,3];
for(k in c(5,7,9)){
	idx=which(is.na(bx[,k])==TRUE);
	bx[idx,k]=10000;
	bx$bpmin_max=pmin(bx$bpmin_max,bx[,k]);
}
idx=which(bx$bpmin_max==10000);
if (length(idx)>0){bx[idx,]$bpmin_max=NA;}
length(which(is.na(bx$bpmin_max)==TRUE))
