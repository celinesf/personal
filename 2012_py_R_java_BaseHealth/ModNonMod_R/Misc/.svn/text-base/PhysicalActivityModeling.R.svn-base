##############################
# Load Physical Activity Data
##############################
Waveid='d';
PAQ <-read.csv (paste('paq_', Waveid, '.csv',sep=""), header=T);
PhysAct=cbind(PAQ$SEQN,as.numeric(as.character(PAQ$PAQ180)));
VigorAct=cbind(PAQ$SEQN,as.numeric(as.character(PAQ$PAD200)));
ModAct=cbind(PAQ$SEQN,as.numeric(as.character(PAQ$PAD320)));
Waveid='c';
PAQ <-read.csv (paste('paq_', Waveid, '.csv',sep=""), header=T);
PhysAct =rbind(PhysAct,cbind(PAQ$SEQN,as.numeric(as.character(PAQ$PAQ180))));
VigorAct =rbind(VigorAct,cbind(PAQ$SEQN,as.numeric(as.character(PAQ$PAD200))));
ModAct =rbind(ModAct,cbind(PAQ$SEQN,as.numeric(as.character(PAQ$PAD320))));
Waveid='b';
PAQ <-read.csv (paste('paq_', Waveid, '.csv',sep=""), header=T);
PhysAct =rbind(PhysAct,cbind(PAQ$SEQN,as.numeric(as.character(PAQ$PAQ180))));
VigorAct =rbind(VigorAct,cbind(PAQ$SEQN,as.numeric(as.character(PAQ$PAD200))));
ModAct =rbind(ModAct,cbind(PAQ$SEQN,as.numeric(as.character(PAQ$PAD320))));
Waveid='a';
PAQ <-read.csv (paste('paq_', Waveid, '.csv',sep=""), header=T);
PhysAct =rbind(PhysAct,cbind(PAQ$SEQN,as.numeric(as.character(PAQ$PAQ180))));
VigorAct =rbind(VigorAct,cbind(PAQ$SEQN,as.numeric(as.character(PAQ$PAD200))));
ModAct =rbind(ModAct,cbind(PAQ$SEQN,as.numeric(as.character(PAQ$PAD320))));

VigorAct[which(VigorAct[,2]>2),2]='.';
VigorAct[which(VigorAct[,2]==2),2]=0;
ModAct[which(ModAct[,2]>2),2]='.';
ModAct[which(ModAct[,2]==2),2]=0;

PhysActC = PhysAct;
PhysAct[which(PhysActC[,2]==1),2] = '0LowestActivity';
PhysAct[which(PhysActC[,2]==2),2] = '1LowActivity';
PhysAct[which(PhysActC[,2]==3),2] = '2HighActivity';
PhysAct[which(PhysActC[,2]==4),2] = '3HighestActivity';
PhysAct[which(PhysActC[,2]>4),2] = NA;


loc=match(Data$SEQN, PhysActC[,1]);
Data$PhysActC =as.numeric(PhysActC[loc,2]);
loc=match(Data$SEQN, PhysAct[,1]);
Data$PhysAct =as.numeric(PhysAct[loc,2]);
loc=match(Data$SEQN, VigorAct[,1]);
Data$VigorAct =as.numeric(VigorAct[loc,2]);
loc=match(Data$SEQN, ModAct[,1]);
Data$ModAct =as.numeric(ModAct[loc,2]);

