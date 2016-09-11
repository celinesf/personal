################################
# Load CVD data (MCQ)
################################
MCQ <- read.csv ("mcq_d.csv", header=T);
CHF=cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160B)));
CAD=cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160C)));
ANG=cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160D)));
MI=cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160E)));
STRK=cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160F)));
CVD=cbind(MCQ$SEQN,rep(0,length(MCQ$SEQN)));

MCQ <- read.csv ("mcq_c.csv", header=T);
CHF=rbind(CHF,cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160B))));
CAD=rbind(CAD,cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160C))));
ANG=rbind(ANG,cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160D))));
MI=rbind(MI,cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160E))));
STRK=rbind(STRK,cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160F))));
CVD=rbind(CVD, cbind(MCQ$SEQN,rep(0,length(MCQ$SEQN))));

MCQ <- read.csv ("mcq_b.csv", header=T);
CHF=rbind(CHF,cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160B))));
CAD=rbind(CAD,cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160C))));
ANG=rbind(ANG,cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160D))));
MI=rbind(MI,cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160E))));
STRK=rbind(STRK,cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160F))));
CVD=rbind(CVD, cbind(MCQ$SEQN,rep(0,length(MCQ$SEQN))))

MCQ <- read.csv ("mcq_a.csv", header=T);
CHF=rbind(CHF,cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160B))));
CAD=rbind(CAD,cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160C))));
ANG=rbind(ANG,cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160D))));
MI=rbind(MI,cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160E))));
STRK=rbind(STRK,cbind(MCQ$SEQN,as.numeric(as.character(MCQ$MCQ160F))));
CVD=rbind(CVD, cbind(MCQ$SEQN,rep(0,length(MCQ$SEQN))))

################################
# Processing
################################
CHF[which(CHF[,2]>2),2]='.';
CHF[which(CHF[,2]==2),2]=0;
CAD[which(CAD[,2]>2),2]='.';
CAD[which(CAD[,2]==2),2]=0;
ANG[which(ANG[,2]>2),2]='.';
ANG[which(ANG[,2]==2),2]=0;
MI[which(MI[,2]>2),2]='.';
MI[which(MI[,2]==2),2]=0;
STRK[which(STRK[,2]>2),2]='.';
STRK[which(STRK[,2]==2),2]=0;

# CVD being defined as MI + Angina + CAD
CVD[which(CAD[,2]==1), 2] = 1
CVD[which(ANG[,2]==1), 2] = 1
CVD[which(MI[,2]==1), 2] = 1

# Not including CHF & stroke as part of CVD
#CVD[which(CHF[,2]==1), 2] = 1
#CVD[which(STRK[,2]==1), 2] = 1

################################
# Matching
################################
loc=match(Data$SEQN,CVD[,1]);
Data$CVD=as.numeric(as.character(CVD[loc,2]));

#loc=match(Data$SEQN,CHF[,1]);
#Data$CHF=as.numeric(as.character(CHF[loc,2]));
#loc=match(Data$SEQN,CAD[,1]);
#Data$CAD=as.numeric(as.character(CAD[loc,2]));
#loc=match(Data$SEQN,ANG[,1]);
#Data$ANG=as.numeric(as.character(ANG[loc,2]));
#loc=match(Data$SEQN,MI[,1]);
#Data$MI=as.numeric(as.character(MI[loc,2]));
#loc=match(Data$SEQN,STRK[,1]);
#Data$STRK=as.numeric(as.character(STRK[loc,2]));



