#################################
#Smoking
#################################
Waveid='d';
SMQ<-read.csv (paste('smq_', Waveid, '_smokingIndex','.csv',sep=""), header=T);
SMQ040=cbind(SMQ$SEQN,as.numeric(as.character(SMQ$SMQ040)));
SMQ050Q=cbind(SMQ$SEQN,as.numeric(as.character(SMQ$SMQ050Q)));
SMQ050U=cbind(SMQ$SEQN,as.numeric(as.character(SMQ$SMQ050U)));
SmokingInd=cbind(SMQ$SEQN,as.numeric(as.character(SMQ$smokingIndex)));
Waveid='c';
SMQ<-read.csv (paste('smq_', Waveid, '_smokingIndex','.csv',sep=""), header=T);
SMQ040 =rbind(SMQ040,cbind(SMQ$SEQN,as.numeric(as.character(SMQ$SMQ040))));
SMQ050Q =rbind(SMQ050Q,cbind(SMQ$SEQN,as.numeric(as.character(SMQ$SMQ050Q))));
SMQ050U =rbind(SMQ050U,cbind(SMQ$SEQN,as.numeric(as.character(SMQ$SMQ050U))));
SmokingInd=rbind(SmokingInd,cbind(SMQ$SEQN,as.numeric(as.character(SMQ$smokingIndex))));
Waveid='b';
SMQ<-read.csv (paste('smq_', Waveid, '_smokingIndex','.csv',sep=""), header=T);
SMQ040 =rbind(SMQ040,cbind(SMQ$SEQN,as.numeric(as.character(SMQ$SMQ040))));
SMQ050Q =rbind(SMQ050Q,cbind(SMQ$SEQN,as.numeric(as.character(SMQ$SMQ050Q))));
SMQ050U =rbind(SMQ050U,cbind(SMQ$SEQN,as.numeric(as.character(SMQ$SMQ050U))));
SmokingInd=rbind(SmokingInd,cbind(SMQ$SEQN,as.numeric(as.character(SMQ$smokingIndex))));
Waveid='a';
SMQ<-read.csv (paste('smq_', Waveid, '_smokingIndex','.csv',sep=""), header=T);
SMQ040 =rbind(SMQ040,cbind(SMQ$SEQN,as.numeric(as.character(SMQ$SMQ040))));
SMQ050Q =rbind(SMQ050Q,cbind(SMQ$SEQN,as.numeric(as.character(SMQ$SMQ050Q))));
SMQ050U =rbind(SMQ050U,cbind(SMQ$SEQN,as.numeric(as.character(SMQ$SMQ050U))));
SmokingInd=rbind(SmokingInd,cbind(SMQ$SEQN,as.numeric(as.character(SMQ$smokingIndex))));


# Convert SmokingIndex into pack-years (20 cigarettes/pack)
SmokingIndPackYrs = cbind(SmokingInd[,1], round(SmokingInd[,2]/20));
# Categorical variable
#SmokingIndPackYrsCat = SmokingIndPackYrs;
#SmokingIndPackYrsCat[which(SmokingIndPackYrs[,2]==0), 2] = '0NoSmoking';
#SmokingIndPackYrsCat[which(SmokingIndPackYrs[,2]>0 & SmokingIndPackYrs[,2]<=   ), 2] = '1';
#SmokingIndPackYrsCat[which(SmokingIndPackYrs[,2]>    & SmokingIndPackYrs[,2]<=  ), 2] = '2';
#SmokingIndPackYrsCat[which(SmokingIndPackYrs[,2]>   ), 2] = '3';



#### Modeling variables for time since subject quit smoking
FormerSmkr = SMQ040;
FormerSmkr[which(SMQ040[,2]==1),2]='0CurrentSmoker';
FormerSmkr[which(SMQ040[,2]==2),2]='0CurrentSmoker';
FormerSmkr[which(SMQ040[,2]==3),2]='1FormerSmoker';
FormerSmkr[which(SMQ040[,2]>3),2]='.';

T_QuitSmking = SMQ050Q;
T_QuitSmking[which(SMQ050U[,2] < 4),2]= 0;
T_QuitSmking[which(SMQ050U[,2] > 4),2]= '.';
T_QuitSmking[which(SMQ050Q[,2] > 99),2]= '.';

T_QuitSmkingCat = cbind(SMQ050Q[,1], round(SMQ050Q[,2]/10));
T_QuitSmkingCat[which(SMQ050U[,2] < 4),2]= 0;
T_QuitSmkingCat[which(SMQ050U[,2] > 4),2]= '.';
T_QuitSmkingCat[which(SMQ050Q[,2] > 99),2]= '.';
####

loc=match(Data$SEQN, SmokingInd[,1]);
Data$SmokingInd = as.numeric(SmokingInd[loc,2]);
loc=match(Data$SEQN, SmokingIndPackYrs[,1]);
Data$SmokingIndPackYrs = as.numeric(SmokingIndPackYrs[loc,2]);
loc=match(Data$SEQN, FormerSmkr[,1]);
Data$FormerSmkr =as.factor(FormerSmkr[loc,2]);
loc=match(Data$SEQN, T_QuitSmking[,1]);
Data$T_QuitSmking =as.numeric(T_QuitSmking[loc,2]);
#loc=match(Data$SEQN, T_QuitSmkingCat[,1]);
#Data$T_QuitSmkingCat =as.factor(T_QuitSmkingCat[loc,2]);