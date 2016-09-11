###### Number of drinks over last yr & drink_score (0, 1, 2, or 3) ######
Waveid='d';
ALQ<-read.csv (paste('alq_', Waveid, '_alcIndex','.csv',sep=""), header=T);
nDrinks=cbind(ALQ$SEQN,as.numeric(as.character(ALQ$n_drinks)));
DRNK=cbind(ALQ$SEQN,as.numeric(as.character(ALQ$drink_score)));
Waveid='c';
ALQ<-read.csv (paste('alq_', Waveid, '_alcIndex','.csv',sep=""), header=T);
nDrinks=rbind(nDrinks,cbind(ALQ$SEQN,as.numeric(as.character(ALQ$n_drinks))));
DRNK =rbind(DRNK,cbind(ALQ$SEQN,as.numeric(as.character(ALQ$drink_score))));
Waveid='b';
ALQ<-read.csv (paste('alq_', Waveid, '_alcIndex','.csv',sep=""), header=T);
nDrinks=rbind(nDrinks,cbind(ALQ$SEQN,as.numeric(as.character(ALQ$n_drinks))));
DRNK =rbind(DRNK,cbind(ALQ$SEQN,as.numeric(as.character(ALQ$drink_score))));
Waveid='a';
ALQ<-read.csv (paste('alq_', Waveid, '_alcIndex','.csv',sep=""), header=T);
nDrinks=rbind(nDrinks,cbind(ALQ$SEQN,as.numeric(as.character(ALQ$n_drinks))));
DRNK =rbind(DRNK,cbind(ALQ$SEQN,as.numeric(as.character(ALQ$drink_score))));

nDrinksPerDay = cbind(nDrinks[,1], round(nDrinks[,2]/365))
DrinksPerDayCat = nDrinksPerDay
DrinksPerDayCat[which(nDrinksPerDay[,2]==0), 2] = '0NoDrinks';
DrinksPerDayCat[which(nDrinksPerDay[,2]>0 & nDrinksPerDay[,2]<=5), 2] = '1UpTo5Drinks';
DrinksPerDayCat[which(nDrinksPerDay[,2]>5 & nDrinksPerDay[,2]<=10), 2] = '2Betw6and10Drinks';
DrinksPerDayCat[which(nDrinksPerDay[,2]>10), 2] = '3MoreThan10Drinks';


loc=match(Data$SEQN,nDrinksPerDay[,1]);
Data$numDrinksPerDay=as.numeric(nDrinksPerDay[loc,2]);
loc=match(Data$SEQN,DRNK[,1]);
Data$DrinkScore=as.numeric(DRNK[loc,2]);
loc=match(Data$SEQN, DrinksPerDayCat[,1]);
Data$DrinksPerDayCat =as.factor(DrinksPerDayCat[loc,2]);

###### Number of binge drinking incidents over last year & Binge Drinking Score (1= above median, 0=below median) ######
Waveid='d';
ALQ<-read.csv (paste('alq_', Waveid, '_alcIndex','.csv',sep=""), header=T);
nBinges=cbind(ALQ$SEQN,as.numeric(as.character(ALQ$n_binges)));
BingeScore=cbind(ALQ$SEQN,as.numeric(as.character(ALQ$binge_score)));
Waveid='c';
ALQ<-read.csv (paste('alq_', Waveid, '_alcIndex','.csv',sep=""), header=T);
nBinges =rbind(nBinges,cbind(ALQ$SEQN,as.numeric(as.character(ALQ$n_binges))));
BingeScore =rbind(BingeScore,cbind(ALQ$SEQN,as.numeric(as.character(ALQ$binge_score))));
Waveid='b';
ALQ<-read.csv (paste('alq_', Waveid, '_alcIndex','.csv',sep=""), header=T);
nBinges =rbind(nBinges,cbind(ALQ$SEQN,as.numeric(as.character(ALQ$n_binges))));
BingeScore =rbind(BingeScore,cbind(ALQ$SEQN,as.numeric(as.character(ALQ$binge_score))));
Waveid='a';
ALQ<-read.csv (paste('alq_', Waveid, '_alcIndex','.csv',sep=""), header=T);
nBinges =rbind(nBinges,cbind(ALQ$SEQN,as.numeric(as.character(ALQ$n_binges))));
BingeScore =rbind(BingeScore,cbind(ALQ$SEQN,as.numeric(as.character(ALQ$binge_score))));

#BingesCat = nBinges;
#BingesCat[which(nBinges[,2]==0), 2] = '0NoBinges';
#BingesCat[which(nBinges[,2]>0 & nBinges[,2]<=24), 2] = '1LessThan24Binges';
#BingesCat[which(nBinges[,2]>24), 2] = '2MoreThan24Binges';
#loc=match(Data$SEQN, BingesCat[,1]);
#Data$BingesCat =as.factor(BingesCat[loc,2]);

loc=match(Data$SEQN, nBinges[,1]);
Data$numBinges =as.numeric(nBinges[loc,2]);
loc=match(Data$SEQN, BingeScore[,1]);
Data$BingeScore =as.numeric(BingeScore[loc,2]);

###### History of 'Alcoholism' (5+ drinks per day every day) ######
Waveid='d';
ALQ<-read.csv (paste('alq_', Waveid, '_alcIndex','.csv',sep=""), header=T);
HxAlcoholism=cbind(ALQ$SEQN,as.numeric(as.character(ALQ$hx_alc)));
Waveid='c';
ALQ<-read.csv (paste('alq_', Waveid, '_alcIndex','.csv',sep=""), header=T);
HxAlcoholism =rbind(HxAlcoholism,cbind(ALQ$SEQN,as.numeric(as.character(ALQ$hx_alc))));
Waveid='b';
ALQ<-read.csv (paste('alq_', Waveid, '_alcIndex','.csv',sep=""), header=T);
HxAlcoholism =rbind(HxAlcoholism,cbind(ALQ$SEQN,as.numeric(as.character(ALQ$hx_alc))));
Waveid='a';
ALQ<-read.csv (paste('alq_', Waveid, '_alcIndex','.csv',sep=""), header=T);
HxAlcoholism =rbind(HxAlcoholism,cbind(ALQ$SEQN,as.numeric(as.character(ALQ$hx_alc))));

loc=match(Data$SEQN, HxAlcoholism[,1]);
Data$HxAlcoholism =as.numeric(HxAlcoholism[loc,2]);
