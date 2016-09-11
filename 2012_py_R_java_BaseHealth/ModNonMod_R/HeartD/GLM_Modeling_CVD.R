setwd("/home/pouria3/workspace/GenStat/trunk/CSV")
#rm(list=ls())

source('../HeartD/1_Load_Demo_BMI_Diab_Chol_Trigly.R')
source('../HeartD/2_Load_Smoking.R')
source('../HeartD/3_Load_BP.R')
source('../HeartD/4_Load_CVD.R')
setwd("/home/pouria3/workspace/GenStat/trunk/HeartD")

Data$HBP=as.factor(Data$HBP)
Data0=Data[which(Data$AgeC> 18),];
idx<-which(is.na(Data0$BMIC2)==FALSE & is.na(Data0$WaistAdBMI2)==FALSE);
Data1=Data0[idx,];
Mbase=glm(CVD~BMIC2+Age+Gender+WaistAdBMI2+Race,data=Data1,family=binomial);
summary(Mbase)
Data1$BasePred=predict(Mbase);
T2DM=glm(CVD~T2D,data=Data1,family=binomial,offset=BasePred)
summary(T2DM)
HBPM=glm(CVD~HBP,data=Data1,family=binomial,offset=BasePred);
summary(HBPM)
TRIGM=glm(CVD~TRIG2,data=Data1,family=binomial,offset=BasePred);
summary(TRIGM)

SMKM=glm(CVD~SmokingIndPackYrs,data=Data1,family=binomial,offset=BasePred)
summary(SMKM)
AllButTrig=glm(CVD~BMIC2+HBP+T2D+SmokingIndPackYrs+Age+Gender+WaistAdBMI2+Race,data=Data1,family=binomial);
summary(AllButTrig)
idx<-which(is.na(Data1$HBP)==FALSE & is.na(Data1$T2D)==FALSE& is.na(Data1$SmokingIndPackYrs)==FALSE);
x=predict(AllButTrig, newdata=Data1[idx,]);
names(x)=NULL;
Data1$BasePredAllButTrig2=1;
Data1[idx,]$BasePredAllButTrig2=x;
TRIGM2=glm(CVD~TRIG2,data=Data1[idx,],family=binomial,offset=BasePredAllButTrig2);

k=1;
DataP=Data1[1,];
##########################################################
#######Test
##########################################################
DataP[k,]$Age=4;
DataP[k,]$AgeC=40;
DataP[k,]$Gender='M';
DataP[k,]$Race='1NHWhite';
DataP[k,]$BMIC2=35;
DataP[k,]$WaistAdBMI2=10;
DataP[k,]$HBP='1';
DataP[k,]$TRIG2=500;
DataP[k,]$T2D=1;
DataP[k,]$SmokingIndPackYrs=120;
x=predict(Mbase,newdata = DataP[k,]);
names(x)=NULL;
DataP[k,]$BasePred=x;
x=predict(AllButTrig,newdata = DataP[k,]);
names(x)=NULL;
DataP[k,]$BasePredAllButTrig2=x;
#predict(Mbase,newdata = DataP[1,],type = "response")
#predict(HBPM,newdata=DataP,type="response")
#predict(T2DM,newdata=DataP,type="response")
#predict(SMKM,newdata=DataP,type="response")
#predict(TRIGM,newdata=DataP,type="response")
x=(predict(HBPM,newdata=DataP)+predict(T2DM,newdata=DataP)+predict(SMKM,newdata=DataP)+predict(TRIGM,newdata=DataP)-3*DataP$BasePred)
current_model_pred=exp(x)/(exp(x)+1)
current_model_pred
#predict(AllButTrig,newdata=DataP,type="response")
x=(predict(TRIGM2,newdata=DataP))
AllBut_model_pred=exp(x)/(exp(x)+1)
AllBut_model_pred
