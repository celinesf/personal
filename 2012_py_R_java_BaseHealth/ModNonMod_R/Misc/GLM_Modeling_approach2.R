rm(list=ls())
#source('R_code/Load_Demo_BMI_Chol_Diab.R')
#source('R_code/BP.R')
source('R_code/Load_Demo_BMI_Chol_Diab_BP.R')
source('R_code/BPModeling.R')
source('R_code/Proc_Match_Corr2.R')
source('R_code/PhysicalActivityModeling.R')
source('R_code/SmokingModeling.R')
source('R_code/AlcoholModeling.R')
idx=which((is.na(Data$T2D)==FALSE)&(Data$AgeC>=15));
Data=Data[idx,];
Data2=data.frame(list(Data$Age,Data$BMI,Data$Waist,Data$T2D,Data$CHR2,Data$TRIG,Data$HBP,Data$SmokingIndPackYrs,Data$BingeScore))
names(Data2)[1:9]=list("Age","BMI","Waist","T2D","CHR2","TRIG","HBP","Smoke","Alcohol");
#idx=which(((is.na(Data$CHR2)==FALSE&is.na(Data$BMIC2)==FALSE)&(is.na(Data$T2D)==FALSE
#&is.na(Data$Race)==FALSE))&((is.na(Data$Age)==FALSE&is.na(Data$Gender)==FALSE)&is.na(Data$WaistAdBMI)==FALSE));
#Data=Data[idx,];
#idx=which(Data$AgeC>15);
#Data=Data[idx,];
summary(Model1 <- glm(T2D ~ Age + BMIC2+Gender+Race+WaistAdBMI, data=Data,family=binomial));
#41000-7627
Data$Pred=predict(Model1,Data)
summary(Model_CHR2<-glm(T2D~CHR2, offset=Data$Pred,data=Data,family=binomial))
summary(Model_TRIG<-glm(T2D~TRIG2, offset=Data$Pred,data=Data,family=binomial))
summary(Model_SmokingIndPackYrs<-glm(T2D~SmokingIndPackYrs, offset=Data$Pred,data=Data,family=binomial))
summary(Model_Alc<-glm(T2D~DrinksPerDayCat + BingeScore + HxAlcoholism, offset=Data$Pred,data=Data,family=binomial))
summary(Model_HBP<-glm(T2D~HBP, offset=Data$Pred,data=Data,family=binomial))
summary(Model_Physc<-glm(T2D~PhysActC + VigorAct , offset=Data$Pred,data=Data,family=binomial))










summary(Test <- glm(T2D ~ Age + BMIC2+Gender+Race+WaistAdBMI  + PhysActC + VigorAct + 
SmokingIndPackYrs+CHR2+HBP+DrinksPerDayCat + BingeScore + HxAlcoholism, data=Data,family=binomial));
Data$Pred2=predict(Test,Data)
Data$Pred2=pmin(0,abs(Data$Pred2))
summary(Model1_2 <- glm(T2D ~ Age + BMIC2+Gender+Race+WaistAdBMI   ,offset=Data$Pred2, data=Data,family=binomial));


Test <- glm(T2D ~ Age + BMIC2+Gender+Race+WaistAdBMI  
+ DrinksPerDayCat + BingeScore + HxAlcoholism + PhysActC + VigorAct + 
SmokingIndPackYrs+CHR2+HBP, data=Data,family	=binomial);










Model1P <- glm(T2D ~ Age + BMIC2+Gender+Race+WaistAdBMI, offset=Data$Pred, data=Data,family=binomial);
Model2 <- glm(T2D ~ Age + BMIC2+Gender+Race+WaistAdBMI+CHR2, data=Data,family=binomial);
Test <- glm(T2D ~ Age + BMIC2+Gender+Race+WaistAdBMI  
+ DrinksPerDayCat + BingeScore + HxAlcoholism + PhysActC + VigorAct + 
SmokingIndPackYrs+CHR2+HBP+TRIG2, data=Data,family	=binomial);
p1=predict(Model1,Data);
p2=predict(Model2,Data);
p3=predict(Model_CHR2,Data);
p4=predict(Model1P,Data);


Model_TRIG=glm(T2D~TRIG, offset=-Data$Pred,data=Data);
Model_SmokingIndPackYrs=glm(T2D~SmokingIndPackYrs, offset=-Data$Pred,data=Data);
Model_PhysActC_VigorAct =glm(T2D~PhysActC+VigorAct, offset=-Data$Pred,data=Data);



summary(Model1)


















#Model_Alc <- glm(T2D ~ Age + BMI+Gender+Race+WaistAdBMI + DrinksPerDayCat + BingeScore + HxAlcoholism, data=Data,family=binomial);
#Model_PhysActivity <- glm(T2D ~ Age + BMI+Gender+Race+WaistAdBMI + PhysActC + VigorAct, data=Data,family=binomial);
#Model_Alc_PhysActivity <- glm(T2D ~ Age + BMI+Gender+Race+WaistAdBMI + DrinksPerDayCat + BingeScore + HxAlcoholism + PhysActC + VigorAct, data=Data,family=binomial);
#Model_Smk <- glm(T2D ~ Age + BMI+Gender+Race+WaistAdBMI + SmokingIndPackYrs, data=Data,family=binomial);
#Model_Alc_PhysActivity_Smk <- glm(T2D ~ Age + BMI+Gender+Race+WaistAdBMI + DrinksPerDayCat + BingeScore + HxAlcoholism + PhysActC + VigorAct + );
Test <- glm(T2D ~ Age + BMIC2+Gender+Race+WaistAdBMI  
+ DrinksPerDayCat + BingeScore + HxAlcoholism + PhysActC + VigorAct + 
SmokingIndPackYrs+CHR2+HBP, data=Data,family	=binomial);

Test2 <- glm(T2D ~ Age + BMIC2+Gender+Race+WaistAdBMI+CHR2+HBP, data=Data,family	=binomial);


BaseModel <- glm(T2D ~ Age + Gender+Race,data=Data,family	=binomial);

Test3 <- glm(TRIG ~ Age + BMIC2+Gender+Race+WaistAdBMI  
+ DrinksPerDayCat + BingeScore + HxAlcoholism + PhysActC + VigorAct + 
SmokingIndPackYrs+CHR2+HBP, data=Data,family	=Gamma);
BaseModel <- glm(T2D ~ Age + Gender+Race,data=Data,family	=binomial);