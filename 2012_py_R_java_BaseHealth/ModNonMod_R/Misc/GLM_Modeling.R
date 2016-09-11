rm(list=ls())
#source('R_code/Load_Demo_BMI_Chol_Diab.R')
#source('R_code/BP.R')
source('R_code/Load_Demo_BMI_Chol_Diab_BP.R')
source('R_code/BPModeling.R')
source('R_code/Proc_Match_Corr2.R')
source('R_code/PhysicalActivityModeling.R')
source('R_code/SmokingModeling.R')
source('R_code/AlcoholModeling.R')
#Model_Alc <- glm(T2D ~ Age + BMI+Gender+Race+WaistAdBMI + DrinksPerDayCat + BingeScore + HxAlcoholism, data=Data,family=binomial);
#Model_PhysActivity <- glm(T2D ~ Age + BMI+Gender+Race+WaistAdBMI + PhysActC + VigorAct, data=Data,family=binomial);
#Model_Alc_PhysActivity <- glm(T2D ~ Age + BMI+Gender+Race+WaistAdBMI + DrinksPerDayCat + BingeScore + HxAlcoholism + PhysActC + VigorAct, data=Data,family=binomial);
#Model_Smk <- glm(T2D ~ Age + BMI+Gender+Race+WaistAdBMI + SmokingIndPackYrs, data=Data,family=binomial);
#Model_Alc_PhysActivity_Smk <- (T2D ~ Age + BMI+Gender+Race+WaistAdBMI + DrinksPerDayCat + BingeScore + HxAlcoholism + PhysActC + VigorAct + );
Test <- glm(T2D ~ Age + BMIC2+Gender+Race+WaistAdBMI  
+ DrinksPerDayCat + BingeScore + HxAlcoholism + PhysActC + VigorAct + 
SmokingIndPackYrs+CHR2+HBP, data=Data,family	=binomial);
Test <- glm(T2D ~ Age + BMIC2+Gender+Race+WaistAdBMI  
+CHR2, data=Data,family	=binomial);


BaseModel <- glm(T2D ~ Age + Gender+Race,data=Data,family	=binomial);