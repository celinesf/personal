rm(list=ls())
source("GLM_Modeling_CVD.R")

fitSBP <- function()
{
	
	RiskFactorValue <- c(110:190)
	comment         <- rep("SBP", length(RiskFactorValue))
	RiskFactorID    <- rep(100, length(RiskFactorValue))
	x               <- c(115, 125, 135, 150, 170, 190)
	y               <- c(log(1), log(1.43), log(1.85), log(2.82), log(3.91), log(6.58) )
	FactorValue     <- predict(lm(y~x), data.frame(x=RiskFactorValue) )
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
	
	
}
 
fitDBP <- function()
{
	RiskFactorValue <- c(70:120)
	comment         <- rep("DBP", length(RiskFactorValue))
	RiskFactorID    <- rep(101, length(RiskFactorValue))
	x               <- c(78, 82, 87, 95, 105, 115)
	y               <- c(log(1), log(1.43), log(1.85), log(2.82), log(3.91), log(6.58) )
	FactorValue     <- predict(lm(y~x), data.frame(x=RiskFactorValue) )
	FactorValue[RiskFactorValue<76] <- FactorValue[RiskFactorValue==76] # 
	FactorValue[RiskFactorValue>115] <- FactorValue[RiskFactorValue==115] # 
		
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)

}
		 
fitLDL <- function()
{
	RiskFactorValue <- c(100:250)
	comment         <- rep("LDL", length(RiskFactorValue))
	RiskFactorID    <- rep(109, length(RiskFactorValue))
	FactorValue     <- c(-30:120) * log(1.07)/10
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
}
	

fitEthnicity <- function()
{
	#Ethnicity
	RiskFactorValue <- c(0, 1, 2, 3, 4, 5, 6) #c(Asian, "Hispanic", "NonHispWhite", "Native American", "NonHisp Black", "Othermixed", "PacificIslander")
	comment         <- rep("Ethnicity", length(RiskFactorValue))
	RiskFactorID    <- rep(111, length(RiskFactorValue))
	FactorValue     <- c(-0.65392647, 0.412109651, 0, 0, 0.693147181, 0, 0)
	Data     <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
	
}

fitDiet <- function()
{
	FactorValue <- c(0.10237290870956, 0.10237290870956, 0.1010482489419, 0.09972358917424, 0.09839892940658, 0.09707426963892, 0.09574960987126, 
	0.0944249501036, 0.09310029033594, 0.09177563056827, 0.09045097080061, 0.08912631103295, 0.08780165126529, 0.08647699149763, 0.08515233172997,
	0.08382767196231, 0.08250301219465, 0.08117835242699, 0.07985369265933, 0.07852903289167, 0.07720437312401, 0.07587971335635, 0.07455505358869,
	0.07323039382103, 0.07190573405337, 0.07058107428571, 0.06049806367346, 0.05041505306122, 0.04033204244898, 0.03024903183673, 0.02016602122449,
	0.01008301061224, 0, -0.01431045771792, -0.02862091543585, -0.04293137315377, -0.05724183087169, -0.07155228858962, -0.07664540955087, -0.08173853051213,
	-0.08683165147339, -0.09192477243464, -0.0970178933959, -0.10211101435716, -0.10720413531842, -0.11229725627967)
	RiskFactorValue <- c(0:45)
	comment         <- rep("DIET(AHEI)", length(RiskFactorValue))
	RiskFactorID    <- rep(114, length(RiskFactorValue))
	Data     <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
}

fitPA  <- function()
{
	FactorValue <- c(0.046, 0.041, 0.035, 0.030, 0.025, 0.019, 0.014, 0.012, 0.011, 0.009, 0.007, 0.005, 0.004, 0.002, 
	0.000, -0.001, -0.003, -0.004, -0.005, -0.007, -0.008, -0.009, -0.011, -0.012, -0.013, -0.015, -0.019, -0.023, -0.027, 
	-0.031, -0.035, -0.040, -0.044, -0.048, -0.052, -0.056, -0.060, -0.064, -0.068, -0.073, -0.077, -0.081, -0.085)
	RiskFactorValue <- c(0:42)
	comment         <- rep("PA", length(RiskFactorValue))
	RiskFactorID    <- rep(115, length(RiskFactorValue))
	Data     <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)


}
fitIntercept <- function()
{
	RiskFactorValue <- c(1)
	comment         <- rep("Intercept", length(RiskFactorValue) )
	RiskFactorID    <- rep(116, length(RiskFactorValue) )
	FactorValue     <- c(-6.349)
	
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
	
}

fitAlcohol <- function()
{
	RiskFactorValue <- c(1, 2, 3, 4)
	comment         <- rep("Alcohol", length(RiskFactorValue) )
	RiskFactorID    <- rep(117, length(RiskFactorValue) )
	#FactorValue     <- c(0, -0.28369005, -0.28369005, 0.032088315)
	FactorValue     <- c(0, 0, 0, (0.28369005 + 0.032088315) )
	
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
	
}

fitT2D <- function()
{
	RiskFactorValue <- c(0, 1)
	comment         <- rep("T2D", length(RiskFactorValue))
	RiskFactorID    <- rep(18, length(RiskFactorValue))
	FactorValue     <- c(0, log(2.78))
	
	Data            <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
}

fitSleep <- function()
{
	RiskFactorValue <- c(0, 1, 2) # sleep Apnea, Sleep disorder
	comment         <- rep("Sleep", length(RiskFactorValue))
	RiskFactorID    <- rep(19, length(RiskFactorValue))
	FactorValue     <- c(0, log(2.89), log(2.89))
	
	Data            <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
}
	
fitPastCVD <- function()
{
	RiskFactorValue <- c(0, 1, 2, 3) #atrail fibrilation, myocardial infaction, heartfailure or strok
	comment         <- rep("PastCVD", length(RiskFactorValue))
	RiskFactorID    <- rep(20, length(RiskFactorValue))
	FactorValue     <- c(0, log(2.86), log(2.09), log(2.2) )
	
	Data            <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
}

	#modelExtraction <- function()
	
	# Mbase is the base model
	#BMI
	bmi_min <- 10*min(Data1$BMIC2); bmi_max <- 10*max(Data1$BMIC2)
	bmi_cof <- Mbase$coefficient[["BMIC2"]]
	RiskFactorValue <- c(bmi_min:bmi_max)/10
	comment         <- rep("BMI", length(RiskFactorValue))
	RiskFactorID    <- rep(103, length(RiskFactorValue))
	FactorValue     <- c(0, c(1:(length(RiskFactorValue)-1))*bmi_cof/10)
	BMI_Data        <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	BMI             <- data.frame(RiskFactorID=103, name="BMI"         , Min=bmi_min/10  , Max=bmi_max/10  , normal=FactorValue[RiskFactorValue==25])
	
	#AGE
	age_min <- min(Data1$AgeC); age_max <- max(Data1$AgeC)
	age_cof <- c(Mbase$coefficient[["Age3"]], 
		     Mbase$coefficient[["Age4"]] - Mbase$coefficient[["Age3"]], 
		     Mbase$coefficient[["Age5"]] - Mbase$coefficient[["Age4"]], 
		     Mbase$coefficient[["Age6"]] - Mbase$coefficient[["Age5"]], 
		     Mbase$coefficient[["Age7"]] - Mbase$coefficient[["Age6"]], 
		     Mbase$coefficient[["Age8"]] - Mbase$coefficient[["Age7"]])
	RiskFactorValue <- c(age_min:age_max)
	comment         <- rep("Age", length(RiskFactorValue))
	RiskFactorID    <- rep(104, length(RiskFactorValue))
	FactorValue     <- rep(0, length(RiskFactorValue))
	for (i in c((age_min+1): age_max) )
	{
		FactorValue[i-(age_min-1)] <- FactorValue[i-age_min] + age_cof[floor(i/10)-1]/10
	}
	FactorValue[RiskFactorValue>79] <- FactorValue[RiskFactorValue==79]
	Age_Data        <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	Age             <- data.frame(RiskFactorID=104, name="Age"         , Min=age_min  , Max=age_max  , normal=NA)
	
	
	#Gender
	RiskFactorValue <- c(0, 1) #c("M", "F")
	comment         <- rep("Gender", length(RiskFactorValue))
	RiskFactorID    <- rep(105, length(RiskFactorValue))
	FactorValue     <- c(Mbase$coefficient[["GenderM"]], 0)
	Gender_Data     <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	Gender          <- data.frame(RiskFactorID=105, name="Gender"         , Min=0  , Max=1  , normal=NA)
	
	
	#Waist
	waist_min = min(Data1$WaistAdBMI2); waist_max=max(Data1$WaistAdBMI2)
	RiskFactorValue <- c(waist_min : waist_max)
	comment         <- rep("Waist", length(RiskFactorValue))
	RiskFactorID    <- rep(106, length(RiskFactorValue))
	FactorValue     <- RiskFactorValue * Mbase$coefficient[["WaistAdBMI2"]]
	FactorValue[RiskFactorValue < -2] <- FactorValue[RiskFactorValue == -2]
	FactorValue[RiskFactorValue > 10] <- FactorValue[RiskFactorValue == 10]
	Waist_Data      <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	Waist           <- data.frame(RiskFactorID=106, name="Waist"         , Min=waist_min  , Max=waist_max  , normal=0)
	
	
	#T2D (was decided to use the paper infor instead of CDC dataset
	#RiskFactorValue <- c(0, 1)
	#comment         <- rep("T2D", length(RiskFactorValue))
	#RiskFactorID    <- rep(107, length(RiskFactorValue))
	#FactorValue     <- c(T2DM$coefficient[[1]], T2DM$coefficient[[1]] + 1*T2DM$coefficient[[2]] )
	#T2D_Data        <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	#T2D             <- data.frame(RiskFactorID=107, name="T2D"         , Min=0  , Max=1  , normal=NA)
	
	
	#TRIG
	trig_min <- min(Data1$TRIG2[!is.na(Data1$TRIG2)]) ; trig_max <- max(Data1$TRIG2[!is.na(Data1$TRIG2)])
	RiskFactorValue <- c(trig_min : trig_max )
	comment         <- rep("TRIG", length(RiskFactorValue))
	RiskFactorID    <- rep(108, length(RiskFactorValue))
	FactorValue     <- TRIGM$coefficient[[1]] + (RiskFactorValue-trig_min) * TRIGM$coefficient[[2]]
	Trig_Data       <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	Trig            <- data.frame(RiskFactorID=108, name="Trig"         , Min=trig_min  , Max=trig_max  , normal=FactorValue[RiskFactorValue==150])
	
	#Smoking
	smk_min <- min(Data1$SmokingIndPackYrs[!is.na(Data1$SmokingIndPackYrs)]) ; smk_max <- 120 # max pack per year, refer to CVD doc
	RiskFactorValue <- c(smk_min : smk_max )
	comment         <- rep("Smoking", length(RiskFactorValue))
	RiskFactorID    <- rep(112, length(RiskFactorValue))
	FactorValue     <- RiskFactorValue * SMKM$coefficient[[2]]*365 #because the model is smokingyear(not per year)
	Smk_Data       <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	Smk            <- data.frame(RiskFactorID=112, name="Smoking"         , Min=smk_min  , Max=smk_max  , normal=NA)
	
	#Future Smoking
	comment         <- rep("FutSmoking", length(RiskFactorValue))
	RiskFactorID    <- rep(113, length(RiskFactorValue))
	FactorValue     <- RiskFactorValue * SMKM$coefficient[[2]]*365 #because the model is smokingyear(not per year)
	FutSmk_Data       <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	FutSmk            <- data.frame(RiskFactorID=113, name="FutSmoking"         , Min=smk_min  , Max=smk_max  , normal=NA)
	
	
	Data <- rbind(fitSBP(), fitDBP(), fitPastCVD(), fitSleep(), fitT2D(), fitLDL(), fitDiet(),fitPA(), fitIntercept(),fitEthnicity(), fitAlcohol(),
		      Age_Data, Gender_Data, Waist_Data, Trig_Data, BMI_Data, Smk_Data, FutSmk_Data)
	
	
	
	
	
	
	
	
	
	
	
	
	




