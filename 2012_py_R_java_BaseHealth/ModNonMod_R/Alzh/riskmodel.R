library(SiZer)
# To Generate linear models for Risk Factors x
fitIntercept <- function()
{
	RiskFactorValue <- c(1)
	comment         <- rep("Intercept", length(RiskFactorValue) )
	RiskFactorID    <- rep(112, length(RiskFactorValue) )
	FactorValue     <- c(-4.5)
	
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
	
}



fitInheritance <- function()
{
	RiskFactorValue <- c(50:100)
	comment         <- rep("Family", length(RiskFactorValue)+1 )
	RiskFactorID    <- rep(111, length(RiskFactorValue)+1 )
	x               <- c(59, 79)
	y               <- c(log(4), log(2.3))
	FactorValue     <- predict(lm(y~x), data.frame(x=RiskFactorValue) )
	FactorValue[RiskFactorValue <59] <- FactorValue[RiskFactorValue==59]
	FactorValue[RiskFactorValue>80] <- FactorValue[RiskFactorValue==80]
	
	FactorValue[RiskFactorValue==50] <- 0 # if you don't have family history
	# if use selects +100000, map to minimum
	RiskFactorValue <- c(+100000, RiskFactorValue)
	FactorValue     <- c(0, FactorValue) # this is an exception. if you don't have the family history
	
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
		
}


fitDownSynd  <- function()
{
	RiskFactorValue <- c(0, 1)
	comment         <- rep("DownSyndrome", length(RiskFactorValue)+1 )
	RiskFactorID    <- rep(110, length(RiskFactorValue)+1 )
	FactorValue     <- c(0, log(4))
	
	# if use selects +100000, map to minimum
	RiskFactorValue <- c(+100000, RiskFactorValue)
	FactorValue     <- c(min(FactorValue), FactorValue)
	
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
	
}

fitEthnicity <- function()
{
	RiskFactorValue <- c(0:6) #nonHispanic white=0, Asian=1, Hispanic=2, NonHispanic Black=3, Native American=4, Mixed=5, Pacific Islander=6
	comment         <- rep("Ethnicity", length(RiskFactorValue)+1 )
	RiskFactorID    <- rep(109, length(RiskFactorValue)+1 )
	FactorValue     <- c(0, 0, log(2.3), log(2.6), 0, 0, 0)
	
	# if use selects +100000, map to minimum
	RiskFactorValue <- c(+100000, RiskFactorValue)
	FactorValue     <- c(min(FactorValue), FactorValue)
	
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
	
}


fitAge <- function()
{
	RiskFactorValue <- c(50:100)
	comment         <- rep("Age", length(RiskFactorValue)+1 )
	RiskFactorID    <- rep(108, length(RiskFactorValue)+1 )
	x               <- c(77, 82, 87, 97)
	y               <- c(0, log(2.20), log(2.33), log(2.92))
	FactorValue     <- predict(lm(y~x), data.frame(x=RiskFactorValue) )
	FactorValue[RiskFactorValue <77] <- 0 # Age below 77 is normal
	FactorValue[RiskFactorValue>97] <- FactorValue[RiskFactorValue==97]
	
	# if use selects +100000, map to minimum
	RiskFactorValue <- c(+100000, RiskFactorValue)
	FactorValue     <- c(min(FactorValue), FactorValue)
	
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
	
	
}



fitGender <- function()
{
	RiskFactorValue <- c(0, 1)
	comment         <- rep("Gender", length(RiskFactorValue)+1 )
	RiskFactorID    <- rep(107, length(RiskFactorValue)+1 )
	FactorValue     <- c(0, log(1.56))
	
	# if use selects +100000, map to minimum
	RiskFactorValue <- c(+100000, RiskFactorValue)
	FactorValue     <- c(min(FactorValue), FactorValue)
	
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
	
	
	
}



fitSBP <- function()
{
	# Midlife blood pressure
	# was classified according to the following groups: normal
	# (140 mm Hg), borderline (140 to 159 mm Hg), or high
	# (160 mm Hg) for systolic blood pressure
	RiskFactorValue <- c(120:180)
	comment         <- rep("SBP", length(RiskFactorValue)+1 )
	RiskFactorID    <- rep(106, length(RiskFactorValue)+1 )
	x               <- c(130, 150, 170)
	y               <- c(0, log(2.3), log(3.6))
	FactorValue     <- predict(lm(y~x), data.frame(x=RiskFactorValue) )
	
	FactorValue[RiskFactorValue<140] <- 0 # BP below 140 is normal
	FactorValue[RiskFactorValue>160] <- FactorValue[RiskFactorValue==160]
	
	# if use selects NA, map to minimum
	RiskFactorValue <- c(+100000, RiskFactorValue)
	FactorValue     <- c(min(FactorValue), FactorValue)
	
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
	
	
}



fitT2D <- function()
{
	RiskFactorValue <- c(0, 1)
	comment         <- rep("T2D", length(RiskFactorValue)+1 )
	RiskFactorID    <- rep(105, length(RiskFactorValue)+1 )
	FactorValue     <- c(0, log(1.75))
	
	# if use selects +100000, map to minimum
	RiskFactorValue <- c(+100000, RiskFactorValue)
	FactorValue     <- c(min(FactorValue), FactorValue)
	
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
	
}



fitGlucose <- function()
{
	RiskFactorValue <- c(198:500)
	comment         <- rep("Cholesterol", length(RiskFactorValue))
	RiskFactorID    <- rep(104, length(RiskFactorValue) )
		
}


fitChol <- function()
{
	RiskFactorValue <- c(198:500)
	comment         <- rep("TCholesterol", length(RiskFactorValue)+1)
	RiskFactorID    <- rep(104, length(RiskFactorValue)+1 )
	x               <- c(210, 234, 275)
	y               <- c(0, log(1.31), log(1.58))
	FactorValue     <- predict(lm(y~x), data.frame(x=RiskFactorValue) )
	
	FactorValue[FactorValue < 0] <- 0
	FactorValue[RiskFactorValue > 300] <- FactorValue[RiskFactorValue == 300]
	
	
	# if use selects NA, map to minimum
	RiskFactorValue <- c(+100000, RiskFactorValue)
	FactorValue     <- c(min(FactorValue), FactorValue)
	
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
	
}



fitHearDisease <- function()
{
	RiskFactorValue <- c(0, 1)
	comment         <- rep("HeartDisease", length(RiskFactorValue)+1 )
	RiskFactorID    <- rep(103, length(RiskFactorValue)+1 )
	FactorValue     <- c(0, log(1.84))
	
	# if use selects NA, map to minimum
	RiskFactorValue <- c(+100000, RiskFactorValue)
	FactorValue     <- c(min(FactorValue), FactorValue)
	
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
}


fitStroke <- function()
{
	RiskFactorValue <- c(0, 1)
	comment         <- rep("Stroke", length(RiskFactorValue)+1 )
	RiskFactorID    <- rep(102, length(RiskFactorValue)+1 )
	FactorValue     <- c(0, log(4))
	
	# if use selects NA, map to minimum
	RiskFactorValue <- c(+100000, RiskFactorValue)
	FactorValue     <- c(min(FactorValue), FactorValue)
	
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
}



fitTSH <- function()
{
	# refer to Alz Odds ratio document for the details
	# Turns out, this is only applicable to women, men's result wasn't statistically significant
	# TSH level is devided into three levels, [0-1]:OR=2.03, [1.01-2.10]:OR=1 and [2.20-10]:OR=1.86
	
	RiskFactorValue       <- c(0:1000)/100
	comment               <- rep("TSH", length(RiskFactorValue) +1)
	RiskFactorID          <- rep(101, length(RiskFactorValue) +1 )
	FactorValue           <- as.numeric(matrix(0, 1, length(RiskFactorValue)))
	#plot
	# fit a piece wise linear
	x1 <- c(0.5, 1.2)
	y1 <- c(log(2.03), 0)
	m1 <- lm(y1~x1)
	
	x2 <- c(1.7, 6.1)
	y2 <- c(0, log(1.86) )
	m2 <- lm(y2~x2)
	
	for (i in c(1:1001) )
	{
		if (i <= 120)
		{
			FactorValue[i] <- m1$coefficient[[1]] + m1$coefficient[[2]]* i/100
	#		FactorValue[i] <- log(T1)
		}
		if ( 120 <i && i<170)
		{
			FactorValue[i] <- 0
		}
		if ( i >= 170 && i<10000 )
		{
			FactorValue[i] <- m2$coefficient[[1]] + m2$coefficient[[2]]*i/100
	#		FactorValue[i] <- log(T3)
		}
	}
	FactorValue[FactorValue < 0 ] <- 0
	# if use selects NA, map to minimum
	RiskFactorValue <- c(+100000, RiskFactorValue)
	FactorValue     <- c(min(FactorValue), FactorValue)
	
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	return(Data)
}



fitThcy <- function()
{
	# refer to the Alz odds ratio documents for the details
	# Consts:
	T3 = 2.61;
	
	RiskFactorValue        <- c(5:20) # range of Thcy reported in the paper
	FactorValue            <- as.numeric(matrix(0, 1, length(RiskFactorValue)))
	RiskFactorID           <- rep(100, length(RiskFactorValue)+1)
	comment                <- rep("Thcy", length(RiskFactorValue)+1)
	
	#piece wise linear 
	x1 = c(5,15 )
	y1 = c(0, log(T3))
	m1 <- lm(y1~x1)
	
	for (i in RiskFactorValue)
	{
	
		FactorValue[i-4] <- m1$coefficients[[1]] + m1$coefficients[[2]]*i
	
	}
	
	FactorValue[1:5] <- 0 # level [5-10] is normal			
	FactorValue[11:16] <- FactorValue[11]
	
	# if use selects NA, map to minimum
	RiskFactorValue <- c(+100000, RiskFactorValue)
	FactorValue     <- c(min(FactorValue), FactorValue)
	
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, comment=comment)
	
	return(Data)
}

