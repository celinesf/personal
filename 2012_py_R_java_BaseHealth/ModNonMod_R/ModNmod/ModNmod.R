rm(list=ls())
library(RMySQL) # for Database stuff
library(SiZer) # for piecewise linear regression





getmap <- function()
{
	# DB stuff; query IDmap table to map ids to names for GUI
	con    <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
	qur    <- ('select * from IDmap where (upper(IsRiskFactorID)=\'YES\') ')
	res    <- dbGetQuery(con, qur)
	dbDisconnect(con)
	
	return(res)
}
	
	
chkConf <- function(ors, cis, ranges)
{
	# Funtion to check confidence intervals, if it corsses 1, set OR back to 1
	#cis refers to Confidence intervals, ors refere to odds ratios
	newor = c(); k <- 0; newrange <- c()
	cat("CIs : ", cis, "\n")
	cat("Original ORs: ", ors, "\n")
	if ( (length(ors)-1) != length(cis) )
		stop('CIs and ORs dont have the same length')
	
	for (j in c(1:length(ors)) )
	{
		if (as.numeric(ors[j]) == 1 && k==0)
		{
			newor <- c(newor, ors[j])
			newrange <- c(newrange, ranges[j])
			k <-1 # set a flag since there is no ci for reference or
		}else
		{
			ci  <- as.numeric(unlist(strsplit(cis[j-k], '-')) )
			or  <- as.numeric(ors[j])
			if ( (ci[1] <= 1) & (ci[2] >= 1) & !is.na(ci))
			{
				#newor <- c(newor, 1)
				cat("Non significant CI was observed :", ci[1], "-", ci[2], " Excluding this point\n")
			}else
			{
				newor <- c(newor, or)
				newrange <- c(newrange, ranges[j])
			}
		}
	}
	cat("ORs after inspection", newor, "\n", "New ranges : ", newrange, "\n") 
	return(data.frame(newor=newor, newrange=newrange))
		
}


categoricals <- function(map, did)
{
	# model Categorical riskFactors
	# DB stuff
	con    <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
	#qur1    <- ('select * from ModNmod where (upper(Genophen)=\'YES\' and DiseaseID=111 and upper(isFactor)=\'YES\') order by PMID DESC')
	qur1    <- paste('select * from ModNmod where (DiseaseID=', did, ' and upper(isFactor)=\'YES\' and Genophen=\'yes\') order by PMID DESC' )
	res1    <- dbGetQuery(con, qur1)
	dbDisconnect(con)
	
	# initialize Data object
	Data <- data.frame(RiskFactorID=c(), RiskFactorValue=c(), FactorValue=c(), Min=c(), Max=c(), normal=c(), name=c(), comment=c(), Gender=c())
	
	par(lty=2, lwd=2, mfrow=c(ceiling(nrow(res1)/3), 3))
	
	for (i in c(1:nrow(res1)) )
	{
		cats  <- unlist(strsplit(res1[i,]$Categories, ',')) 
		ors   <- log(as.numeric(unlist(strsplit(res1[i,]$ORs, ',')) ))
		if (length(cats) != length(ors) )
			stop('categories and ORs have different lengths')
		
		Data <- rbind(Data, data.frame(RiskFactorID=res1[i,]$RiskID, RiskFactorValue=c(0:(length(cats)-1)), FactorValue=ors, 
		Min=0, Max=length(cats)-1, normal=res1[i, ]$NormalVal, name=map$IDName[map$IDNo==res1[i,]$RiskID], comment=cats , Gender=res1[i,]$Gender) )
		
		#plot
		barplot(exp(ors), names.arg=cats, col='blue', main=map$IDName[map$IDNo==res1[i,]$RiskID])
	}
	Data$RiskFactorValue <- as.character(Data$RiskFactorValue)
	return(Data)
	
}




continuos <- function(map, did)
{
	# model continuous riskFactors
	con    <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
	qur2     <- paste('select * from ModNmod where (upper(Genophen)=\'YES\' and DiseaseID=', did, ' and upper(isFactor)=\'NO\') order by PMID DESC')
	res2     <- dbGetQuery(con, qur2)
	dbDisconnect(con)
	
	# initialize Data object
	Data <- data.frame(RiskFactorID=c(), RiskFactorValue=c(), FactorValue=c(), Min=c(), Max=c(), normal=c(), name=c(), comment=c(), Gender=c())
	
	dev.new()
	par(lty=2, lwd=2, mfrow=c(ceiling(nrow(res2)/3), 3))
	
	for (i in c(1:nrow(res2) ) )
	{
		x = c(); y=c(); xx<- c()
		ranges <- unlist(strsplit(res2[i,]$Ranges, ','))
		for (j in c(1:length(ranges)) ) # to keep original points
		{
			r    <- as.numeric(unlist(strsplit(ranges[j], '-')) )
			xx  <- c(xx, 0.5*(r[1] + r[2]) )
		}
		
		ors <- unlist(strsplit(res2[i,]$ORs, ','))
		cis <- unlist(strsplit(res2[i,]$CI, ','))
				
		# check confidence intevals
		cat ("\n")
		cat("----------------------\n")
		cat(map$IDName[map$IDNo==res2[i,]$RiskID], "\n")
		cat("Ranges : ",ranges, "\n")
		df <- chkConf(ors, cis, ranges)
		newor <- as.character(df$newor)
		ranges <- as.character(df$newrange)

		
		
		# find our points for linear regression
		for (j in c(1:length(ranges)) )
		{
			range <- as.numeric(unlist(strsplit(ranges[j], '-')) )
			or    <- as.numeric(newor[j])
			
			x  <- c(x, 0.5*(range[1] + range[2]) )
			y  <- c(y, log(or))
			
		}
			
		# regression is either linear or piecewise linear
		
		if (res2[i,]$RiskID == 3)
		{
			RiskFactorValue <- c((10*res2[i, ]$MinVal) : (10*res2[i,]$MaxVal) )/10
		}else
		{
			RiskFactorValue <- c(res2[i, ]$MinVal: res2[i,]$MaxVal)
		}
		#if (all( diff(y) >= 0 ) == TRUE | all(diff(y) <= 0 )==TRUE ) 
		if (res2[i,]$RegressionModel == 'L' )
		{
			FactorValue     <- predict(lm(y~x), data.frame(x=RiskFactorValue) )
			print(paste("Went to linear fit...", map$IDName[map$IDNo==res2[i,]$RiskID]))
		}else if (res2[i,]$RegressionModel == 'PL' )
		{
			# this is piecewise linear
			FactorValue     <- predict(piecewise.linear(x, y, middle=1, bootstrap.samples=10000), RiskFactorValue)
			print(paste("Went to piece wise...", map$IDName[map$IDNo==res2[i,]$RiskID]))
		}
		
		
		# capping if needed
		FactorValue[FactorValue < min(y) ]            <- min(y)
		FactorValue[FactorValue > max(y) ]            <- max(y)

		#plot
		plot(RiskFactorValue, FactorValue, col='blue', , type='l', main=paste(map$IDName[map$IDNo==res2[i,]$RiskID], res2[i,]$Gender, sep="  " ))
		points(x, y, col='red')
		points(xx, log(as.numeric(ors)), col='green')
		grid(NULL)
		
		Data <- rbind(Data, data.frame(RiskFactorID=res2[i,]$RiskID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, 
		Min=res2[i, ]$MinVal, Max=res2[i, ]$MaxVal, normal=res2[i, ]$NormalVal, name=map$IDName[map$IDNo==res2[i,]$RiskID], comment=paste(map$IDName[map$IDNo==res2[i,]$RiskID], res2[i,]$Gender, sep="_" ) ,
		Gender=res2[i,]$Gender )
		)
		
	}
	return(Data)
}

fitIntercept <- function(did)
{
	con    <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
	qur2     <- paste('select * from ModNmod where (upper(Genophen)=\'YES\' and DiseaseID=', did, ' and RiskID=11) order by PMID DESC')
	res2     <- dbGetQuery(con, qur2)
	dbDisconnect(con)
	
	incident        <- as.numeric(res2$ORs)
	RiskFactorValue <- c(1)
	comment         <- rep("Intercept", length(RiskFactorValue) )
	RiskFactorID    <- rep(11, length(RiskFactorValue) )
	FactorValue     <- c(log(incident/(1 - incident) ) )
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue,Min=NA, Max=NA, normal=NA, name=comment, comment=comment, Gender=res2$Gender)
	return(Data)
}

# prevalence of Disease
prevDisease <- function(did)
{
	PrevD <- data.frame(
	'9'   = 0.335, #0.1  ,
	'8'   = 0.2  ,#0.525  ,
	'7'   = 0.025,#0.001745443,
	'6'   = 0.06 ,#0.001745443,
	'5'   = 0.07 ,#0.016 ,#Alz
	'4'   = 0.06 ,#0.001743701, #CRC
	'1'   = 0.1  ,#0.002814927 #T2D
	'101' = 0.30 ,#Anx
	'102' = 0.70 ,#Bcp
	'103' = 0.13 ,#breast cancer
	'104' = 0.0030,#crohn
	'105' = 0.02 ,#Hrl
	'106' = 0.0008, #lung cancer
	'107' = 0.2  , #migrane
	'108' = 0.32 , #adult obesity
	'109' = 0.05 , #osteo
	'110' = 0.17 , #prostate cancer
	'111' = 0.03 , #apnea
	'112' = 0.001, #Ulcerative_Colitis
	'113' = 0.14  #Edy 
	)
	names(PrevD)  <- c("9", "8", "7", "6", "5", "4", "1", "101", "102", "103", "104", "105", "106","107","108", "109", "110", "111", "112", "113")
	
	
	return(as.numeric(PrevD[[as.character(did)]]))
}
	
#----------Action begins here ---------#

# Mapping
map <- getmap()

# Model RiskFactors
did   <- 111 #DiseaseID
Data1 <- categoricals(map, did) 
Data2 <- continuos(map, did)    
pd    <- prevDisease(did)
# aggregation
Data     <- rbind(Data1, Data2, fitIntercept(did))

