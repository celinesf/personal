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
	
	
chkConf <- function(ors, cis)
{
	# Funtion to check confidence intervals, if it corsses 1, set OR back to 1
	#cis refers to Confidence intervals, ors refere to odds ratios
	newor = c(); k <- 0
	print(cis)
	print(ors)
	if ( (length(ors)-1) != length(cis) )
		stop('CIs and ORs dont have the same length')
	
	for (j in c(1:length(ors)) )
	{
		if (as.numeric(ors[j]) == 1)
		{
			newor <- c(newor, ors[j])
			k <-1 # set a flag since there is no ci for reference or
		}else
		{
			ci  <- as.numeric(unlist(strsplit(cis[j-k], '-')) )
			or  <- as.numeric(ors[j])
			if ( (ci[1] <= 1) & (ci[2] >= 1) )
			{
				newor <- c(newor, 1)
			}else
			{
				newor <- c(newor, or)
			}
		}
	}
	print(newor) 
	return(newor)
		
}


categoricals <- function(map)
{
	# model Categorical riskFactors
	# DB stuff
	con    <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
	qur1    <- ('select * from ModNmod where (upper(Genophen)=\'YES\' and DiseaseID=9 and upper(isFactor)=\'YES\') order by PMID DESC')
	res1    <- dbGetQuery(con, qur1)
	dbDisconnect(con)
	
	# initialize Data object
	Data <- data.frame(RiskFactorID=c(), RiskFactorValue=c(), FactorValue=c(), Min=c(), Max=c(), normal=c(), name=c(), comment=c())
	
	par(lty=2, lwd=2, mfrow=c(ceiling(nrow(res1)/3), 3))
	
	for (i in c(1:nrow(res1)) )
	{
		cats  <- unlist(strsplit(res1[i,]$Categories, ',')) 
		ors   <- log(as.numeric(unlist(strsplit(res1[i,]$ORs, ',')) ))
		if (length(cats) != length(ors) )
			stop('categories and ORs have different lengths')
		
		Data <- rbind(Data, data.frame(RiskFactorID=res1[i,]$RiskID, RiskFactorValue=c(0:(length(cats)-1)), FactorValue=ors, 
		Min=0, Max=length(cats)-1, normal=res1[i, ]$NormalVal, name=map$IDName[map$IDNo==res1[i,]$RiskID], comment=cats ) )
		
		#plot
		barplot(ors, names.arg=cats, col='blue', main=map$IDName[map$IDNo==res1[i,]$RiskID])
	}
	Data$RiskFactorValue <- as.character(Data$RiskFactorValue)
	return(Data)
	
}




continuos <- function(map)
{
	# model continuous riskFactors
	con    <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
	qur2     <- ('select * from ModNmod where (upper(Genophen)=\'YES\' and DiseaseID=9 and upper(isFactor)=\'NO\') order by PMID DESC')
	res2     <- dbGetQuery(con, qur2)
	dbDisconnect(con)
	
	# initialize Data object
	Data <- data.frame(RiskFactorID=c(), RiskFactorValue=c(), FactorValue=c(), Min=c(), Max=c(), normal=c(), name=c(), comment=c())
	
	dev.new()
	par(lty=2, lwd=2, mfrow=c(ceiling(nrow(res2)/3), 3))
	
	for (i in c(1:nrow(res2) ) )
	{
		x = c(); y=c()
		ranges <- unlist(strsplit(res2[i,]$Ranges, ','))
		ors <- unlist(strsplit(res2[i,]$ORs, ','))
		cis <- unlist(strsplit(res2[i,]$CI, ','))
				
		# check confidence intevals
		print ("")
		print("----------------------")
		print(map$IDName[map$IDNo==res2[i,]$RiskID])
		newor <- chkConf(ors, cis)
		print(ranges)
		
		if (length(ranges) != length(newor) )
			stop('ranges and ORs have different lengths')
		
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
		if (res2[i,]$RiskID%in% c(7, 9, 14, 41) )
		{
			FactorValue     <- predict(lm(y~x), data.frame(x=RiskFactorValue) )
			print(paste("Went to linear fit...", map$IDName[map$IDNo==res2[i,]$RiskID]))
		}else
		{
			# this is piecewise linear
			FactorValue     <- predict(piecewise.linear(x, y, middle=1, bootstrap.samples=10000), RiskFactorValue)
			print(paste("Went to piece wise...", map$IDName[map$IDNo==res2[i,]$RiskID]))
		}
		
		
		# capping if needed
		#FactorValue[RiskFactorValue > floor(max(x))  ]      <- y[x==max(x)]
		#FactorValue[RiskFactorValue <= ceiling(min(x))]     <- y[x==min(x)]
		FactorValue[FactorValue < min(y) ]            <- min(y)
		FactorValue[FactorValue > max(y) ]            <- max(y)
		#plot
		plot(RiskFactorValue, FactorValue, col='blue', , type='l', main=paste(map$IDName[map$IDNo==res2[i,]$RiskID], res2[i,]$Gender, sep="  " ))
		points(x, y, col='red')
		points(x, log(as.numeric(ors)), col='green')
		grid(NULL)
		
		Data <- rbind(Data, data.frame(RiskFactorID=res2[i,]$RiskID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, 
		Min=res2[i, ]$MinVal, Max=res2[i, ]$MaxVal, normal=res2[i, ]$NormalVal, name=map$IDName[map$IDNo==res2[i,]$RiskID], comment=paste(map$IDName[map$IDNo==res2[i,]$RiskID], res2[i,]$Gender, sep="_" )) )
		
	}
	return(Data)
}

fitIntercept <- function()
{
	#pouria found this paper, PMID=16533126
	RiskFactorValue <- c(1)
	comment         <- rep("Intercept", length(RiskFactorValue) )
	RiskFactorID    <- rep(116, length(RiskFactorValue) )
	FactorValue     <- c(log(10/100 /(1 - 10/100) ) )
	Data <- data.frame(RiskFactorID=RiskFactorID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue,Min=NA, Max=NA, normal=NA, name=comment, comment=comment)
	return(Data)
}
# intercept : look at Hasin_2005_16203955
	
#----------Action begins here ---------#

# Mapping
map <- getmap()

# Model RiskFactors
Data1 <- categoricals(map) 
Data2 <- continuos(map)    

# aggregation
Data  <- rbind(Data1, Data2, fitIntercept())