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
	
	

categoricals <- function(map)
{
	# model Categorical riskFactors
	# DB stuff
	con    <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
	qur1    <- ('select * from ModNmod where (upper(Genophen)=\'YES\' and DiseaseID=8 and upper(isFactor)=\'YES\') order by PMID DESC')
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
	qur2     <- ('select * from ModNmod where (upper(Genophen)=\'YES\' and DiseaseID=8 and upper(isFactor)=\'NO\') order by PMID DESC')
	res2     <- dbGetQuery(con, qur2)
	dbDisconnect(con)
	
	# initialize Data object
	Data <- data.frame(RiskFactorID=c(), RiskFactorValue=c(), FactorValue=c(), Min=c(), Max=c(), normal=c(), name=c(), comment=c())
	
	dev.new()
	par(lty=2, lwd=2, mfrow=c(1, 3))
	
	for (i in c(1:nrow(res2) ) )
	{
		x = c(); y=c()
		ranges <- unlist(strsplit(res2[i,]$Ranges, ','))
		ors <- unlist(strsplit(res2[i,]$ORs, ','))
		if (length(ranges) != length(ors) )
			stop('ranges and ORs have different lengths')
		
		# find our points for linear regression
		for (j in c(1:length(ranges)) )
		{
			range <- as.numeric(unlist(strsplit(ranges[j], '-')) )
			or    <- as.numeric(ors[j])
			
			x <- c(0.5*(range[1] + range[2]),x )
			y <- c(log(or), y)
		}
			
		# regression is either linear or piecewise linear
		
		if (res2[i,]$RiskID == 3)
		{
			RiskFactorValue <- c((10*res2[i, ]$MinVal) : (10*res2[i,]$MaxVal) )/10
		}else
		{
			RiskFactorValue <- c(res2[i, ]$MinVal: res2[i,]$MaxVal)
		}
		if (any( (diff(y) * diff(x) > 0 ) == FALSE ) )
		{
			# this is piecewise linear, sine the sign of 1st derivative has changed
			FactorValue     <- predict(piecewise.linear(x, y, middle=1), RiskFactorValue)
		}else
		{
			FactorValue     <- predict(lm(y~x), data.frame(x=RiskFactorValue) )
		}
		
		
		# capping if needed
		FactorValue[FactorValue < min(y) ]            <- min(y)
		FactorValue[FactorValue > max(y) ]            <- max(y)

		#plot
		plot(RiskFactorValue, FactorValue, col='blue', main=map$IDName[map$IDNo==res2[i,]$RiskID])
		points(x, y, col='red')
		grid(NULL)
		
		Data <- rbind(Data, data.frame(RiskFactorID=res2[i,]$RiskID, RiskFactorValue=RiskFactorValue, FactorValue=FactorValue, 
		Min=res2[i, ]$MinVal, Max=res2[i, ]$MaxVal, normal=res2[i, ]$NormalVal, name=map$IDName[map$IDNo==res2[i,]$RiskID], comment="") )
		
	}
	return(Data)
}

fitIntercept <- function()
{
	RiskFactorValue <- c(1)
	comment         <- rep("Intercept", length(RiskFactorValue) )
	RiskFactorID    <- rep(116, length(RiskFactorValue) )
	FactorValue     <- c(log(5.25/100 /(1 - 5.25/100) ) )
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