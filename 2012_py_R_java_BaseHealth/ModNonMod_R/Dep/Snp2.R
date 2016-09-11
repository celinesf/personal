#rm(list=ls())
#setwd("/home/pouria3/workspace/GenStat/trunk/Hyp")
source("Snp_functions.R")


# function to check confidence interval before sending for calculations
chkConf <- function(res)
{	
	# Funtion to check confidence intervals, if it corsses 1, set OR back to 1
	#cis refers to Confidence intervals, ors refere to odds ratios
	res_cp <- res
	for (i in c(1:nrow(res)) )
	{
		if (!is.na(res[i,]$CI) )
		{
			if (res[i,]$CI != "")
			{
				ci  <- as.numeric(unlist(strsplit(res[i,]$CI, '-')) )
				if (!is.na(ci) )
				{
					if ( (ci[1] <= 1) & (ci[2] >= 1) )
					{
						cat ("Found a statistically insignificant OR : \n")
						print( res[i,] ) 
						res_cp[i,]$ORValue <- 1.000001
					}
				}
			}
		}
	}
	return(res_cp)
}
# where report files go:
Report <- data.frame("9"="/home/pouria3/Dropbox/EngExl/Hyp/",
"8"="/home/pouria3/Dropbox/EngExl/DP/",
"7"="/home/pouria3/Dropbox/EngExl/Strk/",
"6"="/home/pouria3/Dropbox/EngExl/HD/",
"5"="/home/pouria3/Dropbox/EngExl/Alz/",
"4"="/home/pouria3/Dropbox/EngExl/CRC/",
"1"="/home/pouria3/Dropbox/EngExl/T2D/"
)
# prev of Disease
PrevD <- data.frame('9'= 0.335, #0.1  ,
'8'= 0.2  ,#0.525  ,
'7'= 0.025,#0.001745443,
'6'= 0.06 ,#0.001745443,
'5'= 0.07 ,#0.016 ,#Alz
'4'= 0.06 ,#0.001743701, #CRC
'1'= 0.1  #0.002814927 #T2D
)
names(PrevD)  <- c("9", "8", "7", "6", "5", "4", "1")
names(Report) <- c("9", "8", "7", "6", "5", "4", "1")
for (Did in c("8") )
{
	#Redirect output to a file
	#sink(paste(as.character(Report[[Did]]), Did, "_GenCalcReport.txt", sep=""), append=FALSE, split=TRUE)
	
	# Get all the Genophen=yes snps
	cat("\n DiseaseID : "); cat(Did); cat("\n")
	Pd <- as.numeric(PrevD[[Did]] )
	res    <- SnpGen(as.numeric(Did))

	# get hapMap mapping
	EthHap <- EthMapHap()

	# Check CIs 
	res    <- chkConf(res)

	# Lets put Shit into perspective!
	snip <- data.frame(RiskFactorID=c(), RiskFactorValue=c(), FactorValue=c(), Min=c()  , Max=c()  , normal=c(), name=c(), comment=c() )

	for (snp in (unique(res$dbSNP) ) )
	{
		#print(snp)
		id <- snp #randomly assigned ... might change later
		subres  <- res[apply(res, 1, function(x) snp %in% x),] #information related to that snp only
	
		#Ok, subres could contain various ethnicities
		for (eth in unique(subres$BroadEthnicity) )
		{
			sub  <- subres[apply(subres, 1, function(x) eth %in% x),] #information related to that ethnicity only
			#print("")
			cat('**********************************************************************************************************************\n\n')
			#print(sub)
			#now we need to figure out if it is Genotype related or allele related
			if (nrow(sub) == 2 )
			{
				snip <- rbind(snip, AlleleCalc(sub, EthHap, id, as.numeric(Did) ) )
			}else if (nrow(sub) == 3 )
			{
				snip <- rbind(snip, GenotypeCalc(sub, EthHap, id, as.numeric(Did) ) )
			}else
			{
				cat("WTF!!!\n")
				cat("I am guessing same snp, same broadethnicity but different narrow ethnicity ... \n")
				print(sub)
				for (eth1 in unique(sub$NarrowEthnicity) )
				{
					sub1  <- sub[apply(sub, 1, function(x) eth1 %in% x),] #information related to that ethnicity only
					cat("$$$$$$$\n")
					if (nrow(sub1) == 2 )
					{
						snip <- rbind(snip, AlleleCalc(sub1, EthHap, id, as.numeric(Did) ) )
					}else if (nrow(sub1) == 3 )
					{
						snip <- rbind(snip, GenotypeCalc(sub1, EthHap, id) )
					}else
					{
						stop("This is it! I am stoping ...")
					
					}
							
				}		
			}
		
		}
	
	}

}
		
# for demo 
if (nrow(snip) > 20 )#
{
	snip <- snip[1:20,]
}

cat('**********************************************************************************************************************\n\n')
cat(paste("Number of Snps processed is : ", length(unique(res$dbSNP))) )
Data <- rbind(Data, snip[1:length(Data)] )
#ouput back to terminal
#sink()




