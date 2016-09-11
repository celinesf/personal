rm(list=ls())
setwd("/home/pouria3/workspace/algophen/trunk/R/SNProc")
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
con      <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
bible    <- dbGetQuery(con, "select IDNo, IDName from IDmap where IsDiseaseID='yes' ")
dbDisconnect(con)


Report <- data.frame(
"9"   =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==9]  , sep=''),
"8"   =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==8]  , sep=''),
"7"   =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==7]  , sep=''),
"6"   =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==6]  , sep=''),
"5"   =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==5]  , sep=''),
"4"   =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==4]  , sep=''),
"1"   =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==1]  , sep=''),
"101" =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==101], sep=''),
"102" =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==102], sep=''),
"103" =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==103], sep=''),
"104" =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==104], sep=''),
"105" =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==105], sep=''),
"106" =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==106], sep=''),
"107" =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==107], sep=''),
"108" =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==108], sep=''),
"109" =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==109], sep=''),
"110" =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==110], sep=''),
"111" =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==111], sep=''),
"112" =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==112], sep=''),
"113" =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==113], sep=''),
"114" =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==114], sep=''),
"115" =paste("/home/pouria3/Dropbox/EngExl/"  , bible$IDName[bible$IDNo==115], sep='')
)
# prev of Disease
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
'111' = 0.20 , #apnea
'112' = 0.001, #Ulcerative_Colitis
'113' = 0.14  #Edy 
)

names(PrevD)  <- c("9", "8", "7", "6", "5", "4", "1", "101", "102", "103", "104", "105", "106","107","108", "109", "110", "111", "112", "113")
names(Report) <- c("9", "8", "7", "6", "5", "4", "1", "101", "102", "103", "104", "105", "106","107","108", "109", "110", "111", "112", "113")

for (Did in c("107") )
{
	#Redirect output to a file
	cat('====', paste(as.character(Report[[Did]]), "/", Did, "_GenCalcReport_7.txt", sep=""), '====', '\n')
	sink(paste(as.character(Report[[Did]]), "/", Did, "_GenCalcReport_7.txt", sep=""), append=FALSE, split=TRUE)
	
	# Get all the Genophen=yes snps
	cat("\n DiseaseID : "); cat(Did); cat("\n")
	Pd <- as.numeric(PrevD[[Did]] )
	res    <- SnpGen(as.numeric(Did))

	# get hapMap mapping
	EthHap <- EthMapHap()

	# Check CIs 
	res    <- chkConf(res)

	# Lets put things into perspective!
	snip <- data.frame(RiskFactorID=c(), RiskFactorValue=c(), FactorValue=c(), Min=c()  , Max=c()  , normal=c(), name=c(), comment=c(), Ethn=c() )

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

cat('**********************************************************************************************************************\n\n')
cat(paste("Number of Snps processed is : ", length(unique(res$dbSNP))) )
cat("\n")
#ouput back to terminal
sink()


}



