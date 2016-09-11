library(BB) # for nonlinear equation solving
library(RMySQL) # for Database stuff
library(R.oo)   # for trim function

PushToGenUIDb <- function(data, Did)
{
	
	print(data[3:(length(data)-3)] ) 
	
#	con    <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
#	for (i in c(1:nrow(data)) )
#	{
#		res <- SnpInfo(data[i,]$snpID)
#		
#		qur    <- paste('insert into GenUIMap (DiseaseID, dbSNP, Gene, Genotype, Freq, ORR, NewOR,
#		BroadEthnicity, NarrowEthnicity, A_G, CI, PValue, Chromosome, Position ) values (', Did, ",'", data[i,]$snpID, "',", "'", data[i,]$name, "',", "'", data[i, ]$Genotype, "',", data[i,]$freq, ", ", round(data[i,]$OR,5), ", ", data[i,]$NewOR, 
#			", ", "'", data[i,]$BroadEthnicity, "',", "'", data[i,]$NarrowEthnicity,"', ",  "'", data[i,]$A_G, "',", "'", data[i,]$CI, "',", "'", data[i,]$PValue, "','", res$Chromosome, "','", res$Position, "')", sep='' )
			
		#print(qur)
#		res    <- dbSendQuery(con, qur)
	
#	}
	#dbCommit(con)
#	dbDisconnect(con)
	
}	


SnpGen     <- function(Did)
{
	# DB stuff
	con    <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
	qur    <- (paste('select PMID, dbSNP, BroadEthnicity, NarrowEthnicity, AG_Plus, ORValue, CI, FreqControl, FreqPop, ModelType, 
	Pvalue, IsMinorAllele from SNPOR where (Genophen like \'%es\' and DiseaseID=', Did, ') order by dbSNP, BroadEthnicity', sep='') )
	res    <- dbGetQuery(con, qur)
	dbDisconnect(con)
	
	return(res)
}


SnpInfo     <- function(snp)
{
	# DB stuff
	con    <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
	qur    <- paste("select Chromosome, Position from SNPInfo where (dbSNP = '", snp, "')", sep="")
	res    <- dbGetQuery(con, qur)
	dbDisconnect(con)
	
	return(res)
}
	


	
AveOrCalc <- function()
{
	cat("solving non linear equations ... \n Freqs Are:  ")
	print(Freq)
	params2 <- BBsolve(rep(0.001, 3), nonlifit)
	NewOR  <- rep(NA, 3)
	
	if (params2$convergence == 0 )
	{
		cat(params2$par)
		for (i in seq(1, length(params2$par) ) )
		{
			NewOR[i] <- params2$par[i]/(1 - params2$par[i]) / ( Pd/(1-Pd))
		}
		
	}
	cat('\n')
	return(NewOR)
	
}


nonlifit <- function(p)
{
	k = rep(NA, 3)
	# equations 
	k[1] <- Freq[1]*p[1]  + Freq[2]*p[2] + Freq[3]*p[3] - Pd
	k[2] <- OR[2]*p[1]    + (1-OR[2])*p[1]*p[2] - p[2]
	k[3] <- OR[3]*p[1]    + (1-OR[3])*p[1]*p[3] - p[3]
	k
	
}

EthMapHap  <- function()
{
	# map ethnicities to hapmap population
	ethn    <- c("NORTHERNHANCHINESE", "VARIOUS", "CHINESE_KOREAN", "HONGKONG_CH", "SOUTHASIAN", "JAPANESE_CHINESE_MALAY", "ASIANVARIOUS", "SOUTHEASTASIAN", "CHINESE", "HANCHINESE", "NA", "KOREAN", "SOUTHKOREAN", "JAPANESE", "CAUCASIAN", "HISPANIC", "AFRICANAMERICAN")
	hapmap_ <- c("A", "A", "A", "A", "A", "A", "A", "A", "A"      , "A"         , "A" , "A"     , "A"          , "JPT"     , "CEU"      , "MEX"     , "ASW")
	EthHap  <- data.frame(Ethnicity=ethn, HapMap=hapmap_ )
	
	return(EthHap)
}



#****** Real Deal starts from here!! **********

GenotypeCalc <- function(subres, EthHap, id, Did)
{
	snp <- subres[1,]$dbSNP
	# DB Connection
	con    <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
	
	# Get stuff from HapMAp
	cat("in Genotype sec;quering SNPInfo for frequencies ... \n")
	hapqur <- (paste('select *  from SNPInfo where dbSNP = \'', snp, '\'', sep="") )
	hap    <- dbGetQuery(con, hapqur)
	homo_major <- paste(hap$MajorAllele, hap$MajorAllele, sep="")
	het        <- paste(hap$MajorAllele, hap$MinorAllele, sep="")
	homo_minor <- paste(hap$MinorAllele, hap$MinorAllele, sep="")
	
	#get the gene name!
	qur2<- (paste('select distinct GeneNextBio, dbSNP from SNPOR  where (dbSNP=', '\'', snp, '\'', ')', sep=""))
	genename <- dbGetQuery(con, qur2)
	if (is.na(genename$Gene) )
	{
		genename$Gene <- snp
	}
	dbDisconnect(con)
		
	print(subres)
	
	Freq <- c(); OR <- c(); Geno <- c()
	for (i in c(1:nrow(subres)) )
	{
		if (is.na(subres[i,]$FreqPop) && is.na(subres[i,]$FreqControl))
		{
			#This is where we go to SNPInfo table to get frequencies reported by hapmap
			if ( nrow(hap) == 0 )
				stop("HapMap info is not available")
			
			# figure out what population we are interested in 
			if (toupper(subres[i,]$BroadEthnicity) == "ASIAN")
			{
				pop_ <- EthHap$HapMap[EthHap$Ethnicity== toupper(subres[i,]$NarrowEthnicity)]  
				pop  <- paste("Freq", pop_, sep="_")
			}else
			{
				pop_ <- EthHap$HapMap[EthHap$Ethnicity== toupper(subres[i,]$BroadEthnicity)]  
				pop  <- paste("Freq", pop_, sep="_")
			}
			
			if (pop_ == "A")
			{
				# means average asian populations
				freq <- 0; cnt <- 0
				print("Averaging freqs across all asian populations")
				freq1      <- as.numeric(unlist(strsplit(hap[["Freq_CHD"]], ",")) ) #homomajorFreq, hetFreq, homominorFreq
				freq2      <- as.numeric(unlist(strsplit(hap[["Freq_CHB"]], ",")) ) #homomajorFreq, hetFreq, homominorFreq
				if ( grep("CH", toupper(subres[i,]$NarrowEthnicity) ) == 1 )
				{
					freq3=NA
					cat("Not including JPT for CHINESE \n")
					if ( any( grep("KOREAN", toupper(subres[i,]$NarrowEthnicity) ) ) == 1 )
					{
						freq3      <- as.numeric(unlist(strsplit(hap[["Freq_JPT"]], ",")) ) #homomajorFreq, hetFreq, homominorFreq
						cat("Damn it! Its Chinese and Koreans mixed, including JPT for this mix! \n")
					}
					
				}else
				{
					freq3      <- as.numeric(unlist(strsplit(hap[["Freq_JPT"]], ",")) ) #homomajorFreq, hetFreq, homominorFreq
				}
				if (!is.null(freq1) & any(!is.na(freq1)) )
				{
					freq <- freq + 100 * freq1
					cnt  <- cnt + 100
					#print(freq)
					cat("Fetched Freq for CHD\n")
				}
				if (!is.null(freq2) & any(!is.na(freq2)) )
				{
					freq <- freq + 90 * freq2
					cnt  <- cnt + 90
					#print(freq)
					cat("Fetched Freq for CHB\n")
				}
				if (!is.null(freq3) & any(!is.na(freq3)) )
				{
					freq <- freq + 91 * freq3
					cnt  <- cnt + 91
					#print(freq)
					cat("Fetched Freq for JPT\n")
				}
				
				freq = 1.0/cnt*freq
			}else
			{
				freq       <- as.numeric(unlist(strsplit(hap[[pop]], ",")) ) #homomajorFreq, hetFreq, homominorFreq
			}
						
			if (as.integer(1000*subres[i,]$ORValue) == as.integer(1000) ) #or=1 has to be first for equations to work
			{
				OR   <- c(subres[i,]$ORValue, OR)
				Geno <- c(subres[i,]$AG_Plus, Geno)
				if (subres[i,]$AG_Plus == homo_major)
				{
					Freq <- c(freq[1], Freq)
				}else if (subres[i,]$AG_Plus == het | subres[i,]$AG_Plus == paste(hap$MinorAllele, hap$MajorAllele, sep=""))
				{
					Freq <- c(freq[2], Freq)
				}else if (subres[i,]$AG_Plus == homo_minor )
				{
					Freq <- c(freq[3], Freq)
				}else
				{
					print("Genotypes are not available in HapMap")
				}
					
			}else
			{
				OR   <- c(OR, subres[i,]$ORValue)
				Geno <- c(Geno, subres[i,]$AG_Plus)
				if (subres[i,]$AG_Plus == homo_major)
				{
					Freq <- c(Freq, freq[1] )
				}else if (subres[i,]$AG_Plus == het | subres[i,]$AG_Plus == paste(hap$MinorAllele, hap$MajorAllele, sep=""))
				{
					Freq <- c(Freq, freq[2])
				}else if (subres[i,]$AG_Plus == homo_minor )
				{
					Freq <- c(Freq, freq[3])
				}else
				{
					print("Genotypes are not available in HapMap")
				}

			}

		}else
		{
			print("Frequencies are given, I am not going to hapmap!")
			if (!is.na(subres[i,]$FreqPop) )
			{
				if (as.integer(1000*subres[i,]$ORValue) == as.integer(1000) ) #or=1 has to be first for equations to work
				{
					OR   <- c(subres[i,]$ORValue, OR)
					Freq <- c(subres[i,]$FreqPop, Freq)
					Geno <- c(subres[i,]$AG_Plus, Geno)
				}else
				{
					OR   <- c(OR, subres[i,]$ORValue)
					Freq <- c(Freq, subres[i,]$FreqPop)
					Geno <- c(Geno, subres[i,]$AG_Plus)
				}
			}else
			{
				if (as.integer(1000*subres[i,]$ORValue) == as.integer(1000) ) #or=1 has to be first for equations to work
				{
					OR   <- c(subres[i,]$ORValue, OR)
					Freq <- c(subres[i,]$FreqControl, Freq)
					Geno <- c(subres[i,]$AG_Plus, Geno)
				}else
				{
					OR   <- c(OR, subres[i,]$ORValue)
					Freq <- c(Freq, subres[i,]$FreqControl)
					Geno <- c(Geno, subres[i,]$AG_Plus)
				}
			}
		}
	}
	
	if (!is.null(Freq) & any(!is.na(Freq)) )
	{
		cat("Data Gathering Completed, Solving the equations ...\n")
		params <- list(OR=OR, Freq=Freq)
		attach(params)
		NewOR  <- AveOrCalc()
		detach(params)

		SNP_frame <- data.frame(RiskFactorID=id, RiskFactorValue=c(0, 1, 2, 3), snpID=snp, Genotype=c(NA, Geno), freq=c(0,Freq), OR=c(0, OR), NewOR=c(0,NewOR), BroadEthnicity=subres$BroadEthnicity[1], NarrowEthnicity=subres$NarrowEthnicity[1], name=genename$Gene, A_G='G', CI=c(NA, subres$CI), PValue=c(NA, subres$Pvalue) )  
		PushToGenUIDb(SNP_frame, Did)
		
		snip <-data.frame(RiskFactorID=rep(id, 4), RiskFactorValue=c(0, 1, 2, 3), FactorValue=c(0, log(NewOR)), Min=rep(0,4)  , Max=rep(3,4)  , normal=rep(0,4), name=rep(genename$Gene,4), comment=c("NA", Geno), Ethn=rep(paste(subres$BroadEthnicity[1], subres$NarrowEthnicity[1], sep="_"),4)) 
		#id <- id + 1
		
		return(snip)
		
	}else
	{
		cat("\n!!ATTENTION!!\n")
		cat("Couldnt continute with calculations ...Check HapMap Freqs\n\n")
		print(hapqur)
		print(hap)
		print(Freq)
		print(OR)
		print(Geno)
	}
	return(NA)

}


AlleleCalc <- function(subres, EthHap, id, Did)
{
	snp <- subres[1,]$dbSNP
	# DB Connection
	con    <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
	cat("in Allele sec;quering SNPInfo for Major Minor Allele ...\n ")
	hapqur <- (paste('select *  from SNPInfo where dbSNP = \'', snp, '\'', sep="") )
	hap    <- dbGetQuery(con, hapqur)
	cat("Chromosome : ", hap$Chromosome, "\n")
	
	#get the gene name!
	qur2<- (paste('select distinct GeneNextBio, dbSNP from SNPOR  where (dbSNP=', '\'', snp, '\'', ')', sep=""))
	genename <- dbGetQuery(con, qur2)
	if (is.na(genename$Gene) )
	{
		genename$Gene <- snp
	}
	dbDisconnect(con)
	
	print(subres)
	Freq <- c(); OR <- c(); Geno <- c()
	
	# Lets get busy!
	RefAllele   <- subres$AG_Plus[subres$ORValue == 1] 
	RefOR       <- subres$ORValue[subres$ORValue == 1] 
	RiskAllele  <- subres$AG_Plus[subres$ORValue != 1] 
	RiskOR      <- subres$ORValue[subres$ORValue != 1] 
	
	
	if (nrow(hap) > 0 )
	{
		# Special XY chromosome
		AdjustOR_XY(subres, hap, id, EthHap, genename$Gene )
		
		
		homo_major <- paste(hap$MajorAllele, hap$MajorAllele, sep="")
		het        <- paste(hap$MajorAllele, hap$MinorAllele, sep="")
		homo_minor <- paste(hap$MinorAllele, hap$MinorAllele, sep="")
		
		# See if paper has provided Freq for alleles
		Hardy      <- HardyFreq(subres)
		if(!is.na(Hardy) )
		{
					
			freq       <- c(Hardy[[homo_major]], Hardy$Het, Hardy[[homo_minor]] )
		}
		
		
		# figure out what population we are interested in 
		if (toupper(subres[1,]$BroadEthnicity) == "ASIAN")
		{
			pop_ <- EthHap$HapMap[EthHap$Ethnicity== toupper(subres[1,]$NarrowEthnicity)]  
			pop  <- paste("Freq", pop_, sep="_")
		}else if (subres[1,]$BroadEthnicity == "Native_American" || subres[1,]$BroadEthnicity == "NativeAmerican" )
		{
			Hardy <- HardyFreq_Ctr(subres)
			if (is.na(Hardy))
			{
				cat("!!ATTENTION!!\n")
				print("Couldnt continue with calculations... hapmap empty")
				snip      <- data.frame(RiskFactorID=c(), RiskFactorValue=c(), FactorValue=c(), Min=c()  , Max=c()  , normal=c(), name=c(), comment=c(), Ethn=c() )  
				
				return(snip)
			}else
			{
				freq  <- c(Hardy[[homo_major]], Hardy$Het, Hardy[[homo_minor]] )
				pop   <- subres[1,]$BroadEthnicity  
			}
			
			
		}else
		{
			pop_ <- EthHap$HapMap[EthHap$Ethnicity== toupper(subres[1,]$BroadEthnicity)]  
			pop  <- paste("Freq", pop_, sep="_")
		}
		
		if (!exists("freq") )
		{	
			if (pop_ == "A")
			{
				# means average asian populations
				freq <- 0; cnt <- 0
				print("Averaging freqs across all asian populations")
				freq1      <- as.numeric(unlist(strsplit(hap[["Freq_CHD"]], ",")) ) #homomajorFreq, hetFreq, homominorFreq
				freq2      <- as.numeric(unlist(strsplit(hap[["Freq_CHB"]], ",")) ) #homomajorFreq, hetFreq, homominorFreq
				if ( any( grep("CH", toupper(subres[1,]$NarrowEthnicity) )  ) )
				{
					freq3=NA
					cat("Not including JPT for CHINESE \n")
					if ( any( grep("KOREAN", toupper(subres[1,]$NarrowEthnicity) )) == 1 )
					{
						freq3      <- as.numeric(unlist(strsplit(hap[["Freq_JPT"]], ",")) ) #homomajorFreq, hetFreq, homominorFreq
						cat("Damn it! Its Chinese and Koreans mixed, including JPT for this mix! \n")
					}
						
				}else
				{
					freq3      <- as.numeric(unlist(strsplit(hap[["Freq_JPT"]], ",")) ) #homomajorFreq, hetFreq, homominorFreq
				}
				if (!is.null(freq1) & any(!is.na(freq1)) )
				{
					freq <- freq + 100 * freq1
					cnt  <- cnt + 100
					cat("Fetched Freq for CHD\n")
					#print(freq)
				}
				if (!is.null(freq2) & any(!is.na(freq2)) )
				{
					freq <- freq + 90 * freq2
					cnt  <- cnt + 90
					#print(freq)
					cat("Fetched Freq for CHB \n")
				}
				if (!is.null(freq3) & any(!is.na(freq3)) )
				{
					freq <- freq + 91 * freq3
					cnt  <- cnt + 91
					#print(freq)
					cat("Fetched Freq for JPT\n")
				}

				freq = 1.0/cnt *(freq)
						
			}else
			{
				freq       <- as.numeric(unlist(strsplit(hap[[pop]], ",")) ) #homomajorFreq, hetFreq, homominorFreq
			}
		}
		hapframe   <- data.frame(SNP=snp, Freq=freq,  Geno=c(homo_major, het, homo_minor), Pop=pop )
		
		if ( (RefAllele%in%c(hap$MajorAllele, hap$MinorAllele) ) & (RiskAllele%in%c(hap$MajorAllele, hap$MinorAllele) ) )
		{
			mtype <- subres$ModelType[1]
			if (is.na(mtype) )
			{
				mtype = "additive"
			}
			if (toupper(mtype) == "DOMINANT")
			{
				refGen     <- paste(RefAllele, RefAllele, sep="")
				refGenOr   <- 1
				refGenFreq <- hapframe$Freq[hapframe$Geno==refGen]
				
				hetGen     <- paste(RefAllele, RiskAllele, sep="")
				hetOr      <- RiskOR
				hetFreq    <- hapframe$Freq[hapframe$Geno==het]
				
				rskGen     <- paste(RiskAllele, RiskAllele, sep="")
				rskOr      <- RiskOR
				rskFreq    <- hapframe$Freq[hapframe$Geno==rskGen]
				
				# putting it all together
				OR   <- c(refGenOr, hetOr, rskOr)
				Freq <- c(refGenFreq, hetFreq, rskFreq)
				Geno <- c(refGen, het, rskGen)
			
			}else if(toupper(mtype) == "RECESSIVE")
			{
				refGen     <- paste(RefAllele, RefAllele, sep="")
				refGenOr   <- 1
				refGenFreq <- hapframe$Freq[hapframe$Geno==refGen]
				
				hetGen     <- paste(RefAllele, RiskAllele, sep="")
				hetOr      <- 1
				hetFreq    <- hapframe$Freq[hapframe$Geno==het]
				
				rskGen     <- paste(RiskAllele, RiskAllele, sep="")
				rskOr      <- RiskOR
				rskFreq    <- hapframe$Freq[hapframe$Geno==rskGen]
				
				# putting it all together
				OR   <- c(refGenOr, hetOr, rskOr)
				Freq <- c(refGenFreq, hetFreq, rskFreq)
				Geno <- c(refGen, het, rskGen)

			}else if(toupper(mtype) == "ADDITIVE")
			{
				refGen     <- paste(RefAllele, RefAllele, sep="")
				refGenOr   <- 1
				refGenFreq <- hapframe$Freq[hapframe$Geno==refGen]
				
				hetGen     <- paste(RefAllele, RiskAllele, sep="")
				hetOr      <- RiskOR
				hetFreq    <- hapframe$Freq[hapframe$Geno==het]
				
				rskGen     <- paste(RiskAllele, RiskAllele, sep="")
				rskOr      <- RiskOR^2
				rskFreq    <- hapframe$Freq[hapframe$Geno==rskGen]
				
				# putting it all together
				OR   <- c(refGenOr, hetOr, rskOr)
				Freq <- c(refGenFreq, hetFreq, rskFreq)
				Geno <- c(refGen, het, rskGen)
				
				
			}
			# solve it!
			if (!is.null(Freq) & any(!is.na(Freq)) )
			{
				print("Data Gathering Completed, Solving the equations ...")
				params <- data.frame(OR=OR, Freq=Freq)
				attach(params)
				NewOR  <- AveOrCalc()
				detach(params)
				
				# capture it
				SNP_frame <- data.frame(RiskFactorID=id, RiskFactorValue=c(0, 1, 2, 3), snpID=snp, Genotype=c(NA, Geno), freq=c(0,Freq), OR=c(0, OR), NewOR=c(0,NewOR), BroadEthnicity=subres$BroadEthnicity, NarrowEthnicity=subres$NarrowEthnicity, name=genename$Gene, A_G='A', CI=NA, PValue=NA )  
				PushToGenUIDb(SNP_frame, Did)
				
				snip      <- data.frame(RiskFactorID=rep(id, 4), RiskFactorValue=c(0, 1, 2, 3), FactorValue=c(0, log(NewOR)), Min=rep(0,4)  , Max=rep(3,4)  , normal=rep(0,4), name=rep(genename$Gene,4), comment=c("NA", Geno), Ethn=rep(paste(subres$BroadEthnicity[1], subres$NarrowEthnicity[1], sep='_'),4) )  
				
				return(snip)
			}else
			{
				cat("\n!!ATTENTION!!\n")
				cat("Couldnt continute with calculations ...Check HapMap Freq\n\n")
				print(hapqur)
				print(hap)
				print(paste("Ref Allele : ", RefAllele) )
				print(paste("Risk Allele : ", RiskAllele) )
				print(Freq)
				print(OR)
				print(Geno)
			}
				
		}else
		{
			print("Major/Minor Allele doesn't match  HapMap")
			
		}
	}else
	{
		cat("!!ATTENTION!!\n")
		print("Couldnt continue with calculations... hapmap empty")
		print(hapqur)
		print(hap)
	}
	
	return(NA)
}



HardyFreq <- function(subres)
{
	cat ("\nCalculate Genotype Frequencies based on Hardy-Weinberg principle \n")
	
	FreqPop <- subres$FreqPop
	#FreqCtr <- subres$FreqControl
	
	if (any(!is.na(FreqPop) ))
	{
		cat ("Using FreqPop ... \n")
		A1 <- subres[1,]$AG_Plus; F1 <- subres[1,]$FreqPop
		A2 <- subres[2,]$AG_Plus; F2 <- subres[2,]$FreqPop
		if (!any(c(F1, F2) ) )
			stop("freqpop not there")
		if (is.na(F1) )
			F1 <- 1 - F2
		if (is.na(F2) )
			F2 <- 1 - F1
		
		Hardy <- data.frame(A1A1=F1^2, Het=2*F1*F2, A2A2=F2^2 )
		names(Hardy) <- c(paste(A1, A1, sep=''), "Het", paste(A2, A2, sep='') )
	}#else if(any(!is.na(FreqCtr) ) )
	#{
	#	cat ("Using FreqControl ... \n")
	#	A1 <- subres[1,]$AG_Plus; F1 <- subres[1,]$FreqControl
	#	A2 <- subres[2,]$AG_Plus; F2 <- subres[2,]$FreqControl
	#	if (!any(c(F1, F2) ) )
	#	stop("freqctr not there")
	#	if (is.na(F1) )
	#		F1 <- 1 - F2
	#	if (is.na(F2) )
	#		F2 <- 1 - F1
		
	#	Hardy <- data.frame(A1A1=F1^2, Het=2*F1*F2, A2A2=F2^2 )
	#	names(Hardy) <- c(paste(A1, A1, sep=''), "Het", paste(A2, A2, sep='') )
	#}
	
	if(exists("Hardy") )
	{
		cat("Freqencies are : \n")
		print(Hardy)
		return(Hardy)
	}else
	{
		cat("No Freq Found... Going to HapMap\n")
		return(NA)
	}
}

HardyFreq_Ctr <- function(subres)
{
	cat ("\nBack to Hardy-Weinberg principle \n")
	
	#FreqPop <- subres$FreqPop
	FreqCtr <- subres$FreqControl
	
	#if (any(!is.na(FreqPop) ))
	#{
	#	cat ("Using FreqPop ... \n")
	#	A1 <- subres[1,]$AG_Plus; F1 <- subres[1,]$FreqPop
	#	A2 <- subres[2,]$AG_Plus; F2 <- subres[2,]$FreqPop
	#	if (!any(c(F1, F2) ) )
	#		stop("freqpop not there")
	#	if (is.na(F1) )
	#		F1 <- 1 - F2
	#	if (is.na(F2) )
	#		F2 <- 1 - F1
	#	
	#	Hardy <- data.frame(A1A1=F1^2, Het=2*F1*F2, A2A2=F2^2 )
	#	names(Hardy) <- c(paste(A1, A1, sep=''), "Het", paste(A2, A2, sep='') )
	#}#
	if(any(!is.na(FreqCtr) ) )
	{
		cat ("Using FreqControl ... \n")
		A1 <- subres[1,]$AG_Plus; F1 <- subres[1,]$FreqControl
		A2 <- subres[2,]$AG_Plus; F2 <- subres[2,]$FreqControl
		if (!any(c(F1, F2) ) )
			stop("freqctr not there")
		if (is.na(F1) )
			F1 <- 1 - F2
		if (is.na(F2) )
			F2 <- 1 - F1
	
		Hardy <- data.frame(A1A1=F1^2, Het=2*F1*F2, A2A2=F2^2 )
		names(Hardy) <- c(paste(A1, A1, sep=''), "Het", paste(A2, A2, sep='') )
	}
	
	if(exists("Hardy") )
	{
		cat("Freqencies are : \n")
		print(Hardy)
		return(Hardy)
	}else
	{
		cat("No Freq Found... Going to HapMap\n")
		return(NA)
	}
}



#**************************************This for the special XY chromosome*********************************
AveOrCalc_XY <- function()
{
	cat("solving non linear equations ... \n Freqs Are:  ")
	print(Freq)
	cat("ORs are ", OR, "\n")
	params2 <- BBsolve(rep(0.001, 2), nonlifit_XY)
	NewOR  <- rep(NA, 2)
	
	if (params2$convergence == 0 )
	{
		cat(params2$par)
		for (i in seq(1, length(params2$par) ) )
		{
			NewOR[i] <- params2$par[i]/(1 - params2$par[i]) / ( Pd/(1-Pd))
		}
		
	}
	cat('\n')
	return(NewOR)
	
}


nonlifit_XY <- function(p)
{
	k = rep(NA, 2)
	# equations 
	k[1] <- Freq[1]*p[1]  + Freq[2]*p[2] - Pd
	k[2] <- OR[2]*p[1]    + (1-OR[2])*p[1]*p[2] - p[2]
	k
	
}


#******** special 'XY' chromosome situatioin ******************
AdjustOR_XY <- function(subres, hap, id, EthHap, gene )
{
	snp <- subres$dbSNP[1]
		
	
	if (hap$Chromosome == 'X' )
	{
		cat("\n*********You probably have some balls if you have come this far!!!********\n")
		A1 <- subres$AG_Plus[subres$ORValue==1]; F1 <- subres$FreqPop[subres$ORValue==1]
		A2 <- subres$AG_Plus[subres$ORValue!=1]; F2 <- subres$FreqPop[subres$ORValue!=1]
		
		if (any(is.na(subres$FreqPop) ))
		{
			cat("*********Going to HapMap for X chromosome********\n")
			pop_ <- EthHap$HapMap[EthHap$Ethnicity== toupper(subres[1,]$BroadEthnicity)]  
			pop  <- paste("MajMin", pop_, sep="_")
			
			AlleleFreqs <- unlist(strsplit(hap[[pop]], ",")) 
			MajAllele   <- unlist(strsplit(AlleleFreqs[1], ":" )) 
			MinAllele   <- unlist(strsplit(AlleleFreqs[2], ":" )) 
			cat("Major Allele : ", MajAllele, "   MinorAllele : ", MinAllele, "\n")
			
			if (!subres$AG_Plus[subres$IsMinorAllele=='yes'] %in% trim(MinAllele[1]) )
				stop("Major Mainor don't match HapMap")
			
			RefAllele   <- subres$AG_Plus[subres$ORValue == 1] 
			RefOR       <- subres$ORValue[subres$ORValue == 1]
			RiskAllele  <- subres$AG_Plus[subres$ORValue != 1] 
			RiskOR      <- subres$ORValue[subres$ORValue != 1] 
			if (RefAllele == MajAllele[1])
			{
				RefFreq  <- as.numeric(MajAllele[2] )
				RiskFreq <- as.numeric(MinAllele[2] )
			}else
			{
				RefFreq  <- as.numeric(MinAllele[2] )
				RiskFreq <- as.numeric(MajAllele[2] )
			}
			
			Freq   <- c(RefFreq, RiskFreq)
			allele <- c(RefAllele, RiskAllele)
			OR     <- c(1, RiskOR)
			
		}else
		{
			cat("*********Using FreqPop for X chromosome********\n")
			riskor <- subres$ORValue[subres$ORValue!=1]
			if (is.na(F1) )
				F1 <- 1 - F2
			if (is.na(F2) )
				F2 <- 1 - F1
			
			Freq   <- c(F1, F2)
			allele <- c(A1, A2)
			OR     <- c(1, riskor) 
		}
		
		
		cat("*****Adjusted Allele OR for chromosome 'X' ********\n Solving the equations ...\n")
		params <- data.frame(OR=OR, Freq=Freq)
		attach(params)
		NewOR  <- AveOrCalc_XY()
		detach(params)
		
		# capture it
		SNP_frame <- data.frame(RiskFactorID=id, RiskFactorValue=c(0, 1), snpID=snp, Genotype=c(allele), freq=c(Freq), OR=c(OR), NewOR=c(NewOR), BroadEthnicity=subres$BroadEthnicity[1], NarrowEthnicity=subres$NarrowEthnicity[1], name=gene, A_G='A', CI=NA, PValue=NA )  
		PushToGenUIDb(SNP_frame, Did)
		
		
	}
		
		#******** End of special 'XY' chromosome situatioin ***********
		
		
		
}