library(BB) # for nonlinear equation solving
library(RMySQL) # for Database stuff
#setwd("/home/pouria3/workspace/GenStat/trunk/Stroke")

AveOrCalc <- function(gene)
{
	params <- BBsolve(rep(0.001, 3), nonlifit)
	NewOR  <- rep(NA, 3)
	
	if (params$convergence == 0 )
	{
		print(params$par)
		for (i in seq(1, length(params$par) ) )
		{
			NewOR[i] <- params$par[i]/(1 - params$par[i]) / ( Pd/(1-Pd))
		}
		
	}
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

# map ethnicities to hapmap population
ethn    <- c("SOUTHASIAN", "CHINESE", "SOUTHKOREAN", "JAPANESE", "CAUCASIAN", "HISPANIC", "AFRICANAMERICAN")
hapmap_ <- c("CHB", "CHB", "CHB", "JPT", "CEU", "MEX", "ASW")
EthHap  <- data.frame(Ethnicity=ethn, HapMap=hapmap_ )

# DB stuff
con    <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
qur    <- ('select PMID, dbSNP, BroadEthnicity, NarrowEthnicity, Allele_Genotype, ORValue, FreqPop, ModelType from SNPOR where (Genophen like \'%es\' and DiseaseID=6)
		order by dbSNP, BroadEthnicity')
snpqur <- ('select Distinct dbSNP from SNPOR where (Genophen like \'%es\' and DiseaseID=6)')
res    <- dbGetQuery(con, qur)
snps   <- dbGetQuery(con, snpqur)
Pd     <- exp(Data$FactorValue[Data$comment=="Intercept"])/ (1 + exp(Data$FactorValue[Data$comment=="Intercept"]) )  # prevalence of disease
SNP_frame <- data.frame(RiskFactorID=c(), RiskFactorValue=c(), snpID=c()  , Genotype=c()   ,freq=c()    , OR=c()         , NewOR=c(), ethnicity=c() )
snip <- data.frame(RiskFactorID=c(), RiskFactorValue=c(), FactorValue=c(), comment=c())
id <- 200
#GUI related
genetics_gui <- data.frame(RiskFactorID=c() , name=c()         , Min=c()  , Max=c()  , normal=c())
snpnames <- c(); genocount <- 0; alcount <- 0
for (snp in snps$dbSNP[1:1])
{
	print (snp)
	erflag <- 1
	subres  <- res[apply(res, 1, function(x) snp %in% x),] #information related to that snp only
	#print( subres )
	
	#Lets deal with sub now
	#caseI) Genotype OR OR/AND Freq is specified
	# Note that same snp could be associated with different ethnicities
	popno <- length(unique(subres$BroadEthnicity)) # this should help us figure out how many ethnicites we are dealing with
	if (nrow(subres) == ( popno*3 ) )
	{
		totalsub <- subres
		case <- 1
		while (case <= (3*popno))
		{
			subres <- totalsub[case:(case+2),]
			
			print(subres)
			genocount <- genocount + 1
			Freq <- c(); OR <- c(); Geno <- c()
			for (i in c(1:nrow(subres)) )
			{
				if (is.na(subres[i,]$FreqPop) )
				{
					#This is where we go to SNPInfo table to get frequencies reported by hapmap
					erflag <- 0
					print("quering SNPInfo for frequencies ... ")
					hapqur <- (paste('select *  from SNPInfo where dbSNP = \'', snp, '\'', sep="") )
					hap    <- dbGetQuery(con, hapqur)
					homo_major <- paste(hap$MajorAllele, hap$MajorAllele, sep="")
					het        <- paste(hap$MajorAllele, hap$MinorAllele, sep="")
					homo_minor <- paste(hap$MinorAllele, hap$MinorAllele, sep="")
					
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
						
					freq       <- as.numeric(unlist(strsplit(hap[[pop]], ",")) ) #homomajorFreq, hetFreq, homominorFreq
					
					if (as.integer(1000*subres[i,]$ORValue) == as.integer(1000) ) #or=1 has to be first for equations to work
					{
						OR   <- c(subres[i,]$ORValue, OR)
						Geno <- c(subres[i,]$Allele_Genotype, Geno)
						if (subres[i,]$Allele_Genotype == homo_major)
						{
							Freq <- c(freq[1], Freq)
							}else if (subres[i,]$Allele_Genotype == het | subres[i,]$Allele_Genotype == paste(hap$MinorAllele, hap$MajorAllele, sep=""))
						{
							Freq <- c(freq[2], Freq)
						}else if (subres[i,]$Allele_Genotype == homo_minor )
						{
							Freq <- c(freq[3], Freq)
						}else
						{
							print("Genotypes are not available in HapMap")
						}
							
					}else
					{
						
						OR   <- c(OR, subres[i,]$ORValue)
						Geno <- c(Geno, subres[i,]$Allele_Genotype)
						if (subres[i,]$Allele_Genotype == homo_major)
						{
							Freq <- c(Freq, freq[1] )
						}else if (subres[i,]$Allele_Genotype == het | subres[i,]$Allele_Genotype == paste(hap$MinorAllele, hap$MajorAllele, sep=""))
						{
							Freq <- c(Freq, freq[2])
						}else if (subres[i,]$Allele_Genotype == homo_minor )
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
					if (as.integer(1000*subres[i,]$ORValue) == as.integer(1000) ) #or=1 has to be first for equations to work
					{
						OR   <- c(subres[i,]$ORValue, OR)
						Freq <- c(subres[i,]$FreqPop, Freq)
						Geno <- c(subres[i,]$Allele_Genotype, Geno)
					}else
					{
						OR   <- c(OR, subres[i,]$ORValue)
						Freq <- c(Freq, subres[i,]$FreqPop)
						Geno <- c(Geno, subres[i,]$Allele_Genotype)
					}
				}
			}
			if (!is.null(Freq) & any(!is.na(Freq)) )
			{
				print("Data Gathering Completed, Solving the equations ...")
				NewOR  <- AveOrCalc()
				SNP_frame <- rbind(SNP_frame, data.frame(RiskFactorID=id, RiskFactorValue=c(1, 2, 3), snpID=snp, Genotype=Geno, freq=Freq, OR=OR, NewOR=NewOR, ethnicity=subres$BroadEthnicity) )
				snip <-rbind(snip, data.frame(RiskFactorID=rep(id, 4), RiskFactorValue=c(0, 1, 2, 3), FactorValue=c(0, log(NewOR)), comment=c("NA", Geno)) )
				#id <- id + 1
				
			}
			case <- case + 3
		}		
	}else if(nrow(subres) == 2*popno) #case  II) allele OR is specified
	{
		totalsub <- subres
		case <- 1
		while (case <= (2*popno))
		{
			subres <- totalsub[case:(case+1),]
			print(subres)
			Freq <- c(); OR <- c(); Geno <- c()
			
			# Lets get busy!
			RefAllele   <- subres$Allele_Genotype[subres$ORValue == 1] 
			RefOR       <- subres$ORValue[subres$ORValue == 1] 
			RiskAllele  <- subres$Allele_Genotype[subres$ORValue != 1] 
			RiskOR      <-  subres$ORValue[subres$ORValue != 1] 
			
			print("quering SNPInfo for Major Minor Allele ... ")
			hapqur <- (paste('select *  from SNPInfo where dbSNP = \'', snp, '\'', sep="") )
			hap    <- dbGetQuery(con, hapqur)
			if (nrow(hap) > 0 )
			{
				
				homo_major <- paste(hap$MajorAllele, hap$MajorAllele, sep="")
				het        <- paste(hap$MajorAllele, hap$MinorAllele, sep="")
				homo_minor <- paste(hap$MinorAllele, hap$MinorAllele, sep="")
				
				# figure out what population we are interested in 
				if (toupper(subres[1,]$BroadEthnicity) == "ASIAN")
				{
					pop_ <- EthHap$HapMap[EthHap$Ethnicity== toupper(subres[1,]$NarrowEthnicity)]  
					pop  <- paste("Freq", pop_, sep="_")
				}else
				{
					pop_ <- EthHap$HapMap[EthHap$Ethnicity== toupper(subres[1,]$BroadEthnicity)]  
					pop  <- paste("Freq", pop_, sep="_")
				}
				
				freq       <- as.numeric(unlist(strsplit(hap[[pop]], ",")) ) #homomajorFreq, hetFreq, homominorFreq
				hapframe   <- data.frame(SNP=snp, Freq=freq,  Geno=c(homo_major, het, homo_minor), Pop=pop )
				
				if ( (RefAllele%in%c(hap$MajorAllele, hap$MinorAllele) ) & (RiskAllele%in%c(hap$MajorAllele, hap$MinorAllele) ) )
				{
					mtype <- subres$ModelType[1]
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
						Geno <- c(refGen, hetGen, rskGen)
						
						# solve it!
						print("Data Gathering Completed, Solving the equations ...")
						NewOR  <- AveOrCalc()
						
						# capture it
						SNP_frame <- rbind(SNP_frame, data.frame(RiskFactorID=id, RiskFactorValue=c(1, 2, 3), snpID=snp, Genotype=Geno, freq=Freq, OR=OR, NewOR=NewOR, ethnicity=subres$BroadEthnicity) )
						snip <-rbind(snip, data.frame(RiskFactorID=rep(id, 4), RiskFactorValue=c(0, 1, 2, 3), FactorValue=c(0, log(NewOR)), comment=c("NA", Geno)) )
						#id <- id + 1
						
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
						Geno <- c(refGen, hetGen, rskGen)
						
						# solve it!
						print("Data Gathering Completed, Solving the equations ...")
						NewOR  <- AveOrCalc()
						
						# capture it
						SNP_frame <- rbind(SNP_frame, data.frame(RiskFactorID=id, RiskFactorValue=c(1, 2, 3), snpID=snp, Genotype=Geno, freq=Freq, OR=OR, NewOR=NewOR, ethnicity=subres$BroadEthnicity) )
						snip <-rbind(snip, data.frame(RiskFactorID=rep(id, 4), RiskFactorValue=c(0, 1, 2, 3), FactorValue=c(0, log(NewOR)), comment=c("NA", Geno)) )
						#id <- id + 1
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
						Geno <- c(refGen, hetGen, rskGen)
						
						# solve it!
						print("Data Gathering Completed, Solving the equations ...")
						NewOR  <- AveOrCalc()
						
						# capture it
						SNP_frame <- rbind(SNP_frame, data.frame(RiskFactorID=id, RiskFactorValue=c(1, 2, 3), snpID=snp, Genotype=Geno, freq=Freq, OR=OR, NewOR=NewOR, ethnicity=rep(subres$BroadEthnicity[1], 3)) )
						snip <-rbind(snip, data.frame(RiskFactorID=rep(id, 4), RiskFactorValue=c(0, 1, 2, 3), FactorValue=c(0, log(NewOR)), comment=c("NA", Geno)) )
						#id <- id + 1
					}
					alcount <- alcount + 1
						#print(subres)
						#print("Observed an Allele; Still under construction!")
				}else
				{
					print("Major/Minor Allele doesn't match  HapMap")
				}
			}
			case <- case + 2
		}
	}	
id <- id + 1
}	


# GUI related
genetics_gui <- data.frame(RiskFactorID=c() , name=c()         , Min=c()  , Max=c()  , normal=c())
snpnames <- c()
snpunits <- c()
snpnamecounter <- 2
for (snp in unique(SNP_frame$snpID) )
{
	qur2<- (paste('select dbSNP from SNPInfo  where (dbSNP=', '\'', snp, '\'', ')', sep=""))
	genename <- dbGetQuery(con, qur2)
	genename$Gene = "NA"
	if (is.na(genename$Gene) )
	{
		genename$Gene <- "NA"
	}
	if (genename$Gene%in%snpnames)
	{
		genename$Gene <- paste(genename$Gene, snpnamecounter, sep="_")
		snpnamecounter <- snpnamecounter + 1
	}
	snpnames <- c(snpnames, genename$Gene )
	tmp <- as.character(snip$comment[snip$RiskFactorID==unique(SNP_frame$RiskFactorID[SNP_frame$snpID==snp]) ] )
	tmp1 <- paste(tmp[1], tmp[2], tmp[3], tmp[4])
	snpunits <- c(snpunits, tmp1 )
	
	genetics_gui <- rbind(genetics_gui, data.frame(RiskFactorID=unique(SNP_frame$RiskFactorID[SNP_frame$snpID==snp]) , 
						       name=genename$Gene         , Min=0  , Max=3  , normal=0) )
	
}

# for demo 
if (nrow(genetics_gui) > 4 )
{
	genetics_gui <- genetics_gui[nrow(genetics_gui):(nrow(genetics_gui)-3),]
	snpnames     <- snpnames[length(snpnames):(length(snpnames)-3)]
	snpunits     <- snpunits[length(snpunits):(length(snpunits)-3)]
}
print(paste("Number of Snps processed is : ", length(snps$dbSNP)))
dbDisconnect(con)
Data <- rbind(Data, snip)


