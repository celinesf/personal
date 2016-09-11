# implement haplotides
setwd("/home/pouria3/workspace/GenStat/trunk/Alzh")
#source("Snp.R")

AveOrCalc <- function()
{
	params <- BBsolve(rep(0.001, 6), nonlifit)
	NewOR  <- rep(NA, 6)
	
	if (params$convergence == 0 )
	{
		#print(params$par)
		for (i in seq(1, length(params$par) ) )
		{
			NewOR[i] <- params$par[i]/(1 - params$par[i]) / ( Pd/(1-Pd))
		}
		return(NewOR)
	}
	else
	{
		retrun(NULL)
	}

}
			
nonlifit <- function(p)
			
{
	k = rep(NA, 6)
	# equations 
	k[1] <- freq[1]*p[1]  + freq[2]*p[2] + freq[3]*p[3] + freq[4]*p[4] + freq[5]*p[5] + freq[6]*p[6] - Pd
	k[2] <- OR[2]*p[1]    + (1-OR[2])*p[1]*p[2] - p[2]
	k[3] <- OR[3]*p[1]    + (1-OR[3])*p[1]*p[3] - p[3]
	k[4] <- OR[4]*p[1]    + (1-OR[4])*p[1]*p[4] - p[4]
	k[5] <- OR[5]*p[1]    + (1-OR[5])*p[1]*p[5] - p[5]
	k[6] <- OR[6]*p[1]    + (1-OR[6])*p[1]*p[6] - p[6]
	k
				
}

HapComb <- function()
{
	comb <- c()
	#SNPs involved are RS429358 & RS7412
	hapal <- list(Variant=c("e2", "e3", "e4"), 
	      rs429358=list(allel=c("T", "T", "C"), genotype=c("TT", "TC", "CC")), 
	      rs7412  =list(allel=c("T", "C", "C"), genotype=c("TT", "TC", "CC")) )

	for (f in hapal$rs429358$genotype)
	{
		for (s in hapal$rs7412$genotype)
		{
			ff = unlist(strsplit(f, ""))
			ss = unlist(strsplit(s, ""))
			#print(paste(f, s))
			apoe=c()
			if (any( (hapal$rs429358$allel == ff[1]) & (hapal$rs7412$allel == ss[1]) )  & any( (hapal$rs429358$allel == ff[2]) & (hapal$rs7412$allel == ss[2]) ) )
			{
				apoe <- hapal$Variant[(hapal$rs429358$allel == ff[1]) & (hapal$rs7412$allel == ss[1])] 
				apoe <- paste(apoe, hapal$Variant[(hapal$rs429358$allel == ff[2]) & (hapal$rs7412$allel == ss[2])], sep="/" )
				comb <- c(comb, apoe)
				print(paste(apoe, paste(f, s, sep=" + "), sep=" => "))
			
			}
		
		}
	}
	return(factor(comb))
}

comb      <- HapComb()		
hap       <- read.csv("Hap.csv", header=TRUE)
haplotype <- data.frame(RiskFactorID=c(), RiskFactorValue=c(), FactorValue=c(), comment=c())
ORtable   <- data.frame(gen=c(), OR=c(), NewOR = c(), Ethnicity=c() )
id <- 115
ethmap <- list(Caucasian=0, Japanese=1, Hispanics=2, AfricanAmerican=3, NativeAmerican=4, Mixed=5, PacificIslander=6)
ethnicities <- c("Caucasian", "AfricanAmerican", "Hispanics", "Japanese") 
for (i in c(1:(length(ethnicities))) ) # fix OR for bad CIs
{
	ci <- paste("CI", i, sep="")
	hap[[ethnicities[i]]][hap[ci]==1] <- 1

	OR   <- hap[[ethnicities[i]]]
	freq <- hap[[paste(ethnicities[i], "Freq", sep="")]]
	Pd = 0.016
	NewOR <- AveOrCalc()
	haplotype <-rbind(haplotype, data.frame(RiskFactorID=id, RiskFactorValue=c(0:length(comb)), FactorValue=c(0, log(NewOR) ), ethnicity=ethmap[[ethnicities[i]]], comment=c(NA, as.character(hap$OR)) ) )
	ORtable   <-rbind(ORtable, data.frame(gen=hap$OR, OR=OR, NewOR=NewOR, Ethnicity=ethnicities[i]))
	#id <- id + 1
}
	

