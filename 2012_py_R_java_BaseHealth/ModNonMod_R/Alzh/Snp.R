library(BB) # for nonlinear equation solving
setwd("/home/pouria3/workspace/GenStat/trunk/Alzh")


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
	k[1] <- freq[1]*p[1]  + freq[2]*p[2] + freq[3]*p[3] - Pd
	k[2] <- OR[2]*p[1]    + (1-OR[2])*p[1]*p[2] - p[2]
	k[3] <- OR[3]*p[1]    + (1-OR[3])*p[1]*p[3] - p[3]
	k
	
}



SNP <- read.csv("SNP.csv", header=TRUE)
genotypes <- c("Genotype1", "Genotype2", "Genotype3")
SNP_frame <- list()
snip <- data.frame(RiskFactorID=c(), RiskFactorValue=c(), FactorValue=c(), comment=c())
id <- 113
Pd     <- 0.016  # prevalence of disease

for (r in c(1:2) )
{
	freq <- c();OR   <- c(); geno <- c()
	for (gen in genotypes)
	{
		if (SNP[r, ][paste(gen, "OR", sep="")] == 1 )
		{
			OR   <- c(as.numeric(SNP[r,][paste(gen, "OR",   sep="")]), OR) 
			freq <- c(as.numeric(SNP[r,][paste(gen, "Freq", sep="")]), freq)
			geno <- c(as.character(SNP[r,][gen][[1]]), geno )
		}
		else
		{
			OR   <- c(as.numeric(SNP[r,][paste(gen, "OR",   sep="")]), OR )
			freq <- c(as.numeric(SNP[r,][paste(gen, "Freq", sep="")]), freq )
			geno <- c(as.character(SNP[r,][gen][[1]]), geno )
		}

	}
	NewOR  <- AveOrCalc()
	SNP_frame <- rbind(list(snpID=SNP[r,]$SNPID, Genotype=geno, OR=OR, NewOR=NewOR), SNP_frame )
	
	snip <-rbind(snip, data.frame(RiskFactorID=rep(id, 4), RiskFactorValue=c(0, 1, 2, 3), FactorValue=c(0, log(NewOR)), comment=c("NA", geno)) )
	id <- id + 1
	
	
}






