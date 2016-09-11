# Experimening with Celine's scoring system
library(RMySQL) # for Database stuff


con    <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
qur1   <- ("select * from SNPOR where Genophen='yes' order by dbSNP")
res    <- dbGetQuery(con, qur1)
dbDisconnect(con)

# Diseases in the system so far:
Did <- c("9", "8", "7", "6", "5", "4", "1", "101", "102", "103", "104", "105", "106","107","108", "109", "110", "111", "112", "113")

# Lets just look at Snps and Scores
snpsc       <- data.frame()


for (did in Did)
{
	snpor_disease <- res[res$DiseaseID==did, ]
	
	for (snp in (unique(snpor_disease$dbSNP) ) )
	{
		mt<-c(); re<-c(); pm<-c();rep<-c(); zee<-c();
		subres  <- snpor_disease[apply(snpor_disease, 1, function(x) snp %in% x),] #information related to that snp 
		row     <- subres[subres$ORValue==max(subres$ORValue),]
		if (nrow(row) > 1)
			row <- row[1,]
		#print(row)
			
		st     <- as.numeric(unlist(strsplit(unlist(strsplit(row$Score, "st="))[2], " "))[1])
		if (!is.na(st) )
		{
			if (as.numeric(st) >= 2)
			{
				mt     <- as.numeric(unlist(strsplit(unlist(strsplit(row$Score, "mt="))[2], "/"))[1])
				re     <- as.numeric(unlist(strsplit(unlist(strsplit(row$Score, "re="))[2], "-"))[1])
				pm     <- as.numeric(unlist(strsplit(unlist(strsplit(row$Score, "pm="))[2], "/"))[1])
				if (is.na(mt))
					mt <- 0
				if (is.na(re) )
					re <- 0
				if (is.na(pm))
					pm <- 0
				if (mt==0 | re==0 )
				{
					rep <- mt + re + pm
				}else
				{
					rep <- mt + re
				}
				
			}else if( as.numeric(st) <= 1 )
			{
				mt     <- as.numeric(unlist(strsplit(unlist(strsplit(row$Score, "mt="))[2], "/"))[1])
				re     <- as.numeric(unlist(strsplit(unlist(strsplit(row$Score, "re="))[2], "-"))[1])
				pm     <- as.numeric(unlist(strsplit(unlist(strsplit(row$Score, "pm="))[2], "/"))[1])
				if (is.na(mt))
				mt <- 0
				if (is.na(re) )
				re <- 0
				if (is.na(pm))
				pm <- 0
				
				rep <- mt + re + pm
				
				
			}
			
		}else
		{
			cat("no st for ", "Disease ", did, "\n")
			print(row)
		}
		
		#Confidence interval or Pvalue
		if (!is.na(row$CI) )
		{
			# estimate(OR) divide by standard err (Z)
			ci  <- as.numeric(unlist(strsplit(as.character(row$CI), "-")))
			sc  <- log(row$ORValue)/((ci[2] - ci[1])/(2 * 1.96)) # CI ~ 4*std
			
		}else
		{
			pval <- row$Pvalue
			zee  <- abs(qnorm(pval/2))
			if (is.na(zee))
				zee <- 0
			
			smpl <- row$SampleSize
			s    <- row$ORValue * sqrt(smpl) / zee
			
			sc   <- log(row$ORValue)/s
		}
			
		#cat(row$dbSNP, " <> ", zee, " <> ", rep, "\n")
			
		snpsc  <- rbind(snpsc, data.frame(Did=did, PMID=row$PMID, dbSNP=row$dbSNP, BE=row$BroadEthnicity, OR=row$ORValue, Rep=rep, CI=row$CI, Pval=row$Pvalue, Point=sc*sqrt(rep) ) )
	}
	snpsc_nonzero <- snpsc[snpsc$Did==did & snpsc$Point != 0,]
	cat("Average Score for DiseaseID ", did, " is : ", sum(snpsc_nonzero$Point)/nrow(snpsc_nonzero) , "\n")
	
}
	


	
