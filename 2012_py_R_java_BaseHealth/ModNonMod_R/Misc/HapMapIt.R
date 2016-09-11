library(biomaRt)
library(RMySQL)

# Phase1, get the rows from SNPOR table 
con    <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
snpqur <- ('select Distinct dbSNP from SNPOR order by dbSNP DESC')
DBres    <- dbGetQuery(con, snpqur)


# Phase2, for each dbSNP find chromosome, position, Major/Minor allele from total population
hapmap <- useMart("HapMap_rel27", dataset="hm27_variation") 
res  <- getBM(attributes=c("chrom", "start", "strand", "alleles", "other_allele", "refallele", "ref_allele", "ref_allele_freq", 
			   "other_allele_freq"), filters=c("marker_name"), values=list(dbsnp), mart=hapmap)



# Phase3, find frequency for each population
datasets <- c("hm27_variation_yri", "hm27_variation_tsi", "hm27_variation_chd", "hm27_variation_chb", "hm27_variation_asw",
	      "hm27_variation_ceu", "hm27_variation_mkk", "hm27_variation_lwk", "hm27_variation_jpt", "hm27_variation_mex", "hm27_variation_gih"
	      )
	      
freq     <- data.frame(yri=c(), tsi=c(), chd=c(), chb=c(), asw=c(), ceu=c(), 
		       mkk=c(), lwk=c(), jpt=c(), mex=c(), gih=c()
		      )

for (dataset in datasets)
{	
	ext      <- unlist(strsplit(dataset, split="_"))[2]
	hapmap   <- useMart("HapMap_rel27", dataset)



# Phase4, populate the SNPInfo database with all the information
snpqur <- (paste('update SNPInfo set Chromosome=', chromosome, ', Position=', position, ', HS_version=HapMap27', 
		 ', MajorAllele=', majorallele, ', MinorAllele=', minorallele, 
		 ', Freq_ASW=', freq$asw, ', Freq_CEU=', freq$ceu, ', Freq_CHB=', freq$chb, ', Freq_CHD=', freq$chd, 
		 ', Freq_GIH=', freq$gih, ', Freq_JPT=', freq_jpt, ', Freq_LWK=', freq$lwk, ', Freq_MEX=', freq$mex,
		 ', Freq_MKK=', freq$mkk, ', Freq_TSI=', freq_tsi, ', Freq_YRI=', freq$yri, 
		 ' where dbSNP=', res$dbSNP, sep="")
	  )