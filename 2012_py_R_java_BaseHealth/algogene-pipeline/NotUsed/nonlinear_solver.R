library(rjson)
library(BB)
library(quadprog)

#setwd("/home/pouria3/Dropbox")



prepareInputs <- function(did)
{
	msg = sprintf(' R Function: prepareInputs -- disease: %s', did)
	print(msg)
	equations <- function(x)
	{
		n <- length(x)
		f <- rep(NA, n)
		
		eq   <- 2:n
		
		f[1]  <- freq %*% x - pd
		f[eq] <- or[eq]*x[1]  + (1-or[eq])*x[1] *x[eq]  - x[eq]
		f
	}
	# first read the json file and convert to R object
	jsn <- fromJSON(file=paste(did, "_Hap_Geno.json", sep='') ) 
	fin <- c() #final result
	
	# now loop through to find our reference (equations need to be ordered, x[1] is the reference
	for ( hap in names(jsn) )
	{
		for (ethn in names(jsn[[hap]]) )
		{
			msg = sprintf(' disease: %s, hap: %s, ethn: %s (prepareInputs)', did, hap, ethn)
			print(msg)
			or <- c(); freq <- c(); gent <- c()
			for (gt in names(jsn[[hap]][[ethn]]) )
			{
				gtinfo <- jsn[[hap]][[ethn]][[gt]] 
				# find the reference Genotype. Genotype with OR=1 has to be first
				if (gtinfo$CIissue == "No" & gtinfo$OR==1)
				{
					or   <- c(gtinfo$OR      , or)
					freq <- c(gtinfo$FControl, freq)
					gent <- c(gt             , gent)
				}else
				{
					or   <- c(or  , gtinfo$OR)
					freq <- c(freq, gtinfo$FControl)
					gent <- c(gent, gt)
				}
			}
			
			# find the prevalence of disease
			pd <- prevDisease(did)
			
			# solve the equations
			#sol <- BBsolve(par=rep(0.001,6), method=1, control=list(M=500, NM=FALSE), fn=equations)
			sol <- BBsolve(par=or/(1+or), method=1, control=list(M=500, NM=FALSE), fn=equations)
			
			if (sol$convergence == 0)
			{
				# calculate adjusted OR
				newor <- rep(NA, length(or) )
				for (i in seq(1, length(sol$par) ) )
				{
					newor[i] <- sol$par[i]/(1 - sol$par[i]) / ( pd/(1-pd))
				}
				# aggregate results into a data frame
				fin <- rbind(fin, data.frame(Hap=hap, Ethn=ethn, GT=gent, OR=or, NewOR=newor, Freq=freq) )
			}else
			{
				cat('\nEquations Didnt merge\n', hap, '\n', or, '\n', freq, '\n', gent, '\n')
			}
			
		}
	}
	# write the results to file
	prepareOutput(jsn, fin)	
}



prepareOutput <- function(jsn, fin)
{
	msg = sprintf(' R Function: prepareOutput')
	print(msg)
	jsn_out <- jsn
	outfile <- paste(did, "_Hap_Geno_out.json", sep='') 
	
	for ( hap in names(jsn) )
	{
		for (ethn in names(jsn[[hap]]) )
		{
			or <- c(); freq <- c(); gent <- c()
			for (gt in names(jsn[[hap]][[ethn]]) )
			{
				gtinfo <- jsn[[hap]][[ethn]][[gt]] 
				# find the reference Genotype. Genotype with OR=1 has to be first
				f      <- fin[fin$OR==gtinfo$OR & fin$Ethn==ethn & fin$Freq==gtinfo$FControl & fin$GT==gt & fin$Hap==hap,]
				newor  <- fin$NewOR[fin$OR==gtinfo$OR & fin$Ethn==ethn & fin$Freq==gtinfo$FControl & fin$GT==gt & fin$Hap==hap]
				#print (gt)
				#print(f)
				#print(newor)
				gtinfo$NewOR <- newor
				jsn_out[[hap]][[ethn]][[gt]] <- gtinfo
			}
			#NA case:
			gtinfo$OR       <-  1
			gtinfo$NewOR    <- 1
			gtinfo$Freq     <- 0
			gtinfo$CIissue  <- ""
			gtinfo$FCase    <- 0
			gtinfo$FControl <- 0
			jsn_out[[hap]][[ethn]][["NA"]] <- gtinfo
			
		}
	}
	write(toJSON(jsn_out), outfile)
}

prevDisease <- function(did)
{
	msg = sprintf(' R Function: prevDisease -- disease: %s',did)
	print(msg)
					
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
	return(as.numeric(PrevD[[did]]))
}

##### Main of R script
cmd_args = commandArgs()

did = cmd_args[6] # disease name

sprintf('*** I am calculating adjusted OR for haplotypic genotype of disease %s',did)
prepareInputs(did)
