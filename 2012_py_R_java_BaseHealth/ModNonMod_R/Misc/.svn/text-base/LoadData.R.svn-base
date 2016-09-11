Main <- function()
{
	### Load Demographic
	DEMO <- LoadDemoData("demo", c("RIDAGEYR", "RIAGENDR", "RIDRETH1"), c("Age", "Gender", "Race") )
	
	### Load BMI-Waist
	BMX <- LoadDemoData("bmx", c("BMXBMI", "BMXWAIST"), c("BMIC", "Waist") )
	
	### Load Diabetes
	T2D <- LoadDiab()
	
	### Load Cholestrol (Total,HDL,Total/HDL)
	TCH <- LoadChol("l13", "TCH")
	
	
	### Load Triglycride
	TRIG <- LoadTrig("l13am", c("LBXTR"), c("TRIG"))
	
}


LoadDemoData <- function(nam, headers, returns)
{
		lout <- list()
	cntflag <- 0
	for (file in c("_d.csv", "_c.csv", "_b.csv", "_a.csv"))
	{
		csv <- read.csv (paste(nam, file, sep=""), header=T); # read the csv file in and assign it to csv variable
		tmp <- unlist(strsplit(file, NULL))
		
		# marking what number comes from what source
		if (cntflag == 0)
		{
			Data=data.frame(list(SEQN=csv$SEQN,Wave=rep(tmp[2],nrow(csv))));
		}
		else
		{
			Data <- rbind(Data, data.frame(list(SEQN=csv$SEQN,Wave=rep(tmp[2],nrow(csv)))))
		}
		# combine data from variouse sources
		for (i in 1:length(returns))
		{
			if (cntflag != 0)
			{
				assign(returns[i], rbind(get(returns[i]), cbind(csv$SEQN, csv[headers[i]]) ) )
			}
			else
			{
				assign(returns[i], cbind(csv$SEQN, csv[headers[i]]))
			}
		}
		cntflag <-1
	}
	
	#return some stuff
	for (j in 1:length(returns))
	{
		lout[[returns[j] ]] <- get(returns[j])
	} 
	lout[["Data"]] <- Data
	return(lout)
}


LoadBMIWaist <- function(nam, headers, returns)
{
	lout <- list()
	BMX   <- (read.csv (paste(nam, "_d.csv", sep=""), header=T)); # read the first file to initialize
	for (i in 1:length(returns))
	{
		assign(returns[i], cbind(BMX$SEQN, as.numeric(as.character(BMX[headers[i]]))) ) 
	}
	

	for (file in c("_c.csv", "_b.csv", "_a.csv"))
	{
		for (i in 1:length(returns))
		{
			assign(returns[i], rbind(get(returns[i]), cbind(BMX$SEQN, as.numeric(as.character(BMX[headers[i]]))) ) ) 
		}
		
	}
	for (j in 1:length(returns))
	{
		lout[[returns[j] ]] <- get(returns[j])
	} 
	return(lout)
}


LoadDiab <- function()
{
	DIQ<-read.csv ("diq_d.csv", header=T);
	T2D=cbind(DIQ$SEQN,DIQ$DIQ010);
	DIQ <- read.csv ("diq_c.csv", header=T);
	T2D=rbind(T2D,cbind(DIQ$SEQN,DIQ$DIQ010));
	DIQ <- read.csv ("diq_b.csv", header=T);
	T2D=rbind(T2D,cbind(DIQ$SEQN,as.numeric(as.character(DIQ$DIQ010))));
	DIQ <- read.csv ("diq_a.csv", header=T);
	T2D=rbind(T2D,cbind(DIQ$SEQN,as.numeric(as.character(DIQ$DIQ010))));
	
	lout <- list("T2D"= T2D)
	return(lout)
	
}


LoadChol <- function(nam, returns)
{
	### Load Cholestrol (Total,HDL,Total/HDL)
	TCHOL<-read.csv ("tchol_d.csv", header=T);
	TCHOL2<-read.csv ("hdl_d.csv", header=T);
	TCH=cbind(TCHOL$SEQN,as.numeric(as.character(TCHOL$LBXTC)), as.numeric(as.character(TCHOL2$LBDHDD)), as.numeric(as.character(TCHOL$LBXTC))/as.numeric(as.character(TCHOL2$LBDHDD)));
	TCHOL<-read.csv ("l13_c.csv", header=T);
	TCH=rbind(TCH,cbind(TCHOL$SEQN,as.numeric(as.character(TCHOL$LBXTC)), as.numeric(as.character(TCHOL$LBXHDD)), as.numeric(as.character(TCHOL$LBXTC))/as.numeric(as.character(TCHOL$LBXHDD))));
	TCHOL<-read.csv ("l13_b.csv", header=T);
	TCH=rbind(TCH,cbind(TCHOL$SEQN,as.numeric(as.character(TCHOL$LBXTC)), as.numeric(as.character(TCHOL$LBDHDL)), as.numeric(as.character(TCHOL$LBXTC))/as.numeric(as.character(TCHOL$LBDHDL))));
	TCHOL<-read.csv ("lab13_a.csv", header=T);
	TCH=rbind(TCH,cbind(TCHOL$SEQN,as.numeric(as.character(TCHOL$LBXTC)), as.numeric(as.character(TCHOL$LBDHDL)), as.numeric(as.character(TCHOL$LBXTC))/as.numeric(as.character(TCHOL$LBDHDL))));
	lout <- list("TCH"= TCH)
	return(lout)	

}



LoadTrig <- function(nam, headers, returns)
{
	TRIGLY<-read.csv ("trigly_d.csv", header=T);
	TRIG=cbind(TRIGLY$SEQN,as.numeric(as.character(TRIGLY$LBXTR)));
	for (file in c("_c.csv", "_b.csv"))
	{
		TRIGLY<-read.csv (paste(nam, file, sep=""), header=T);
		TRIG=rbind(TRIG,cbind(TRIGLY$SEQN,as.numeric(as.character(TRIGLY$LBXTR))));
		
	}
	TRIGLY<-read.csv ("lab13am_a.csv", header=T);
	TRIG=rbind(TRIG,cbind(TRIGLY$SEQN,as.numeric(as.character(TRIGLY$LBXTR))));
	
	lout <- list("TRIG"= TRIG)
	return(lout)	
}

