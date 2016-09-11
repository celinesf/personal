rm(list=ls())
#library(lattice)
setwd("/home/pouria3/workspace/GenStat/trunk/Alzh")
source("riskmodel.R")
source("Snp.R")
source("Hap.R")

# Min & Max for every RiskFactor
RiskFactorMinMax <- function()
{
	
	HAP          <- data.frame(RiskFactorID=115, name="Apoe"        , Min=0    , Max=6    , normal=0)
	
	TNF          <- data.frame(RiskFactorID=114, name="_TNF_rs4647198"   , Min=0    , Max=3    , normal=0)
	CST3         <- data.frame(RiskFactorID=113, name="_CST3_rs1064039"   , Min=0    , Max=3    , normal=0)
	
	Family       <- data.frame(RiskFactorID=111, name="Family"      , Min=50   , Max=100  , normal=NA)
	DownSynd     <- data.frame(RiskFactorID=110, name="DownSyndrome", Min=0    , Max=1    , normal=NA)
	Ethnicity    <- data.frame(RiskFactorID=109, name="Ethnicity"   , Min=0    , Max=6    , normal=NA)
	Age          <- data.frame(RiskFactorID=108, name="Age"         , Min=50   , Max=100  , normal=NA)
	GNDR         <- data.frame(RiskFactorID=107, name="Gender"      , Min=0    , Max=1    , normal=NA)
	SBP          <- data.frame(RiskFactorID=106, name="SBP"         , Min=120  , Max=180  , normal=NA)
	T2D          <- data.frame(RiskFactorID=105, name="T2D"         , Min=0    , Max=1    , normal=NA)
	Chol         <- data.frame(RiskFactorID=104, name="TChol"       , Min=198  , Max=500  , normal=NA)
	HeartDisease <- data.frame(RiskFactorID=103, name="HeartDisease", Min=0    , Max=1    , normal=NA)
	Stroke       <- data.frame(RiskFactorID=102, name="Stroke"      , Min=0    , Max=1    , normal=NA)
	TSH          <- data.frame(RiskFactorID=101, name="TSH"         , Min=0    , Max=10   , normal=NA)  #female only
	Thcy         <- data.frame(RiskFactorID=100, name="Thcy"        , Min=5    , Max=20   , normal=NA)
	
	FacLimits <- rbind(SBP, T2D, Chol, HeartDisease, Stroke, TSH, Thcy, GNDR, Age, Ethnicity, DownSynd, Family, CST3, TNF, HAP )
	return(FacLimits)
	
}





#Gather Data
Data <- rbind(fitThcy(), fitTSH(), fitStroke(), fitHearDisease(), 
	fitChol(), fitSBP(), fitT2D(), fitGender(), fitAge(), 
	fitEthnicity(), fitDownSynd(), fitInheritance(), fitIntercept())

#Gather Genes
Data <- rbind(Data, snip)

limits <- RiskFactorMinMax()






# Lay out the GUI
units <- c("mm/Hg", "", "mg/dl", "", "", "mU/L", "umol/L", "", "years", "", "", "", "NA GG AG AA", "NA TT CT CC", "NA e3/e3 e2/e2 e2/e3 e2/e4 e3/e4 e4/e4")
require(tcltk)
tt <- tktoplevel()
frameOverall <- tkframe(tt)
frameUpper <- tkframe(frameOverall,relief="groove",borderwidth=2)
tkgrid(tklabel(frameUpper,text="RiskFactors"))
frameLower <- tkframe(frameOverall,relief="groove",borderwidth=2)
tkgrid(tklabel(frameLower,text="Action"))

tkgrid(frameUpper)
tkgrid(tklabel(frameOverall,text="Hit run when you are done with sliders!"))
tkgrid(frameLower)
tkgrid(frameOverall)

tt <- frameUpper
counter=1
for (fac in c("SBP", "T2D", "TChol", "HeartDisease", "Stroke", "TSH", "Thcy", "Gender", "Age", "Ethnicity", "DownSyndrome", "Family", "_TNF_rs4647198", "_CST3_rs1064039", "Apoe") )
{
	varnm <- paste(fac)
	slval <- paste("sl_", fac)
	assign(slval, tclVar("0") )
	assign(varnm,  tklabel(tt,text=as.character(tclvalue(get(slval))) ) )
	tkgrid(get(varnm), tklabel(tt,text=fac), tklabel(tt, text=units[counter]) )
	tkconfigure(get(varnm),textvariable=get(slval))
	if (fac=="TSH")
	{
		slider <- tkscale(tt, to=limits$Max[limits$name==fac], from=limits$Min[limits$name==fac],
		showvalue=F, variable=get(slval),
		resolution=0.01, orient="horizontal")
	}else
	{
		slider <- tkscale(tt, to=limits$Max[limits$name==fac], from=limits$Min[limits$name==fac],
		showvalue=F, variable=get(slval),
		resolution=1, orient="horizontal")
	}
	tkgrid(slider)
	counter <- counter + 1
}
# Run botton
# button function call
PressedRun <- function()
{
	par(lty=2, lwd=2)
	asses <- data.frame(Name=c(), RiskFactorID = c(), FactorValue = c(), DR=c())
	for (fac in c("SBP", "TChol", "TSH", "Thcy") )
	{
		limits     <- RiskFactorMinMax()
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac]
		normal     <- limits$normal[limits$name==fac]
		riskvalue  <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
		tmp <- data.frame(Name=fac, RiskFactorID = riskid, FactorValue = riskvalue, DR = riskvalue - min(Data$FactorValue[Data$RiskFactorID==riskid]) )
		asses <- rbind(asses, tmp)
	}
		
	# Ready for charts
	slices <- asses$DR
	lbls   <- c("SBP", "TChol", "TSH", "Thcy")
	blbls  <- lbls
	pct <- round(slices/sum(slices)*100)
	lbls <- paste(lbls, pct) # add percents to labels 
	lbls <- paste(lbls,"%",sep="") # ad % to labels 
	par(mfrow=c(2, 2))
	
	if (any(slices) )
	{
		pie(slices,labels = lbls, col=rainbow(length(lbls)),
		main="Modifiable Risk Factors")
	}
	

	
	barplot(slices, names.arg=blbls, col=rainbow(length(lbls)),width=0.5, horiz=TRUE, beside=TRUE, legend=TRUE, axes=FALSE)
	axis(1, at=seq(log(1),log(16),by=log(2)), labels=c("base","2","4","8","16"))
		
	
	# Non Modifiables, we are comparing the person to her/his self so "Gender", "Age", "Ethnicity", should match
	nonmod <- data.frame(Name=c(), RiskFactorID = c(), FactorValue = c(), DR=c())
	for (fac in c("T2D", "HeartDisease", "Stroke", "DownSyndrome", "Family") )
	{
		limits     <- RiskFactorMinMax()
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac]
		normal     <- limits$normal[limits$name==fac]
		riskvalue  <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
		tmp <- data.frame(Name=fac, RiskFactorID = riskid, FactorValue = riskvalue, DR = riskvalue - min(Data$FactorValue[Data$RiskFactorID==riskid]) )
		nonmod <- rbind(nonmod, tmp)
		fh_risk <- riskvalue
	}
	genetic_risk=0
	for (fac in c("_TNF_rs4647198", "_CST3_rs1064039") )
	{
		limits     <- RiskFactorMinMax()
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac]
		normal     <- limits$normal[limits$name==fac]
		riskvalue  <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
		genetic_risk<- genetic_risk + riskvalue
		tmp <- data.frame(Name=fac, RiskFactorID = riskid, FactorValue = riskvalue, DR = riskvalue - normal )
		nonmod <- rbind(nonmod, tmp)
	}
	for (fac in c("Apoe") )
	{
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		eth        <- as.numeric(tclvalue(get(paste("sl_", "Ethnicity"))) )
		riskid     <- limits$RiskFactorID[limits$name==fac]
		normal     <- limits$normal[limits$name==fac]
		riskvalue  <- haplotype$FactorValue[haplotype$RiskFactorID==riskid & haplotype$RiskFactorValue==var & haplotype$ethnicity==eth]
		genetic_risk <- genetic_risk + riskvalue
		tmp <- data.frame(Name=fac, RiskFactorID = riskid, FactorValue = riskvalue, DR = riskvalue - normal )
		nonmod <- rbind(nonmod, tmp)
	}
		
		
	slices <- nonmod$DR
	slices[slices<0]=0
	lbls   <- c("T2D", "HeartD", "Str", "Dsynd", "FH", "TNF", "CST3", "Apoe")
	blbls  <- lbls
	pct <- round(slices/sum(slices)*100)
	lbls <- paste(lbls, pct) # add percents to labels 
	lbls <- paste(lbls,"%",sep="") # ad % to labels 
	
	if (any(slices) )
	{
		pie(slices,labels = lbls, col=rainbow(length(lbls)),
		main="NonModifiable Risk Factors")
	}
	
	#now, "Gender", "Age", "Ethnicity",  might add to ur risk, so :
	GAE <- data.frame(Name=c(), RiskFactorID = c(), FactorValue = c(), DR=c())
	for (fac in c("Gender", "Age", "Ethnicity") )
	{
		limits     <- RiskFactorMinMax()
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac]
		normal     <- limits$normal[limits$name==fac]
		riskvalue  <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
		tmp <- data.frame(Name=fac, RiskFactorID = riskid, FactorValue = riskvalue, DR = riskvalue - min(Data$FactorValue[Data$RiskFactorID==riskid]) )
		GAE <- rbind(GAE, tmp)
	}
	# Age, Gender and Ethnicity shouldn't increase the lifetime risk , sum(GAE$FactorValue)  
		
	# from Genetic Factors and Family history, pick the max
	if (fh_risk > 0 )
	{
		fh_gen    <- max(genetic_risk, fh_risk)
	}else
	{
		fh_gen    <- genetic_risk
	}
	
	otherrisks <- 0
	for (fac in c("T2D", "HeartDisease", "Stroke", "DownSyndrome") )
	{
		limits     <- RiskFactorMinMax()
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac]
		normal     <- limits$normal[limits$name==fac]
		riskvalue  <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
		otherrisks <- otherrisks + riskvalue
	}
		
	lifetimerisk <- round(100*exp(sum(asses$FactorValue) + fh_gen + otherrisks ))/100
	achievable   <- round(100*exp(sum(fh_gen + otherrisks) ))/100
	
	barplot(nonmod$DR, names.arg=blbls, col=rainbow(length(lbls)),xlim=c(log(1/4), log(16)),
	,main=paste("LifeTime Risk : ", lifetimerisk, "Fold", "\n", "Achievable Risk : ", achievable, "Fold", "\n", "You Could reduce your risk",
	round(100*lifetimerisk/achievable)/100, "times!", sep=" "), 
	width=0.5, horiz=TRUE, beside=TRUE, legend.text=blbls, axes=FALSE)
	axis(1, at=seq(log(1/4),log(16),by=log(4)), labels=c("4", "base","4","16"))
	
					
					
	
}

PressedVis <- function()
{
	factors <- c("SBP", "T2D", "TChol", "HeartDisease", "Stroke", "TSH", "Thcy", "DownSyndrome", "Family", "_TNF_rs4647198", "_CST3_rs1064039")
	mt <- matrix(0, nrow=length(factors), ncol=1,  dimnames=list(factors, c("Risk")) )
	cnt <- 1
	for ( fac in factors)
	{
		limits     <- RiskFactorMinMax()
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac]
		normal     <- limits$normal[limits$name==fac]
		riskvalue  <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
		if (riskvalue==0)
		{
			mt[cnt, 1] <- 0
		}
		else
		{
			mt[cnt, 1] <- round(100*exp(riskvalue))/100
		}
		cnt        <- cnt + 1
	}	
	dev.new()
	par(xpd=T, mar=c(5, 4, 4, 2) + c(0, 0, 0, 14) )
	bp <- barplot(mt, axes=FALSE, col=rainbow(length(factors)), legend=FALSE, ylab="Fold", space=1)
	legend(2, 4, paste(factors, ":", as.numeric(mt[,1]) ), fill=rainbow(length(factors)), cex=0.8 )
	#text(bp, as.numeric(mt[,1]), labels=as.numeric(mt[,1]), col="black" )
	
}

Pressedlife <- function()
{
	paf <- function(x)
	{
		x2 = x
		intervals = c(0, 0.04, 0.2, 0.4, 0.6, 0.8, 1 )
		for (j in c(1:length(x)) )
		{
			for ( i in seq(1, length(intervals)-1 ) )
			{
				if ( intervals[i]<x[j] & x[j]<intervals[i+1] )
				{
					x2[j] = 1/(intervals[i+1] - intervals[i]) * (x[j]-intervals[i]) + (i-1)
				}
				#if (x[j]>0.99)
				#{
			#		x2[j] = 5
		#		}
			}
		}
		return(x2)
	}
		
	
	x=c(); y=c(); n=c();
	factors <- c("SBP", "T2D", "TChol", "HeartDisease", "Stroke", "TSH", "Thcy", "Gender", "Ethnicity", "DownSyndrome") #notice the Age is not listed
	asses <- data.frame(Name=c(), RiskFactorID = c(), FactorValue = c())
	for ( fac in factors)
	{
		limits     <- RiskFactorMinMax()
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac]
		normal     <- limits$normal[limits$name==fac]
		riskvalue  <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
		tmp        <- data.frame(Name=fac, RiskFactorID = riskid, FactorValue = riskvalue ) 
		asses      <- rbind(asses, tmp)
		fh_risk <- riskvalue
		
	}
	
	limits     <- RiskFactorMinMax()
	var        <- as.numeric(tclvalue(get(paste("sl_", "Family"))) )
	riskid     <- limits$RiskFactorID[limits$name=="Family"]
	normal     <- limits$normal[limits$name=="Family"]
	riskvalue  <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
	fh_risk    <- riskvalue
	
	
	healthy <- data.frame(Name=c(), RiskFactorID = c(), FactorValue = c(), DR=c())
	for (fac in c("SBP", "TChol", "TSH", "Thcy") )
	{
		limits     <- RiskFactorMinMax()
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac]
		normal     <- limits$normal[limits$name==fac]
		riskvalue  <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
		tmp        <- data.frame(Name=fac, RiskFactorID = riskid, FactorValue = min(Data$FactorValue[Data$RiskFactorID==riskid]) )
		healthy    <- rbind(healthy, tmp)
	}
	for (fac in c("Gender", "Ethnicity", "HeartDisease", "Stroke", "DownSyndrome", "T2D") )
	{
		limits     <- RiskFactorMinMax()
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac]
		normal     <- limits$normal[limits$name==fac]
		riskvalue  <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
		tmp        <- data.frame(Name=fac, RiskFactorID = riskid, FactorValue = riskvalue )
		healthy    <- rbind(healthy, tmp)
	}
	
	genetic_risk=0
	for (fac in c("_TNF_rs4647198", "_CST3_rs1064039") )
	{
		limits     <- RiskFactorMinMax()
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac]
		normal     <- limits$normal[limits$name==fac]
		riskvalue  <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
		genetic_risk<- genetic_risk + riskvalue
		#tmp <- data.frame(Name=fac, RiskFactorID = riskid, FactorValue = riskvalue, DR = riskvalue - normal )
		#nonmod <- rbind(nonmod, tmp)
	}
		
	for (fac in c("Apoe") )
	{
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		eth        <- as.numeric(tclvalue(get(paste("sl_", "Ethnicity"))) )
		riskid     <- limits$RiskFactorID[limits$name==fac]
		normal     <- limits$normal[limits$name==fac]
		riskvalue  <- haplotype$FactorValue[haplotype$RiskFactorID==riskid & haplotype$RiskFactorValue==var & haplotype$ethnicity==eth]
		genetic_risk <- genetic_risk + riskvalue
		#tmp <- data.frame(Name=fac, RiskFactorID = riskid, FactorValue = riskvalue, DR = riskvalue - normal )
		#healthy <- rbind(healthy, tmp)
		#asses      <- rbind(asses, tmp)
	}
		
	if (fh_risk > 0 )
	{
		fh_gen    <- max(genetic_risk, fh_risk)
	}else
	{
		fh_gen    <- genetic_risk
	}
				
		
	healthy_risk  <- sum(healthy$FactorValue) + fh_gen + Data$FactorValue[Data$comment == "Intercept"] 
	current_risk  <- sum(asses$FactorValue)   + fh_gen + Data$FactorValue[Data$comment == "Intercept"] 
	
	riskid_age   <- limits$RiskFactorID[limits$name=="Age"]
	
			
	for ( age in seq(limits$Min[limits$name=="Age"], limits$Max[limits$name=="Age"], by=5 ) )
	{
		x <- c(x, age)
		y <- c(y, current_risk + Data$FactorValue[Data$RiskFactorID==riskid_age & Data$RiskFactorValue==age] )
		n <- c(n, healthy_risk + Data$FactorValue[Data$RiskFactorID==riskid_age & Data$RiskFactorValue==age] )
	}
	
	current_age_risk    = current_risk + Data$FactorValue[Data$RiskFactorID==riskid_age & Data$RiskFactorValue==as.numeric(tclvalue(get(paste("sl_", "Age"))))]
	c_a_r_exp           = exp(current_age_risk)/(exp(current_age_risk) + 1 )
	current_age_healthy = healthy_risk + Data$FactorValue[Data$RiskFactorID==riskid_age & Data$RiskFactorValue==as.numeric(tclvalue(get(paste("sl_", "Age"))))]
	c_a_h_exp           = exp(current_age_healthy)/(exp(current_age_healthy) + 1 )
	reduction_ratio     = round(100*(c_a_r_exp - c_a_h_exp)/c_a_r_exp)  
	
	prob  <- exp(y)/(exp(y) + 1) 
	nprob <- exp(n)/(exp(n) + 1) 
	dev.new()
	par(lwd=2)
	plot  (x, paf(prob), type="l", col="blue", xlab="Age", ylab="Risk ", ylim=c(0, 7), axes=FALSE, xlim=c(50, 100) )
	axis(2, at=c(0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5 ), labels=c("", "VL", "L", "M", "H", "VH", "XH", ""))
	axis(1, at=x )
	par(lwd=1)
	grid(NULL, NULL, col="black")
	par(lty=2, lwd=2)
	points(x, paf(nprob), type="l", col="red")
	abline(v=as.numeric(tclvalue(get(paste("sl_", "Age"))) ),col=3,lty=3)
	



}

tt <- frameLower
run.but <- tkbutton(tt,text="Run",command=PressedRun)
#vis.but <- tkbutton(tt,text="Risks exposed!", command=PressedVis)
life.but<- tkbutton(tt,text="Life Time Risk", command=Pressedlife)

tkgrid(run.but)
#tkgrid(vis.but)
tkgrid(life.but)

					

	