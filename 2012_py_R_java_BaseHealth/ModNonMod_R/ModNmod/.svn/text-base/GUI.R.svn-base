setwd("/home/pouria3/workspace/R/ModNmod")
source("ModNmod.R")
#source("Snp2.R")

# DB stuff; 
con     <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
qur1    <- ('select IDName, Unit from IDmap where (upper(isRiskFactorID)=\'YES\' and upper(IsModifiable)=\'YES\') ')
qur2    <- ('select IDNo, IDName, Unit, IDParent from IDmap where upper(isRiskFactorID)=\'YES\' ')

mods    <- dbGetQuery(con, qur1)
allfacs <- dbGetQuery(con, qur2)
dbDisconnect(con)


#Gather Data
limits <- Data

# Lay out the GUI
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
Modifiable <- c()
NonModifiable <- c()
for (ID in unique(Data$RiskFactorID))
{
	fac    <- allfacs$IDName[allfacs$IDNo==ID]	
	if (ID != 11)
	{
			
		if (fac%in%mods$IDName )
		{
			Modifiable <- c(Modifiable, fac)
		}else 
		{
			NonModifiable <- c(NonModifiable, fac)
		}
		varnm  <- paste(fac)
		slval  <- paste("sl_", fac)
		normal <- limits$normal[limits$name==fac][1]
		riskid <- limits$RiskFactorID[limits$name==fac][1]
		
		if ( is.na(normal) )
		{
			normal <- Data$RiskFactorValue[Data$FactorValue== min(Data$FactorValue[Data$RiskFactorID==riskid])  & Data$RiskFactorID==riskid ][1]
		}
		
		print (fac)
		print(normal)
		assign(slval, tclVar(as.character(normal)) )
		assign(varnm,  tklabel(tt,text=as.character(tclvalue(get(slval))) ) )
		unit_    <- allfacs$Unit[allfacs$IDNo==ID]
		if (is.na(unit_) )
		{
			unit_ <- unique(as.character(Data$comment[Data$RiskFactorID==riskid ] ) )
			if ( length(unit_) > 1 )
			{
				u <- c()
				for (un in unit_)
				{
					u <- paste(u, un)
				}
				unit_ <- u
			}
		}
		print(unit_)
		tkgrid(get(varnm), tklabel(tt,text=fac), tklabel(tt, text=unit_) )
		tkconfigure(get(varnm),textvariable=get(slval))
		if (fac == "BMI")
		{
			slider <- tkscale(tt, to=limits$Max[limits$name==fac][1], from=limits$Min[limits$name==fac][1],
					showvalue=F, variable=get(slval),
					resolution=0.1, orient="horizontal")
		}
		else
		{
			slider <- tkscale(tt, to=limits$Max[limits$name==fac][1], from=limits$Min[limits$name==fac][1],
			showvalue=F, variable=get(slval),
			resolution=1, orient="horizontal")
		}
		tkgrid(slider)
	}
}

GE               <- c("Gender", "Ethnicity")
NonModifiable    <- NonModifiable[!NonModifiable%in%c("Age", "Gender", "Ethnicity", "intercept")]



#----------------------Functions ------------------#

getgender <- function()
{
	if (did == 103)
	{
		return("F")
	}else if (did == 105 | did==110)
	{
		return("M")
	}else
	{
		var        <- as.numeric(tclvalue(get(paste("sl_", 'Gender'))) )
		if (var==0)
			return("M")
		else
			return("F")
	}
	
}

getrisk   <- function(fac)
{
	not_ur_gender <- 0 # risk-gender flag
	gender = getgender()
	
	var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
	riskid     <- Data$RiskFactorID[Data$name==fac][1]
	genders    <- as.character(unique(Data$Gender[Data$RiskFactorID==riskid]))
	if (length(genders) == 1 ) 
	{
		if (genders == "MF")
		{
			gender <- "MF" #risk-factor applies to both genders
		}else if (genders != gender)
		{
			not_ur_gender <- 1 # doesnt apply to ur gender
		}
	}
	
	if (not_ur_gender == 1)
	{
		riskvalue <- 0
	}else
	{
		riskvalue    <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var & gender==Data$Gender]
	}
	
	return(riskvalue)

}

getminrisk   <- function(fac)
{
	not_ur_gender <- 0
	gender = getgender()
	
	var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
	riskid     <- Data$RiskFactorID[Data$name==fac][1]
	genders    <- as.character(unique(Data$Gender[Data$RiskFactorID==riskid]))
	if (length(genders) == 1 ) 
	{
		if (genders == "MF")
		{
			gender <- "MF" #risk-factor applies to both genders
		}else if (genders != gender)
		{
			not_ur_gender <- 1 # doesnt apply to ur gender
		}
	}
	
	if (not_ur_gender == 1)
	{
		minrisk <- 0
	}else
	{
		minrisk      <- min(Data$FactorValue[Data$RiskFactorID==riskid & gender==Data$Gender])
	}
	

	return(minrisk)	

}

getnormrisk   <- function(fac)
{
	not_ur_gender <- 0
	gender <- getgender()
	riskid     <- Data$RiskFactorID[Data$name==fac][1]
	genders    <- as.character(unique(Data$Gender[Data$RiskFactorID==riskid]))
	if (length(genders) == 1 ) 
	{
		if (genders == "MF")
		{
			gender <- "MF" #risk-factor applies to both genders
		}else if (genders != gender)
		{
			not_ur_gender <- 1 # doesnt apply to ur gender
		}
	}
	
	if (not_ur_gender == 1)
	{
		normalrisk <- 0
	}else
	{
		isgen  <- gender==Data$Gender
		normalvalue  <- Data$normal[Data$name==fac & isgen][1]
		if ( is.na(normalvalue) )
		{
			normalrisk     <- min(Data$FactorValue[Data$RiskFactorID==riskid & isgen])
		}else
		{
			normalrisk     <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==normalvalue]
		}
	}
	
	
	return(normalrisk)
	
}

# Run botton
# button function call
PressedRun <- function()
{
	par(lty=2, lwd=2, mfrow=c(1, 2))
	
	modrisk <- 0; minmodrisk <- 0; normmodrisk <- 0; mods <- c()
	for (fac in Modifiable )
	{
		modrisk     <- modrisk     + getrisk(fac)
		minmodrisk  <- minmodrisk  + getminrisk(fac)
		normmodrisk <- normmodrisk + getnormrisk(fac)
		
		mods        <- rbind(mods, data.frame(RiskName=fac, RiskVal=getrisk(fac)-getnormrisk(fac), CurrentVal = getrisk(fac), NormRisk=getnormrisk(fac), MinRisk=getminrisk(fac) ) )
	}
	
	# Non Modifiables, we are comparing the person to her/his self so "Gender", "Age", "Ethnicity", should match
	nmodrisk <- 0; nmods <- c();
	for (fac in NonModifiable )
	{
		nmodrisk <- nmodrisk + getrisk(fac)	
		
		nmods    <- rbind(nmods, data.frame(RiskName=fac, RiskVal=getrisk(fac) ) )
	}
	

	
	lifetimerisk <- round(100*exp(abs( nmodrisk + modrisk    - normmodrisk )  ))/100
	achievable   <- round(100*exp(abs( nmodrisk + minmodrisk - normmodrisk )  ))/100
	
	# gui setting for protective vs causative
	if (( nmodrisk + modrisk - normmodrisk) < 0 )
	{
		main = paste("LifeTime Risk : ", lifetimerisk, "Fold", "-", "\n")
	}else
	{
		main = paste("LifeTime Risk : ", lifetimerisk, "Fold", "+", "\n")
	}
	if ( (nmodrisk + minmodrisk - normmodrisk)  < 0 )
	{
		main= paste(main, "Achievable Risk : ", achievable, "Fold", "-", "\n",sep=" ")
	}else
	{
		main= paste(main, "Achievable Risk : ", achievable, "Fold", "+", "\n",sep=" ")
	}
	
	
	blbls  <- as.character(nmods$RiskName)
	bp1    <- barplot(nmods$RiskVal, names.arg=blbls, col=rainbow(length(blbls)),xlim=c(log(1/4), log(16)),
		,main=main, 
		width=0.5, horiz=TRUE, beside=TRUE, axes=FALSE)
	axis(1, at=seq(log(1/4),log(16),by=log(4)), labels=c("4", "base","4","16"))
	mtext(round(100*exp(abs(nmods$RiskVal)))/100, side=4, at=bp1,line=0.5,  cex=1.5 )
	
	
	blbls <- Modifiable
	bp <- barplot(mods$RiskVal, names.arg=blbls, col=rainbow(length(blbls)), xlim=c(log(1/4), log(16)), width=0.5, horiz=TRUE,  beside=TRUE, legend=TRUE, axes=FALSE)
	axis(1, at=seq(log(1/4),log(16),by=log(4)), labels=c("4", "base","4","16"))
	mtext(round(100*exp(abs(mods$RiskVal)))/100, side=4, at=bp,line=0.5,  cex=1.5 )
	
	
	
}

Pressedlife <- function()
{
		modrisk <- 0; minmodrisk <- 0; normmodrisk <- 0; mods <- c()
	for (fac in Modifiable )
	{
		modrisk     <- modrisk     + getrisk(fac)
		minmodrisk  <- minmodrisk  + getminrisk(fac)
		normmodrisk <- normmodrisk + getnormrisk(fac)
		
		mods        <- rbind(mods, data.frame(RiskName=fac, RiskVal=getrisk(fac)-getnormrisk(fac), CurrentVal = getrisk(fac), NormRisk=getnormrisk(fac), MinRisk=getminrisk(fac) ) )
	}
	
	# Non Modifiables, we are comparing the person to her/his self so "Gender", "Age", "Ethnicity", should match
	nmodrisk <- 0; nmods <- c();
	for (fac in NonModifiable )
	{
		nmodrisk <- nmodrisk + getrisk(fac)	
		
		nmods    <- rbind(nmods, data.frame(RiskName=fac, RiskVal=getrisk(fac) ) )
	}

	riskid_age   <- Data$RiskFactorID[limits$name=="Age"][1] 
	max_age      <- limits$Max[limits$name=="Age"][1]
	max_agerisk <- Data$FactorValue[Data$RiskFactorID==riskid_age & Data$RiskFactorValue==max_age]
	
	lifetime   <- modrisk + nmodrisk - normmodrisk 
	achievable <- nmodrisk + minmodrisk - normmodrisk 
	
	paf <- function(x)
	{
		inc = exp(Data$FactorValue[Data$RiskFactorID==11]) / (1 + exp(Data$FactorValue[Data$RiskFactorID==11]) ) #incidence 
		vlow = (1-inc)/4
		low  = vlow + 1*(1-inc)/4
		med  = vlow + 2*(1-inc)/4
		high = vlow + 3*(1-inc)/4
		vhigh= 0.9
		x2 = x
		if (did == 111 )
		{
			intervals <- c(0, 0.03, 0.25, 0.5, 0.7, 0.8, 0.99 )
		}else 
		{
			intervals <- c(0, 0.01, 0.03, 0.05, 0.1, 0.4, 0.99 )
		}
		
		for (j in c(1:length(x)) )
		{
			for ( i in seq(1, length(intervals)-1 ) )
			{
				if ( intervals[i]<=x[j] & x[j]<intervals[i+1] )
				{
					x2[j] = 1/(intervals[i+1] - intervals[i]) * (x[j]-intervals[i]) + (i-1)
				}
				if (x[j]>=0.99)
				{
					#x2[j] = intervals[length(intervals)]
					x2[j] = 5
				}
			}
		}
		return(x2)
	}
		
	modrisk <- 0; minmodrisk <- 0; 
	for (fac in Modifiable )
	{
		modrisk     <- modrisk     + getrisk(fac)
		minmodrisk  <- minmodrisk  + getminrisk(fac)
	}
	
	nmodrisk <- 0
	for (fac in NonModifiable )
	{
		nmodrisk <- nmodrisk + getrisk(fac)	
	}

	
	if (did %in% c(103, 105, 110))
	{
		GE = c("Ethnicity")  # (gender doesn't exist for breascancer or hairloss or prostae cancer
		
	}else
	{
		GE = c("Ethnicity", "Gender") # for others gender contributes to risk
	}
	
	ge_risk <- 0
	for (fac in GE )
	{
		ge_risk    <- ge_risk    + getrisk(fac)
	}
	
	intercept <- Data$FactorValue[Data$RiskFactorID==11]		
	
	riskid_age   <- Data$RiskFactorID[limits$name=="Age"][1] ; x<- c() ; y <- c() ; n <- c(); newy<-c(); ach<- c()
	for ( age in seq(Data$Min[limits$name=="Age"][1], limits$Max[limits$name=="Age"][1], by=5 ) )
	{
		age_risk <- Data$FactorValue[Data$RiskFactorID==riskid_age & Data$RiskFactorValue==age]
		x        <- c(x, age)
		y        <- c(y, nmodrisk + modrisk    + intercept + ge_risk + age_risk  )
		n        <- c(n, nmodrisk + minmodrisk + intercept + ge_risk + age_risk  )
		
		newy     <- c(newy, lifetime + age_risk - max_agerisk)
		ach      <- c(ach, achievable + age_risk - max_agerisk)
	}
	
	# protective effect?
	for (i in c(1:length(newy)) )
	{
		
		if (newy[i] < 0 )
		{
			newy[i] <- -exp(abs(newy[i]) )
		}else
		{
			newy[i] <- exp(abs(newy[i]) )
		}
		if (ach[i] < 0 )
		{
			ach[i] <- -exp(abs(ach[i]) )
		}else
		{
			ach[i] <- exp(abs(ach[i]) )
		}
		
	}
			
	prob  <- exp(y)/(exp(y) + 1) 
	nprob <- exp(n)/(exp(n) + 1) 
	#dev.new()
	par(lwd=2)
	plot  (x, paf(prob), type="l", col="blue", xlab="Age", ylab="Risk ", ylim=c(0, 7), axes=FALSE, xlim=c(limits$Min[limits$name=="Age"][1], limits$Max[limits$name=="Age"][1]) )
	axis(2, at=c(0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5 ), labels=c("", "VL", "L", "M", "H", "VH", "XH", ""))
	axis(1, at=x )
	par(lwd=1)
	grid(NULL, NULL, col="black")
	par(lty=2, lwd=2)
	points(x, paf(nprob), type="l", col="red")
	abline(v=as.numeric(tclvalue(get(paste("sl_", "Age"))) ),col=3,lty=3)
	
	#dev.new()
	par(lwd=2)
	plot(x,newy, type="l", col='blue', xlab="Age", ylab="Risk", axes=FALSE, ylim=c(-20,30) )
	axis(2, at=c(-20, -10, 0, 10, 20, 30, 40, 50) )
	axis(1, at=x )
	axis(4, at=c(-20, -10, 0, 10, 20, 30, 40, 50 ), labels=c("WOW!", "WOW!", "0", "L", "M", "H", "VH", "XH"))
	abline(v=as.numeric(tclvalue(get(paste("sl_", "Age"))) ),col=3,lty=3)
	par(lty=2, lwd=2)
	points(x, ach, type="l", col="red")
	par(lwd=1)
	grid(NULL, NULL, col="black")
	

}

tt <- frameLower
run.but <- tkbutton(tt,text="Run",command=PressedRun)
#vis.but <- tkbutton(tt,text="Risks exposed!", command=PressedVis)
life.but<- tkbutton(tt,text="Life Time Risk", command=Pressedlife)

tkgrid(run.but)
#tkgrid(vis.but)
tkgrid(life.but)

					

	