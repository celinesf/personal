setwd("/home/pouria3/workspace/GenStat/trunk/Hyp")
source("ModNmod.R")
source("Snp2.R")

# DB stuff; 
con     <- dbConnect(MySQL(), user='root', password='genophen', host="localhost", dbname="AlgoPhen")
qur1    <- ('select IDName, Unit from IDmap where (upper(isRiskFactorID)=\'YES\' and upper(IsModifiable)=\'YES\') ')
qur2    <- ('select IDNo, IDName, Unit from IDmap where upper(isRiskFactorID)=\'YES\' ')

mods    <- dbGetQuery(con, qur1)
allfacs <- dbGetQuery(con, qur2)
# add the snp portion
unique_snps    <- as.character(unique(snip$RiskFactorID) ) 
for (uni in unique_snps) 
{
	snip_  <- snip[apply(snip, 1, function(x) uni %in% x),] 
	if ( nrow(snip_) > 4 )
	{
		print(snip_)
		snip_          <- snip_[1:4,]
	}
	id             <- as.character(unique(snip_$RiskFactorID) )   
	nm             <- as.character(unique(snip_$name))        
	unit           <- as.character(snip_$comment)
	allfacs        <- rbind(allfacs, data.frame(IDNo=id, IDName=nm, Unit=paste(unit[1], unit[2], unit[3], unit[4], sep=" " )  ))
}

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
	if (ID != 116)
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
GAE           <- c("Age", "Gender", "Ethnicity")
GE            <- c("Gender", "Ethnicity")
FH            <- c("FH")
Genetics      <- as.character(unique(snip$name))
NonModifiable    <- NonModifiable[!NonModifiable%in%c(GAE, FH, Genetics)]


#Genetics <- c()

AllFactors    <- c(Modifiable, NonModifiable, GAE, FH, Genetics)


# Run botton
# button function call
PressedRun <- function()
{
	par(lty=2, lwd=2, mfrow=c(1, 2))
	
	asses <- data.frame(Name=c(), RiskFactorID = c(), FactorValue = c(), DR=c(), minFacVal=c())
	for (fac in Modifiable )
	{
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac][1]
		normal     <- limits$normal[limits$name==fac][1]
		riskvalue  <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
		if ( is.na(normal) )
		{
			normal <- min(Data$FactorValue[Data$RiskFactorID==riskid])
		}else
		{
			normal     <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==normal]
		}
		tmp <- data.frame(Name=fac, RiskFactorID = riskid, FactorValue = riskvalue, DR = riskvalue - normal, minFacVal=min(Data$FactorValue[Data$RiskFactorID==riskid])-normal )
		asses <- rbind(asses, tmp)
	}
	
	# Non Modifiables, we are comparing the person to her/his self so "Gender", "Age", "Ethnicity", should match
	nonmod <- data.frame(Name=c(), RiskFactorID = c(), FactorValue = c(), DR=c())
	for (fac in NonModifiable )
	{
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac][1]
		normal     <- limits$normal[limits$name==fac][1]
		riskvalue  <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
		if (length(riskvalue) > 1 )
		{
			riskvalue <- riskvalue[1]
		}
		if ( is.na(normal) )
		{
			normal <- min(Data$FactorValue[Data$RiskFactorID==riskid])
		}else
		{
			normal     <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==normal]
		}
		
		tmp <- data.frame(Name=fac, RiskFactorID = riskid, FactorValue = riskvalue, DR = riskvalue - normal )
		nonmod <- rbind(nonmod, tmp)
	}
	
	#Blood pressure, max of DBP and SBP and noth both
	#BP  <- min(asses$DR[asses$Name=="DBP"], asses$DR[asses$Name=="SBP"] )
	#BP_ <- min(asses$minFacVal[asses$Name=="DBP"], asses$minFacVal[asses$Name=="SBP"] )
	BP=0; BP_=0
	
	#Family History and Genetics
	family_h <- 0; genetic <- 0
	genefam  <- data.frame(Name=c(), RiskFactorID = c(), FactorValue = c(), DR=c())
	genetic_detail  <- data.frame(Name=c(), RiskFactorID = c(), FactorValue = c(), DR=c())
	for (fac in FH)
	{
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac][1]
		normal     <- limits$normal[limits$name==fac][1]
		if (is.na(normal) )
		{
			normal <- 0
		}
		riskvalue  <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
		family_h   <- family_h + riskvalue
		tmp <- data.frame(Name=fac, RiskFactorID = riskid, FactorValue = riskvalue, DR = riskvalue - normal )
		genefam <- rbind(genefam, tmp)
	}
	for (fac in Genetics)
	{
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac][1]
		normal     <- limits$normal[limits$name==fac][1]
		riskvalue  <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
		if (length(riskvalue) > 1 )
		{
			riskvalue <- riskvalue[1]
		} # revise later
		genetic    <- genetic + riskvalue
		tmp  <- data.frame(Name="Genetics", RiskFactorID = riskid, FactorValue = genetic, DR = genetic - normal )
		tmp2 <- data.frame(Name=fac, RiskFactorID = riskid, FactorValue = riskvalue, DR = riskvalue - normal )
		genetic_detail <- rbind(genetic_detail, tmp2)
	}
	genefam <- rbind(genefam, tmp)
	
	lifetimerisk <- round(100*exp(abs(sum(asses$DR) + sum(nonmod$DR) - BP + max(family_h, genetic) )))/100
	achievable   <- round(100*exp(abs(sum(nonmod$DR) + sum(asses$minFacVal) - BP_ + max(family_h, genetic) )  ))/100
	
	# gui setting for protective vs causative
	if ((sum(asses$DR) + sum(nonmod$DR) - BP + max(family_h, genetic) ) < 0 )
	{
		main = paste("LifeTime Risk : ", lifetimerisk, "Fold", "-", "\n")
	}else
	{
		main = paste("LifeTime Risk : ", lifetimerisk, "Fold", "+", "\n")
	}
	if (sum(nonmod$DR) + sum(asses$minFacVal)+ - BP_ + max(family_h, genetic) < 0 )
	{
		main= paste(main, "Achievable Risk : ", achievable, "Fold", "-", "\n",sep=" ")
	}else
	{
		main= paste(main, "Achievable Risk : ", achievable, "Fold", "+", "\n",sep=" ")
	}
	
	nonmod <- rbind(nonmod, genefam)
	blbls  <- as.character(nonmod$Name)
	bp1    <- barplot(nonmod$DR, names.arg=blbls, col=rainbow(length(blbls)),xlim=c(log(1/4), log(16)),
		,main=main, 
		width=0.5, horiz=TRUE, beside=TRUE, axes=FALSE)
	axis(1, at=seq(log(1/4),log(16),by=log(4)), labels=c("4", "base","4","16"))
	mtext(round(100*exp(abs(nonmod$DR)))/100, side=4, at=bp1,line=0.5,  cex=1.5 )
	
	
	blbls <- Modifiable
	bp <- barplot(asses$DR, names.arg=blbls, col=rainbow(length(blbls)), xlim=c(log(1/4), log(16)), width=0.5, horiz=TRUE, 
		      beside=TRUE, legend=TRUE, axes=FALSE)
	axis(1, at=seq(log(1/4),log(16),by=log(4)), labels=c("4", "base","4","16"))
	mtext(round(100*exp(abs(asses$DR)))/100, side=4, at=bp,line=0.5,  cex=1.5 )
	
	# details of Genetics:
	dev.new()
	par(lty=2, lwd=2)
	nonmod <- genetic_detail
	blbls  <- Genetics
	bp1    <- barplot(nonmod$DR, names.arg=blbls, col=rainbow(length(blbls)),xlim=c(log(1/4), log(16)),
	,main="Genetics Expanded", 
	width=0.5, horiz=TRUE, beside=TRUE, axes=FALSE)
	axis(1, at=seq(log(1/4),log(16),by=log(4)), labels=c("4", "base","4","16"))
	mtext(round(100*exp(abs(nonmod$DR)))/100, side=4, at=bp1,line=0.5,  cex=1.5 )
	
	
}

Pressedlife <- function()
{
	paf <- function(x)
	{
		x2 = x
		intervals = c(0, 0.3, 0.6, 0.8, 0.85, 0.9, 0.99 )
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
					x2[j] = 6
				}
			}
		}
		return(x2)
	}
		
	AllFactors    <- c(Modifiable, NonModifiable, GE)
	x=c(); y=c(); n=c();
	asses <- data.frame(Name=c(), RiskFactorID = c(), FactorValue = c())
	healthy <- data.frame(Name=c(), RiskFactorID = c(), FactorValue = c())
	
	for ( fac in Modifiable)
	{
					
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac][1]
		normal     <- limits$normal[limits$name==fac][1]
		normal     <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==normal]
		riskvalue  <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
		tmp        <- data.frame(Name=fac, RiskFactorID = riskid, FactorValue = riskvalue ) 
		asses      <- rbind(asses, tmp)
		tmp        <- data.frame(Name=fac, RiskFactorID = riskid, FactorValue = min(Data$FactorValue[Data$RiskFactorID==riskid]) )
		healthy    <- rbind(healthy, tmp)
			
	}
	
	for (fac in c(GE, NonModifiable) )
	{
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac][1]
		normal     <- limits$normal[limits$name==fac][1]
		normal     <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==normal]
		riskvalue  <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
		tmp        <- data.frame(Name=fac, RiskFactorID = riskid, FactorValue = riskvalue )
		healthy    <- rbind(healthy, tmp)
		asses      <- rbind(asses, tmp)
	}
	#Family history and Genetics
	family_h <- 0; genetic <- 0
	for (fac in FH)
	{
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac][1]
		normal     <- limits$normal[limits$name==fac][1]
		normal     <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==normal]
		family_h   <- family_h + Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var]
	}
	for (fac in Genetics)
	{
		var        <- as.numeric(tclvalue(get(paste("sl_", fac))) )
		riskid     <- limits$RiskFactorID[limits$name==fac][1]
		normal     <- limits$normal[limits$name==fac][1]
		normal     <- Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==normal]
		genetic    <- genetic + Data$FactorValue[Data$RiskFactorID==riskid & Data$RiskFactorValue==var][1]
	}
			
			
	healthy_risk  <- sum(healthy$FactorValue) + Data$FactorValue[Data$comment == "Intercept"] + max(family_h, genetic)
	current_risk  <- sum(asses$FactorValue)   + Data$FactorValue[Data$comment == "Intercept"] + max(family_h, genetic)
	
	riskid_age   <- limits$RiskFactorID[limits$name=="Age"][1]
	for ( age in seq(limits$Min[limits$name=="Age"][1], limits$Max[limits$name=="Age"][1], by=5 ) )
	{
		x <- c(x, age)
		y <- c(y, current_risk + Data$FactorValue[Data$RiskFactorID==riskid_age & Data$RiskFactorValue==age] )
		n <- c(n, healthy_risk + Data$FactorValue[Data$RiskFactorID==riskid_age & Data$RiskFactorValue==age] )
	}
			
	prob  <- exp(y)/(exp(y) + 1) 
	nprob <- exp(n)/(exp(n) + 1) 
	dev.new()
	par(lwd=2)
	plot  (x, paf(prob), type="l", col="blue", xlab="Age", ylab="Risk ", ylim=c(0, 7), axes=FALSE, xlim=c(limits$Min[limits$name=="Age"][1], limits$Max[limits$name=="Age"][1]) )
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

					

	