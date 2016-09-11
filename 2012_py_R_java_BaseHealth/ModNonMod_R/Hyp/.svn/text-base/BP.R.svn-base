
dbp <- Data[Data$RiskFactorID==41, ]
sbp <- Data[Data$RiskFactorID==40, ]

BP <- c()

for (i in c(1:nrow(dbp)) )
	for (j in c(1:nrow(sbp) ) )
	{
		BP <- rbind(BP, data.frame(ID=-10, dp_id=dbp[i,]$RiskFactorValue, sp_id=sbp[j,]$RiskFactorValue, FactorID=(as.numeric(dbp[i,]$RiskFactorValue) * 1000+ as.numeric(sbp[j,]$RiskFactorValue) ), sp_fac=sbp[j,]$FactorValue, bp_fac=dbp[i,]$FactorValue, FactorValue=max(sbp[j,]$FactorValue, dbp[i,]$FactorValue) ) )
	}
write.csv(BP, 'BP.csv')
print (BP)
