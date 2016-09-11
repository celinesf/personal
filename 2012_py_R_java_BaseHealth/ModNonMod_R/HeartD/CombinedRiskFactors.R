# Blood pressure is a maximum risk you get from Dbp and Sbp
# for documentation we need to combine these together

sbp <- Data[Data$RiskFactorID==40, ]
dbp <- Data[Data$RiskFactorID==41, ]

for (i in nrow(sbp) )
{
	row <- sbp[i,]
	for (j in nrow(dbp))
		{
			row2 <- dbp[j,]
			c    <- data.frame(ID=(row$RiskFactorID + 1000*row2$RiskFactorID), SBP=row$RiskFactorValue, 
			SBP_risk=row$RiskFactorFactor, DBP=row2$RiskFactorValue, 
			DBP_risk=row2$RiskFactorFactor )
			print(c)
			
		}
}
