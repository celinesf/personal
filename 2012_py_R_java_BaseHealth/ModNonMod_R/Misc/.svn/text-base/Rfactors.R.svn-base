# to generate riskfactors values for riskfactors 
# model :  Test <- glm(T2D ~ Age + BMIC2+Gender+Race+WaistAdBMI  +CHR2, data=Data,family	=binomial);

factors <- c("Age", "BMIC2", "Gender", "Race", "WaistAdBMI", "CHR2") # model predictors
RF_id=c(); RF_val=c(); F_val=c(); comment = c();
cats  <- Test$xlevels #categorical variables in our model

#==========================================
#Age 0-20 has zero risk
RF_id      <- c(RF_id, rep(1, 21))
RF_val     <- c(RF_val, c(0:20))
F_val      <- c(F_val, rep(0, 21))
comment    <- c(comment, rep("Age", 21))
#==========================================
id <- 1
for (fac in factors)
{
	if (fac%in%names(cats))
	{
		print(fac)
		levels <- as.character(cats[[fac]])
		print(levels)
		if (fac == "Age")
		{
			levels <- levels[3:length(levels)] # under Age of 20, risks are zero
			prev <- Test$coefficients[["Age2"]]
			base <- 0
		}
			
		for (lev in levels) 
		{
			var <- paste(fac, lev, sep="")
			if (var%in%names(Test$coefficients))
			{
				if (fac == "Age")
				{
					F_val = c(F_val, base + c(1:10)*(Test$coefficients[[var]]-prev) )
					RF_id   <- c(RF_id, rep(id, 10))
					RF_val  <- c(RF_val, c(((as.numeric(lev)-1)*10+1):((as.numeric(lev)-1)*10 +10)) )
					comment <- c(comment, rep(fac, 10))
					base    <- 10 *(Test$coefficients[[var]]-prev ) + base
					prev    <- Test$coefficients[[var]]
					
				}
				else
				{
					F_val   <- c(F_val,Test$coefficients[[var]])
					RF_id   <- c(RF_id, id)
					RF_val  <- c(RF_val, lev)
					comment <- c(comment, fac)
				}				
			}
			else
			{
				F_val   <- c(F_val, 0)
				RF_id   <- c(RF_id, id)
				RF_val  <- c(RF_val, lev)
				comment <- c(comment, fac)
			}
		}
	}
	else
	{
		
	}
	id <- id + 1
		
			
}
table <- data.frame(RiskFactorID=RF_id, RiskFactorValue=RF_val, FactorValue=F_val, comment = comment)		
