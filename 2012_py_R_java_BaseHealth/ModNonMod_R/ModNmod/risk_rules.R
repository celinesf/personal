library(SiZer)
library(RMySQL) 

getall <- function()
{
	# DB stuff; query IDmap table to map ids to names for GUI
	con    <- dbConnect(MySQL(), user='pouria', password='mojabi', host="localhost", dbname="AlgoPhen")
	qur    <- ('select * Reza_all')
	res    <- dbGetQuery(con, qur)
	dbDisconnect(con)
	
	return(res)
}

hdlGen <- function()
{
	did             <- c(1:5)
	riskfactorvalue <- c(30:80)
	hdlvar          <- data.frame(did=c(1:5), mean=rep(1.7, length(did)), std= abs(rnorm(length(did) ) ) )
	hdldata         <- c()
	
	for (d in did)
	{
		m = hdlvar$mean[hdlvar$did==d]
		s = hdlvar$std[hdlvar$did==d]
		
		factorvalue     <- predict(piecewise.linear(riskfactorvalue[1:5], sort(rnorm(5, m, s )), middle=1, bootstrap.samples=10000), riskfactorvalue )

		plot(riskfactorvalue, factorvalue, col="blue" )
		
		hdldata <- rbind(hdldata, data.frame(did=d,factorv=factorvalue ) )
	}
	
	return(hdldata)

}
