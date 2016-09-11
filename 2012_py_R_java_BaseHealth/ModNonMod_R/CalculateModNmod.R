##### functions to optain mod non-mod and intercepts for a disease
### TODO addo logs
### TODO change GUI to use data_modNmod as input
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
print(args)
#library(RMySQL) # for Database stuff
library(SiZer) # for piecewise linear regression
library(MASS)
library(rmongodb)

### -------------- checkConfidence --------- ###
# Funtion to check confidence intervals, if it includes 1, set OR back to 1
checkConfidence <- function(odd_ratios, confidence, ranges)
{
	cat("In function checkConfidence for:\t ORs: ",odd_ratios, '\t CI: ',confidence,'\n')
	new_or = c(); 
	new_range <- c()
	
	### check CI and OR number is the same
	if ( (length(odd_ratios)) != length(confidence) )
		stop('CIs and ORs dont have the same length')
	
	for (nodd in c(1:length(odd_ratios)) )
	{
		### reference 
		if (as.numeric(odd_ratios[nodd]) == 1)
		{
			if (!(confidence[nodd]=='na'))
				stop('Reference CI is not NA')
			new_or <- c(new_or, odd_ratios[nodd])
			new_range <- c(new_range, ranges[nodd])
		}else
		{
			ci_values  <- as.numeric(unlist(strsplit(confidence[nodd], '-')) )
			or_value  <- as.numeric(odd_ratios[nodd])
			### ignore non-significant data point
			if ( (ci_values[1] <= 1) & (ci_values[2] >= 1) & length(ci_values)==2)
			{
				cat("CI Non-significant (point excluded):\tOR: ",or_value, "\t CI", ci_values[1], "-", ci_values[2], "\n")
			}else
			{
				new_or <- c(new_or, or_value)
				new_range <- c(new_range, ranges[nodd])
			}
		}
	}
	cat("DONE checkConfidence:\t ORs: ",new_or, '\n ranges: ',new_range,'\n')
	
	return(data.frame(new_or=new_or, new_range=new_range))
}

### ---------- capPrediction -------------- ###
# cap predicted values to min/max OR
capPrediction <- function(x_4_pred,prediction, x_val, y_val, model,data, mid, ORs)
{
	cat("In function capPrediction ", model,'\n')
	prediction[prediction < min(y_val) ]<- min(y_val)
	prediction[prediction > max(y_val) ]<- max(y_val)
	# plot for record
	plot(x_4_pred, prediction, col='blue', , type='l', main=paste(model,data$risk_id, data$gender, sep="  " ))
	points(x_val, y_val, col='red')
	points(mid, log(as.numeric(ORs)), col='green')
	grid(NULL)
	return(prediction)
}

### ---------- chooseRegressionModel -------------- ###
# choose between linear and piecewise gregression models
chooseRegressionModel <- function(x, y, value_int,data, mid, ORs)
{
	cat("In function chooseRegressionModel for ", data$risk_id,"\n")
	cat(x,'\n')
	cat(value_int,'\n')
	cat(mid,'\n')
	### --- Linear Model --- ###
	linear_or <- predict(lm(y~x), data.frame(x=value_int) )
	linear_or = capPrediction(value_int, linear_or, x, y, "Linear",data, mid, ORs)
	
	### --- PieceWise Linear Model --- ###
	piece_or  <- predict(piecewise.linear(x, y, middle=1, bootstrap.samples=10000), value_int)
	piece_or = capPrediction(value_int,piece_or, x, y, "PieceWise",data, mid, ORs)
	
	predictions = matrix(ncol=2,nrow=length(value_int))
	predictions[,1]=as.matrix(as.vector(linear_or))[,1]
	predictions[,2]=as.matrix(as.vector(piece_or))[,1]
	
	gof = c(0,0)
	for (v in c(1:length(x)) ) 
	{
		val = x[v] # midpoint value
		### cap to min/max
		if(val<min(value_int))
			val=min(value_int)+0.5
		if(val>max(value_int))
			val=max(value_int)-0.5
		### get integer values around midpint value
		if(as.integer(val) == val){
			index = which(value_int==as.integer(val), arr.ind=TRUE)
			flanking=c(index,index)
		}else{
			min_val = which(value_int==as.integer(val-0.5), arr.ind=TRUE)
			max_val = which(value_int==as.integer(val+0.5), arr.ind=TRUE)
			flanking = c(min_val,max_val)
#			index = which(value_int==as.integer(val), arr.ind=TRUE)
#			flanking=c(index,index)
		}
		error=c(0,0)
		# calculate chi^2 error
		for(m in c(1:length(flanking))) 
			for(e in c(1:2)){
				error[e]=error[e]+(predictions[flanking[m],e]- y[v])^2
			}
		# mean of chi^2 error between values around midpoint
		for(e in c(1:2)){
			error[e]=error[e]/2
			gof[e] = gof[e]+error[e]
		}
	}
	### return fitest model
	print((data$is_step == 'no'))
	print(data$is_step)
	best_model = which(gof==min(gof), arr.ind=TRUE)
	
	if(is.na(data$is_step)){
		if(best_model == 1){
			cat('Chose Linear for ',data$risk_id,gof,'\n')
		}else{
			cat('Chose Piecewise Linear for ',data$risk_id,gof,'\n')
		}
	}else if(data$is_step == 'no'){
		best_model = 1
		cat('Chose Linear for ',data$risk_id,gof,'\n')
	}else{
		best_model = 2
		cat('Chose Piecewise Linear for ',data$risk_id,gof,'\n')
	}
	return(predictions[,best_model])
}

### ---------- getCapIndex -------------- ###
#  find min or max index that dont have min or max or
getCapIndex <- function(ref, start, end, data){
	cat("In function getCapIndex ", ref,'\n')
	index	= start
	for(v in c(start:end)){
		if(data[v] != ref){
			index = v 
			break
		}
	}
	return(index)
}

### ---------- capContinousData -------------- ###
# cap the cooked data to not repeat min and max OR values
capContinousData <- function(real_val, value_int,ORs){
	cat("In function capContinousData \n")
	new_log_or = c(); new_real_value =c(); new_value_as_int = c()
	### cap bigining
	min_index = getCapIndex(ORs[1], 2, length(ORs), ORs) -1 
	max_index = getCapIndex(ORs[length(ORs)],length(ORs)-1,1 , ORs) +1
	cat('new min:',min_index,',old max:',length(ORs), ',new max:',max_index ,'\n')
	no = 1
	
	for(o in c(min_index:max_index)){
		new_log_or[no] = ORs[o]
		new_real_value[no] = real_val[o]
		new_value_as_int[no] = value_int[o]
		no = no + 1
	}
	return(data.frame(new_log_or=new_log_or, new_real_value=new_real_value, new_value_as_int=new_value_as_int))
}

### ---------- getContinuouslRiskFactors -------------- ###
getContinuouslRiskFactors <- function(data_cooked, data_nrisk, disease_id)
{
	cat("In function getContinuouslRiskFactors for ",disease_id,"\n")
	x_axis = c(); y_axis=c(); mid_points<- c()
	
	### record ranges mid point
	ranges <- unlist(strsplit(data_nrisk$ranges, ','))
	for (j in c(1:length(ranges)) ) 
	{
		min_max    <- as.numeric(unlist(strsplit(ranges[j], '-')) )
		mid_points  <- c(mid_points, 0.5*(min_max[1] + min_max[2]) )
	}
	
	odd_ratios <- unlist(strsplit(data_nrisk$odds_ratio, ','))
	confidence <- unlist(strsplit(data_nrisk$ci, ','))
	
	### check confidence intevals significance
	cat ("\n")
	cat("----------------------\n")
	cat('risk factor: ',data_nrisk$risk_id, "\tRanges : ",ranges, "\n")
	df <- checkConfidence(odd_ratios, confidence, ranges)
	newor <- as.character(df$new_or)
	ranges <- as.character(df$new_range)
	
	### find  points for linear regression
	for (j in c(1:length(ranges)) )
	{
		range <- as.numeric(unlist(strsplit(ranges[j], '-')) )
		or    <- as.numeric(newor[j])
		x_axis  <- c(x_axis, 0.5*(range[1] + range[2]) )
		y_axis  <- c(y_axis, log(or))	
	}
	
	### ---------- TODO change values as integer for weird risk factors here --------- ###
	if (data_nrisk$risk_id == 'bmi')
	{
		x_axis = x_axis*10
		mid_points = mid_points*10
		value_as_int <- c((10*as.numeric(data_nrisk$min_val)) : (10*as.numeric(data_nrisk$max_val)) )
		real_value = c((10*as.numeric(data_nrisk$min_val)) : (10*as.numeric(data_nrisk$max_val)) )/10
		data_nrisk$normal_val = 10*as.numeric(data_nrisk$normal_val)
		cat(value_as_int ,'\n')
	}
	else if (data_nrisk$risk_id == 'diesel_exhaust' & data_nrisk$risk_cat == 'high')
	{
		x_axis = x_axis*10
		mid_points = mid_points*10
		value_as_int <- c((10*as.numeric(data_nrisk$min_val)) : (10*as.numeric(data_nrisk$max_val)) )
		real_value = c((10*as.numeric(data_nrisk$min_val)) : (10*as.numeric(data_nrisk$max_val)) )/10
		data_nrisk$normal_val = 10*as.numeric(data_nrisk$normal_val)
		cat(value_as_int ,'\n')
	}
	else if (!is.na(data_nrisk$risk_cat) & (data_nrisk$risk_cat == 'female' ))
	{
		x_axis = x_axis+100
		mid_points = mid_points+100
		value_as_int <- c((100+as.numeric(data_nrisk$min_val)) : (100+as.numeric(data_nrisk$max_val)) )
		real_value = c((100+as.numeric(data_nrisk$min_val)) : (100+as.numeric(data_nrisk$max_val)) )
		data_nrisk$normal_val = 100+as.numeric(data_nrisk$normal_val)
		cat(value_as_int ,'\n')
	}
	else
	{
		value_as_int <- c(data_nrisk$min_val: data_nrisk$max_val)
		real_value  = value_as_int
	}
	### regression is either linear or piecewise linear
	log_or = chooseRegressionModel(x_axis, y_axis, value_as_int, data_nrisk, mid_points, odd_ratios)
	### add capping of value_as_int/real_value/log_or
	capped_data = capContinousData(real_value, value_as_int,log_or)
	new_log_or = capped_data$new_log_or
	new_real_value = capped_data$new_real_value
	new_value_as_int = capped_data$new_value_as_int
	### add NA for this risk factor
	reference = which(new_log_or==min(new_log_or), arr.ind=TRUE)
	data_cooked = rbind(data_cooked,
			data.frame(
					risk_id=data_nrisk$risk_id,
					value = new_real_value[reference[1]],
					value_as_int=NA, 
					log_or= min(new_log_or), 
					min= min(new_value_as_int), 
					max=max(new_value_as_int),  
					normal=new_value_as_int[reference[1]],
					modifiable=data_nrisk$risk_type, 
					protective=data_nrisk$protective,
					gender=data_nrisk$gender,
					risk_cat =data_nrisk$risk_cat,
					is_step = data_nrisk$is_step))
	### fill data object with new continuous risk factor
	data_cooked = rbind(data_cooked,
			data.frame(
					risk_id=data_nrisk$risk_id,
					value= new_real_value,
					value_as_int=new_value_as_int, 
					log_or=new_log_or, 
					min= min(new_value_as_int), 
					max=max(new_value_as_int),
					normal=data_nrisk$normal_val,
					modifiable=data_nrisk$risk_type, 
					protective=data_nrisk$protective,
					gender=data_nrisk$gender,
					risk_cat =data_nrisk$risk_cat,
					is_step = data_nrisk$is_step))
	return(data_cooked)
}

### ---------- getCategoricalRikFactors -------------- ###
getCategoricalRikFactors <- function(data_cooked, data_nrisk, disease_id)
{
	cat("In function getCategoricalRikFactors for ",disease_id, ' risk factor: ',data_nrisk$risk_id,"\n")
	categories  <- unlist(strsplit(data_nrisk$cat, ',')) 
	odd_ratios   <- log(as.numeric(unlist(strsplit(data_nrisk$odds_ratio, ',')) ))
	### sanity check # categories # ORS
	if (length(categories) != length(odd_ratios) )
		stop(paste0('categories and ORs have different lengths for risk: ',data_nrisk$risk_id))
	### add NA for this risk factor
	reference = which(odd_ratios==min(odd_ratios), arr.ind=TRUE)
	data_cooked = rbind(data_cooked,
			data.frame(
					risk_id=data_nrisk$risk_id,
					value= categories[reference],
					value_as_int=NA, 
					log_or= min(odd_ratios), 
					min= "0", 
					max=as.character(length(categories)-1), 
					normal=as.character(reference-1),
					modifiable=data_nrisk$risk_type, 
					protective=data_nrisk$protective,
					gender=data_nrisk$gender,
					risk_cat =data_nrisk$risk_cat,
					is_step = data_nrisk$is_step))
	### fill data object with new categorical risk factor
	data_cooked = rbind(data_cooked,
			data.frame(
					risk_id=data_nrisk$risk_id,
					value= categories,
					value_as_int=c(0:(length(categories)-1)), 
					log_or=odd_ratios, 
					min= "0", 
					max=as.character(length(categories)-1), 
					normal=data_nrisk$normal_val,
					modifiable=data_nrisk$risk_type, 
					protective=data_nrisk$protective,
					gender=data_nrisk$gender,
					risk_cat = data_nrisk$risk_cat,
					is_step = data_nrisk$is_step))
	### plot categories + OR
	barplot(exp(odd_ratios), names.arg=categories, col='blue', main=data_nrisk$risk_id)
	
	return(data_cooked)
}

### -------------- checkConfidence --------- ###
# Funtion to check confidence intervals, if it includes 1, set OR back to 1
getIntercept <- function(data_cooked, data_nrisk, disease_id)
{
	incident       <- as.numeric(data_nrisk$odds_ratio)
	if(incident>1 | incident< 0)
		stop(paste0("INCIDENCE WRONG ",incidence))
	odd_ratios     <- c(log(incident/(1 - incident) ) )
	### fill data object with new categorical risk factor
	data_cooked = rbind(data_cooked,
			data.frame(
					risk_id="intercepts",
					value= data_nrisk$risk_id,
					value_as_int=1, 
					log_or=odd_ratios, 
					min= NA, 
					max=NA, 
					normal=NA,
					modifiable=NA, 
					protective=NA,
					gender=NA,
					risk_cat =NA,
					is_step = NA))
	return(data_cooked)
}

### ---------- getModNonModRiskFactors -------------- ###
getModNonModRiskFactors <- function(data_table, disease_id)
{
	cat("In function getModNonModRiskFactors for ",disease_id,"\n")
	### initialize data_object object
	data_object <- data.frame(risk_id=c(), value =c(),value_as_int=c(), log_or=c(), min=c(), max=c(), normal=c(),modifiable=c(), protective=c(), gender=c() )
	
	### graphs of risk factor saved in pdf
	pdf(file=paste0(disease_id,'_GraphsModNmod','.pdf'))
	total_risk = length(which(data_table$genophen=='yes', arr.ind=TRUE))
	par(lty=2, lwd=2, mfrow=c(ceiling(total_risk/3), 3), mar = rep(2, 4)) 
	for (nrisk in c(1:nrow(data_table)) )
	{
		if (data_table[nrisk,]$genophen == "yes"){
			cat(nrisk, ' ',data_table[nrisk,]$risk_id,"\n")
			### consider only categorical data
			if (is.na(data_table[nrisk,]$is_cat)){
				if(data_table[nrisk,]$risk_id == "incidence"){
					data_object = getIntercept(data_object, data_table[nrisk,], disease_id)
				}else{
					cat('NOT USED RISK FACTOR:',data_table[nrisk,]$risk_id,'\n')
				}
			}
			else if (data_table[nrisk,]$is_cat == "yes")
				data_object = getCategoricalRikFactors(data_object, data_table[nrisk,], disease_id)
			else if(data_table[nrisk,]$is_cat == "no")
				data_object = getContinuouslRiskFactors(data_object, data_table[nrisk,], disease_id)
			else{
				cat('ELSE', data_table[nrisk,]$is_cat,'\n')
			} 
		}
	}
	dev.off()
	data_object$risk_id = as.character(data_object$risk_id)
	cat('DONE getModNonModRiskFactors for ',disease_id,'\n')
	return(data_object)
}

### ---------- getModNonModRiskFactors -------------- ###
getDataFromMongo <- function(host, disease_id)
{
	cat("In function getDataFromMongo for ",disease_id,"\n")
	### TODO change if needed
	db = mongo.create(host=host,username="adm1n", password="adm1npw0rd", db="admin", timeout=0L)
	mod_collection = "comprehensive.mod.n.mod"
	raw_data = mongo.find(db, mod_collection, query =  mongo.bson.from.list(list("disease_id" = disease_id,"genophen" ="yes" )),fields = list("_id"=0))
	ndoc = 1
	while (mongo.cursor.next(raw_data)){
		doc = mongo.bson.to.list(mongo.cursor.value(raw_data))
		if(ndoc == 1){
			keys = names(doc)
			data <-data.frame(matrix(matrix(rep(1,length(keys)),1),1)) 
			colnames(data)<-keys
		}
		values = c()
		for(k in c(1:length(keys))){
			if (is.null(doc[[keys[k]]]))
				doc[[keys[k]]] = NA
			values = c(values,doc[[keys[k]]])
#			cat(keys[k],'-',doc[[keys[k]]],'\n')
		}
		data[ndoc,] = values
		ndoc = ndoc + 1
	}
	return(data)
}

#------------ MAIN ----------#
disease_id = args[1]
mongo_host = args[2]
wd = args[3]
############### TODO TO RUN FROM R ############
#disease_id = "anxiety"
#mongo_host = "54.245.239.40"
#wd = "/Users/celine/Box\ Documents/Celine/ScriptInOut/ModNonMod/"

### change working directory
setwd(wd)
getwd()

cat(disease_id, '\n')

mydata = getDataFromMongo( mongo_host,disease_id)
#print(mydata)

data_modNmod <- getModNonModRiskFactors(mydata, disease_id)

write.matrix(data_modNmod, file = paste0(disease_id,"_data_modNmod",".xls"), sep = "\t")

