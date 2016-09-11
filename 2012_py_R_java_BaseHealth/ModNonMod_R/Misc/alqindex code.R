# Alcohol Consumption Score
# Vars used: ALQ100/101, ALQ110, ALQ120Q, ALQ120U, ALQ130, ALQ140Q, ALQ140U, ALQ150
# All variables in alq_[abcd].csv


create.alqIndex <- function(alq.dataset) {

alq <- alq.dataset

### number of drinks in last year
ALQ100_101 <- as.numeric(as.vector(alq[,2]))
ALQ110 <- as.numeric(as.vector(alq[,3]))
ALQ120Q <- as.numeric(as.vector(alq[,4]))
ALQ120U <- as.numeric(as.vector(alq[,5]))
ALQ130 <- as.numeric(as.vector(alq[,6]))

n_drinks <- vector(mode="numeric", length=dim(alq)[1])
n_drinks[1:dim(alq)[1]] <- NA



### days of 5+ drinks in last year
ALQ140Q <- as.numeric(as.vector(alq[,7]))
ALQ140U <- as.numeric(as.vector(alq[,8]))

n_binges <- vector(mode="numeric", length=dim(alq)[1])
n_binges[1:dim(alq)[1]] <- NA


### 5+ drinks every day during period of life (indicating history of alcoholism)
ALQ150 <- as.numeric(as.vector(alq[,9]))

hx_alc <- vector(mode="numeric", length=dim(alq)[1])
hx_alc[1:dim(alq)[1]] <- NA






for( i in 1:dim(alq)[1] ) {

	# people who have never drank >12 drinksin a year (ALQ100/ALQ101)
	if( !is.na(ALQ100_101[i]) & ALQ100_101[i] == 2) {
		n_drinks[i] <- 0
		n_binges[i] <- 0
		
		if( ALQ110[i] == 2 ) {
			hx_alc[i] <- 0
		}
	}

	
	# people who have had 0 drinks in the last year (ALQ120Q)
	if( !is.na(ALQ120Q[i]) & ALQ120Q[i]  == 0 ) {
		n_drinks[i] <- 0
		n_binges[i] <- 0
	}

	# calculate the number of drinks subjects had in the last year (ALQ120Q, ALQ120U, ALQ130)
	if( ALQ100_101[i] == 1 & !is.na(ALQ120Q[i]) & !is.na(ALQ120U[i]) & !is.na(ALQ130[i]) & ALQ120Q[i] != 999 & ALQ120U[i] != 999 & ALQ130[i] != 999 & ALQ130[i] != 99) {
		
		# check frequency unit of drinking (ALQ120U)
		if( ALQ120U[i] == 1) { # per week
			n_drinks[i] <- ALQ120Q[i] * 52 * ALQ130[i]
		}

		if( ALQ120U[i] == 2) { # per month
			n_drinks[i] <- ALQ120Q[i] * 12 * ALQ130[i]
		}

		if( ALQ120U[i] == 3) { # per year
			n_drinks[i] <- ALQ120Q[i] * 1 * ALQ130[i]
		}
	}


	# calculate the number of binge drinking episodes in the last year (ALQ140Q, ALQ140U)
	if( !is.na(ALQ140Q[i]) & ALQ140Q[i] == 0) {
		n_binges[i] <- 0
	}
	if( !is.na(ALQ140Q[i]) & ALQ140Q[i] > 0 & ALQ140Q[i] != 9999 & !is.na(ALQ140U[i]) ) {
		
		# check frequency unit of binge drinking (ALQ140U)
		if( ALQ140U[i] == 1) { # per week
			n_binges[i] <- ALQ140Q[i] * 52
		}

		if( ALQ140U[i] == 2) { # per month
			n_binges[i] <- ALQ140Q[i] * 12
		}

		if( ALQ140U[i] == 3) { # per year
			n_binges[i] <- ALQ140Q[i] * 1
		}

	}

	# assess history of alcoholism
	if( !is.na(ALQ150[i]) ) {
		if( ALQ150[i] == 1 ) { # Yes
			hx_alc[i] <- 1
		}
		if( ALQ150[i] == 2 ) { # No
			hx_alc[i] <- 0

		}
	}

}



drink_score <- vector(mode="numeric", length=dim(alq)[1])
drink_score[1:dim(alq)[1]] <- NA

drink_score[n_drinks == 0] <- 0
drink_score[n_drinks > 0 & n_drinks <= 365] <- 1
drink_score[n_drinks > 365 & n_drinks <= 1000] <- 2
drink_score[n_drinks > 1000] <- 3


# median number of binge drinking incidents is 12, make this the cutoff
binge_score_cutoff <- median(n_binges[n_binges>0], na.rm=T)

binge_score <- vector(mode="numeric", length=dim(alq)[1])
binge_score[1:dim(alq)[1]] <- NA

binge_score[n_binges <= binge_score_cutoff] <- 0
binge_score[n_binges > binge_score_cutoff] <- 1


# total score
total_alcohol_score <- drink_score + binge_score + hx_alc

alq_index <- cbind(alq, n_drinks, drink_score, n_binges, binge_score, hx_alc, total_alcohol_score)
write.table(alq_index, file="alq_index.csv", sep=",")

}