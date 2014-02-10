library(plotrix); library(lattice); library(pso); library(car)
# George Kachergis July 7, 2012
# fit associative model and Fazly model using MLE to the 
# frequency/CD experiments, subject-by-subject

make_exp3_ind_graphs <- function() {
	source("paper_graphs_ggplot2.R")
	load("humans/freqCD_data_exps8_9_10_12.Rdata") # raw 159 Ss
	load("freqCD_tilles_model_MLE_PSOfits.RData")
	freq_graph(raw, tilfreq, "tilles")
}

make_exp3_group_graphs <- function() {
	source("paper_graphs_ggplot2.R")
	load("humans/freqCD_data_exps8_9_10_12.Rdata") # raw 159 Ss
	load("freq_tilles_model_MLE_fitPSO_group.RData")
	freq_graph_group(raw, tilfreq, "tilles_group")
}


#assfreq <- fit_by_subj("model", PSO=F)
#save(assfreq, file="freqCD_assoc_model_MLE_fits.RData") # 32.96

# Exp 3: CD and Frequency
#best <- fit(c(.2146, 6.931, .982)) # Grow M SSE=.064 1/100 startval

binomial_likelihood <- function(cdat, M) {
	est_prob = diag(M) / rowSums(M) # prob correct
	lik = 0
	for(i in 1:dim(M)[1]) { 
		resp = cdat[which(cdat$CorrectAns==i),]$Response
		if(resp==i) {
			lik = lik + log(est_prob[i])
		} else {
			lik = lik + log(1-est_prob[i])
		}
	}
	return(lik) # 18*log(1/18) = -52.02669 for guessing
}

multinomial_likelihood <- function(cdat, M) {
	M = M / rowSums(M) # p(o|w)
	lik = 0
	for(i in 1:dim(M)[1]) { 
		resp = cdat[which(cdat$CorrectAns==i),]$Response
		lik = lik + log(M[i,resp]) 
	}
	return(lik) # 18*log(1/18) = -52.02669 for guessing
}

fit_subj <- function(par, orders, sdat) {
	# ords is all of the orders
	# sdat has order file names and responses
	tot_lik = 0
	for(f in unique(sdat$file)) {
		o <- orders[[f]] # $trials ?
		#print(paste("Sub",sdat$Subject[1]," File",f))
		M <- model(par, ord=o)
		# need: for each word (CorrectAns), what object they chose (Response)
		tot_lik = tot_lik + binomial_likelihood(subset(sdat, file==f), M)
	}
	mlik = tot_lik / length(unique(sdat$file))
	#print(par)
	#print(mlik)
	return(-mlik)
}

fit_all <- function(par, orders, dat) {
	# ords is all of the orders
	# dat has order file names and responses
	tot_lik = 0
	for(s in unique(dat$Subject)) {
		sdat <- subset(dat, Subject==s)
		for(f in unique(sdat$file)) {
      #print(paste(s,f))
			o <- orders[[f]]
			#print(paste("Sub",sdat$Subject[1]," File",f))
			M <- model(par, ord=o)
			# need: for each word (CorrectAns), what object they chose (Response)
			tot_lik = tot_lik + binomial_likelihood(subset(sdat, file==f), M)
		}
	}
	mlik = tot_lik #/ length(unique(dat$Subject))
	#print(par)
	#print(mlik)
	return(-mlik)
}

fit_all_MLE <- function(modeln, PSO=FALSE) {
	source(paste(modeln,".R",sep="")) 
	load("master_orders.RData") # orders
	load("humans/freqCD_data_exps8_9_10_12.Rdata") # raw 159 Ss
	agg <- subset(raw, Subject!=1224) # block2_369-3x3hiCD somehow had 24 test trials??
	agg <- subset(agg, Exp==12) # 
	agg <- with(agg, aggregate(Correct, list(CorrectAns=CorrectAns, Response=Response, Subject=Subject, file=file), mean))
	controls = list(maxit=200, max.restart=2, reltol=.01)
	#agg <- with(agg, aggregate(Correct, list(Subject=Subject, file=file), mean))
	
	agg$X <- NA
	agg$B <- NA
	agg$A <- NA
	parcols = c("X","B","A")
	initp = c(.6, .8, .85) 
	lowerp = c(0,0,.27) # parm 3 > .268 for "block2_369-3x3hiCD" ...see tilles.R
	upperp = c(1,1,1)

	if(PSO) {
		best <- psoptim(initp, fit_all, orders=orders, dat=agg, lower=lowerp, upper=upperp, control=controls) 
		# control=list(maxit=500, reltol=1e-4) # something like "max.reps=2", too?
	} else {
		best <- optim(initp, fit_all, orders=orders, dat=agg, lower=lowerp, upper=upperp, method="L-BFGS-B") # , control=list(parscale=c(10,10,10), maxit=20) 
	}
	ret <- list()
	for(f in unique(agg$file)) {
		mp = model(best$par, orders[[f]])
		ret[[f]] <- diag(mp) / rowSums(mp)
	}
	ret$par = best$par
	ret$ML = best$value
	return(ret)	
}

fit_by_subj <- function(modeln, PSO=FALSE) {
	source(paste(modeln,".R",sep="")) 
	load("master_orders.RData") # orders
	load("humans/freqCD_data_exps8_9_10_12.Rdata") # raw 159 Ss
	agg <- subset(raw, Subject!=1224) # block2_369-3x3hiCD somehow had 24 test trials??
	agg <- subset(agg, Exp==12) # 
	agg <- with(agg, aggregate(Correct, list(CorrectAns=CorrectAns, Response=Response, Subject=Subject, file=file), mean))
	controls = list(maxit=100, max.restart=2) # reltol=.1
	#agg <- with(agg, aggregate(Correct, list(Subject=Subject, file=file), mean))
	
	agg$X <- NA
	agg$B <- NA
	agg$A <- NA
	parcols = c("X","B","A")
	initp = c(.6, .8, .85) 
	lowerp = c(0.0001,0.0001,0.0001)
	upperp = c(1,1,1)
	
	agg$ML <- NA
	agg$Model <- NA
	
	for(s in unique(agg$Subject)) {
		sdat <- subset(agg, Subject==s)
		if(PSO) {
			best <- psoptim(initp, fit_subj, orders=orders, sdat=sdat, lower=lowerp, upper=upperp, control=controls) 
		} else {
			best <- optim(initp, fit_subj, orders=orders, sdat=sdat, lower=lowerp, upper=upperp, method="L-BFGS-B") # , control=list(parscale=c(10,10,10), maxit=20) 
		}
		for(f in unique(sdat$file)) {
			rows <- which(agg$Subject==s & agg$file==f)
			mp = model(best$par, orders[[f]])
			agg[rows,parcols] = matrix(rep(best$par, length(rows)), nrow=length(rows), byrow=T)
			agg[rows,]$ML = best$value
			agg[rows,]$Model <- diag(mp) / rowSums(mp)
		}
		
		print(agg[rows[1],])
	}
	return(agg)
}


#assfreq <- fit_by_subj("model", PSO=T)
#save(assfreq, file="freqCD_assoc_model_MLE_PSOfits.RData")
#print(paste("assoc PSO:",mean(assfreq$ML))) # 

#fazfreq <- fit_by_subj("fazly", PSO=T)
#save(fazfreq, file="freqCD_fazly_model_MLE_PSOfits.RData") 
#print(paste("fazly PSO:",mean(fazfreq$ML))) # 


#tilfreq <- fit_all_MLE("tilles", PSO=T)
#save(tilfreq, file="freq_tilles_model_MLE_fitPSO_group.RData")
#print(paste("tilles all PSO:",tilfreq$ML)) 
# "tilles all PSO: 2450.6612888634"


tilfreq <- fit_by_subj("fazly", PSO=T)
save(tilfreq, file="freqCD_tilles_model_MLE_PSOfits.RData") 
print(paste("tilles PSO:",mean(tilfreq$ML)))
# "tilles PSO: 12.725226234322"

make_exp3_group_graphs()



