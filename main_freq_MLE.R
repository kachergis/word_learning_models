library(plotrix); library(lattice); library(pso); library(car)
# George Kachergis July 7, 2012
# fit associative model and Fazly model using MLE to the 
# frequency/CD experiments, subject-by-subject

make_exp3_ind_graphs <- function() {
	source("paper_graphs_ggplot2.R")
	load("humans/freqCD_data_exps8_9_10_12.Rdata") # raw 159 Ss
	#load("all_fazly_model_MLE_PSOfits.RData")
	load("freqCD_fazly_model_MLE_PSOfits.RData")
	freq_graph(raw, fazfreq, "fazly")
	load("freqCD_assoc_model_MLE_PSOfits.RData")
	freq_graph(raw, assfreq, "assoc")
}

make_exp3_group_graphs <- function() {
	source("paper_graphs_ggplot2.R")
	load("humans/freqCD_data_exps8_9_10_12.Rdata") # raw 159 Ss
	
	load("freq_fazly_model_MLE_fitPSO_group.RData")
	freq_graph_group(raw, fazfreq, "fazly_group")
	load("freq_assoc_model_MLE_fitPSO_group.RData")
	freq_graph_group(raw, assfreq, "assoc_group")
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
	controls = list(maxit=200, max.restart=2, reltol=.05)
	#agg <- with(agg, aggregate(Correct, list(Subject=Subject, file=file), mean))
	if(modeln=="fazly") {
		agg$lambda <- NA
		agg$beta <- NA
		agg$theta <- NA
		parcols = c("lambda","beta","theta")
		initp = c(.00001,1000,.7)
		#lowerp = c(.0000001,100,.69) # from their paper
		#upperp = c(.5,10000, .71)
		lowerp = c(.0000001,10,.01) # let's be more lax
		upperp = c(.5,10000, 1)
	} else if(modeln=="model" | modeln=="modelBigM") {
		agg$X <- NA
		agg$B <- NA
		agg$C <- NA
		parcols = c("X","B","C")
		initp = c(.2146, 6.931, .982) # Exp 3: Grow M SSE=.064 1/100 startval
		lowerp = c(0,.01,.7)
		upperp = c(40,18, 1)
	}
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
	controls = list(maxit=100, max.restart=2, reltol=.1)
	#agg <- with(agg, aggregate(Correct, list(Subject=Subject, file=file), mean))
	if(modeln=="fazly") {
		agg$lambda <- NA
		agg$beta <- NA
		agg$theta <- NA
		parcols = c("lambda","beta","theta")
		initp = c(.00001,1000,.7)
		lowerp = c(.0000001,100,.69)
		upperp = c(.5,10000, .71)
	} else if(modeln=="model" | modeln=="modelBigM") {
		agg$X <- NA
		agg$B <- NA
		agg$C <- NA
		parcols = c("X","B","C")
		initp = c(.2146, 6.931, .982) # Exp 3: Grow M SSE=.064 1/100 startval
		lowerp = c(0,.01,.7)
		upperp = c(40,15, 1)
	}
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


fazfreq <- fit_all_MLE("fazly", PSO=T)
save(fazfreq, file="freq_fazly_model_MLE_fitPSO_group.RData")
assfreq <- fit_all_MLE("model", PSO=T)
save(assfreq, file="freq_assoc_model_MLE_fitPSO_group.RData")
print(paste("assoc all PSO:",assfreq$ML)) 
print(paste("fazly all PSO:",fazfreq$ML))
# "assoc all PSO: 2445.0 .33 18 .99
# "fazly all PSO: 2507.57 (2627.3 OLD) .03 1822 .29

make_exp3_group_graphs()

# model comparison for dissertation:
n = 3690 # total data points (responses)
# BIC:
2*2445.0 + 3*log(n) # 4914.64 # wins! (BIC is stricter)
2*fazfreq$ML + 2*log(n) # 5031.57

#assfreq <- fit_by_subj("model", PSO=F)
#save(assfreq, file="all_assoc_model_MLE_fits.RData") # 32.97
#print(paste("assoc no PSO:",mean(assfreq$ML))) # 35.83

#fazfreq <- fit_by_subj("fazly", PSO=F)
#save(fazfreq, file="all_fazly_model_MLE_fits.RData") # 
#print(paste("fazly no PSO:",mean(fazfreq$ML))) # 46.47

#assfreq <- fit_by_subj("modelBigM", PSO=T)
#save(assfreq, file="all_assBigM_model_MLE_PSOfits.RData")
#print(paste("assBigM PSO:",mean(assfreq$ML))) # no PSO: 35.73


