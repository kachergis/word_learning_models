library(plotrix); library(lattice); library(pso); library(car)
# George Kachergis July 7, 2012
# fit associative model using MLE to 
# asymmetric experiments, subject-by-subject

binomial_likelihood <- function(cdat, M) {
	est_prob = diag(M) / rowSums(M) # prob correct
	lik = 0
  #print(cdat)
  #print(est_prob)
	for(i in 1:length(cdat)) { # dim(M)[1]
		resp = cdat[i]
		if(resp==i) {
			lik = lik + log(est_prob[i])
		} else {
			lik = lik + log(1-est_prob[i])
		}
	}
	return(lik) # 18*log(1/18) = -52.02669 for guessing
}


fit_subj <- function(par, ord, sdat) {
	tot_lik = 0
	M <- model(par, ord=ord)
	# need: for each word (CorrectAns), what object they chose (Response)
	mlik = binomial_likelihood(sdat, M)
	#print(par)
	#print(mlik)
	return(-mlik)
}

fit_all <- function(par, ord, dat) {
	tot_lik = 0
	for(s in dim(dat)[1]) {
		sdat <- unlist(dat[s,])
		M <- model(par, ord=ord)
		# need: for each word (CorrectAns), what object they chose (Response)
		tot_lik = tot_lik + binomial_likelihood(sdat, M)
	}
	mlik = tot_lik #/ length(unique(dat$Subject))
	#print(par)
	#print(mlik)
	return(-mlik)
}

fit_all_MLE <- function(modeln, ord, aggm) {
	source(paste(modeln,".R",sep="")) 
	agg = data.frame(cbind(Subject=1:dim(aggm)[1]))
	controls = list(maxit=200, max.restart=2, reltol=.05)
	if(modeln=="fazly") {
		parcols = c("lambda","beta","theta")
		initp = c(.00001,1000,.7)
		#lowerp = c(.0000001,100,.69) # from their paper
		#upperp = c(.5,10000, .71)
		lowerp = c(.0000001,10,.01) # let's be more lax
		upperp = c(.5,10000, 1)
	} else if(modeln=="model" | modeln=="modelBigM") {
		parcols = c("X","B","C")
		initp = c(.2146, 6.931, .982) # Exp 3: Grow M SSE=.064 1/100 startval
		lowerp = c(0,.01,.7)
		upperp = c(40,18, 1)
	}
	best <- psoptim(initp, fit_all, ord=ord, dat=aggm, lower=lowerp, upper=upperp, control=controls) 
	ret <- list()
	mp = model(best$par, ord)
	ret[["model_perf"]] <- diag(mp) / rowSums(mp)
	ret$par = best$par
	ret$ML = best$value
	return(ret)	
}

fit_by_subj <- function(modeln, ord, aggm) {
	source(paste(modeln,".R",sep="")) 
	controls = list(maxit=100, max.restart=2, reltol=.1)
  agg = data.frame(cbind(Subject=1:dim(aggm)[1]))
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
	agg$Human <- NA
	for(s in 1:dim(aggm)[1]) {
		sdat <- aggm[s,] # cols = responses (objs) to words 1:18
    agg[s,]$Human = mean(sdat==1:18)
		best <- psoptim(initp, fit_subj, ord=ord, sdat=sdat, lower=lowerp, upper=upperp, control=controls) 
		mp = model(best$par, ord)
		agg[s,parcols] = best$par
		agg[s,]$ML = best$value
		agg[s,]$Model <- mean( diag(mp) / rowSums(mp) )
	}
	return(agg)
}

fit_asym_conds <- function(modeln, bySubj=T) {
  dir = "exp_asym/"
  load("asym_conds.RData") # conds[[num]]$train, $resps, $test
  afc18 = c(201,202,203, 204, 205,206,207,215:225) # 204,  Error in if (resp == i) { : missing value where TRUE/FALSE needed
  all_fits = list()
  for(n in afc18) {
    ord = conds[[n]]$train #$trials
    aggm = read.table(paste(dir,n,"_data.txt",sep=''))
    if(bySubj) {
      by = "byS"
      fit = fit_by_subj(modeln, ord, aggm)
    } else {
      by = "byGroup"
      fit = fit_all_MLE(modeln, ord, aggm)
    }
    save(fit, file=paste("asym_cond_",n,"_fit_",by,".RData",sep=''))
    all_fits[[n]] = fit
    print(fit)
    print(paste("done fitting cond",n))
  }
  save(all_fits, file=paste("asym_conds_fit_",by,".RData",sep=''))
}

fit_asym_conds("model", bySubj=F)
#fit_asym_conds("model", bySubj=T)

#fit_asym_conds("fazly", bySubj=F)
#fit_asym_conds("fazly", bySubj=T)

