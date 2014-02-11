library(plotrix); library(lattice); library(pso); library(car)
# George Kachergis Dec 16, 2012
# fit associative model using MLE by subject
setwd("~/Dropbox/projects/priming/NEWanalyses")

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


fit_subj <- function(par, orders, sdat) {
	# ords is all of the orders
	# sdat has order file names and responses
	tot_lik = 0
	for(f in unique(sdat$file)) {
		o <- orders[[f]] # $trials ?
		#print(paste("Sub",sdat$Subject[1]," File",f))
    M=NA
    for(rep in 1:max(sdat$Block)) {
      sdatb = subset(sdat, Block==rep)
		  M <- model(par, ord=o, m=M)
		  # need: for each word (CorrectAns), what object they chose (Response)
		  tot_lik = tot_lik + binomial_likelihood(subset(sdatb, file==f), M)
    }
	}
	#mlik = tot_lik / (length(unique(sdat$file))*max(sdat$Block))
  mlik = tot_lik
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
			M=NA
			for(rep in 1:max(sdat$Block)) {
			  sdatb = subset(sdat, Block==rep)
			  M <- model(par, ord=o, m=M)
			  # need: for each word (CorrectAns), what object they chose (Response)
			  tot_lik = tot_lik + binomial_likelihood(subset(sdatb, file==f), M)
			}
			# need: for each word (CorrectAns), what object they chose (Response)
			#tot_lik = tot_lik / (length(unique(sdat$file))*max(sdat$Block))
		}
	}
	mlik = tot_lik / length(unique(dat$Subject))
	#print(par)
	#print(mlik)
	return(-mlik)
}

fit_all_MLE <- function(modeln, PSO=FALSE) {
	source(paste(modeln,".R",sep="")) 
	load("all_trajectory.RData")
  load("master_orders.RData")
	agg = all
	controls = list(maxit=200, max.restart=2, reltol=.05)
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
	load("all_trajectory.RData")
  #all$file = with(all, ifelse(Exp=="4 Pairs/Trial, 6x", "orig_4x4", ifelse(Exp=="3 Pairs/Trial 3,6,9x", "block1_369-3x3loCD", ifelse(Exp=="4 Pairs/Trial 3,6,9x", "3_x8_369_4x4", NA))))
	agg <- all # subset() # just one exp?
	agg$Subject = with(all, ifelse(Exp=="4 Pairs/Trial, 6x", Subject+100, ifelse(Exp=="3 Pairs/Trial 3,6,9x", Subject+200, ifelse(Exp=="4 Pairs/Trial 3,6,9x", Subject+300, NA))))
  load("master_orders.RData")
	controls = list(maxit=250, max.restart=2, reltol=.05)
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
			agg[rows,parcols] = matrix(rep(best$par, length(rows)), nrow=length(rows), byrow=T)
			agg[rows,]$ML = best$value
      M = NA
			for(rep in 1:max(sdat$Block)) {
			  blrows <- with(agg, which(Subject==s & file==f & Block==rep))
			  M = model(best$par, orders[[f]], m=M)
			  agg[blrows,]$Model <- diag(M) / rowSums(M)
			}
		}
		
		print(agg[rows[1],])
	}
	return(agg)
}


traj <- fit_by_subj("model", PSO=T)
save(traj, file="traj_assoc_MLE_ind_fitPSO.RData")
print(paste("assoc ind PSO:",mean(traj$ML))) # 36.9 / 4


# all might not work yet -- check it
#traj <- fit_all_MLE("model", PSO=T)
#save(traj, file="traj_assoc_MLE_group_fitPSO.RData")
#print(paste("assoc all PSO:",full$ML)) # 894.03 / 77 = 11.61

