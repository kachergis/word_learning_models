library(pso)

source("fit.R")
#source("test.R")
#source("graphs.R")

ord_folder = "orderings/"
ord = read.table(paste(ord_folder,"orig_4x4.txt",sep=''))

load("master_orders.RData") # orders
load("asym_master.RData")
load("asym_conds.RData")


fit_by_subj_hyp <- function(modeln, PSO=FALSE) {
	source(paste(modeln,".R",sep="")) 
	ord = read.table(paste(folder,"orig_order_study_9s-corrected.txt",sep=''))
	load("humans/priming_all_trajectory.RData") # all
	agg <- subset(all, Exp=="4 Pairs/Trial, 6x")
	agg <- with(agg, aggregate(Correct, list(Block=Block, Subject=Subject), mean))
	agg$forget <- NA
	agg$store <- NA
	agg$SSE <- NA
	agg$Model <- NA
	for(s in unique(agg$Subject)) {
		sdat <- subset(agg, Subject==s)
		if(PSO) {
			best <- psoptim(c(.05), fit_subj_traj_samp, ord=ord, perf=sdat$x, lower=c(0), upper=c(1)) 
		} else {
			best <- optim(c(.05,.9), fit_subj_traj_samp, ord=ord, perf=sdat$x, lower=c(0,0), upper=c(1,1), method="L-BFGS-B", control=list(parscale=c(10,10), maxit=30)) 
		}
		rows <- with(sdat, which(agg$Subject==Subject))
		mp = model(best$par, ord, reps=length(rows))
		agg[rows,c("forget","store")] = matrix(rep(best$par, length(rows)), nrow=length(rows), byrow=T)
		agg[rows,]$SSE = best$value
		agg[rows,]$Model <- mp
		print(agg[rows,])
	}
	return(agg)
}


fit_by_subj <- function(modeln, PSO=FALSE, K=1) {
	source(paste(modeln,".R",sep="")) 
	ord = read.table(paste(ord_folder,"orig_order_study_9s-corrected.txt",sep=''))
	load("humans/priming_all_trajectory.RData") # all
	agg <- subset(all, Exp=="4 Pairs/Trial, 6x")
	agg <- with(agg, aggregate(Correct, list(Block=Block, Subject=Subject), mean))
	agg$X <- NA
	agg$B <- NA
	agg$C <- NA
	agg$SSE <- NA
	agg$Model <- NA
	for(s in unique(agg$Subject)) {
		sdat <- subset(agg, Subject==s)
		if(PSO) {
			best <- psoptim(c(.8,2.25,.95), fit_subj_traj_samp, ord=ord, perf=sdat$x, K=K, lower=c(0,1,.8), upper=c(30,7,1)) 
		} else {
			best <- optim(c(.03,1.08,.964), fit_subj_traj_samp, ord=ord, perf=sdat$x, K=K, lower=c(0,1,.8), upper=c(30,7,1), method="L-BFGS-B", control=list(parscale=c(10,10,10), maxit=20)) 
		}
		rows <- with(sdat, which(agg$Subject==Subject))
		mp = model(best$par, ord, reps=length(rows))
		agg[rows,c("X","B","C")] = matrix(rep(best$par, length(rows)), nrow=length(rows), byrow=T)
		agg[rows,]$SSE = best$value
		agg[rows,]$Model <- mp
		print(agg[rows,])
	}
	return(agg)
}

fit_sampling_model <- function(exp, files, Nsubj, PSO=FALSE) {
	m = "model_sampling"
	source(paste(m,".R",sep="")) 
	ords = list()
	for(f in files) {
		ords[[f]] <- read.table(paste(folder,f,".txt",sep=""), header=T)
	}
	if(PSO) {
		best <- psoptim(c(.8,2.25,.95), fit_generic, ords=ords, Nsubj=Nsubj, exp=exp, lower=c(0,1,.8), upper=c(30,7,1), control=list(parscale=c(10,10,10))) 
	} else {
		best <- optim(c(.8,2.25,.95), fit_generic, ords=ords, Nsubj=Nsubj, exp=exp, lower=c(0,1,.8), upper=c(30,7,1), method="L-BFGS-B", control=list(parscale=c(10,10,10))) 
	}
	return(best)
}

