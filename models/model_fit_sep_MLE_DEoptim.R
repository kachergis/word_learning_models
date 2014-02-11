# try model on Umay's child XSL - fit via Maximum Likelihood Estimation (by ind) 
# August 13, 2013
library(DEoptim)
save_ords <- function() {
	ords = list()
	ords$TR1 = as.matrix(read.table("orderings/TR1-hiCD.txt")) 
	ords$TR2 = as.matrix(read.table("orderings/TR2-hiCD.txt")) 
	ords$TR3 = as.matrix(read.table("orderings/TR3-medCD.txt")) 
	ords$TR4 = as.matrix(read.table("orderings/TR4-medCD.txt"))
	ords$TR5 = as.matrix(read.table("orderings/TR5-medCD.txt")) 
	ords$TR6 = as.matrix(read.table("orderings/TR6-medCD.txt")) 
	ords$TR7 = as.matrix(read.table("orderings/TR7-lowCD.txt"))
	ords$TR8 = as.matrix(read.table("orderings/TR8-lowCD.txt"))
	ords$T1noCD = as.matrix(read.table("orderings/T1_noCD.txt"))
	ords$T2noCD = as.matrix(read.table("orderings/T2_noCD.txt"))
	save(ords, file="child_xsl_2x2_orders.RData")
}
#save_ords()
load("child_xsl_2x2_orders.RData")
source("model.R") 
INIT_PAR = c(0.227, 1.176, 0.963) # filtering: SSE=0.177
# from kids data averages: INIT_PAR = c(.179, 5.16, 0.647)

dat1 = read.table("data/Exp1_data.txt", header=TRUE) # Exp1 has TR1-TR8, but 4AFC tests included only non-cooc objects
dat2 = read.table("data/Exp2_data.txt", header=TRUE) # Exp2 has only TR3 & TR4 (med CD), and tests included cooc competitors
dat2$TrainOrd <- with(dat2, ifelse(TrainOrd=="TR1", "TR3", "TR4")) # Umay says TR1 & TR2 from Exp2 are TR3 & TR4 from Exp1
dat1$CD <- with(dat1, ifelse(TrainOrd=="TR1" | TrainOrd=="TR2", "High", ifelse(TrainOrd=="TR7" | TrainOrd=="TR8", "Low", "Medium")))


# look at performance for different orderings with the same parameters (but different test trials--importantly?)
compare_orders_same_parms <- function(pars, ords, acc, order_names) {
  for(o in order_names) {
    op <- which(acc$TrainOrd==o)
    mp = model(fit$par, ord=ords[[o]])
    for(i in 1:length(op)) {
      wordi = acc[op[i],]$Word
      response_probs = mp[wordi,unlist(acc[op[i],c("Obj1","Obj2","Obj3","Obj4")])]
      acc[op[i],]$Model = mp[wordi,wordi] / sum(response_probs)
    }
  }
}

#compare_orders_same_parms(INIT_PAR, ords, c("TR3","TR4","TR5","TR6"))

evaluate_model_perf <- function(pars, ord) {
  mp = model(pars, ord=ord)
  perf = diag(mp) / rowSums(mp)
  print(format(c(perf, mean(perf)), digits=3))
}

# avg of best fitting parms for dat1: PAR = c(0.1766065, 6.2574313, 0.6717503)

#evaluate_model_perf(INIT_PAR, ords[["T1noCD"]]) # .4
#evaluate_model_perf(INIT_PAR, ords[["T2noCD"]]) # .4
#evaluate_model_perf(INIT_PAR, ords[["TR3"]]) # .63
#evaluate_model_perf(INIT_PAR, ords[["TR4"]]) # .64
#evaluate_model_perf(INIT_PAR, ords[["TR5"]]) # .59
#evaluate_model_perf(INIT_PAR, ords[["TR6"]]) # .59


#model(INIT_PAR, ord=ord1)

fit_group <- function(acc) {
	acc = acc[order(acc$Subject, acc$Word),]
	acc$Model <- NA
	acc$MLL <- NA
	acc$X <- NA
	acc$B <- NA
	acc$alpha <- NA
	lower = c(.01, .01, .3) # was .001 -- but that's lower than the 'new' item value; lots of low-perf Ss got .001
	upper = c(10,    15,  1)
	#controls = list(maxit=500, max.restart=5, reltol=.0001) # mean(fit1$MLL) = -124.81   mean(fit2$MLL) = -155.30
	controls = list(reltol=.00001, steptol=150, itermax=5000, trace=250) 
	set.seed(1234)
	start <- Sys.time()
	for(o in unique(acc$TrainOrd)) {
		op <- which(acc$TrainOrd==o)
		ord = ords[[o]]
		perf = acc[op,c("TrainOrd","Trial","Correct","Word","Obj1","Obj2","Obj3","Obj4","Response")]
		fit <- DEoptim(fn=multinomial_likelihood, lower=lower, upper=upper, ord=ord, cdat=perf, control=controls)
		MLL = fit$optim$bestval
		bestpar = fit$optim$bestmem
		acc[op,]$MLL = -MLL
		acc[op,]$X = bestpar[1]
		acc[op,]$B = bestpar[2]
		acc[op,]$alpha = bestpar[3]
		
		mp = model(bestpar, ord=ord)
		for(i in 1:length(op)) {
			wordi = acc[op[i],]$Word
			response_probs = mp[wordi,unlist(acc[op[i],c("Obj1","Obj2","Obj3","Obj4")])]
			acc[op[i],]$Model = mp[wordi,wordi] / sum(response_probs)
		}
	}
	#print(acc)
	end <- Sys.time()
	print(end-start)
	return(acc)
}

fit_by_subject <- function(acc_s) {
	acc_s = acc_s[order(acc_s$Subject, acc_s$Word),] # acc_s$Trial - if modeling test-trial order
	Ss <- unique(acc_s$Subject)
	acc_s$Model <- NA
	acc_s$MLL <- NA
	acc_s$X <- NA
	acc_s$B <- NA
	acc_s$alpha <- NA
	set.seed(1234)
	start <- Sys.time()
	for(s in Ss) {
		subp <- which(acc_s$Subject==s)
		ord = ords[[unlist(acc_s[subp,]$TrainOrd[1])]]
		# fit the subject
		lower = c(.01, .01, .3) # was .001 -- but that's lower than the 'new' item value; lots of low-perf Ss got .001
		upper = c(10,  15,  1)
		np = 30
		initpop <- cbind(runif(np,lower[1],upper[1]),runif(np,lower[2],upper[2]), runif(np,lower[3],upper[3]))
		controls = list(reltol=.00001, steptol=150, itermax=5000, trace=250, initial=initpop) 
		perf = acc_s[subp,c("TrainOrd","Trial","Correct","Word","Obj1","Obj2","Obj3","Obj4","Response")]
		fit <- DEoptim(fn=multinomial_likelihood, lower=lower, upper=upper, ord=ord, cdat=perf, control=controls)
		MLL = fit$optim$bestval
		bestpar = fit$optim$bestmem
		acc_s[subp,]$MLL = -MLL
		acc_s[subp,]$X = bestpar[1]
		acc_s[subp,]$B = bestpar[2]
		acc_s[subp,]$alpha = bestpar[3]
		
		mp = model(bestpar, ord=ord)
		for(i in 1:length(subp)) {
			response_probs = mp[i,unlist(acc_s[subp[i],c("Obj1","Obj2","Obj3","Obj4")])]
			acc_s[subp[i],]$Model = mp[i,i] / sum(response_probs)
		}
		print(acc_s[subp,])
	}
	end <- Sys.time()
	print(end-start)
	return(acc_s)
}

multinomial_likelihood <- function(par, ord, cdat) {
	M <- model(unlist(par), ord=ord)
	#M = M / rowSums(M) # p(o|w)
	lik = 0
	for(i in 1:dim(cdat)[1]) { 
		wordi = cdat[i,"Word"]
		response_probs = M[wordi,unlist(cdat[i,c("Obj1","Obj2","Obj3","Obj4")])] # strengths of test objects
		resp = cdat[wordi,]$Response  #resp = cdat[which(cdat$Word==i),]$Response
		lik = lik + log(M[wordi,resp] / sum(response_probs)) 
	}
	return(-lik) # minimize -loglik
}


binomial_likelihood <- function(par, ord, cdat) {
	M <- model(unlist(par), ord=ord)
	est_prob = diag(M)
	lik = 0
	for(i in 1:dim(cdat)[1]) { 
		resp = cdat[which(cdat$Word==i),]$Response
		response_probs = M[i,unlist(cdat[i,c("Obj1","Obj2","Obj3","Obj4")])] # strengths of test objects
		if(resp==i) {
			lik = lik + log(est_prob[i] / sum(response_probs))
		} else {
			lik = lik + log(1-(est_prob[i] / sum(response_probs)))
		}
	}
	return(-lik) # 18*log(1/18) = -52.02669 for guessing
}


# if want to do group binomial fit, gotta fix binomial_likelihood for that (resp for many Ss)
#sfit1 <- fit_by_subject(dat1)
#save(sfit1, file="sfit_Exp1_binom_DEoptim.RData")
#sfit2 <- fit_by_subject(dat2)
#save(sfit2, file="sfit_Exp2_binom_DEoptim.RData")

#fit1 <- fit_group(dat1)
#save(fit1, file="groupfit_Exp1_multinom_DEoptim.RData")
#fit2 <- fit_group(dat2)
#save(fit2, file="groupfit_Exp2_multinom_DEoptim.RData")

sfit1 <- fit_by_subject(dat1)
save(sfit1, file="sfit_Exp1_multinom_DEoptim.RData")
sfit2 <- fit_by_subject(dat2)
save(sfit2, file="sfit_Exp2_multinom_DEoptim.RData")

# re-scale parameters:
# [-inf,inf] -> [0,1] : 1/(1+exp(-alpha))
# x>=0: exp(x) or x^2 or log(1+exp(x))

INIT_PAR = colMeans(fit1[,c("X","B","alpha")])
colMeans(fit2[,c("X","B","alpha")])