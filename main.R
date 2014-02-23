library(pso)

#source("fit.R")
source("graphics.R")

order_dir = "orderings/"
model_dir = "models/"

# load lists of many trial orderings (and some means of human performance):
load(paste(order_dir,"master_orders.RData",sep='')) # orders
load(paste(order_dir,"asymmetric_conditions.RData",sep='')) # conds



run_model <- function(cond, model_name, parameters, print_perf=F) {
	require(pso) # or require(DEoptim)
	source(paste(model_dir,model_name,".R",sep=''))
	mod = model(parameters, ord=cond$train)
	if(print_perf) print(mean(mod$perf))
	return(mod)
}

mod = run_model(conds[["201"]], "fazly", c(.0001,8000,.7), print_perf=T)
animate_trajectory(mod)

mod = run_model(conds[["201"]], "kachergis", c(1,3,.97), print_perf=T)
mod = run_model(conds[["201"]], "strength", c(1,.97), print_perf=T)
mod = run_model(conds[["201"]], "uncertainty", c(1,3,.97), print_perf=T)
mod = run_model(conds[["201"]], "novelty", c(1,3,.97), print_perf=T)
mod = run_model(conds[["201"]], "Bayesian_decay", c(.7,1.7,1), print_perf=T)
mod = run_model(conds[["201"]], "rescorla-wagner", c(.7,1.7,1), print_perf=T)

coocs3x4 = make_cooccurrence_matrix(conds[["201"]], print_matrix=T, heatmap_filename="201")
filt3e6l = make_cooccurrence_matrix(orders[["filt3E_6L"]])

fit_model <- function(model_name, orders, par_lower, par_upper) {
	require(pso)
	source(paste(model_dir,model_name,".R",sep=''))
	fits = list()
	startt = Sys.time()
	cat("Order\tSSE\tParameters\n")
	for(i in 1:length(names(orders))) {
		best <- psoptim(c(.1,1,.97), meanSSE, ord=orders[[i]]$train, human_perf=unlist(orders[[i]]$hum_perf), lower=par_lower, upper=par_upper) 
		cat(names(orders)[i],'\t',best$value,'\t',best$par,'\n')
		mod = model(best$par, ord=orders[[i]]$train)
		fits[[names(orders)[i]]] = list(SSE=best$value, par=best$par, perf=mod$perf)
	}
	stopt = Sys.time()
	print(stopt-startt)
	return(fits)
}


meanSSE <- function(par, order, human_perf) {
	mod = model(par, order)
	# if there is an item grouping factor, can first aggregate item perf by the factor
	return((mean(mod$perf)-human_perf)^2)
}


fit_model("kachergis", orders[1], c(.001,.1,.5), c(5,10,1))
fit_model("kachergis", orders[c(1,5,9)], c(.001,.1,.5), c(5,10,1))
fit_model("fazly", orders[c("orig_4x4","orig_3x3")], c(1e-10,5,.1), c(.5,20000,1))
fit_model("Bayesian_decay", orders[c("orig_4x4","orig_3x3")], c(1e-5,1e-5,1e-5), c(10,10,10))

multinomial_likelihood_perfect <- function(par, ord) {
	M = model(par, ord=ord)
	pOgW = diag(M) / rowSums(M) # p(o|w)
	lik = sum(log(pOgW))
	return(-lik) # 18*log(1/18) = -52.02669 for AFC guessing
}

multinomial_likelihood <- function(cdat, M) {
	M = M / rowSums(M) # p(o|w)
	lik = 0
	for(i in 1:dim(cdat)[1]) { 
		wordi = cdat[i,"Word"]
		response_probs = M[wordi,unlist(cdat[i,c("Obj1","Obj2","Obj3","Obj4")])] # strengths of test objects
		resp = cdat[wordi,]$Response  #resp = cdat[which(cdat$Word==i),]$Response
		lik = lik + log(M[wordi,resp] / sum(response_probs)) 
	}
	return(-lik) # minimize -loglik
}


binomial_likelihood <- function(cdat, M) {
	est_prob = diag(M) / rowSums(M) # prob correct
	lik = 0
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
	return(-mlik)
}