library(pso)

#source("fit.R")
#source("test.R")
#source("graphs.R")

order_dir = "orderings/"
model_dir = "models/"

# can load individual trial orderings like so:
ord = read.table(paste(order_dir,"orig_4x4.txt",sep=''))

# or lists of many trial orderings (and some means of human performance):
load(paste(order_dir,"master_orders.RData",sep='')) # orders
load(paste(order_dir,"asymmetric_conditions.RData",sep='')) # conds

make_cooccurrence_matrix <- function(cond, print_matrix=F, heatmap_filename=c()) {
	# makes word x object co-occurrence matrix from training list of words and objects by trial
	# prints a heatmap if filename is specified
	words = cond$train$words
	objs = cond$train$objs
	m = matrix(0, nrow=max(words), ncol=max(objs))
	for(t in 1:nrow(words)) {
		m[words[t,], objs[t,]] = m[words[t,], objs[t,]] + 1
	}
	
	if(print_matrix==T) print(m)
	if(length(heatmap_filename>0)) {
		pdf(paste(heatmap_filename,".pdf",sep="")) 
		heatmap(m, Rowv = NA, Colv = "Rowv",  scale="none", margin=c(3,3), xlab="Object", ylab="Word", col=heat.colors(10)) 
		# labRow=NA, labCol=NA
		dev.off()
	}
	return(m)
}

run_model <- function(cond, model_name, parameters) {
	source(paste(model_dir,model_name,".R",sep=''))
	mod = model(parameters, ord=cond$train)
	return(mod)
}

mod = run_model(conds[["201"]], "kachergis", c(.03,1.08,.964))
mod = run_model(conds[["201"]], "fazly", c(.03,1.08,.964))


coocs3x4 = make_cooccurrence_matrix(conds[["201"]], print_matrix=T, heatmap_filename="201")
filt3e6l = make_cooccurrence_matrix(orders[["filt3E_6L"]])

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
	ord = read.table(paste(order_dir,"orig_order_study_9s-corrected.txt",sep=''))
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

