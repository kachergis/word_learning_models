rm(list=ls()) #clear memory
library(plotrix); library(lattice); library(pso); library(car)
# non-sampling version: no need to simulate many subjects 

source("fit.R")
#source("test.R")
#source("graphs.R")

folder = "orderings/"
exp2 = c("orig_order_study_9s-corrected", "3x3_orig_order", "cont_div6-12") # Exp 2
exp3 = c("block1_3_6_9-3x3-lo_cd", "block2_3_6_9-3x3", "block3_3_6_9_lo_medCD", "block4_369_39mx")
NSUBJECTS = 1000

# sampling version of associative model
fit_sampling_model("Exp2", exp2, NSUBJECTS)
fit_sampling_model("Exp3", exp3, NSUBJECTS, PSO=TRUE)

traj6fits <- fit_by_subj("model")
traj6fits_sampK1 <- fit_by_subj("model_sampling", K=1)
save(traj6fits_sampK1, file="traj6fits_sampK1.RData")
traj6fits_sampK2 <- fit_by_subj("model_sampling", K=2)
save(traj6fits_sampK2, file="traj6fits_sampK2.RData")
traj6fits_sampK3 <- fit_by_subj("model_sampling", K=3)
save(traj6fits_sampK3, file="traj6fits_sampK3.RData")

# incomplete mess of a function (run internally!) 5-9-12
graph_trajs <- function(fits) {
	fits$SSEbl = (fits$x-fits$Model)^2
	cor.test(fits$x, fits$Model)
	#with(fits, aggregate(SSEbl, list(Block=Block), mean)) # fits block 1 best, but all are good
	#parms = with(fits, aggregate(SSE, list(Subject=Subject, X=X, B=B, C=C), mean))
	parms = aggregate(cbind(SSE, Model) ~ X + B + C, data=fits, mean)
	#parms = aggregate(cbind(SSE, Model) ~ Subject + f, data=fits, mean) # for hyp
	cor(parms)
	plot(parms)
	scatterplotMatrix(~f+SSE+Model, data=parms, smooth=F)
	cor.test(parms$X, parms$B) # .69
	cor.test(parms$X, parms$C) # nope
	cor.test(parms$C, parms$B) # -.58 
	# no sig correlation of parms with 
	cor.test(parms$X, parms$SSE)
	cor.test(parms$X, parms$Model)
	
	avgM = aggregate(Model ~ Block, data=fits, mean)
	sdM = aggregate(Model ~ Block, data=fits, sd)
	avgM$SE = sdM$Model / sqrt(length(unique(parms$Subject))-1)
	names(avgM) = c("Block", "x", "SE")
	avgH = aggregate(x ~ Block, data=fits, mean)
	sdH = aggregate(x ~ Block, data=fits, sd)
	avgH$SE = sdH$x / sqrt(length(unique(parms$Subject))-1)
	avgH$Type = "Human"
	avgM$Type = "Model"
	avg = rbind(avgM, avgH)
	avg$Block <- as.factor(as.character(avg$Block))
	xyplot(x ~ Block, groups=Type, auto.key=T, data=avg, type='l', ylab="Accuracy")
	
	fits$M <- "Model"
	fits$H <- "Human"
	gr <- with(fits, rbind(cbind(Subject, Block, H, x), cbind(Subject, Block, M, Model)))
	names(gr) <- c("Subject","Block","Type","Performance")
	gr <- data.frame(gr)
	gr$x <- as.numeric(as.character(gr$x))
	gr$Subject <- as.numeric(gr$Subject)
	
	xyplot(x ~ Block|H, groups=Subject, data=gr, type='l', ylab="Accuracy")
	cor(fits$x, fits$Model) # .977
	
	hist(fits$SSEbl, main='', xlab="Block SSE")
	
	fits1 <- traj6fits_sampK1
	fits1$SSEbl = (fits1$x-fits1$Model)^2
	parms1 = aggregate(cbind(SSE, Model) ~ X + B + C, data=fits1, mean)
	plot(parms1)
	cor(fits1$x, fits1$Model) # .577
}

#traj6fits_hyp <- fit_by_subj_hyp("hypoth_model")
load("traj_4x4_fits_Medina_hyp1000Ss.RData")

traj6fits_hyp2par <- fit_by_subj_hyp("hypoth_model")
save(traj6fits_hyp2par, file="traj_4x4_fits_Medina_v2_hyp1000Ss.RData")

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
	#load("humans/freqCD_data_exps8_9_10_12.Rdata") # add ord file name per Cond??
	# table(raw$Cond, raw$Exp)
	#uni = subset(raw, Cond=="3x3" | Cond=="4x4")
	#tri = subset(raw, Cond=="Varied Freq High CD" | "Varied CD" | "3 Subsets" | "Varied Freq Low CD")
	#by_s <- with(raw, aggregate(Correct, list(Exp=Exp, Cond=Cond, Freq=Freq, CD=CD, PairsPerTrial=PairsPerTrial, Subject=Subject), mean))
	#by_s <- with(uni, aggregate(Correct, list(Exp=Exp, Cond=Cond, PairsPerTrial=PairsPerTrial, Subject=Subject), mean))
	
	source(paste(modeln,".R",sep="")) 
	#load("master_orders.RData") # orders
	ord = read.table(paste(folder,"orig_order_study_9s-corrected.txt",sep=''))
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

dist_by_k = show_perf_distro("orig_order_study_9s-corrected", c(0.80, 2.24, 0.95), Nsubj=10000)
dist_by_k$Samples = ifelse(dist_by_k$K==1, paste(dist_by_k$K,"Sample"), paste(dist_by_k$K,"Samples"))
histogram(~Performance | Samples, data=dist_by_k, layout=c(4,1), strip=TRUE, nint=15)

ks.test(subset(test10000, K==1)$Performance, subset(test10000, K==2)$Performance)

show_perf_distro <- function(fname, parms, Nsubj=1000) {
	source("model_sampling.R")
	ord = read.table(paste(folder,fname,".txt",sep=""), header=T)
	
	krange = c(1,2,3,4)
	ans <- c() # rows = Ss, cols = item correctness
	mans <- data.frame(K=character(0), Performance=numeric(0))
	for(k in krange) {
		for (s in 1:Nsubj) {
			s <- model(parms, ord, K=k)
			ans = rbind(ans, s)
			mans = rbind(mans, data.frame(K=k, Performance=mean(s)))
		}
	}
	print(histogram(~Performance | K, data=mans, layout=c(4,1)))
	return(mans)
}


m = "model"
source(paste(m,".R",sep="")) 
ords = list()
for(f in files) {
	ords[[f]] <- read.table(paste(folder,f,".txt",sep=""), header=T)
}
best <- psoptim(c(.8,2.25,.95), fit, ords=ords, lower=c(0,1,.8), upper=c(30,7,1))


# Exp 3: CD and Frequency
best <- fit(c(.2146, 6.931, .982)) # Grow M SSE=.064 1/100 startval
#best <- fit(c(.2146, 5.0135, .974)) # Grow M SSE=.10 1/100 startval


#write.table(best, file=paste("models/",m,"_params.txt", sep=""), append=T)
	
# what to graph
#results <- read.table("tmp_model_output.txt")
#graph(results, m) 
#write.table(results, paste("models/",m,"_output.txt", sep=""))
	
