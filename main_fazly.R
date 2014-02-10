rm(list=ls()) #clear memory
library(plotrix); library(lattice); library(pso); library(car)
# non-sampling version: no need to simulate many subjects 

source("fit.R")
#source("test.R")
#source("graphs.R")

folder = "orderings/"
exp2 = c("orig_order_study_9s-corrected", "3x3_orig_order", "cont_div6-12") # Exp 2
exp3 = c("block1_3_6_9-3x3-lo_cd", "block2_3_6_9-3x3", "block3_3_6_9_lo_medCD", "block4_369_39mx")

graph_fits <- function(fits, modeln) {
	print(cor.test(fits$x, fits$Model)) # .79
	fits$SSEbl = (fits$x-fits$Model)^2
	avgM = aggregate(Model ~ Block, data=fits, mean)
	sdM = aggregate(Model ~ Block, data=fits,sd)
	avgM$SE = sdM$Model / sqrt(length(unique(fits$Subject))-1)
	names(avgM) = c("Block", "x", "SE")
	avgH = aggregate(x ~ Block, data=fits, mean)
	sdH = aggregate(x ~ Block, data=fits, sd)
	avgH$SE = sdH$x / sqrt(length(unique(fits$Subject))-1)
	avgH$Type = "Human"
	avgM$Type = "Model"
	avg = rbind(avgM, avgH)
	avg$Block <- as.factor(as.character(avg$Block))
	png(paste(modeln,"_traj6fits_average.png",sep=''), width=400, height=400, pointsize=18)
	print(xyplot(x ~ Block, groups=Type, auto.key=T, data=avg, type='l', ylab="Accuracy"))
	dev.off()
	
	fits$M <- "Model"
	fits$H <- "Human"
	gr <- with(fits, rbind(cbind(Subject, Block, H, x), cbind(Subject, Block, M, Model)))
	names(gr) <- c("Subject","Block","Type","Performance")
	gr <- data.frame(gr)
	gr$x <- as.numeric(as.character(gr$x))
	gr$Subject <- as.numeric(gr$Subject)
	png(paste(modeln,"_traj6fits_bysubj.png",sep=''), width=500, height=400, pointsize=18)
	print(xyplot(x ~ Block|H, groups=Subject, data=gr, type='l', ylab="Accuracy"))
	dev.off()
	
	parms = aggregate(cbind(SSE, Model) ~ lambda + beta + theta, data=fits, mean)
	png(paste(modeln,"_traj6fits_parms.png",sep=''), width=600, height=600, pointsize=18)
	print(scatterplotMatrix(~lambda+beta+theta+SSE+Model, data=parms, smooth=F))
	dev.off()
}

fit_by_subj <- function(modeln, PSO=FALSE) {
	source(paste(modeln,".R",sep="")) 
	#load("master_orders.RData") # orders
	ord = read.table(paste(folder,"orig_order_study_9s-corrected.txt",sep=''))
	load("humans/priming_all_trajectory.RData") # all
	agg <- subset(all, Exp=="4 Pairs/Trial, 6x")
	agg <- with(agg, aggregate(Correct, list(Block=Block, Subject=Subject), mean))
	agg$lambda <- NA
	agg$beta <- NA
	agg$theta <- NA
	agg$SSE <- NA
	agg$Model <- NA
	initp = c(.00001,1000,.7)
	lowerp = c(.0000001,100,.69)
	upperp = c(.5,10000, .71)
	
	for(s in unique(agg$Subject)) {
		sdat <- subset(agg, Subject==s)
		if(PSO) {
			best <- psoptim(initp, fit_subj_traj, ord=ord, perf=sdat$x, lower=lowerp, upper=upperp) 
		} else {
			best <- optim(initp, fit_subj_traj, ord=ord, perf=sdat$x, lower=lowerp, upper=upperp, method="L-BFGS-B") # , control=list(parscale=c(10,10,10), maxit=20) 
		}
		rows <- with(sdat, which(agg$Subject==Subject))
		mp = model(best$par, ord, reps=length(rows))
		agg[rows,c("lambda","beta","theta")] = matrix(rep(best$par, length(rows)), nrow=length(rows), byrow=T)
		agg[rows,]$SSE = best$value
		agg[rows,]$Model <- mp
		print(agg[rows,])
	}
	return(agg)
}

traj6fits <- fit_by_subj("fazly", PSO=T) # ,PSO=T
graph_fits(traj6fits, "fazly")

model(c(.00001,8000,.7), ord=ord, reps=4, print_matrix=T)
