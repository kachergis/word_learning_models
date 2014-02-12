# 2.26.2009 - Bayesian model

# cross-situational statistical word learning models
# 
# probabilistic choice: sample(vector, 1, replace=TRUE, prob=vector)
# optim(par, fn)

folder = "orderings/"
files = c("18-temp_cont_1olap2tr", "19-temp_cont_1olap3tr", "20-temp_cont_2olap2tr", "orig_order_study_9s-corrected")
NSUBJECTS = 200
VOCAB_SZ = 18
initK = 1.11800808949661; initA = .5; initAp = .08
INIT_PAR = c(initA, initAp)

norm <- function(x) {
	tot = sum(x)
	y <- c()
	for(i in x){ y=c(y,i/tot) } 
	}

fit <- function(par) {
	sse = 0 # sum of squared error
	for(f in files) {
		o <- read.table(paste(folder,f,".txt",sep=""), header=T) # ordering
		ans <- c() # rows = Ss, cols = item correctness
		for (s in 1:NSUBJECTS) {
			s <- model(o, par, VOCAB_SZ)
			ans = rbind(ans, s)
		}
		# compare simulation avg to each condition's human avg
		if(f=="18-temp_cont_1olap2tr") { 
			contig <- c( mean(ans[,1]), mean(ans[,2]), mean(ans[,3]), mean(ans[,7]), mean(ans[,8]), mean(ans[,11]), mean(ans[,14]), mean(ans[,15]), mean(ans[,18]))
			discontig <- c( mean(ans[,4]), mean(ans[,5]), mean(ans[,6]), mean(ans[,9]), mean(ans[,10]), mean(ans[,12]), mean(ans[,13]), mean(ans[,16]), mean(ans[,17]))
			sse = sse + ( mean(contig)-.47 )^2 + ( mean(discontig)-.27 )^2
			print(paste(f,"contig:",mean(contig),", discontig:",mean(discontig)))
		} else if(f=="19-temp_cont_1olap3tr") {
			contig <- c( mean(ans[,1]), mean(ans[,3]), mean(ans[,5]), mean(ans[,6]), mean(ans[,7]), mean(ans[,11]), mean(ans[,13]), mean(ans[,14]), mean(ans[,16]))
			discontig <- c( mean(ans[,2]), mean(ans[,4]), mean(ans[,8]), mean(ans[,9]), mean(ans[,10]), mean(ans[,12]), mean(ans[,15]), mean(ans[,17]), mean(ans[,18]))
			sse = sse + ( mean(contig)-.40 )^2 + ( mean(discontig)-.36 )^2
			print(paste(f,"contig:",mean(contig),", discontig:",mean(discontig)))
		} else if(f=="20-temp_cont_2olap2tr") {
			sse = sse + ( mean(ans)-.33 )^2
			print(paste(f,"contig:",mean(ans)))
		} else { # original
			sse = sse + ( mean(ans)-.29 )^2
			print(paste(f,"discontig:",mean(ans)))
		}
	}
	print(paste("params: ",par))
	print(paste("SSE: ",sse))
	return(sse)
}

model <- function(ord, params, voc_sz) {
	K <- initK # assoc units per trial
	m <- matrix(1/voc_sz, voc_sz, voc_sz)
	A <-  params[1]
	Ap <- params[2]
	
	# training
	for(t in 1:dim(ord)[1]) { 
		# temporal contiguity: separate cont & discont
		if(t>0) { 
			contig = intersect(ord[t,], ord[t-1,]) 
			discontig = setdiff(ord[t,], ord[t-1,])
			}
		else { contig = c(); discontig = ord[t,] }
		
		
		# apportion K to each possible contig & discontig association
		if(length(contig)>0) { 
			Kc = (thisA * K) / length(contig) 
			Kd = ((1-thisA) * K) / length(discontig)
		} else { Kd = K / length(discontig) }
		
		# shuffle words
		if(length(contig)>0) { contig <- sample(contig, length(contig)) }
		discontig <- sample(discontig, length(discontig))
		
		for(w in contig) {
			strengths = c()
			for (p in contig) { strengths = c(strengths, m[w,p]) } # get assocs
			x <- sample(contig, 1, replace=FALSE, prob=strengths) # 
			for (obj in x) { m[w,obj] = m[w,obj] + Kc }
		}
		
		for(w in discontig) {
			strengths = c()
			for (p in discontig) { strengths = c(strengths, m[w,p]) } # get assocs
			x <- sample(discontig, 1, replace=FALSE, prob=strengths) # 
			for (obj in x) { m[w,obj] = m[w,obj] + Kd } 
		}
		
	}	
	
	# testing
	responses = c()
	correct = c()
	for(w in 1:voc_sz) {
		r <- sample(1:18, 1, prob=m[w,])
		responses = c(responses, r)
		correct = c(correct, r==w)
		}
	return(correct)
	}


#optim(INIT_PAR, fit)
optim(INIT_PAR, fit, method="L-BFGS-B", lower=c(0,0), upper=c(1,1), control=list(ndeps=c(.1,.1), lmm=50))
#optimize(fit, interval=c(.01,100)) # for one parameter
#optgrid(fit, INIT_PAR, incr, lower, upper, verbose=1, ...)