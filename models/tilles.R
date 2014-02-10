# Tilles and Fontanari's (2013) cross-situational word learning model
# George Kachergis  george.kachergis@gmail.com  September 16, 2013

shannon.entropy <- function(p) {
	if (min(p) < 0 || sum(p) <= 0)
		return(NA)
	p.norm <- p[p>0]/sum(p)
	-sum(log2(p.norm)*p.norm)
	}


model <- function(params, ord=c(), verbose=FALSE) {
	X <- params[1] # reinforcement parameter for current stimuli
	B <- params[2] # inference parameter, regulating ME and prior info integration
	alpha_0 <- params[3] # baseline efficiency corresponding to maximal uncertainty about referent of target word
	# all parms in [0,1]
	# B=1 enforces ME: new word is associated with equal prob to any available new objects
	# B=0 forces new word to be associated to all objects already seen (including current)
	
	voc_sz = max(unlist(ord), na.rm=TRUE) # vocabulary size
	novel = rep(TRUE,voc_sz) # set novel[i]=0 once i has appeared
	m <- matrix(0, voc_sz, voc_sz) # association matrix - probabilities
	
	# training
	tr = as.integer(ord[1,])
	ppt = length(tr) # pairs per trial
	m[tr,tr] = 1/ppt # eq 1
	novel[tr] = FALSE
	for(t in 2:dim(ord)[1]) { 
    tr = as.integer(ord[t,])
		
		Nprev = length(which(novel==FALSE)) # N_{t-1} for alpha_{t-1}
		
		cur_novel = tr[which(novel[tr])] # ~ = current novel stimuli
		cur_old = tr[which(!novel[tr])] # current trial's old stimuli
		other_old = setdiff(which(!novel), cur_old) # stimuli appearing before, but not currently
		novel[tr] = FALSE
		
		N = length(which(novel==FALSE)) # number of distinct stimuli observed so far: N_t
		
		# preprocessing of novel stimuli
		if(length(cur_novel)!=0) {
			m[cur_novel,cur_novel] = B/length(cur_novel) + (1-B)/N
			m[cur_novel,cur_old] = (1-B) / N
			m[cur_novel,other_old] = (1-B) / N
		}
		
		# alpha_{t-1} -- must be previous trial!
		alpha = rep(0,voc_sz) # length(tr)
		for(w in 1:voc_sz) {
			alpha[w] = alpha_0 + (1-alpha_0) * (1-shannon.entropy(m[w,])/log(Nprev)) # eq 7
		}
		
		# now for all stimuli on the current trial...
		# net gain of confidence for association w_i and o_j:
		
		flux = rep(0,voc_sz)
		r <- matrix(0, voc_sz, voc_sz)
		
		# slow iterative method; rowSums(m) = 1
		for(w in tr) {
			flux[w] = sum(m[w,other_old]) / sum(m[w,tr]) # eq 5.1
			r[w,tr] = X*m[w,tr]*flux[w] # eq 5.2
			m[w,tr] = m[w,tr] + alpha[w]*r[w,tr] + (1-alpha[w])*(1/N - m[w,tr]) # eq 8
			m[w,other_old] = m[w,other_old] - alpha[w]*X*m[w,other_old] + (1-alpha[w])*(1/N - m[w,other_old]) # eq 9
		}
    	
    if(verbose) {
    	print(paste("Trial",t))
		  print(rowSums(m))
    }

		# X=1 transfers confidences for absent objects to the currently ones
		# model also assumes that the more certain an association is, the more
		# efficiently that information should be used for reinforcement
		mt = m
		# now update the other_old words
    # if params[3] < .268 on "block2_369-3x3hiCD" m[18,18] becomes negative at trial 9
		flux = rep(0,voc_sz)
		r <- matrix(0, voc_sz, voc_sz)
		flux[other_old] = rowSums(m[other_old,tr]) / rowSums(m[other_old,other_old]) # eq 11.1
		r[other_old,other_old] = B*mt[other_old,other_old]*flux[other_old] # eq 11.2
		m[other_old,other_old] = mt[other_old,other_old] + alpha[other_old]*r[other_old,other_old] + (1-alpha[other_old])*(1/N - mt[other_old,other_old]) # eq 12
		m[other_old,tr] = mt[other_old,tr] - alpha[other_old]*B*mt[other_old,tr] + (1-alpha[other_old])*(1/N - mt[other_old,tr]) # eq 13
		
		if(verbose) {
			print(rowSums(m)) # should be 1...
			print(m)
		}
		
	}
	
	return(m)
	}

ord = read.table("block2_3_6_9-3x3.txt", header=F)

print(model(c(.6, .8, .85), ord)) # works

print(model(c(.6, .8, .23), ord, verbose=TRUE)) # doesn't work
print(model(c(1, 1, .25), ord, verbose=TRUE)) # doesn't work