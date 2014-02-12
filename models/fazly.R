

model <- function(params, ord=c(), reps=1) {
	lambda <- params[1] # small smoothing factor (1e-5)
	beta <- params[2] # upper bound on number of symbol types to expect? (8500)
	theta <- params[3] # threshold for knowledge (.7)
	
	traj = list()
	voc_sz = max(unlist(ord$words), na.rm=TRUE) # vocabulary size
	ref_sz = max(unlist(ord$objs), na.rm=TRUE) # number of objects
	dummy = voc_sz + 1
	ppt = dim(ord)[2] # pairs per trial
	mean_ent = c()
	m <- matrix(0, voc_sz+1, ref_sz) # association matrix + 1 dummy/nonreferring word
	assm <- matrix(0, voc_sz+1, ref_sz) # track assoc scores SEPARATELY
	m[dummy,] = lambda # maybe?? unspecified in description
	# training
	for(rep in 1:reps) { # for trajectory experiments, train multiple times
	  for(t in 1:nrow(ord$words)) {
	  	#print(format(m, digits=3))
		
		tr_w = as.integer(ord$words[t,])
		tr_w = tr_w[!is.na(tr_w)]
		tr_o = as.integer(ord$objs[t,])
		tr_o = tr_o[!is.na(tr_o)]
		
		# calculate word-referent alignment probabilities
		utt <- c(tr_w,dummy) # Fazly paper doesn't actually specify adding dummy to U(t)
		# except for in Eqn 1's denominator, but it must be done everywhere...
		align = m[utt,tr_o] / colSums(m[utt,tr_o]) # Eqn 1
		#assoc = m[utt,tr_o] + align # Eqn 2
		assm[utt,tr_o] = assm[utt,tr_o] + align # new 11/1/12
		#m[utt,tr_o] = (assoc+lambda) / (rowSums(m[utt,]) + lambda*beta) # Eqn 3
		m[tr_w,tr_o] = (assm[tr_w,tr_o]+lambda) / (rowSums(assm[tr_w,]) + lambda*beta) # # new 11/1/12
		
		#comp = diag(m) / rowSums(m) # comprehension score p(oi|wi)
		# assume learned if comp score > theta
		#known = sum(comp>theta)
		index = (rep-1)*length(ord$trials) + t # index for learning trajectory
		traj[[index]] = m
	  }
	  comp = diag(m) / rowSums(m)[1:voc_sz]
	}
	perf = diag(m) / rowSums(m)[1:voc_sz] # not using theta threshold, but it only converts to binary know/not (can't help)
	want = list(perf=perf, matrix=m[1:voc_sz,], traj=traj)
	return(want)
} 
