

model <- function(params, ord=c(), reps=1) {
	lambda <- params[1] # small smoothing factor (1e-5)
	beta <- params[2] # upper bound on number of symbol types to expect? (8500)
	theta <- params[3] # threshold
  
	traj = list()
	voc_sz = max(unlist(ord$words), na.rm=TRUE) # vocabulary size
	ref_sz = max(unlist(ord$objs), na.rm=TRUE) # number of objects
	dummy = voc_sz + 1
	ppt = dim(ord)[2] # pairs per trial
	mean_ent = c()
	m <- matrix(0, voc_sz+1, ref_sz) # association matrix + 1 dummy/nonreferring word -- normalized (every trial)
	assm <- matrix(0, voc_sz+1, ref_sz) # track assoc scores SEPARATELY -- non-normalized
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
    
    align = matrix(0, voc_sz+1, ref_sz)
    for(w in utt) {
      for(f in tr_o) {
        align[w,f] = ((m[w,f]+1e-10) / sum(m[w,]+1e-10)) / sum((m[utt,f]+1e-10) / (m[utt,]+1e-10))
      }
    }
    align = align[utt,tr_o]
		# except for in Eqn 1's denominator, but it must be done everywhere...
		align2 = m[utt,tr_o] / colSums(m[utt,tr_o]) # Eqn 1
		#print(align)
    #print(align2)
		assm[utt,tr_o] = assm[utt,tr_o] + align # new 11/1/12 Eqn 2
		#m[tr_w,tr_o] = (assm[tr_w,tr_o]+lambda) / (rowSums(assm[tr_w,]) + lambda*beta) # # new 11/1/12 Eqn 3
    for(w in tr_w) {
      m[w,tr_o] = (assm[w,tr_o] + lambda) / (sum(assm[w,]) + beta*lambda) # beta*lambda inside sum or not?? prob not...
    }
		
		#comp = diag(m) / rowSums(m) # comprehension score p(oi|wi)
		# assume learned if comp score > theta
		#known = sum(comp>theta)
		index = (rep-1)*length(ord$trials) + t # index for learning trajectory
		traj[[index]] = m
	  }
	  comp = diag(m)[1:voc_sz] / rowSums(m)[1:voc_sz]
	}
	#perf = diag(m)[1:voc_sz] / rowSums(m)[1:voc_sz] # not using theta threshold, but it only converts to binary know/not (can't help)
	mt = m
  mt[which(mt<theta)] = 0
  perf = diag(mt)[1:voc_sz] / (rowSums(mt)[1:voc_sz]+1e-12)
  want = list(perf=perf, matrix=m[1:voc_sz,], traj=traj)
	return(want)
} 
