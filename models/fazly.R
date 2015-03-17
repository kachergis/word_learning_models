# reimplemented following Fazly et al.'s code 
# George Kachergis March 17, 2015

model <- function(params, ord=c(), reps=1) {
	lambda <- params[1] # small smoothing factor (1e-5)
	beta <- params[2] # upper bound on number of symbol types to expect? (8500)
	theta <- params[3] # threshold for knowledge (.7)
	
	# for alignment prob calculation Fazly et al use two extra fixed parameters for smoothing:
  # (lines 174-5 learn.py)
	alpha = 10
	epsilon = 0.001
  
	traj = list()
	voc_sz = max(unlist(ord$words), na.rm=TRUE) # vocabulary size
	ref_sz = max(unlist(ord$objs), na.rm=TRUE) # number of objects
	dummy = voc_sz + 1
	m <- matrix(0, voc_sz+1, ref_sz) # association matrix
  t_unseen = rep(1/beta,voc_sz+1) # in transp - used for updating probs of unseen objs
	alignSum = matrix(0, voc_sz+1, ref_sz)
  
	#assm <- matrix(0, voc_sz+1, ref_sz) # track assoc scores SEPARATELY
	# training
	for(rep in 1:reps) { # for trajectory experiments, train multiple times
	  for(t in 1:nrow(ord$words)) {
	  	#print(format(m, digits=3))
		
		tr_w = as.integer(ord$words[t,])
		tr_w = tr_w[!is.na(tr_w)]
		tr_o = as.integer(ord$objs[t,])
		tr_o = tr_o[!is.na(tr_o)]
		
		# calculate word-referent alignment probabilities
		utt <- c(tr_w,dummy) # add dummy to utterance U(t)
		# following Fazly's updateMappingTables code (learn.py line 168)
		align = matrix(0, voc_sz+1, ref_sz) # align is their 'alignp' (they have alignp.sum and alignp.value)
    for(f in tr_o) {
      sumT = 0
      for(w in utt) {
        if(m[w,f]==0) {  # transp.getValue(w,f): unseen or m[w,f]
          sumT = sumT + t_unseen[w]
        } else {
          sumT = sumT + m[w,f]
        }
      }
      for(w in utt) {
        # m[w,f] = transp.getValue: unseen or m[m,f]
        if(m[w,f]==0) {
          tmp = t_unseen[w]
        } else {
          tmp = m[w,f]
        }
        align[w,f] = tmp + epsilon / (sumT + alpha*epsilon) # m is their 'transp'
        alignSum[w,f] = alignSum[w,f] + align[w,f] # alignp.setValue and alignp.setSum
      }
    }
		# from Fazly code (learn.py lines 193-204)
    # "update translation probs for words in current sentence,
		# taking into account all semantic primitives seen so far"
    for(w in utt) {
      sumA = sum(alignSum[w,]) # for all seen (allprimsD) -- tho alignp.prob was cleared
      denom = sumA + beta*lambda
      for(f in 1:ref_sz) { # again for all observed referents
        m[w,f] = (alignSum[w,f] + lambda) / denom
      }
      t_unseen[w] = lambda / denom
    }
    
		index = (rep-1)*length(ord$trials) + t # index for learning trajectory
		traj[[index]] = m
	  }
	  comp = diag(m) / rowSums(m)[1:voc_sz]
	}
	minprob = 1/beta
  #m[which(m<minprob)] = 0
	perf = diag(m) / rowSums(m)[1:voc_sz] # not using theta threshold, but it only converts to binary know/not (can't help)
	want = list(perf=perf, matrix=m[1:voc_sz,], traj=traj)
	return(want)
} 
