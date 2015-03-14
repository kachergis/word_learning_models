# Associative Uncertainty- (Novelty) & Familiarity-Biased Model
# George Kachergis  george.kachergis@gmail.com


update_known <- function(m, tr_w, tr_o) {
  startval = .01
  
  for(i in tr_w) {
    for(c in 1:dim(m)[2]) {
      if(sum(m[,c]>0) & m[i,c]==0) {
        m[i,c] = startval
        if(c<nrow(m)) m[c,i] = startval
      }
    }
    for(j in tr_o) {
      if(m[i,j]==0) m[i,j] = startval
      if(j<nrow(m) && m[j,i]==0) m[j,i] = startval
    }
  }
  return(m)
}


model <- function(params, ord=c(), reps=1, test_noise=0) {
	X <- params[1] # associative weight to distribute
	B <- params[2] # weighting of uncertainty vs. familiarity
	C <- params[3] # decay
	
	voc_sz = max(unlist(ord$words), na.rm=TRUE) # vocabulary size
	ref_sz = max(unlist(ord$objs), na.rm=TRUE) # number of objects
	freq_w = rep(0,voc_sz) # freq[i] = times word i has appeared
	freq_o = rep(0,ref_sz)
	traj = list()
	m <- matrix(0, voc_sz, ref_sz) # association matrix
	perf = matrix(0, reps, voc_sz) # a row for each block
	# training
	for(rep in 1:reps) { # for trajectory experiments, train multiple times
	  for(t in 1:nrow(ord$words)) { 
		#print(format(m, digits=3))
		
		tr_w = as.integer(ord$words[t,])
		tr_w = tr_w[!is.na(tr_w)]
		tr_o = as.integer(ord$objs[t,])
		tr_o = tr_o[!is.na(tr_o)]
		m = update_known(m, tr_w, tr_o) # what's been seen so far?
		freq_w[tr_w] = freq_w[tr_w] + 1
		freq_o[tr_o] = freq_o[tr_o] + 1
		novelty_w = 1/(1+freq_o[tr_w])
		novelty_o = 1/(1+freq_o[tr_o])
		
		novelty_w = exp(B*novelty_w)
		novelty_o = exp(B*novelty_o)
		nov = (novelty_w %*% t(novelty_o))

		# get all current w,o strengths and normalize to distr X
		assocs = m[tr_w,tr_o]
		denom = sum(assocs * nov)
		m = m*C # decay everything
		# update associations on this trial
		m[tr_w,tr_o] = m[tr_w,tr_o] + (X * assocs * (novelty_w %*% t(novelty_o))) / denom 

		index = (rep-1)*length(ord$trials) + t # index for learning trajectory
		traj[[index]] = m
	  }
	m_test = m+test_noise # test noise constant k
	perf[rep,] = diag(m_test) / rowSums(m_test)
	}
	
	want = list(perf=perf, matrix=m, traj=traj)
	return(want)
	}

