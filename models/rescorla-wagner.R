# Associative Uncertainty- (Entropy) & Familiarity-Biased Model
# George Kachergis  george.kachergis@gmail.com


update_known <- function(m, tr_w, tr_o) {
  startval = .01
  
  for(i in tr_w) {
    for(c in 1:dim(m)[2]) {
      if(sum(m[,c]>0) & m[i,c]==0) {
        m[i,c] = startval
        m[c,i] = startval
      }
    }
    for(j in tr_o) {
      if(m[i,j]==0) m[i,j] = startval
      if(m[j,i]==0) m[j,i] = startval
    }
  }
  return(m)
}


model <- function(params, ord=c(), reps=1, test_noise=0) {
	lambda <- params[1] # learning rate
	C <- params[2] # decay
	
	voc_sz = max(unlist(ord$words), na.rm=TRUE) # vocabulary size
	ref_sz = max(unlist(ord$objs), na.rm=TRUE) # number of objects
	ppt = length(ord$trials[[1]]$words) # pairs per trial
	traj = list()
	m <- matrix(0, voc_sz, ref_sz) # association matrix
	# training
	for(rep in 1:reps) { # for trajectory experiments, train multiple times
	  for(t in 1:nrow(ord$words)) { 
		#print(format(m, digits=3))
		
		tr_w = as.integer(ord$words[t,])
		tr_w = tr_w[!is.na(tr_w)]
		tr_o = as.integer(ord$objs[t,])
		tr_o = tr_o[!is.na(tr_o)]
		m = update_known(m, tr_w, tr_o) # what's been seen so far?
		
		outcome = 1 
		# objects=cues and words=outcomes, or vice-versa?
		pred = rowSums(m[tr_w,tr_o]) # or colSums...and not only for trial stimuli, but all
		delta = lambda*(outcome - pred)
		m[tr_w,tr_o] = m[tr_w,tr_o] + delta 
		
		m = m*C
		
		index = (rep-1)*length(ord$trials) + t # index for learning trajectory
		traj[[index]] = m
	  }
	}
	m = m + test_noise 
	perf = diag(m) / rowSums(m)
	want = list(perf=perf, matrix=m, traj=traj)
	return(want)
	}

