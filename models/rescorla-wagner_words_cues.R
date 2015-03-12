# Rescorla-Wagner (1972) associative learning model 
# adapted for cross-situational word learning by
# George Kachergis  george.kachergis@gmail.com


model <- function(params, ord=c(), reps=1, test_noise=0) {
	beta = params[1] # learning rate
	C = params[2] # decay
	lambda = params[3] # maximum associative value that a CS can achieve
	alpha = 1 # salience (fix at 1 unless manipulated)
	
	voc_sz = max(unlist(ord$words), na.rm=TRUE) # vocabulary size
	ref_sz = max(unlist(ord$objs), na.rm=TRUE) # number of objects
	ppt = length(ord$trials[[1]]$words) # pairs per trial
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
	  
		# if words are cues, use rowSums--but only of the currently-presented stimuli
    if(length(tr_o)>1 & length(tr_w)>1) {
		  pred = rowSums(m[tr_w,tr_o])
    } else {
      pred = m[tr_w,tr_o]
    }
		delta = alpha*beta*(lambda - pred)
		m[tr_w,tr_o] = m[tr_w,tr_o] + delta 
		
		m = m*C
		
		index = (rep-1)*length(ord$trials) + t # index for learning trajectory
		traj[[index]] = m
	  }
	m_test = m+test_noise # test noise constant k
	perf[rep,] = diag(m_test) / rowSums(m_test)
	}
	want = list(perf=perf, matrix=m, traj=traj)
	return(want)
	}

