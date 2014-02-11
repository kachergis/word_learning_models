

model <- function(params, ord=c(), ord_name="", reps=1, name="fazly", save_traj=FALSE, print_matrix=FALSE) {
	lambda <- params[1] # small smoothing factor (1e-5)
	beta <- params[2] # upper bound on number of symbol types to expect? (8500)
	theta <- params[3] # threshold for knowledge (.7)
	
	# learning trajectory data
	if(save_traj) {
		traj <- data.frame(Trial=numeric(0), testtype=character(0), x=numeric(0.0))
		ord <- read.table(paste(ord_dir,ord_name,".txt",sep=""), header=T) # ordering
		}
	
	voc_sz = max(unlist(ord), na.rm=TRUE) # vocabulary size
	dummy = voc_sz + 1
	ppt = dim(ord)[2] # pairs per trial
	mean_ent = c()
	m <- matrix(0, voc_sz+1, voc_sz) # association matrix + 1 dummy/nonreferring word
	assm <- matrix(0, voc_sz+1, voc_sz) # new 11/1/12 track assoc scores SEPARATELY
	m[dummy,] = lambda # maybe?? Fazly doesn't specify
	# training
	for(rep in 1:reps) { # for trajectory experiments, train multiple times
	  for(t in 1:length(ord$trials)) {
	  	#print(format(m, digits=3))
		
		tr = as.integer(ord$trials[[t]]$words)
		tr = tr[!is.na(tr)]
		
		# calculate word-referent alignment probabilities
		utt <- c(tr,dummy) # Fazly paper doesn't actually specify adding dummy to U(t)
		# except for in Eqn 1's denominator, but it's gotta be done everywhere...
		align = m[utt,tr] / colSums(m[utt,tr]) # Eqn 1
		#assoc = m[utt,tr] + align # Eqn 2
		assm[utt,tr] = assm[utt,tr] + align # new 11/1/12
		#m[utt,tr] = (assoc+lambda) / (rowSums(m[utt,]) + lambda*beta) # Eqn 3
		m[tr,tr] = (assm[tr,tr]+lambda) / (rowSums(assm[tr,]) + lambda*beta) # # new 11/1/12
		
		#comp = diag(m) / rowSums(m) # comprehension score p(oi|wi)
		# assume learned if comp score > theta
		#known = sum(comp>theta)
		if(print_matrix) print(m)
		if(save_traj & t>1) traj <- rbind(traj, data.frame(Trial=t, testtype="all", x=mean(diag(m))))
	  }
	  comp = diag(m) / rowSums(m)[1:voc_sz]
	  #perf = c(perf, sum(comp>theta)/voc_sz) # mean(comp) # for fazly2
	}
	if(save_traj) save(traj, file=paste(name,"_",ord_name,"_traj.RData",sep=''))
	want = list()
	#want[['tr_ent']] = mean_ent
	#want[['perf']] = perf
	return(m[1:voc_sz,]) # NOT USING THETA RIGHT NOW...BUT IT DOESN'T HELP
} 
