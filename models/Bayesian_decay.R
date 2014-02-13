# Bayesian model of cross situational learning
# originally conceived by Stephen Denton, Apr. 20, 2010

# Define a simple likelihood function that only updates appropriate rows and columns
likelihoodFun = function(words, objs, alpha, delta) {
	## With alpha=0, this enforces a mutual exclusivity constraint...
	mat = outer(words,objs) + outer(!words,!objs)
	## Can relax mutual exclusivity contraint by increasing alpha or allowing for one of the following
	## 1) Objects can have multiple words applied to them
	# mat = outer(words,objs) + outer(!words,ones)
	## 2) Words can refer to multiple objects
	# mat = outer(words,objs) + outer(ones,!objs)
	likelihood = alpha^(1-mat) + (delta-1) * outer(words,objs)
	return(likelihood)
}

model <- function(params, ord=c(), reps=1, verbose=F) {
	# Define noise probability (when alpha = 0, Bayesian model is a deterministic 'ideal observer')
	alpha <- params[1] # 0.1, 0.5, 0.9 decay for word/object non-co-occurrence
	delta <- params[2] # Multiplier for word/object co-occurrence (set to 1 for no increase)
	chDec <- params[3] # Decision parameter (e.g., 1)

	voc_sz = max(unlist(ord$words), na.rm=TRUE) # vocabulary size
	ref_sz = max(unlist(ord$objs), na.rm=TRUE) # number of objects
	traj = list()
	
	ones = rep(1, voc_sz)
	pWgO <- matrix(1/ref_sz, voc_sz, ref_sz) # prob(word|object) matrix
	pW_O <- matrix(1/(ref_sz*voc_sz), voc_sz, ref_sz)
	# training
	for(rep in 1:reps) { # for trajectory experiments, train multiple times
	  for(t in 1:nrow(ord$words)) { 
		tr_w = as.integer(ord$words[t,])
		tr_w = tr_w[!is.na(tr_w)]
		tr_o = as.integer(ord$objs[t,])
		tr_o = tr_o[!is.na(tr_o)]
		
		words_tr = rep(0,voc_sz)
		objects_tr = rep(0,ref_sz)
		words_tr[tr_w] = 1
		objects_tr[tr_o] = 1
		
		likelihood = likelihoodFun(words_tr, objects_tr, alpha, delta)
		pWgO = likelihood * pWgO
		pWgO = pWgO/outer(ones,colSums(pWgO)) # rowSums(pWgO)
		pW_O = likelihood * pW_O
		pW_O = pW_O/sum(pW_O)
		
		if(verbose) {
			cat('Conditional Probability of Words (rows) given Objects (columns)\n')
			print(pWgO)
			cat('Joint Probability of Words (rows) and Objects (columns)\n')
			print(pW_O) 
		}
		
		index = (rep-1)*length(ord$trials) + t # index for learning trajectory
		traj[[index]] = pWgO
	  }
	}
	#perf = diag(pWgO) / rowSums(pWgO) # no, we'll use parameterized choice
	# calculated from joint probs
	#perf = diag(pW_O)^chDec / sum(pW_O^chDec) # rowSums(pW_O^chDec)
	
	# power choice rule
	chProb = pWgO^chDec / outer(ones,colSums(pWgO^chDec))
	perf = diag(chProb) # same as diag(pWgO)^chDec / colSums(pWgO^chDec)

	#cat('\nChoice Probabilities for each word given each object (using exp choice rule):\n')
	#chProb = exp(chDec*pWgO) / outer(ones,colSums(exp(chDec*pWgO)))
	#cat(diag(chProb))
	
	want = list(perf=perf, matrix=chProb, traj=traj)
	return(want)
	}

