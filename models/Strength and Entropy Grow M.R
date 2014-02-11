# cross-situational statistical word learning models
# (was testmodel2 -- works well for freq and CD exp)

shannon.entropy <- function(p) {
	#if (min(p) < 0 || sum(p) <= 0)
	#	return(NA)
	p.norm <- p[p>0]/sum(p)
	-sum(log2(p.norm)*p.norm)
	}

update_known <- function(m, tr) {
	num_known = 0
	for(r in 1:dim(m)[1]) {
		if(sum(m[r,])!=0) num_known = num_known+1
		}
	startval = 1/100
	for(i in tr) {
		for(c in 1:dim(m)[2]) {
			if(sum(m[,c]>0) & m[i,c]==0) {
				m[i,c] = startval
				m[c,i] = startval
			}
		}
		for(j in tr) {
			if(m[i,j]==0) m[i,j] = startval
			if(m[j,i]==0) m[j,i] = startval
			}
		}
		return(m)
	}


model <- function(ord, params, voc_sz) {
	X <- params[1]*10 # associative weight to distribute
	#A <- params[2]
	B <- params[2]
	#recw <- .797 # prob of recognizing each word on a trial
	m <- matrix(0, voc_sz, voc_sz)
	trial_sz = dim(ord)[2]
	# training
	for(t in 1:dim(ord)[1]) { 		
		tr = ord[t,]
		m = update_known(m, tr)
		ent = c() # more entropy = more dispersive
		for(w in tr) { ent = c(ent, shannon.entropy(m[w,])) }
		ent = exp(B*ent)
		ent = ent / sum(ent) # should be row & column entropy
		
		# get all current w,o strengths and normalize to distr X
		assocs = matrix(0, nrow=trial_sz, ncol=trial_sz)
		for(w in 1:trial_sz) {
			for(p in 1:trial_sz) {
				assocs[w,p] = m[as.integer(tr[w]), as.integer(tr[p])]
				}
			}
		assocs = assocs / sum(assocs)
		
		for(w in 1:trial_sz) {
			for(p in 1:trial_sz) {
				m[as.integer(tr[w]), as.integer(tr[p])] = m[as.integer(tr[w]), as.integer(tr[p])] + ent[w]*ent[p]*X*assocs[w,p]
				}
			#m[as.integer(tr[w]),] = m[as.integer(tr[w]),]/sum(m[as.integer(tr[w]),])
			}
		#m = m / sum(m)
		}
	#m = exp(m*A)
	return(test(m)) 
	}

