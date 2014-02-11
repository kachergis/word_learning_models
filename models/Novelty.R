# cross-situational statistical word learning models

shannon.entropy <- function(p) {
	if (min(p) < 0 || sum(p) <= 0)
		return(NA)
	p.norm <- p[p>0]/sum(p)
	-sum(log2(p.norm)*p.norm)
}

model <- function(ord, params, voc_sz) {
	X <- params*10 # associative weight to distribute
	m <- matrix(1/voc_sz, voc_sz, voc_sz)
	freq = rep(0,voc_sz) # freq[i] = times pair i has appeared
	trial_sz = dim(ord)[2]
	# training
	for(t in 1:dim(ord)[1]) { 
		tr = ord[t,]
		# update pair counts
		for(p in tr) { freq[p] = freq[p] + 1 }
		nov = c() # more entropy = more dispersive
		for(i in 1:length(tr)) { 
			nov = c(nov, 1/(1+freq[as.integer(tr[i])])) 
			}
		nov = nov / sum(nov) # should be row & column entropy
		
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
				m[as.integer(tr[w]), as.integer(tr[p])] = m[as.integer(tr[w]), as.integer(tr[p])] + nov[w]*nov[p]*X*assocs[w,p]
				}
			
			}
		}
	
	return(test(m))  
	}


