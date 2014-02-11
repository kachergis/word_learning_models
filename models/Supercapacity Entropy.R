# cross-situational statistical word learning models

shannon.entropy <- function(p) {
	if (min(p) < 0 || sum(p) <= 0)
		return(NA)
	p.norm <- p[p>0]/sum(p)
	-sum(log2(p.norm)*p.norm)
}

model <- function(ord, params, voc_sz) {
	R <- params*10 # associative weight to distribute
	m <- matrix(1/voc_sz, voc_sz, voc_sz)
	trial_sz = dim(ord)[2]
	# training
	for(t in 1:dim(ord)[1]) { 		
		tr = ord[t,]
		ent = c() # more entropy = more dispersive
		for(w in tr) { ent = c(ent, shannon.entropy(m[w,])) }
		ent = ent / sum(ent) # should be row & column entropy
		for(w in 1:trial_sz) {
			for(p in 1:trial_sz) {
				prev = m[as.integer(tr[w]), as.integer(tr[p])]
				m[as.integer(tr[w]), as.integer(tr[p])] = prev + ent[w]*ent[p]*R*prev
				}
			}
		}
	
	return(test(m))
	}


