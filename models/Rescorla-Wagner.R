# cross-situational statistical word learning models
# Rescorla-Wagner 1972 model:
# predicted outcome (word) given objects is sum of weights


model <- function(params, ord=c(), ord_name="", name="Rescorla-Wagner") {
	lambda <- params[1] # learning rate
	decay <- params[3]
	voc_sz = max(ord) # vocabulary size
	m <- matrix(0.01, voc_sz, voc_sz)
	trial_sz = dim(ord)[2]
	# training
	for(t in 1:dim(ord)[1]) { 		
		tr = ord[t,]
		#ent = c()
		#for(w in tr) { ent = c(ent, shannon.entropy(m[w,])) }
		#ent = exp(B*ent)
		#ent = ent / sum(ent) # should be row & column entropy
		
		# cues = objects, outcomes = words
		pred = rep(0, trial_sz)
		for(p in 1:trial_sz) {
			for(w in 1:trial_sz) {
				pred[p] = pred[p] + m[as.integer(tr[w]), as.integer(tr[p])]
				}
			}
			
		outcome = 1 # ??
		delta = rep(0, trial_sz)
		for(p in 1:trial_sz) {
			for(w in 1:trial_sz) {
				m[as.integer(tr[w]), as.integer(tr[p])] = lambda*(outcome-pred[p])
				}
			}
		m = m*decay
		}
	
	return(test(m, voc_sz, ord_name))  
	}


