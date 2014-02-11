# try TODAM for cross-situational statistical word learning

# convolve(c(1,2,3,4), c(2,3,4,5), type="circular") # real but slow
#filter(x, rep(1, 3), sides = 1, circular = TRUE)
# convolve associates two items

#correl <- function(cue, mem) {
#	return(convolve(rev(cue), mem, type="circular"))
#	}

correlate <- function(a, b) {
	a = rev(a)
	return(as.real(convo(c(a[length(a)], a[1:length(a)-1]), b)))
	}

convo <- function(a, b) {
	return(as.real(fft(fft(a) * fft(b), inverse=TRUE)/length(a)))
	}

convoOLD <- function(a, b) {
	return( rev( convolve(a, rev(b), type="circular") ) )
	}
	
correlateOLD <- function(cue, mem) {
	D = length(cue)
	v = rep(0, D)
	for(x in 1:D) {
		for(y in 1:length(cue)) {
			v[x] = v[x] + cue[(y-1)%%D + 1]*mem[(x+y-1)%%D + 1]
			}
		}
	# based on Mike's code, but below I roll vector
	# to make it right:
	return(c(v[length(v)], v[1:length(v)-1]))
	}

cos_sim <- function(a,b) {
	# cosine similarity of vectors a and b:
	# a.b / ||a|| ||b|| (dot product and magnitude)
	return( sum(a * b) / (sqrt(sum(a^2))*sqrt(sum(b^2))) )
	# Mike uses Matlab's corrcoef, but that's not right, I think:
	# (takes sqrt of whole denom, instead of mag(a) * mag(b))
	}

generate_vectors <- function(voc_sz, D) {
	vecs <- matrix(0, nrow=voc_sz, ncol=D)
	for(i in 1:voc_sz) { vecs[i,] <- rnorm(D, sd=1/sqrt(D)) }
	return(vecs)
	}

todam_testOLD <- function(M, words, objects) {	
	# testing
	correct = c()
	voc_sz = dim(words)[1] # rows = words
	for(w in 1:voc_sz) {
		#sim = c()
		#for o in 1:voc_sz) {
		#	sim[o] = cos_sim(words[w,], correlate(objects[o,], M))
		#	}
		correct[w] = cos_sim(words[w,], correlate(objects[w,], M))
		}
	return(correct) 
	}

todam_test <- function(M, words, objects, X) {	
	# testing
	correct = c()
	voc_sz = dim(words)[1] # rows = words
	for(w in 1:voc_sz) {
		sim = c()
		for(o in 1:voc_sz) {
			sim[o] = cos_sim(words[w,], correlate(objects[o,], M))
			}
		if(min(sim) < 0) { sim = sim + min(sim) }
		sim = sim/sum(sim)
		sim = exp(sim * X)
		correct[w] = sim[w] / sum(sim)
		}
	return(correct) 
	}

model <- function(ord, params, voc_sz) {
	D = 2000 # dimensionality of vectors (param?)
	X <- params[1] # scalar for cos similarity...
	#alpha <- params[1] # decay for M
	#gamma <- params[2]
	#omega <- params[3]
	M <- rep(0, D)
	words <- generate_vectors(voc_sz, D)
	objects <- generate_vectors(voc_sz, D)
	trial_sz = dim(ord)[2]
	# training
	for(t in 1:dim(ord)[1]) { 		
		tr = ord[t,]
		
		for(w in 1:trial_sz) {
			for(p in 1:trial_sz) {
				w_rep <- words[as.integer(tr[w]),]
				o_rep <- objects[as.integer(tr[p]),]
				M = M + w_rep + o_rep + convo(w_rep, o_rep)
				#M = alpha*M + gamma*(w_rep + o_rep) + omega*convo(w_rep, o_rep)
				}
			}
		
		}
	return(todam_test(M, words, objects, X))
	}


