# try TODAM for cross-situational statistical word learning

#filter(x, rep(1, 3), sides = 1, circular = TRUE)
# convolve associates two items

correlate <- function(a, b) {
	a = rev(a)
	return(as.real(convo(c(a[length(a)], a[1:length(a)-1]), b)))
	}

convo <- function(a, b) {
	return(as.real(fft(fft(a) * fft(b), inverse=TRUE)/length(a)))
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

todam_test <- function(M, words, objects, X) {	
	# testing
	correct = c()
	voc_sz = dim(words)[1] # rows = words
	for(w in 1:voc_sz) {
		sim = c()
		for(o in 1:voc_sz) {
			sim[o] = cos_sim(words[w,], correlate(objects[o,], M))
			}
		sim = sim + min(sim) # necessary for -sims
		sim = sim/sum(sim)
		sim = exp(sim * X)
		correct[w] = sim[w] / sum(sim)
		}
	return(correct) 
	}

model <- function(ord, params, voc_sz) {
	D = 8000 # dimensionality of vectors (param?)
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
		
		w_attn = c()
		o_attn = c()
		for(w in 1:trial_sz) {
			w_rep <- words[as.integer(tr[w]),]
			o_rep <- objects[as.integer(tr[w]),]
			w_attn <- c(w_attn, cos_sim(w_rep, M))
			o_attn <- c(o_attn, cos_sim(o_rep, M))
			}
		w_attn = w_attn + min(w_attn)
		o_attn = o_attn + min(o_attn)
		w_attn = w_attn/sum(w_attn)
		o_attn = o_attn/sum(o_attn)
		
		for(w in 1:trial_sz) {
			w_rep <- words[as.integer(tr[w]),]
			for(p in 1:trial_sz) {
				o_rep <- objects[as.integer(tr[p]),]
				#M = M + w_rep + o_rep + convo(w_rep, o_rep)
				M = M + w_attn[w]*o_attn[p]*convo(w_rep, o_rep)
				#M = alpha*M + gamma*(w_rep + o_rep) + omega*convo(w_rep, o_rep)
				}
			}
		
		}
	return(todam_test(M, words, objects, X))
	}


