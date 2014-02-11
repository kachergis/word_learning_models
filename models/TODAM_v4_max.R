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

# in _max models, picks the best match
todam_test <- function(M, words, objects, X) {	
	# testing
	correct = c()
	voc_sz = dim(words)[1] # rows = words
	for(w in 1:voc_sz) {
		sim = rep(0,voc_sz)
		for(o in 1:voc_sz) {
			sim[o] = cos_sim(words[w,], correlate(objects[o,], M))
			}
		choice = which(sim==max(sim))
		correct[w] = ifelse(choice==w, 1, 0)
		}
	return(correct) 
	}

model <- function(ord, params, voc_sz) {
	D = 8000 # dimensionality of vectors (param?)
	X <- params[1] # scalar for cos similarity...
	#alpha <- params[1] # decay for M
	#gamma <- params[2]
	#omega <- params[3]
	M <- rep(0, D) # memory
	words <- generate_vectors(voc_sz, D)
	objects <- generate_vectors(voc_sz, D)
	trial_sz = dim(ord)[2]
	# training
	for(t in 1:dim(ord)[1]) {
		tr = ord[t,]
		context = rep(0, D)
		for(p in 1:trial_sz) {
			context = context + objects[as.integer(tr[p]),]
			}
		for(w in 1:trial_sz) {
			w_rep <- words[as.integer(tr[w]),]
			M = M + convo(w_rep, context)
			#M = alpha*M + gamma*convo(w_rep, context)
			}
		
		}
	return(todam_test(M, words, objects, X))
	}

