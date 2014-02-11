# George Kachergis  Dec 28, 2010
# adapting MINERVA2 for cross-situational statistical word learning
# MINERVA2 uses uniform random binary vectors (-1 or 1, with 
# missing values=0) -- let's try gaussian centered at 0
# each feature encoded with probability L (else, 0)
# forgetting is complementary: with prob F, a feature reverts to 0
# but Hintzman just uses L=1 (F=0)
# a probe vector comes in, and compared (S=sim) to each trace i:
#  S(i) = normalized dot product of probe and trace
# activation is S(i)^3
# echo intensity: summed activation of all traces when probed (recog)
# echo content: summed content of all memory traces, weighted by resp.
#  activations

sim <- function(a,b) {
	# MINERVA2 sim function: normalized dot product
	return( sum(a * b) / length(a) )
	}

cos_sim <- function(a,b) {
	# cosine similarity of vectors a and b:
	# a.b / ||a|| ||b|| (dot product and magnitude)
	return( sum(a * b) / (sqrt(sum(a^2))*sqrt(sum(b^2))) )
	}

generate_vectors <- function(voc_sz, D) {
	vecs <- matrix(0, nrow=voc_sz, ncol=D)
	for(i in 1:voc_sz) { vecs[i,] <- rnorm(D, sd=1/sqrt(D)) }
	return(vecs)
	}

echo_intensity <- function(M, probe, voc_sz) {
	intensity <- 0
	for(i in 1:voc_sz) {
		intensity = intensity + sim(probe, M[i,])^3
		}
	return(intensity)
	}

echo_content <- function(M, probe, voc_sz) { 
	content <- rep(0, dim(M)[2])
	for(i in 1:voc_sz) {
		content = content + sim(probe, M[i,])^3 * M[i,]
		}
	return(content)
	}

minerva_test <- function(M, words, objects, X, D) {
	# testing
	# retrieved_obj = sum(from i=0 to #_exemps) obect_exemp_i * sim(word_exemp_i, test_word)^3
	correct = c()
	voc_sz = dim(words)[1] # rows = words
	
	for(w in 1:voc_sz) {
		probe = c(words[w,], rep(0,D)) 
		content = echo_content(M, probe, voc_sz)
		sim = c()
		for(o in 1:voc_sz) {
			objprobe = c(rep(0,D), objects[o,])
			sim[o] = sim(content, objects[o,]) # or cos_sim
			}
		choice = which(sim==max(sim))
		correct[w] = ifelse(choice==w, 1, 0)
		}
	return(correct) 
	}

model <- function(ord, params, voc_sz) {
	D = 100 # dimensionality of vectors
	X <- params[1] # scalar for cos similarity...
	num_trials = dim(ord)[1]
	trial_sz = dim(ord)[2]
	
	# simulate 1000 Ss
	num_subs = 1000
	results <- matrix(0, nrow=num_subs, ncol=voc_sz)
	for(sub in 1:num_subs) {
	
	M <- matrix(0, nrow=num_trials, ncol=2*D) # memory for c(words,objects)
	#Mo <- matrix(0, nrow=voc_sz, ncol=D) # memory for objects
	words <- generate_vectors(voc_sz, D)
	objects <- generate_vectors(voc_sz, D)
	# training
	for(t in 1:num_trials) {
		tr = ord[t,]
		trwords = rep(0, D) # words on a trial
		trobjs = rep(0, D) # objects on a trial
		# Brendan says he multiplied these two vectors by a
		# uniform random number to simulate encoding failure (assume (0,1)?)
		for(p in 1:trial_sz) { # add up all words, objs separately
			trobjs = trobjs + objects[as.integer(tr[p]),]*runif(1)
			trwords = trwords + words[as.integer(tr[p]),]*runif(1)
			}

		M[t,] = c(trwords, trobjs) # concatenated wordobj vector
		}
	results[sub,] = minerva_test(M, words, objects, X, D)
	}
	return(colMeans(results))
	}

