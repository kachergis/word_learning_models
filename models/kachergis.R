# Associative Uncertainty- (Entropy) & Familiarity-Biased Model
# George Kachergis  george.kachergis@gmail.com  June 10, 2011

shannon.entropy <- function(p) {
	if (min(p) < 0 || sum(p) <= 0)
		return(NA)
	p.norm <- p[p>0]/sum(p)
	-sum(log2(p.norm)*p.norm)
	}

update_known <- function(m, tr_w, tr_o, startval=.01) {
  for(i in tr_w) {
    for(c in 1:dim(m)[2]) {
      if(sum(m[,c]>0) & m[i,c]==0) {
        m[i,c] = startval
        m[c,i] = startval
      }
    }
    for(j in tr_o) {
      if(m[i,j]==0) m[i,j] = startval
      if(m[j,i]==0) m[j,i] = startval
    }
  }
  return(m)
}


model <- function(params, ord=c(), ord_name="", reps=1, name="model", save_traj=FALSE, print_matrix=FALSE) {
	X <- params[1] # associative weight to distribute
	B <- params[2] # weighting of uncertainty vs. familiarity
	C <- params[3] # decay
	
	# learning trajectory data
	if(save_traj) {
		traj <- data.frame(Trial=numeric(0), testtype=character(0), x=numeric(0.0))
		ord <- read.table(paste(ord_dir,ord_name,".txt",sep=""), header=T) # ordering
		}
	
	voc_sz = max(unlist(ord), na.rm=TRUE) # vocabulary size
	ppt = length(ord$trials[[1]]$words) # pairs per trial
	mean_ent = c()
	m <- matrix(0, voc_sz, voc_sz) # association matrix
	# training
	for(rep in 1:reps) { # for trajectory experiments, train multiple times
	  for(t in 1:length(ord$trials)) { 
		#print(format(m, digits=3))
		
		tr_w = as.integer(ord$trials[[t]]$words)
		tr_w = tr_w[!is.na(tr_w)]
		tr_o = as.integer(ord$trials[[t]]$objs)
		tr_o = tr_o[!is.na(tr_o)]
		m = update_known(m, tr_w, tr_o) # what's been seen so far?
		ent_w = c() # more entropy = more dispersive
		for(w in tr_w) { ent_w = c(ent_w, shannon.entropy(m[w,])) }
		ent_w = exp(B*ent_w)
		
		ent_o = c() # more entropy = more dispersive
		for(o in tr_o) { ent_o = c(ent_o, shannon.entropy(m[,o])) }
		ent_o = exp(B*ent_o)
		
		nent = (ent_w %*% t(ent_o))
		#nent = nent / sum(nent)
		# get all current w,o strengths and normalize to distr X
		assocs = m[tr_w,tr_o]
		denom = sum(assocs * nent)
		m = m*C # decay everything
		
		m[tr_w,tr_o] = m[tr_w,tr_o] + (X * assocs * (ent_w %*% t(ent_o))) / denom # update assocs
		
		if(print_matrix) print(m)
		
		if(save_traj & t>1) { # saving model learning trajectories for graphing
			modelx <- test(m+.01)
			modela <- test_abs(m+.01)
			traj <- rbind(traj, data.frame(Trial=t, x=modelx[1]))
			traja <- rbind(traja, data.frame(Trial=t, x=modela[1]))
		}
	  }
	}
	if(save_traj) {
		save(traj, file=paste(name,"_",ord_name,"_traj.RData",sep=''))
		save(traja, file=paste(name,"_",ord_name,"_abs_traj.RData",sep=''))
		}
	m = m+.01
	want = list()
	#print(plot(1:dim(ord)[1], mean_ent)) # col=perf[perf<.5]
	#want[['tr_ent']] = mean_ent
	#want[['perf']] = perf
	return(m)
	}

