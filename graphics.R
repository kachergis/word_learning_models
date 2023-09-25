# functions to graph word-object co-occurrences from trial orderings and model associations

make_cooccurrence_matrix <- function(cond, print_matrix=F, heatmap_filename=c()) {
	# makes word x object co-occurrence matrix from training list of words and objects by trial
	# prints a heatmap if filename is specified
	words = cond$train$words
	objs = cond$train$objs
	m = matrix(0, nrow=max(words), ncol=max(objs))
	for(t in 1:nrow(words)) {
		m[words[t,], objs[t,]] = m[words[t,], objs[t,]] + 1
	}
	
	if(print_matrix==T) print(m)
	if(length(heatmap_filename>0)) {
		pdf(paste(heatmap_filename,".pdf",sep="")) 
		heatmap(m, Rowv = NA, Colv = "Rowv",  scale="none", margin=c(3,3), xlab="Object", ylab="Word", col=heat.colors(10)) 
		# labRow=NA, labCol=NA
		dev.off()
	}
	return(m)
}


animate_trajectory <- function(mod, modname='', condname='') {
	# creates an animation of a model's word-object associations
	require(animation)
	ani.options(interval=.1) # delay between frames
	col.range <- heat.colors(20)
	breaks = seq(0,1, .05)
	saveGIF({
		#layout(matrix(c(1, rep(2, 5)), 6, 1))
		for(i in 1:length(mod$traj)) {
			mm = mod$traj[[i]] + 1e-7 # to eliminate division by zero
			mm = mm / rowSums(mm) # row normalize
			heatmap(mm, Rowv=NA, Colv="Rowv", scale="none", margin=c(3,3), breaks=breaks, xlab="Object", ylab="Word", col=col.range) 
		}
	}, movie.name=paste(modname,"_model",condname,"_cond_trajectory.gif",sep=''))
}


