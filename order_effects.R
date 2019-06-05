get_forward_reverse_cor <- function(ord, par, modeln="") {
	if(modeln!="") source(paste("models/",modeln,".R",sep=''))
	tr_f = ord$train
	tr_r = ord$train
	tr_r$words = tr_r$words[nrow(tr_r$words):1,]
	tr_r$objs = tr_r$objs[nrow(tr_r$objs):1,]

	m_f = model(par, ord=tr_f)
	m_r = model(par, ord=tr_r)

	print(paste("cor:",round(cor(as.vector(m_f$perf), as.vector(m_r$perf)),3)))
	#plot(as.vector(m_f$perf), as.vector(m_r$perf), xlab="Item Accuracy on Forward Order", ylab="Item Accurac on Reverse Order")
	#print(paste("mean forward perf:",mean(m_f$perf),"mean reverse perf:",mean(m_r$perf)))
	print(paste("mean forward - reverse perf:",round(mean(m_f$perf)-mean(m_r$perf), 3)))
	dat = data.frame(forward=as.vector(m_f$perf), backward=as.vector(m_r$perf))
	dat$Condition = ord$Condition
	return(dat)
}

dat = data.frame()
for(ordn in names(orders)) {
	print(ordn)
	dat = rbind(dat, get_forward_reverse_cor(orders[[ordn]], c(1, 2, .9), modeln="kachergis")) # -.09
	get_forward_reverse_cor(orders[[ordn]], c(1, 2, 1), modeln="kachergis") # -.23
	get_forward_reverse_cor(orders[[ordn]], c(.1, 2, .9), modeln="kachergis") # -.12
	get_forward_reverse_cor(orders[[ordn]], c(.1, .5, .9), modeln="kachergis") # .48
}

#get_forward_reverse_cor(orders[["freq369-3x3hiCD"]], c(1, 2, .9), modeln="kachergis") # -.09
#get_forward_reverse_cor(orders[["freq369-3x3hiCD"]], c(1, 2, 1), modeln="kachergis") # -.23
#get_forward_reverse_cor(orders[["freq369-3x3hiCD"]], c(.1, 2, .9), modeln="kachergis") # -.12
#get_forward_reverse_cor(orders[["freq369-3x3hiCD"]], c(.1, .5, .9), modeln="kachergis") # .44
#get_forward_reverse_cor(orders[["freq369-3x3hiCD"]], c(.1, 4, .9), modeln="kachergis") # -.19


#get_forward_reverse_cor(orders[["freq369-3x3hiCD"]], c(.1, 4, .9), modeln="strength") # -.05
#get_forward_reverse_cor(orders[["freq369-3x3hiCD"]], c(.1, 4, 1), modeln="strength") # -.05
#get_forward_reverse_cor(orders[["freq369-3x3hiCD"]], c(.1, 1, 1), modeln="strength") # .83
#get_forward_reverse_cor(orders[["freq369-3x3hiCD"]], c(.1, 1, .8), modeln="strength") # .83
#get_forward_reverse_cor(orders[["freq369-3x3hiCD"]], c(.1, .5, .8), modeln="strength") # .11


#get_forward_reverse_cor(orders[["freq369-3x3hiCD"]], c(1e-5, 8500), modeln="fazly") # 1
#get_forward_reverse_cor(orders[["freq369-3x3hiCD"]], c(1e-5, 4000), modeln="fazly") # 1 
#get_forward_reverse_cor(orders[["freq369-3x3hiCD"]], c(1e-2, 8500), modeln="fazly") # 1

