# what parameters maximize performance for different training orders?

library(pso)
require(tidyverse)
#ord1 = as.matrix(read.table("orderings/TR1-hiCD.txt"))
#ord2 = as.matrix(read.table("orderings/TR2-hiCD.txt")) 
#ord3 = as.matrix(read.table("orderings/TR3-medCD.txt")) 
#ord4 = as.matrix(read.table("orderings/TR4-medCD.txt"))
#ord7 = as.matrix(read.table("orderings/TR7-lowCD.txt"))
load("orderings/master_orders.Rdata")

source("models/kachergis.R") # base_model

INIT_PAR = c(0.227, 1.176, 0.963) # filtering: SSE=0.177
lower = c(.001, .01, .1)
upper = c(10, 15, 1)
controls = list(maxit=200, max.restart=2, reltol=.001)

multinomial_likelihood_perfect <- function(par, ord) {
	M = model(par, ord=ord)
	#pOgW = diag(M) / rowSums(M) 
	pOgW = M$perf # p(o|w)
	lik = sum(log(pOgW))
	return(-lik) # 18*log(1/18) = -52.02669 for AFC guessing; -16.63 for 8AFC
}


# find paras for each condition that yield highest performance
high_perf = tibble()
for(o in names(orders)) {
  fit <- psoptim(INIT_PAR, multinomial_likelihood_perfect, 
                 ord=orders[[o]]$train, 
                 lower=lower, upper=upper, control=controls)
  M = model(fit$par, ord=orders[[o]]$train) 
  high_perf <- high_perf %>% 
    bind_rows(bind_cols(cond = o, perf = mean(M$perf), 
                        chi=fit$par[1], lambda=fit$par[2], alpha=fit$par[3]))
}

# maximal learning rate (10), lambda of ~2.6, high decay (alpha=.6)

# all for ord1 (high CD)
# with learning rate limited to 1, ~.88 perf is achieved with lambda=.97 and alpha=.82
# chi  lambda	alpha	perf
# .5		.68      .87		.79		(worst perf with high lambda (~12.5))
# 1		.97      .82		.87
# 2	   1.15      .79		.92
# 5	   1.51      .73		.95
# 10	   1.78      .70		.98
# 20	   2.8		.65		.985
# 40	2.12	.63		.991

# medium CD
 # maximal learning rate (60), lambda of ~2.6, high decay (alpha=.6)
# chi	lambda	alpha	perf
# .5	1.59	.85		.77	
# 1		1.61	.81		.85
# 2		1.84	.76		.
# 5		2.44	.70		.
# 10	3.10	.65		.
# 20	3.78	.61		.
# 40	4.42	.57		.99

# low CD: lower perf at same learning rate, lower lambda
# maximal learning rate (60), lambda of ~2.6, high decay (alpha=.6)
# chi	lambda	alpha	perf
# .5	1.31	.77		.67
# 1		1.57	.71		.76

# overall, with increasing learning rate greater decay results in better performance 
# (the faster you can learn from a single trial, the better it is to forget old (less-correct) information)
# lambda also increases with higher learning rates: it is more useful to focus on uncertain stimuli 
# if you can learn a lot (presumably you've already learned the old stuff well)
# as CD increases, max performance increases for the same learning rate, and decay is lower (>alpha). 
# lambda is highest for medium CD, but lowest for low CD (does this hold for adult orderings?)
