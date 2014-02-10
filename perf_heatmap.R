library(plotrix); library(lattice)

ord = read.table("orig_4x4.txt")

INIT_PAR = c(0.0762, 5.861, 0.979)

fit <- function(mod, ord) {
	source(paste(mod,".R",sep=""))
	
	best <- optim(INIT_PAR, freq_sse, method="L-BFGS-B", lower=c(0.001,.01,.8), upper=c(20,30,1), control=list(factr=1e3, maxit=50))
	# control=list(parscale=c(10,10,1))
	# fixed=list(fixedpar=27)
	# for-non BFGS algs: reltol=.001,
	
	results <- read.table("tmp_model_output.txt")
	best$results <- results
	save(best, file=paste(mod,"_freq_SSEfit.RData", sep=""))
	}


eval_perf_3d_grid <- function(mod, ord, p1s=seq(.01,5,.03), p2s=seq(.1,18,.3), p3s=seq(.05,1,.1), fname="") {
	# p1 = assoc weight, p2 = familiarity (<) vs. uncertainty (>), p3 = decay
	source(paste(mod,".R",sep=""))
	its <- length(p1s)*length(p2s)*length(p3s)
	count = 0
	lev = seq(0, 1, by=.05) 
	all <- list()
	for(p3 in 1:length(p3s)) { # memory par is outermost
		grid <- matrix(0, nrow=length(p1s), ncol=length(p2s))
		for(p1 in 1:length(p1s)) {
			for(p2 in 1:length(p2s)) {
				mat = model(c(p1s[p1],p2s[p2],p3s[p3]), ord=ord)
				grid[p1,p2] = mean(diag(mat) / rowSums(mat))
				count = count+1
				if(count%%1000==0) {
					print(paste(count,"/",its))
					}
				}
			}
		if(fname=="") {
			quartz()
		} else {
			pdf(paste(fname,p3s[p3],".pdf",sep=''), width=6, height=6.5, pointsize=11)
			image(p1s, p2s, grid, col=heat.colors(100), axes=F, xlab="Associative Weight (X)", ylab="Familiarity vs. Uncertainty (Lambda)", main=paste("Memory Fidelity (alpha) =",p3s[p3]))
			contour(p1s, p2s, grid, levels=lev, add=TRUE, col="black")
			axis(1, at=p1s)
			axis(2, at=p2s)
			dev.off()
		}
		
		row.names(grid) = p1s
		colnames(grid) = p2s
		
		print(paste(min(grid), max(grid), mean(grid)))
		print(which(grid==min(grid)))
		all[[paste("decay",as.character(p3s[p3]),sep='')]] = grid #
		}
	return(all)
	}

grids <- eval_perf_3d_grid("model", ord, fname="perf_heatmap_chi_vs_lambda_alpha")
save(grids, file="grid_model.RData")
