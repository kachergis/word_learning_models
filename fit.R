

exp <- read.table("humans/cdf_exp2_means.txt", header=T)

fit_subj_traj <- function(par, ord, perf) {
	mp = model(par, ord, reps=length(perf))
	sse = sum((mp-perf)^2)
	return(sse)
}

fit_subj_traj_samp <- function(par, ord, perf, K, Nsubj=1000) {
	mperf = rep(0, length(perf))
	for(s in 1:Nsubj) {
		mperf = mperf + model(par, ord, K=K, reps=length(perf))
	}
	mperf = mperf/Nsubj
	sse = sum((mperf-perf)^2)
	print(par)
	print(mperf)
	print(sse)
	return(sse)
}

fit_subj <- function(par, ords, perf) {
	sse = 0
	files = names(ords)
	for(f in files) {
		o <- ords[[f]]
		cp <- perf[[f]]
	}
}

fit_generic <- function(par, ords, Nsubj, exp) {
	if(exp=="Exp2") return(fit_exp2_samp(par, ords, Nsubj))
	if(exp=="Exp3") return(fit_exp3_samp(par, ords, Nsubj))
}

fit_exp3_samp <- function(par, ords, Nsubj, exp) {
	sse = 0
	results <- data.frame(EarlyRepetitions=numeric(0), Late=numeric(0), testtype=character(0), x=numeric(0))
	for(f in names(ords)) {
		o <- ords[[f]]
		# ans[i] = p(choosing obj i for word i)
		#ans <- model(par, o)
		ans <- c() # rows = Ss, cols = item correctness
		for (s in 1:Nsubj) {
			s <- model(par, o)
			ans = rbind(ans, s)
		}
		# compare simulation avg to each condition's human avg
		fr3 <- mean(ans[,1:6])
		fr6 <- mean(ans[,7:12])
		fr9 <- mean(ans[,13:18])
		
		if(f=="block1_3_6_9-3x3-lo_cd") { 
			cond = "aLow CD"
			sse = sse + (fr3-.25 )^2 + (fr6-.42 )^2 + (fr9-.71 )^2
		} else if(f=="block2_3_6_9-3x3") {
			cond = "High CD"
			sse = sse + ( fr3-.45 )^2 + ( fr6-.64 )^2 + ( fr9-.6 )^2
		} else if(f=="block3_3_6_9_lo_medCD") {
			cond = "Med CD (3/6)"
			sse = sse + (fr3-.38 )^2 + ( fr6-.47 )^2 + ( fr9-.64 )^2
		} else { # block4_369_39mx
			cond = "Med CD (3/9)"
			sse = sse + ( fr3-.63 )^2 + ( fr6-.58 )^2 + ( fr9-.72 )^2
		}
		results <- rbind(results, data.frame(Condition=cond, Freq=3, x=fr3))
		results <- rbind(results, data.frame(Condition=cond, Freq=6, x=fr6))
		results <- rbind(results, data.frame(Condition=cond, Freq=9, x=fr9))
	}
	print(results)
	write.table(results, file="tmp_model_output.txt")
	print("parms:")
	print(par)
	print(paste("SSE:",sse))
	return(sse)
}

fit_exp2_samp <- function(par, ords, Nsubj, exp) {
	# for Exp 2
	sse = 0 # sum of squared error
	files = names(ords)
	for(f in files) {
		o <- ords[[f]]
		ans <- c() # rows = Ss, cols = item correctness
		for (s in 1:Nsubj) {
			s <- model(par, o)
			ans = rbind(ans, s)
		}
		print(paste("SD:",mean(rowMeans(ans))))
		# compare simulation avg to each condition's human avg
		if(f=="orig_order_study_9s-corrected") { 
			sse = sse + ( mean(ans)-.31 )^2
			print(paste(f,"4x4:",mean(ans)))
		} else if(f=="cont_div6-12") {
			it_m <- colMeans(ans)
			loCD <- it_m[1:6]
			hiCD <- it_m[7:18]
			sse = sse + ( mean(loCD)-.333 )^2 + ( mean(hiCD)-.522 )^2
			print(paste(f,"loCD:",mean(loCD)))
			print(paste(f,"hiCD:",mean(hiCD)))
		} else { # 3x3_orig_order
			sse = sse + ( mean(ans)-.43 )^2
			print(paste(f,"3x3:",mean(ans)))
		}
	}
	print("parms:")
	print(par)
	print(paste("SSE:",sse))
	return(sse)
}

fit_exp2_assoc <- function(par, ords) {
	# for Exp 2
	sse = 0 # sum of squared error
	results <- data.frame(Condition=numeric(0), x=numeric(0))
	files = names(ords)
	for(f in files) {
		o <- ords[[f]]
		# rows = Ss, cols = item correctness
		
		# compare simulation avg to each condition's human avg (pairs learned)
		if(f=="orig_order_study_9s-corrected") {
			ans <- model((4/3)*par, o)
			sse = sse + ( mean(ans)*18.0-5.58 )^2
			line <- data.frame(Condition="18 pairs, 4/trial", x=mean(ans)*18)
			results <- rbind(results, line)
			
		} else if(f=="3x3_orig_order") {
			ans <- model(par, o)
			sse = sse + ( mean(ans)*18.0-7.74 )^2
			line <- data.frame(Condition="18 pairs, 3/trial", x=mean(ans)*18)
			results <- rbind(results, line)
			
		} else if(f=="cont_div6-12") {
			ans <- model(par, o)
			loCD <- mean(ans[1:6])
			hiCD <- mean(ans[7:18])
			sse = sse + ( mean(loCD)*6.0-2.04 )^2 + ( mean(hiCD)*12-5.64 )^2
			results <- rbind(results, data.frame(Condition="12 pairs, 3/trial", x=mean(hiCD)*12))
			results <- rbind(results, data.frame(Condition="6 pairs, 3/trial", x=mean(loCD)*6))
			
		} else {
			print("Unrecognized ordering name...")
			}
	}
	print(results)
	#write.table(results, file="tmp_model_output.txt")
	print("parms:")
	print(par)
	print(paste("SSE:",sse))
	return(sse)
}
