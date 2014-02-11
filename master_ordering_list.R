# want to save all trial orderings in a common format:
# a list of orderings, which are in turn lists of trials,
# each with a vector of words and objects

# should also store the item indices of subgroups of interest
# (e.g., different freq), and the human performance on those groups...later


# this saving function works on cross-situational orderings
# where the same words and objects are on each trial
save_master_list <- function(files, orders=list(), folder="orderings", save=FALSE) {
	for(f in files) {
		o <- read.table(paste(folder,"/",f,".txt",sep=""), header=F)
		trials = list()
		for(t in 1:dim(o)[1]) {
			wo <- as.numeric(o[t,])
			trials[[t]] = list(words=wo, objs=wo)
			}
		orders[[f]] = list(trials=trials)
		#orders[[f]]$groups = list(freq3=seq(1,6), freq6=c(7,12))
		#orders[[f]]$hum_perf = list(freq3=.xx, )
		}
	if(save) save(orders, file="master_orders.RData")
	return(orders)
	}

files <- c("filt0E_3L", "filt3E_3L", "filt6E_3L", "filt9E_3L", 
	"filt0E_6L", "filt3E_6L", "filt6E_6L", "filt9E_6L", 
	"filt0E_9L", "filt3E_9L", "filt6E_9L", "filt9E_9L",
	"2_x8_39_4x4", "3_x8_369_4x4", "block1_369-3x3loCD",
	"block2_369-3x3hiCD", "block3_369_36mx", "block4_369_39mx",
	"orig_4x4", "orig_3x3", "cont_div6-12")

orders = save_master_list(files, save=F)

fgr = list(freq3=seq(1,6), freq6=seq(7,12), freq9=seq(13,18))

orders[["orig_3x3"]]$groups = list(all=seq(1,18))
orders[["orig_3x3"]]$hum_perf = list(all=0.43)

orders[["orig_4x4"]]$groups = list(all=seq(1,18))
orders[["orig_4x4"]]$hum_perf = list(all=0.31)

orders[["cont_div6-12"]]$groups = list(loCD=seq(1,6), medCD=seq(7,18))
orders[["cont_div6-12"]]$hum_perf = list(loCD=0.34, medCD=.47)


orders[["block1_369-3x3loCD"]]$groups = fgr
orders[["block1_369-3x3loCD"]]$hum_perf = list(freq3=.25, freq6=.42, freq9=.71)

orders[["block2_369-3x3hiCD"]]$groups = fgr
orders[["block2_369-3x3hiCD"]]$hum_perf = list(freq3=.45, freq6=.64, freq9=.6)

orders[["block3_369_36mx"]]$groups = fgr
orders[["block3_369_36mx"]]$hum_perf = list(freq3=.38, freq6=.47, freq9=.64)

orders[["block4_369_39mx"]]$groups = fgr
orders[["block4_369_39mx"]]$hum_perf = list(freq3=.63, freq6=.58, freq9=.72)


orders[["filt0E_3L"]]$groups = list(all=seq(1,12))
orders[["filt0E_3L"]]$hum_perf = list(all=0.372)

orders[["filt0E_6L"]]$groups = list(all=seq(1,12))
orders[["filt0E_6L"]]$hum_perf = list(all=0.462)

orders[["filt0E_9L"]]$groups = list(all=seq(1,12))
orders[["filt0E_9L"]]$hum_perf = list(all=0.586)

save(orders, file="master_orders.RData")
