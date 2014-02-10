# for dissertation Frequency and Contxtual Diversity chapter figures
# Oct 29, 2012
setwd("~/Dropbox/projects/freqCD/model_for_paper")

folder <- "orderings"

files = c("3_x8_369_4x4", "2_x8_39_4x4") # Exp 1
files = c("cont_div6-12", "orig_3x3", "orig_4x4") # Exp 2
files = c("block1_369-3x3loCD", "block2_369-3x3hiCD", "block3_369_36mx", "block4_369_39mx") # Exp 3

outdir = "heatmaps/"

for (f in files) {
	train <- read.table(paste(folder,"/",f,".txt", sep=""))
	coocs = get_cooc_matrix(train) # for comparison

	#jpeg(paste(outdir,f,"_heatmap.jpg",sep=""))
	pdf(paste(outdir,f,"_heatmap.pdf",sep=""))
	heatmap(coocs, Rowv = NA, Colv = NA, symm=TRUE, scale="none", labRow=NA, margin=c(.2,.2), labCol=NA, xlab="", ylab="", col = heat.colors(15))
	dev.off()
}


get_cooc_matrix <- function(train) {
	m <- matrix(0, nrow=max(train),ncol=max(train))
	for (i in 1:dim(train)[1]) {
		#m[train[i,],train[i,]] = m[train[i,],train[,i]] + 1
		for (wi in 1:dim(train)[2]) {
			w <- train[i,wi]
			for (pi in 1:dim(train)[2]) {
				p <- train[i,pi]
				m[w,19-p] = m[w,19-p] + 1
				}
			}
		}
	return(m)
	}