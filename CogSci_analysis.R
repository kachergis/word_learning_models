# CogSci 2013
# modeling asymmetric and noisy conditions
# 225, 2x2.  Words appear 6 times. objects 6.  18afc.
# 205, 2x4.  extending the original 2x2. Each word appears 6 times.  Each object appears 12 times.  18afc
# 206, 3x3.  extending the original 2x2. 2x2 normal word/referent. plus 1 random word and 1 random object.  18afc
# 207, 4x4.  extending the original 2x2. 2x2 normal word/referent. plus 2 random words and 2 random objects.  18afc

# 201, 3x4 condition. extending the original 2x2 [3x3??]. each word appears 6 times. each object 8 times. Therefore p(w|o) = 80% for all 18 pairs. 18afc
# 202, 3x4 condition.  extending the original 2x2 [3x3??]. Each word appears 6 times.  First 12 objects appear only with their words p(w|o) = p(o|w) = 1. objs 13-18 appear 12 times. p(w|o) = 50%.  18afc
# 203, 3x4.  extending the original 2x2 [3x3??]. All words appear 6 times.  Objects appear with their words p(w|o) = p(o|w) = 1. and objs 3-8 and 13-18 also appear 3 more times p(w|o) = 66%.
# 204, 3x4.  All 18 words appear 6 times. [extending 3x3??]  There are 24 objects. each appears 6 times.  The first 18 are referents. the last 6 are not the referents of any particular word. and therefore noisy items18afc.

w_by_o = c("3x4","3x4","3x4","3x4","2x4","3x3","4x4","1x3","2x3","2x4","3x4","2x3","3x4","1x4","1x3","4x4","3x3","2x2")
o_freq = c(8,12,7.5,6,12,6,6,18,9, 12,8,9,12,24,18,6,6,6)
w_freq = c(6,6,6,6,6,6,6,6,6, 6,6,6,9,6,6,6,6,6)
# 203: objs 3-8,13-18 9x, others 6x
# 204: objs 19-24 are noise

# 215, 1x3.  Words appear 6 times. objects 18.  18afc. 23 Ss
# 222, 1x3.  Words appear 6 times. objects 18.  18afc. 40 Ss
# 216, 2x3.  Words appear 6 times. objects 9.  18afc. 23 Ss
# 219, 2x3.  Words appear 6 times. objects 9.  18afc 32 Ss
# 221, 1x4.  Words appear 6 times. objects 24.  18afc. 40 Ss ~.18 (many bad Ss?)
# 223, 4x4.  Words appear 6 times. objects also 6.  18afc. 77 Ss
# 224, 3x3.  Words appear 6 times. objects 6.  18afc. 36 Ss

# prob. don't use:
# 217, 2x4.  Words appear 6 times. objects 12.  18afc. [but not built from 2x2...nor 4x4?] ~.3 (14 Ss)
# 218, 3x4.  Words appear 6 times. objects 8.  18afc. [not built from 2x2...nor 4x4?] ~.43 (13 Ss)
# 220, 3x4.  Words appear 9 times. objects 12.  18afc - high perf: ~.67 (33 Ss)

#exp3x3 = c(215, 224) # all 3x3...
exp2x2 = c(201,202,203,204,205,206,207,225) # all extending 2x2, no within-cond manipulation (graph together)
exp4x4 = c(215,222,216,219,221,223,224) # all built from 4x4 [um...maybe not]
all$Condition <- with(all, ifelse(CondNum==215 | CondNum==222, "1x3", ifelse(CondNum==216 | CondNum==219, "2x3", 
                           ifelse(CondNum==221, "1x4", ifelse(CondNum==223, "4x4", ifelse(CondNum==224, "3x3", 
                           ifelse(CondNum==225, "2x2", ifelse(CondNum==201, "3x4", ifelse(CondNum==202, "3x4 1/.5", 
                           ifelse(CondNum==203, "3x4 1/.66", ifelse(CondNum==204, "3x4 +6o", ifelse(CondNum==205, "2x4", 
                           ifelse(CondNum==206, "3x3 +1w/o", ifelse(CondNum==207, "4x4 +2w/o", NA))))))))))))))
exp4x4c = c("1x3","2x3","1x4","4x4","3x3") # Exp 2 by Condition name 
exp2x2c = c("2x2","3x4", "3x4 1/.5", "3x4 1/.66", "3x4 +6o", "2x4", "3x3 +1w/o", "4x4 +2w/o")

# for each order, for each word/obj pair print freq, p(w|o) & p(o|w)
ordering_stats <- function(condns) {
  load("asym_conds.RData") # conds[[num]]$train, $resps, $test
  for (n in condns) {
    train = conds[[n]]$train
    coocs = get_cooc_matrix(train) # for comparison
    print(paste("Cond",n,"p(w|o)="))
    print(format(diag(coocs)/rowSums(coocs), digits=2))
    #print("p(o|w)=")
    #print(format(diag(coocs)/colSums(coocs), digits=2))
  }
}

ordering_stats(exp2x2)
ordering_stats(exp4x4)

load("asym_conds_fit_byS.RData") # all_fits[[condn]]
load("asym_conds_fit_byGroup.RData")


stats <- function(n, all_fits, w_by_o, o_freq, w_freq) {
  d = all_fits[[n]]
  d$CondNum = n
  d$Cond = w_by_o
  d$o_freq = o_freq
  d$w_freq = w_freq
  print(paste("Cond",n,"Human:",mean(d$Human),"Model:",mean(d$Model),"MLL:",mean(d$ML)))
  pdf(paste("plots/asym_",n,"_model_vs_human.pdf",sep=''))
  plot(jitter(d$Human), d$Model, xlab="Human Performance", ylab="Model Performance", xlim=c(0,1), ylim=c(0,1))
  abline(0,1, lty=2)
  text(.3,.9, paste(w_by_o,"Human:",format(mean(d$Human), digits=2), "Model:",format(mean(d$Model), digits=2)))
  dev.off()
  return(d)
}

prepare_dataframe <- function(all_fits) {
  afc18 = c(201,202,203,204,205,206,207,215:225)
  all <- data.frame() 
  for(i in 1:length(afc18)) {
    n = afc18[i]
    d <- stats(n, all_fits, w_by_o[i], o_freq[i], w_freq[i])
    all = rbind(all, d)
  }
  return(all)
}

all = prepare_dataframe(all_fits)

table(all$CondNum) # Ss per condition
with(subset(all, Human<.1), table(CondNum)) # below-chance Ss per cond

# old exp definitions:
#exp2x2 = c(201,202,204,205,206,207,225) # all extending 2x2, no within-cond manipulation (graph together)
#exp4x4 = c(217,218,220,221,223) # all built from 4x4 [um...maybe not]


agg <- aggregate(cbind(Human,Model) ~ Condition, data=all, mean)
sda <- aggregate(cbind(Human,Model) ~ Condition, data=all, sd)
agg$humSE <- sda$Human / sqrt(table(all$Condition)-1)
agg$modSE <- sda$Model / sqrt(table(all$Condition)-1)

simple_graph <- function(agg, name) {
  if(name=="exp4x4") {
    wid = 4.5
  } else if(name=="exp2x2") {
    wid = 6.5
  } else {
    wid = 9.3
  }
  require(ggplot2); require(gplots); 
  limits <- with(agg, aes(ymax=Human+humSE, ymin=Human-humSE))
  gg <- ggplot(agg, aes(Human, x=Condition)) + ylab("Proportion Correct") + xlab("Condition") + geom_bar(fill="blue") + geom_errorbar(limits,  width=0.2) + scale_y_continuous(limits=c(0,.85))
  # aes(Model, fill=CondNum, x=Cond)
  gg + geom_hline(aes(yintercept=1/18), linetype="dashed")
  ggsave(paste(name,"_human.pdf", sep=''), width=wid)
  dev.off()
  limits <- with(agg, aes(ymax=Model+modSE, ymin=Model-modSE)) 
  gg <- ggplot(agg, aes(Model, x=Condition)) + ylab("Proportion Correct") + xlab("Condition") + geom_bar(fill="red") + geom_errorbar(limits,  width=0.2) + scale_y_continuous(limits=c(0,.85))
  gg + geom_hline(aes(yintercept=1/18), linetype="dashed")
  ggsave(paste(name,"_model.pdf", sep=''), width=wid)
  dev.off()
}

model_vs_human <- function(d, expn) {
  pdf(paste(expn,"_model_vs_human.pdf",sep=''), width=4.5, height=4.5)
  print(xyplot(jitter(Human) ~ Model, groups=Condition, data=d, ylab="Human Performance", xlab="Model Performance", 
         xlim=c(-.02,1.02), ylim=c(-.02,1.02), panel = function(...) {
           panel.xyplot(...)
           panel.abline(0,1, lty=2) 
           panel.text(.3,.9, paste("R^2 =",format(cor(d$Human,d$Model), digits=2)))
           }))  
  dev.off()
}

agg$Condition <- as.factor(agg$Condition)
e2 = subset(agg, is.element(Condition, exp2x2c))
e4 = subset(agg, is.element(Condition, exp4x4c))
both = rbind(e2, e4)
both$Condition <- factor(both$Condition, levels=c("2x2","2x3","2x4","3x3 +1w/o","4x4 +2w/o","3x3","1x3","3x4","3x4 +6o","3x4 1/.5","3x4 1/.66","4x4","1x4"))
simple_graph(e2,"exp2x2")
simple_graph(e4,"exp4x4")
simple_graph(both,"both")
e2s <- subset(all, is.element(Condition, exp2x2c))
e4s <- subset(all, is.element(Condition, exp4x4c))
boths = rbind(e2s, e4s)
model_vs_human(e2s, "exp2x2")
model_vs_human(e4s, "exp4x4")
model_vs_human(boths, "both")

anova <- function(d, var) {
  d$Condition <- as.factor(d$Condition)
  d$Subject <- as.factor(d$Subject)
  m <- aov(d[,var] ~ Condition + Error(Subject/Condition), data=d)
  print(summary(m))
}

logreg <- function(d, var) {
  require(lme4)
  d$Condition <- as.factor(d$Condition)
  d$Subject <- as.factor(d$Subject)
  m <- glmer(paste(var,"~ Condition + (1|Subject)"), data=d, family=binomial)
  print(summary(m))
}

print_var_comparisons <- function(d) {
  load("asym_item_stats.RData")
  vars = c("Human", "X", "B", "C")
  for(v in vars) {
    print(v)
    logreg(e2s, v)
    logreg(e4s, v)
  }
}

# run once
merge_data <- function(all, modeln) {
  source(paste(modeln,".R",sep="")) 
  load("asym_item_stats.RData")
  load("asym_conds.RData")
  raw$ML <- NA
  raw$Model <- NA
  raw$Condition <- NA
  parcols = c("X","B","C")
  raw[,parcols] = NA
  for(i in 1:dim(all)[1]) {
    ord = conds[[ all[i,]$CondNum ]]$train
    rows = with(raw, which(Subject==all[i,]$Subject & CondNum==all[i,]$CondNum))
    mp = model(as.numeric(all[i,parcols]), ord)
    raw[rows,parcols] = all[i,parcols]
    raw[rows,]$Condition = all[i,]$Condition
    raw[rows,]$ML = all[i,]$ML
    mp = mp[1:18,]
    raw[rows,]$Model <- diag(mp) / rowSums(mp) 
  }
  save(raw, file="asym_master.RData")
  return(raw)
}

raw <- merge_data(all, "model")
#print_var_comparisons

print(glmer("Correct ~ CF + AE + CDit + maxOther + avgOther + testProb + wfreq + ofreq + (1|Subject)", family=binomial, data=raw))
# CDit +.165 avgOther +1.67 testProb 22.80 
print(glmer("Correct ~ wfreq + ofreq + testProb + (1|Subject)", family=binomial, data=raw))

# add trials per condition
tr_per_cond <- function(raw) {
  load("asym_conds.RData")
  afc18 = c(201,202,203, 204, 205,206,207,215:225)
  raw$Ntrials = NA
  for(n in afc18) {
    raw[which(raw$CondNum==n),]$Ntrials = length(conds[[n]]$train$trials)
  }
  return(raw)
}
#raw <- tr_per_cond(raw) # already added

print(aggregate(cbind(Human, X, B, C, ML) ~ Condition, data=e2s, mean))
print(aggregate(cbind(Human, X, B, C, ML) ~ Condition, data=e4s, mean))

print(aggregate(cbind(Human, X, B, C, ML) ~ Condition, data=subset(e2s, Human>.1), mean))
print(aggregate(cbind(Human, X, B, C, ML) ~ Condition, data=subset(e4s, Human>.1), mean))


run_once <- function(all_fits, afc18) {
  all = all_fits[[afc18[1]]]
  all$Cond = afc18[1]
  for(n in 2:length(afc18)) {
    d = all_fits[[afc18[n]]]
    d$Cond = afc18[n]
    all = rbind(all, d)
  }
  return(all)
}


#all <- run_once(all_fits, afc18)
#all = subset(all, Human>.1)
hist(all$Human)
pdf("cor_plot_asym.pdf")
plot(all[2:7], col=all$Cond)
dev.off()
cor.test(all$X, all$B) # -.24***
cor.test(all$B, all$C) # -.05
cor.test(all$C, all$X) # -.03
cor.test(all$C, all$Human) # .05
cor.test(all$B, all$Human) # -.24***
cor.test(all$X, all$Human) # .73***
cor.test(all$ML, all$Human) # -.30*** better fits for high perf?
cor.test(all$ML, all$X) # -.66***
cor.test(all$ML, all$C) # .26***
cor.test(all$ML, all$B) # .07

make_heatmaps <- function(condns) {
  outdir = "heatmaps/"
  load("asym_conds.RData") # conds[[num]]$train, $resps, $test
  for (n in condns) {
    train = conds[[n]]$train
    coocs = get_cooc_matrix(train) # for comparison
    #jpeg(paste(outdir,n,"_heatmap.jpg",sep=""))
    pdf(paste(outdir,n,"_heatmap.pdf",sep=""))
    heatmap(coocs, Rowv = NA, Colv = NA, symm=TRUE, scale="none", labRow=NA, margin=c(.2,.2), labCol=NA, xlab="", ylab="", col = heat.colors(10, alpha=.9))
    dev.off()
  }
}

get_cooc_matrix <- function(ord) {
  voc_sz = max(unlist(ord), na.rm=TRUE) # vocabulary size
  wfreq = rep(0, voc_sz)
  ofreq = rep(0, voc_sz)
  m <- matrix(0, nrow=voc_sz,ncol=voc_sz)
  for(t in 1:length(ord$trials)) { 
    tr_w = as.integer(ord$trials[[t]]$words)
    tr_o = as.integer(ord$trials[[t]]$objs)
    m[tr_w,tr_o] = m[tr_w,tr_o] + 1
    wfreq[tr_w] = wfreq[tr_w] + 1
    ofreq[tr_o] = ofreq[tr_o] + 1
  }
  print("word freq:")
  print(wfreq)
  print("obj freq:")
  print(ofreq)
  return(m)
}


make_heatmaps(afc18)