# CogSci_analysis2 - with item-level data
# after merge_data in CogSci_analysis

load("asym_master.RData")
raw$Condition = ifelse(raw$Condition=="3x4 .75", "3x4", ifelse(raw$Condition=="3x4 +6objs", "3x4 +6o", raw$Condition))

# the two conditions with pair groups:
subset(raw, Condition=="3x4 1/.5") # first 12 words appear 6x, last 6 appear 12x
subset(raw, Condition=="3x4 1/.66") # objs 3-8 and 13-18 appear 9x, others 6x

raw$pOgW = raw$wfreq / raw$ofreq

agg <- aggregate(cbind(Correct,Model) ~ Condition + Item, data=raw, mean)

model_vs_human <- function(d, expn) {
  pdf(paste(expn,"_model_vs_human_item.pdf",sep=''), width=4.5, height=4.5)
  print(xyplot(jitter(Correct) ~ Model, groups=Condition, data=d, ylab="Human Performance", xlab="Model Performance", 
               xlim=c(-.02,1.02), ylim=c(-.02,1.02), panel = function(...) {
                 panel.xyplot(...)
                 panel.abline(0,1, lty=2) 
                 panel.text(.3,.9, paste("R^2 =",format(cor(d$Correct,d$Model), digits=2)))
               }))  
  dev.off()
}

agg$Condition <- as.factor(agg$Condition)
e2 = subset(agg, is.element(Condition, exp2x2c))
e4 = subset(agg, is.element(Condition, exp4x4c))
model_vs_human(e2, "exp2x2")
model_vs_human(e4, "exp4x4")
both = rbind(e2, e4)
model_vs_human(both, "both")

# gotta get subject groups correct *whew*
table(raw$CondNum)
raw$Subject = with(raw, ifelse(CondNum==201 | CondNum==202 | CondNum==203, Subject+100, ifelse(CondNum==215 | CondNum==216, Subject+200, ifelse(CondNum==205 | CondNum==220, Subject+300, ifelse(CondNum==206 | CondNum==207, Subject+400, ifelse(CondNum==221 | CondNum==222, Subject+500, ifelse(CondNum==223, Subject+600, ifelse(CondNum==224, Subject+700, ifelse(CondNum==219, Subject+800, ifelse(CondNum==204, Subject+900, ifelse(CondNum==225, Subject+1000, Subject)))))))))))

e2s = subset(raw, is.element(Condition, exp2x2c))
e4s = subset(raw, is.element(Condition, exp4x4c))
both = rbind(e2s, e4s)

# first report this for Exp 1
#print(glmer("Correct ~ Ntrials + wfreq + ofreq + pOgW + testProb + (1|Subject)", family=binomial, data=e2s))
#print(aggregate(cbind(Correct, Ntrials, wfreq, ofreq, pOgW, testProb) ~ Condition, data=e2s, mean))
# then report this for Exp 2
#print(glmer("Correct ~ Ntrials  + testProb + wfreq + ofreq + pOgW + (1|Subject)", family=binomial, data=e4s))
#print(aggregate(cbind(Correct, testProb, wfreq, ofreq, pOgW) ~ Condition, data=e4s, mean))

#print(aggregate(cbind(Correct, X, B, C, ML) ~ Condition, data=e2s, mean))
#print(aggregate(cbind(Human, X, B, C, ML) ~ Condition, data=e4s, mean))

# now just doing one big Exp...maybe.
#print(glmer("Correct ~ Ntrials + wfreq + ofreq + pOgW + testProb + (1|Subject)", family=binomial, data=both))
# Ntrials***, ofreq*, testProb***
both_sc <- both
both_sc$Ntrials = scale(both$Ntrials)
both_sc$wfreq = scale(both$wfreq)
both_sc$ofreq = scale(both$ofreq)
both_sc$pOgW = scale(both$pOgW)
both_sc$testProb = scale(both$testProb)
both_sc$CF = scale(both$CF)
both_sc$AE = scale(both$AE)
both_sc$CDit = scale(both$CDit)
print(glmer("Correct ~ Ntrials + wfreq + ofreq + pOgW + testProb + CDit + CF + AE + (1|Subject)", family=binomial, data=both_sc))


print(glmer("Correct ~ Ntrials + wfreq + ofreq + pOgW + testProb + CF + AE + CDit + (1|Subject)", family=binomial, data=both))
# Ntrials, ofreq, testProb, CF, AE, CDit
print(format(aggregate(cbind(Correct, wfreq, ofreq, pOgW, testProb, CF, AE, CDit) ~ Condition, data=both, mean), digits=2))

print(glmer("Correct ~ CF + AE + CDit + Ntrials  + testProb + wfreq + ofreq + pOgW + (1|Subject)", family=binomial, data=both))

aggs <- aggregate(cbind(Correct, X, B, C, AE, CF, CDit, testProb) ~ Condition + Subject, data=both, mean)
format(cor(aggs[,c("Correct","X","B","C","AE","CF","CDit","testProb")]), digits=1)

aggs <- aggregate(cbind(Correct, X, B, C, ML) ~ Condition + Subject, data=both, mean)

format(cor(aggs[,c("Correct","X","B","C")]), digits=1)
cor.test(aggs$X, aggs$Correct) # Pearson's r=.72 t(494)=22.74, p<.001
cor.test(aggs$B, aggs$Correct) # Pearson's r=-.22 t(494)=-5.04, p<.001
cor.test(aggs$X, aggs$B)

summary(aov(X ~ Condition + Error(Subject), data=aggs)) # (F(12,482) = 11.63, p<.001)
summary(aov(B ~ Condition + Error(Subject), data=aggs)) # (F(12,482) = 2.13, p=.01)
summary(aov(C ~ Condition + Error(Subject), data=aggs)) # (F(12,482) = 2.70, p<.01)

summary(lm("X~Condition", data=aggs))

print(format(aggregate(cbind(Correct, X, B, C) ~ Condition, data=aggs, mean), digits=2))

print(format(aggregate(cbind(Correct, X, B, C) ~ Condition, data=subset(aggs, Correct>.1), mean), digits=2))



# maxOther + avgOther
#table(both$wfreq, both$Condition) # 6 in all but 3x3 +1w/o (9) and 4x4 +2w/o (12)
table(both$ofreq, both$Condition)
table(both$Ntrials, both$Condition)
table()
print(glmer("Correct ~ CF + AE + CDit + Ntrials + avgOther + testProb + (1|Subject)", family=binomial, data=e4s))


