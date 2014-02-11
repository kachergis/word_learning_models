# Bayesian model of cross situational learning
# Stephen Denton, Apr. 20, 2010
rm(list=ls(all=TRUE))
graphics.off()

# Define noise probability (when alpha = 0, Bayesian model is a deterministic 'ideal observer')
alpha = .5 # 0.1, 0.5, 0.9
# Multiplier for word/object co-occurences (set to 1 for no increase)
delta = 1
chDec = 1

# Define number of words (equals number of objects)
nWords = 10 

ones = rep(1, nWords)
condPrior = rep(1/nWords, nWords)
pWgO = outer(ones,condPrior)
pW_O = outer(condPrior,condPrior)
dimnames(pWgO) = list(LETTERS[1:nWords], letters[1:nWords])
cat('\nPrior Conditional Probability of Words (rows) given Objects (columns)\n')
print(pWgO)

dimnames(pW_O) = dimnames(pWgO)
cat('\nPrior Joint Probability of Words (rows) and Objects (columns)\n')
print(pW_O)

words = c(1,1,1,1,0,0,0,0,0,0)
objs = c(1,1,1,1,0,0,0,0,0,0)

# Define a simple likelihood function that only updates appropriate rows and columns
likelihoodFun = function(words,objs) {
  ## With alpha=0, this enforces a mutual exclusivity constraint...
  mat = outer(words,objs) + outer(!words,!objs)
  ## Can relax mutual exclusivity contraint by increasing alpha or allowing for one of the following
  ## 1) Objects can have multiple words applied to them
  # mat = outer(words,objs) + outer(!words,ones)
  ## 2) Words can refer to multiple objects
  # mat = outer(words,objs) + outer(ones,!objs)
  likelihood = alpha^(1-mat) # + (delta-1) * outer(words,objs)
  return(likelihood)
  }

trnWords = matrix( c(
    c(1,1,1,1,0,0,0,0,0,0),  # A,B,C,D
    c(0,0,0,0,1,1,1,1,0,0),  # E,F,G,H
    c(1,1,0,0,0,0,0,0,1,1),  # A,B,I,J
    c(1,0,0,0,0,0,1,1,1,0)   # A,G,H,I
    ), ncol=10, byrow=T )
    
trnObjs = matrix( c(
    c(1,1,1,1,0,0,0,0,0,0),  # a,b,c,d
    c(0,0,0,0,1,1,1,1,0,0),  # etc.
    c(1,1,0,0,0,0,0,0,1,1),
    c(1,0,0,0,0,0,1,1,1,0)
    ), ncol=10, byrow=T )

# Train the model and update belief probabilities
for (trnIdx in 1:dim(trnWords)[1] ) {
    cat('\nTrial ', trnIdx, ': \n', sep='')
    likelihood = likelihoodFun( trnWords[trnIdx,], trnObjs[trnIdx,] )
    pWgO = likelihood * pWgO
    pWgO = pWgO/outer(ones,colSums(pWgO)) # rowSums(pWgO)
    cat('Conditional Probability of Words (rows) given Objects (columns)\n')
    print(pWgO)
    pW_O = likelihood * pW_O
    pW_O = pW_O/sum(pW_O)
    cat('Joint Probability of Words (rows) and Objects (columns)\n')
    print(pW_O)
}

cat('\nWord Choice Probabilities when object (a) is presented (calculated from joint):\n')
pWg_a = pW_O[,'a']^chDec/sum(pW_O[,'a']^chDec)
print(pWg_a)


cat('\nChoice Probabilities for each word given each object (using power choice rule):\n')
chProb = pWgO^chDec / outer(ones,colSums(pWgO^chDec))
print(chProb)

#cat('\nChoice Probabilities for each word given each object (using exp choice rule):\n')
#chProb = exp(chDec*pWgO) / outer(ones,colSums(exp(chDec*pWgO)))
#print(chProb)