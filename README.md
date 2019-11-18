# Word Learning Models

George Kachergis, December 12, 2013
george.kachergis@gmail.com
http://www.kachergis.com/academic.php

This repository contains a number of learning models that I've implemented for cross-
situational word learning research. I am slowly trying to collect implementations 
from over the years (along with data), so please bear with me as I standardize and 
document functions. Most of these models and data have appeared in at least one of my 
publications, which can all be found on my website: http://www.kachergis.com

Please let me know if you have any questions, or if you have found this collection
useful please cite my dissertation:

Kachergis, G. (2012). Mechanisms for Cross-Situational Learning of Word-Referent Mappings: Empirical and Modeling Evidence. http://kachergis.com/docs/kachergis_thesis.pdf
(précis: http://kachergis.com/docs/kachergis_precis.pdf)


## Overview

main.R - contains functions to load training trial orderings, and load and fit models

graphics.R - contains functions to graph word-object co-occurrence matrices and animate model performance trajectories

orderings/ 
	- contains many .txt files that have trial orders in the format 1 line per trial, and each number on a line represents the appearance of a word-object pair (e.g., 4 means both word 4 and object 4 occurred on that trial)

models/
	- cotains model implementations in R.
	
## Specific Included Conditions 

- asymmetric_conditions.RData 
	- contains the asymmetric (i.e., unequal number of words and objects per trial) trial orderings reported and modeled in: 
	Kachergis, G. & Yu, C. (2013). More naturalistic cross-situational word learning. Proceedings of the 35th Annual Conference of the Cognitive Science Society. http://kachergis.com/docs/kachergis_yu_2013_asym.pdf
	conds[["number"]] (e.g., conds[["225"]] below) will return a list "train" that contains two matrices "words" and "objs" with rows containing the words and objects seen on each training trial. Here are the conditions, and a short description of each:
	225, 2x2 (words x objects per trial).  Words appear 6 times. objects 6. Participant testing was 18-alternative forced choice (18afc).
	205, 2x4.  extending the original 2x2. Each word appears 6 times.  Each object appears 12 times.  18afc
	206, 3x3.  extending the original 2x2. 2x2 normal word/referent. plus 1 random word and 1 random object.  18afc
	207, 4x4.  extending the original 2x2. 2x2 normal word/referent. plus 2 random words and 2 random objects.  18afc
	201, 3x4 condition. extending the original 2x2 [3x3??]. each word appears 6 times. each object 8 times. Therefore p(w|o) = 80% for all 18 pairs. 18afc
	202, 3x4 condition.  extending the original 2x2 [3x3??]. Each word appears 6 times.  First 12 objects appear only with their words p(w|o) = p(o|w) = 1. objs 13-18 appear 12 times. p(w|o) = 50%.  18afc
	203, 3x4.  extending the original 2x2 [3x3??]. All words appear 6 times.  Objects appear with their words p(w|o) = p(o|w) = 1. and objs 3-8 and 13-18 also appear 3 more times p(w|o) = 66%.
	204, 3x4.  All 18 words appear 6 times. [extending 3x3??]  There are 24 objects. each appears 6 times.  The first 18 are referents. the last 6 are not the referents of any particular word. and therefore noisy items18afc.
	203, objs 3-8, and 13-18 appear 9x; others appear 6x
	204, objs 19-24 are noise
	215, 1x3.  Words appear 6 times. objects 18.  18afc. 23 Ss
	222, 1x3.  Words appear 6 times. objects 18.  18afc. 40 Ss
	216, 2x3.  Words appear 6 times. objects 9.  18afc. 23 Ss
	219, 2x3.  Words appear 6 times. objects 9.  18afc 32 Ss
	221, 1x4.  Words appear 6 times. objects 24.  18afc. 40 Ss ~.18
	223, 4x4.  Words appear 6 times. objects also 6.  18afc. 77 Ss
	224, 3x3.  Words appear 6 times. objects 6.  18afc. 36 Ss
	217, 2x4.  Words appear 6 times. objects 12.  18afc. [not built from 2x2] ~.3 (14 Ss)
	218, 3x4.  Words appear 6 times. objects 8.  18afc. [not built from 2x2] ~.43 (13 Ss)
	220, 3x4.  Words appear 9 times. objects 12.  18afc - high perf: ~.67 (33 Ss)

- master_orders.RData - similar structure to "asymmetrical_conditions.RData", but contain training trial orderings for several other papers: 
	Kachergis, G., Yu, C. & Shiffrin, R. M. (submitted). A bootstrapping model of frequency and contextual diversity effects in word learning.
	Kachergis, G., Yu, C., & Shiffrin, R. M. (2012). An Associative Model of Adaptive Inference for Learning Word-Referent Mappings. Psychonomic Bulletin & Review, 19(2), 317-324.

## Models

	- fazly.R - incremental probabilistic model from Fazly, A., Alishahi, A., and Stevenson, S. (2010). A probabilistic computational model of cross-situational word learning, Cognitive Science, 34(6): 1017-1063.
	
	- kachergis.R - Strength- and Uncertainty-biased Model from Kachergis, G., Yu, C., & Shiffrin, R.M. (2012). An Associative Model of Adaptive Inference for Learning Word-Referent Mappings. Psychonomic Bulletin & Review, 19(2), 317-324.

	- novelty.R - like the kachergis.R model, but uses novelty (inverse stimulus frequency) in lieu of associate entropy

	- uncertainty.R - like the kachergis.R model, but without the strength/familiarity bias

	- strength.R - like the kachergis.R model, but without the entropy/uncertainty bias

	- Bayesian_decay.R - a Bayesian model that strengthens all plausible (i.e., co-occurring) pairings and weakens all implausible associations based on strengthening and decay parameters

	- MINERVA2.R - adaptation of the MINERVA 2 episodic memory model from Hintzman, D.L. (1984). MINERVA 2: A simulation model of human memory. Behavior Research Methods, 16(2), 96-101.

	- tilles.R - (possibly flawed) implementation of model from Tilles, P.F.C. & Fontanari, J.F. (2013). Reinforcement and inference in cross-situational word learning. Frontiers in Behavioral Neuroscience. doi: 10.3389/fnbeh.2013.00163

	- TODAM2.R - a few different adaptations of the TODAM 2 model from Murdock, B.B. (1997). Context and mediators in a theory of distributed associative memory (TODAM2). Psychological Review, 104, 839-862.

	- rescorla-wagner.R - Rescorla, R.A. & Wagner, A.R. (1972) A theory of Pavlovian conditioning: Variations in the effectiveness of reinforcement and nonreinforcement, Classical Conditioning II, A.H. Black & W.F. Prokasy, Eds., pp. 64–99. Appleton-Century-Crofts.


