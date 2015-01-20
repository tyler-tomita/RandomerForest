Fig 1: Illustrative Comparisons of RF vs TF vs TF+ vs TF+d
-- we will prob remove TF in the final figs
-- n: # samples, D: # ambient dimensions, sigma: standard deviation (of D=1 when it is decaying)
-- errorbars are standard error of mean, so, run Ntrials large enough that one these errors are sufficiently small (as determined by eye).

for these 3, hold n & sigma fixed
1A) Trunk: Lhat vs. D
1B) XOR: Lhat vs. D
1C) Parity: Lhat vs. D

for these 3, hold n and D fixed
1D) Trunk: Lhat vs sigma
1E) XOR: Lhat vs. sigma
1F) Parity: Lhat vs sigma

Fig 2: same as fig 1, but replace y-axis with time instead of Lhat.

Supp Fig 1: # trees to stabilize per alg.  
why is this fixed? it seems like you must have some threshold of the variance of the outputs or something?  can you plot whatever that metric is as a function of #trees to make the case more clear.

Fig 3: each gaussian is correlated
-- here, do something like generate random means & covariance matrices per class.
-- to generate a random covariance matrix, eg: 
L= randn(D); Sigma=L'*L;
-- to start, let's just do 2 classes, each with 1 gaussian component
-- but then, let's do K classes, each with J gaussian components
-- not yet sure how to visualize when D>2, so we can start with D=2.
-- we might integrate this into Fig 1, we'll see how it turns out

Fig 4:  Robustness to scale & translation & rotation & outliers
﻿-- RobustLoveFest first passes to the marginal pdf for each variable independently prior to doing anything else, ie, as a pre-processing step.
﻿-- so, we can generate data according to some distribution (or use some real data), and (i) rotate randomly, (ii) also randomly shift each variable, (iii) also randomly rescale each variable, (iv) also randomly add points to each class (eg, sampled uniformly from the same).
-- thus, we get 5 panels, the raw data and including each of the transformations, and we plot Lhat vs D for  TF+delta both with & without the pre-processing step.



Fig 5: some real data plots
here, i am thinking for each real data plot, at first, i'd like to see 1 panel per dataset, and i'd like to see walltime vs. Lhat (on held-out data), with errorbars, say averaging over 10-fold cv.  let's run this code on abour 5 different datasets, and see how we do?


after that, i think we can start working on some theorems and/or finish writing stuff and submit??