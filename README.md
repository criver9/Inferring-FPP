# Inferring-FPP

This code inferes phenomenological models of First Passage processes using the multi-path model family introduced in [1]. Completion times are required to implement the software. 

The files of this repository are as follows:

1. The file MultiPathModelDetermMainGit.m finds the best fits (Model parameters and Likelihood) for the first M models in the multi-path model family.

2. The file Hessian Determ Analytics Synthetic Purkinje.nb estimates the Hessian of the posterior distribution of a given model at the optimal fit estimated in the    step 1.
3. ImportanceSamplingGitDic2020.m estimated the log-marginal likelihood for a the selected model.

Note: 1) Steps 2 and 3 need to be evaluated for each model in the hierarchy, until the simplest model that represents the data is found. 2) If optimal values fall at the boundary of the parameter domain Hessian needs to be manually modified as explain in Material and Methods of [1].

# Figure's Generation

1. The file GenerateFig3AutoGit.m generates figures similar to Fig 3 and Fig 4 in [1].
   Note: The file uses entropy estimations for the PDF distributions inferred. The Bayesian entropy estimation was based on the method proposed in [2]. Here the 
   link to the code:
   http://nsb-entropy.sourceforge.net/
   
2. The file GenerateFig3AutoGit.m generates figures similar to Fig5 in [1]. 






# References

1. C. Rivera, D. Hofmann and I. Nemenman. Inferring phenomenological models of first passage processes. Submitted, 2020.
2. I. Nemenman, F. Shafee and W. Bialek. Entropy and inference revisited. Advances in Neural Information Processing Systems, 2002.

