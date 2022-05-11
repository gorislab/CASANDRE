# CASANDRE
Code for CASANDRE ( Confidence-AS-A-Noisy-Decision-Reliability-Estimate) model
CASANDRE package

This package contains example template script code to process psychophysical data, fit the Confidence-AS-A-Noisy-DecisionReliabilityEstimate (CASANDRE) model to an example dataset, plot model fits to an example dataset, illustrate a recovery analysis, and illustrate model behavior. 

Included:

Scripts:
sortAdlerMa2018Expt1.m
fitAdlerMa2018.m
recoveryAnalysis.m
plotModelBehavior.m

Functions:
getLlhChoice.m
getLlhChoiceB.m

Data:
Adler_2018_Expt1-3.mat (Downloaded from the confidence database https://osf.io/s46pr/; Citation: WT Adler, WJ Ma. (2018). Comparing Bayesian and non-Bayesian accounts of human confidence reports. PLOS Computational Biology. 14(11): e1006572. https://doi.org/10.1371/journal. pcbi.1006572)

Adler_2018_Expt1_sorted.mat (Sorted data, original downloaded from the confidence database https://osf.io/s46pr/; Citation: WT Adler, WJ Ma. (2018). Comparing Bayesian and non-Bayesian accounts of human confidence reports. PLOS Computational Biology. 14(11): e1006572. https://doi.org/10.1371/journal. pcbi.1006572) This file contains CASANDRE model fits to data. [Note: Adler_2018_Expt1_sorted_ASYM_CC.mat contains CASANDRE fits with asymmetric confidence criteria.]


File description:
sortAdlerMa2018Expt1.m  
Template code to sort confidence data into one ‘trial’ structure per experiment. This specific example sorts data from Expt 1 of Adler and Ma (2018). Note: this template works well for data from the confidence database saved in .csv formatting.

fitAdlerMa2018.m	
Fit CASANDRE to empirical data. This specific example fits CASANDRE to Expt 1 from Adler & Ma (2018). 

recoveryAnalysis.m 
Generates synthetic data from CASANDRE and fits CASANDRE to these data.  

plotModelBehavior.m 
Illustrates key aspects of CASANDRE for a given set of model parameters including primary decision and confidence variable distributions, choice likelihood-functions, and psychometric and confidence functions split by confidence judgment or primary decision correctness respectively. 

getLlhChoice.m
Function used to calculate likelihood of each choice option in 2-AFC tasks with a single fixed decision criterion for a given stimulus strength. Works for an arbitrary number of confidence rating levels. (Note: this is the key piece of code needed to use model)

getLlhChoiceB.m
Function used to calculate likelihood of each choice option in the 2-AFC “task B” of Adler and Ma (2018) for a given stimulus strength. Works for an arbitrary number of confidence rating levels.



Link to preprint: https://www.biorxiv.org/content/biorxiv/early/2021/12/22/2021.12.17.473249.full.pdf
Read-me written by Zoe Boundy-Singer (zoebsinger@utexas.edu) 04/04/2022
