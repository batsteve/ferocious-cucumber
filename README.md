# ferocious-cucumber
Toy implementation of QRS surface metrics for extreme event prediction quality

Accompanying documentation has been submitted for publication.  

Features:
 - two toy models:  a bimodal distribution and a multivariate Gaussian distribution
 - two dynamical systems:  the Majda-McLaughlin-Tabak model and the Kolmogorov Flow model
 - produces both pdf and QRS plots
 - outputs either analytic results, or simulates sampled results

Instructions:

Run main_sample.m in the latest version of Matlab to produce QRS surface plots for the toy models.

Run /Review Branch/main_optimization.m for the MMT analysis code, and /Review Branch/main_kolm_opt.m for the Kolmogorov analysis code.  These scripts will require large sample data sets, simulate the datasets by running main.m in each corresponding subfolder.
