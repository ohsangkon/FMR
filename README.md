# FMR

A finite mixture of regressions is a valuable technique for situations where heterogeneous subpopulations exhibit distinct relationships between the response and covariates. 
We propose a novel approach that utilizes nonparametric Gaussian scale mixture error distributions to enhance the robustness of the modeling process.

(Oh, S., & Seo, B. (2023). Semiparametric mixture of linear regressions with nonparametric Gaussian scale mixture errors. Advances in Data Analysis and Classification, 1-27.)

# NGSM.R 
Basic code for finite mixture of regressions with nonparametric Gaussian scale mixture errors.
There are 11 functions in this file. 
fmix_reg_scalemix function is the main function, and others are required to save if you want to use fmix_reg_scalemix function.

Arguments

formula: formula between the response and covariates
data: input data
m: the number of components
p: mixing proportion
beta: regression coefficients
con: minimum of scale 
ini_sigma: initial value of sigma
maxitr: maximum iteration
tRatio: trimmed ratio for trimmed likelihood estimator. (trimmed likelihood estimators are used for initial value if initial values are not specified.) 

Values

p: mixing proportion
beta0: initial value for regression coefficients
beta: regression coefficients
Q: support points for sigma
con: minimum of scale 
loglik: log-likelihood
likep: p-deleted log-likelihood
AIC: Akaike information criterion
BIC: Bayesian information criterion
ICL: Integrated Complete-data Likelihood
Conv: convergence or not
z: posterior probabilities
cluster: clustering results based on posterior probabilities

# NGSM(multiple initial values).R
Based on NGSM.R, fmix_reg_scalemix_iter function chooses the solution with the highest d-deleted log-likelihood to prevent singular or spurious solutions.
