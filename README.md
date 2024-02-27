# FMR

A finite mixture of regressions is a valuable technique for situations where heterogeneous subpopulations exhibit distinct relationships between the response and covariates. 
We propose a novel approach that utilizes nonparametric Gaussian scale mixture error distributions to enhance the robustness of the modeling process.

(Oh, S., & Seo, B. (2023). Semiparametric mixture of linear regressions with nonparametric Gaussian scale mixture errors. Advances in Data Analysis and Classification, 1-27.)

# NGSM.R 
Basic codes for finite mixture of regressions with nonparametric Gaussian scale mixture errors.
There are 11 functions in this file. 
fmix_reg_scalemix function is the main function, and others are required to save if you want to use fmix_reg_scalemix function.

Arguments

1. formula: formula between the response and covariates
2. data: input data
3. m: the number of components
4. p: mixing proportion
5. beta: regression coefficients
6. con: minimum of scale
7. ini_sigma: initial value of sigma
8. maxitr: maximum iteration
9. tRatio: trimmed ratio for trimmed likelihood estimator. (trimmed likelihood estimators are used for initial value if initial values are not specified.) 

Values

1. p: mixing proportion
2. beta0: initial value for regression coefficients
3. beta: regression coefficients
4. Q: support points for sigma
5. con: minimum of scale
6. loglik: log-likelihood
7. likep: p-deleted log-likelihood
8. AIC: Akaike information criterion
9. BIC: Bayesian information criterion
10. ICL: Integrated Complete-data Likelihood
11. Conv: convergence or not
12. z: posterior probabilities
13. cluster: clustering results based on posterior probabilities

# NGSM(multiple initial values).R
Based on NGSM.R, fmix_reg_scalemix_iter function chooses the solution with the highest d-deleted log-likelihood to prevent singular or spurious solutions.
