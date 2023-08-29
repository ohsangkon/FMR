# FMR

A finite mixture of regressions is a valuable technique for situations where heterogeneous subpopulations exhibit distinct relationships between the response and covariates. 
We propose a novel approach that utilizes nonparametric Gaussian scale mixture error distributions to enhance the robustness of the modeling process.

# NGSM.R 
Basic code for finite mixture of regressions with nonparametric Gaussian scale mixture errors.

# NGSM(multiple initial values).R
Based on NGSM.R, we choose the solution with the highest d-deleted log-likelihood to prevent singular or spurious solutions.
