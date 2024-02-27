###############################################
### Hitters dataset
###############################################
## Data load
load(".../hitters_names.Rdata")

###############################################
# original dataset
dt = Hitters[,c("PutOuts", "AtBat", "position", "cluster")]
ncol(dt)

plot(PutOuts ~ AtBat, pch = NA_integer_, cex.lab = 1.5, data = dt) 
text(PutOuts ~ AtBat, position, data = Hitters[Hitters$cluster == 0,], cex = 1.2, col = "black")
text(PutOuts ~ AtBat, position, data = Hitters[Hitters$cluster == 1,], cex = 1.2, col = "red")

# dataset added with artificial outliers
dt3 = rbind(dt, c(1500, 150, NA, NA), c(2000, 300, NA,NA))
plot(PutOuts ~ AtBat, data = dt3)

plot(PutOuts ~ AtBat, pch = NA_integer_, cex.lab = 1.5, data = dt3) 
text(PutOuts ~ AtBat, position, data = Hitters[Hitters$cluster == 0,], cex = 1.2, col = "black")
text(PutOuts ~ AtBat, position, data = Hitters[Hitters$cluster == 1,], cex = 1.2, col = "red")
points(150,1500, pch = 8)
text(150, 1450, "Outlier")
points(300, 2000, pch = 8) 
text(300, 1950, "Outlier")

########################################################
### simple linear regression

plot(as.formula("PutOuts ~ AtBat"), data = Hitters, cex.lab = 1.5, cex.axis = 1, cex = 1.3)
abline(lm(PutOuts ~ AtBat, data = Hitters), col = "red")

slr= lm(PutOuts ~ AtBat, data = dt)
summary(slr)
plot(slr)

slr2= lm(PutOuts ~ AtBat, data = dt3)
summary(slr2)
plot(slr2)


#############################################################################
### Fitting NGSM
########################################################
## Original dataset
ngsm1 = list()
ngsm1_BIC = NULL
ngsm1_ICL = NULL

time_ngsm1 = list()

set.seed(81231)
for(m in 1 : 3){
	start_time <- Sys.time()
	ngsm1[[m]] = fmix_reg_scalemix_iter(formula = "PutOuts ~ AtBat", data = Hitters , m = m, iter = 5, p=NULL, beta=NULL,con=NULL,ini_sigma=NULL, maxitr=100, tRatio = 0.05)
	end_time <- Sys.time()
	time_ngsm1[[m]] = as.vector(end_time - start_time)/5 
	ngsm1_BIC = append(ngsm1_BIC, ngsm1[[m]]$BIC)	
	ngsm1_ICL = append(ngsm1_ICL, ngsm1[[m]]$ICL)	

}
ngsm1_BIC  # 2
ngsm1_ICL  # 2
time_ngsm1

ARI(ngsm1[[2]]$cluster, Hitters$cluster)
AMI(ngsm1[[2]]$cluster, Hitters$cluster)

#############
# Dataset added with outliers
ngsm2 = list()
ngsm2_BIC = NULL
ngsm2_ICL = NULL

time_ngsm2 = list()

set.seed(81231)
for(m in 1 : 3){
	start_time <- Sys.time()
	ngsm2[[m]] = fmix_reg_scalemix_iter(formula = "PutOuts ~ AtBat", data = dt3 , m = m, iter = 5, p=NULL, beta=NULL,con=NULL,ini_sigma=NULL, maxitr=100, tRatio = 0.05)
	end_time <- Sys.time()
	time_ngsm2[[m]] = as.vector(end_time - start_time)/5 
	ngsm2_BIC = append(ngsm2_BIC, ngsm2[[m]]$BIC)
	ngsm2_ICL = append(ngsm2_ICL, ngsm2[[m]]$ICL)	
	
}
ngsm2_BIC  # 2
ngsm2_ICL  # 2

ARI(ngsm2[[2]]$cluster[-c(323,324)], Hitters$cluster)
AMI(ngsm2[[2]]$cluster[-c(323,324)], Hitters$cluster)


########################################################

