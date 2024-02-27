#####################################################
#### mtcars dataset
#####################################################
## data load

data(mtcars)
dt = mtcars
##########################################################
# simple linear regression
##########################################################
plot(mpg ~ hp, data = dt, cex.lab = 1.5, cex.axis = 1, cex = 1.3)
abline(lm(mpg ~ hp, data = dt), col = "red")

slr= lm(mpg ~ hp, data = dt)
summary(slr)
plot(slr, cex = 1.3, cex.lab = 1.5 , which = c(1))

##########################################################
####Fitting linear regression for each engine
##########################################################
### vs = 0 

# linear regression under a normal error
first = lm(mpg ~ hp, dt[dt$vs == 0,])
summary(first)

par(mfrow = c(1,2))
plot(first, cex = 1.3, cex.lab = 1.5, which = 2 )
plot(first, cex = 1.3, cex.lab = 1.5, which = 4 )

# linear regression under a laplace error
first_lad = lad(formula = mpg ~ hp, data = dt[dt$vs == 0,], method = "BR" )
summary(first_lad)

### vs = 1
# linear regression under a normal error
second = lm(mpg ~ hp, dt[dt$vs == 1,])
summary(second)

par(mfrow = c(1,2))
plot(second, cex = 1.3, cex.lab = 1.5 , which = 2)
plot(second, cex = 1.3, cex.lab = 1.5 , which = 4)

second_lad = lad(formula = mpg ~ hp, data = dt[dt$vs == 1,], method = "BR" )
summary(second_lad)

##########################################################################
# linear regression under NGSM error
#################################################################################
ngsm1 = list()

ngsm1_BIC = NULL
ngsm1_ICL = NULL

time_ngsm1 = list()

set.seed(91304981)
for(m in 1 : 3){
	start_time <- Sys.time()
	ngsm1[[m]] = fmix_reg_scalemix_iter(formula = "mpg~hp", data = dt , m = m, iter = 5, p=NULL, beta=NULL,con=NULL,ini_sigma=NULL, maxitr=100, tRatio = 0.05)
	end_time <- Sys.time()
	time_ngsm1[[m]] = as.vector(end_time - start_time)/5 
	ngsm1_BIC = append(ngsm1_BIC, ngsm1[[m]]$BIC)	
	ngsm1_ICL = append(ngsm1_ICL, ngsm1[[m]]$ICL)	

}
ngsm1_BIC  # 2
ngsm1_ICL  # 2
time_ngsm1


#####################################################
