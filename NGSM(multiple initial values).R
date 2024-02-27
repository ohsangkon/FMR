###########################################################
### NGSM
###########################################################

fmix_reg_scalemix_iter =function(formula, data , m = 2, iter = 5, p=NULL,beta=NULL,con=NULL,ini_sigma=NULL,maxitr=100, tRatio = 0.05){
	### Find maximum likelihood reg_scalemix  ###
	model = list()
	likelihood <- c()

	if(m != 1){
		for(idx in 1 : iter){
	tryCatch(
			model[[idx]] <- fmix_reg_scalemix(formula, data,m, p, beta, con,ini_sigma,maxitr=100, tRatio = tRatio),
		error = function(e){
			model[[idx]] <- fmix_reg_scalemix(formula, data,m, p, beta, con,ini_sigma,maxitr=100, tRatio = tRatio)
		}

	)
			likelihood[idx] = model[[idx]]$likep #loglik[length(model[[idx]]$loglik)]
			print(idx)
		}
		return(model[[which.max(likelihood)]])
	}else{
  		model = fmix_reg_scalemix(formula, data,m, p, beta, con,ini_sigma,maxitr=100, tRatio = tRatio)
		return(model)
	}

	
}
