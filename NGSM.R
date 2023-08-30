#### 2023/08/29 (new version)

library(pracma)
library(quadprog) 
library(lsei)
library(mixtools)

################## beta_estimation ##########################################

beta_esti <- function(x,v,y){

	x<-model.matrix(~., data=data.frame(x))
	tryCatch(
		return(solve(t(x) %*% v %*% x) %*% t(x) %*% v %*% y),
		error = function(e){
			library(MASS)
			return(ginv(t(x) %*% v %*% x)%*% t(x) %*% v %*% y)
		}
	)
}


################## normal scale mixture density ##########################################
normal.scale.mix=function(y,Q){
  p=nrow(Q)
  n=length(y)
  tmp=matrix(0, nrow = n, ncol = 1)
  for(j in 1:p){
    tmp=tmp+Q[j,2]*(dnorm(y,0,sqrt(Q[j,1]))+exp(-700))
  }
  return(tmp)
}
#########################################################################################
gradient.max=function(y,Q,weight=y/y){
  p=nrow(Q)
  tmp=NULL
  for(i in 1:p){
    tmp[i]=gradient.scale(y,Q[i,1],Q,weight)
  }
  return(max(tmp))
}
########################################################################################
###################Compute Gradient########################################################
gradient.scale=function(y,phi,Q,weight=y/y) {
  n=length(y)
  p=nrow(Q)
  
  tmp=matrix(0, nrow = n, ncol = p)
  for(j in 1:p){
    tmp[,j]=Q[j,2]*(dnorm(y,0,sqrt(Q[j,1]))+exp(-700))
  }
  
  #tmpp=matrix(dnorm(y,0,rep(sqrt(Q[,1]),each=n))+exp(-700),n,p)
  #tmp=tmpp%*%diag(Q[,2])
  
  tmp2=rowSums(tmp)
  grad=sum(weight*((dnorm(y,0,sqrt(phi))+exp(-700))/tmp2))-sum(weight)
  return(grad)
}
###################Compute Gradient########################################################
d_gradient.scale=function(y,phi,Q,weight=y/y) {
  n=length(y)
  p=nrow(Q)
  
  tmp=matrix(0, nrow = n, ncol = p)
  for(j in 1:p){
    tmp[,j]=Q[j,2]*(dnorm(y,0,sqrt(Q[j,1]))+exp(-700))
  }
  
  #tmpp=matrix(dnorm(y,0,rep(sqrt(Q[,1]),each=n))+exp(-700),n,p)
  #tmp=tmpp%*%diag(Q[,2])
  
  tmp2=rowSums(tmp)
  ttemp=(y^2-phi)*dnorm(y,0,sqrt(phi))/2/phi^2
  grad=sum(weight*(ttemp/tmp2))
  return(grad)
}

###################Compute Gradient########################################################
d2_gradient.scale=function(y,phi,Q,weight=y/y) {
  n=length(y)
  p=nrow(Q)
  
  tmp=matrix(0, nrow = n, ncol = p)
  for(j in 1:p){
    tmp[,j]=Q[j,2]*(dnorm(y,0,sqrt(Q[j,1]))+exp(-700))
  }
  
  #tmpp=matrix(dnorm(y,0,rep(sqrt(Q[,1]),each=n))+exp(-700),n,p)
  #tmp=tmpp%*%diag(Q[,2])
  
  tmp2=rowSums(tmp)
  ttemp1=(y^2-phi)*dnorm(y,0,sqrt(phi))/2/phi^2
  grad1=sum(weight*(ttemp1/tmp2))
  ttemp2=(y^4-6*y^2*phi+3*phi^2)*dnorm(y,0,sqrt(phi))/4/phi^4
  grad2=sum(weight*(ttemp2/tmp2))
  return(list(d1=grad1,d2=grad2))
}


###############Find new support############################################################
find.sigma=function(y,Q,con,weight=y/y) {
  #con=IQR(y)/20
  #print(con)
  #####Set grid
  max_range=max(100*IQR(y),2*(diff(range(y)))^2)
  min_range=con*1.01
  grid=seq(min_range,max_range,(max_range-min_range)/100)
  grid=c(grid,seq(grid[1],grid[2],(grid[2]-grid[1])/20))
  grid=sort(unique(c(con,grid,as.vector(Q[,1]))) )
  #######################################
  ## Detect sigma with positive gradient
  #######################################
  ran=NULL;  br=NULL; newx=NULL
  for (j in 1:length(grid)){
    #br[j]=(gradient.scale(y,grid[j]+1.0e-6,Q,weight)-gradient.scale(y,grid[j]-1.0e-6,Q,weight))/2.0e-6
    br[j]=d_gradient.scale(y,grid[j],Q,weight)
    if (j>1 && br[j-1]>0 && br[j]<0) {
      newx=c(newx,mean(grid[c(j-1,j)]))
      ran=rbind(ran,c(grid[c(j-1,j)]))
    }
  }
  #######################################
  # Find the local modes in gradient function
  #######################################
  final.can=NULL
  tmp_fn=function(x){return(-gradient.scale(y,x,Q,weight))}
  tmpp=NULL
  
  for (i in 1:length(newx)) {
    result=optim(newx[i],tmp_fn,lower=ran[i,1],upper=ran[i,2],method = "L-BFGS-B")
    if (result$value<(-0.1) && min(abs(result$par/Q[,1]-1))!=0) 
    {
      tmpp=c(tmpp,-result$value) #Gradient values
      final.can=c(final.can,result$par) #New support
    }
    if (gradient.scale(y,con,Q,weight)>0 && sum(Q[,1]==con)==0) {
      tmpp=c(tmpp,gradient.scale(y,con,Q,weight))
      final.can=c(final.can,con)
    }
  }
  if (length(final.can)>1) {
    grr=rbind(final.can,tmpp)
    grr=grr[,order(-tmpp)]
    return(grr[1,])
  } else { return(final.can) }
}
################ Update Weight ##############################################
update.weight=function(y,Q,weight=y/y){
  eps=1.1e-14
  if (size(Q,1)==1) {
    Q=t(as.matrix(as.vector(Q)))
    return(Q)
  }
  if (size(Q,2)==1) {
    Q=t(Q)
    return(Q)
  }
  
  n=length(y)
  p=nrow(Q)
  tmp=matrix(0, nrow = n, ncol = p)
  tmp11=tmp
  for(j in 1:p){
    tmp[,j]=dnorm(y,0,sqrt(Q[j,1]))+exp(-700)
    tmp11[,j]=Q[j,2]*tmp[,j]
  }
  
  tmp2=rowSums(tmp11)
  S=tmp/tmp2
  #S=S*matrix(rep(weight,p),n,p)
  #print(dim(S))
  S=diag(sqrt(weight))%*%S
  if (kappa(S)<3.0e+3)
  {
    new_weight=lsqlincon(S,2*as.matrix(sqrt(weight)),NULL,NULL,matrix(1,1,p),1,lb=0*Q[,1]+1.0e-14,ub=Q[,1]/Q[,1]-1.0e-14)
  } else {
    gam=sqrt(n*1.0e-6)
    Sc=rbind(gam*S,ones(1,p))
    onec=rbind(2*gam*as.matrix(sqrt(weight)),1)
    new_weight=nnls(Sc,onec)$x
    print("nnls is used instead of lsqlincon")
  }
  Q[,2]=t(new_weight/sum(new_weight))
  Q=Q[Q[,2]>eps,]
  if (size(Q,1)==1) {Q=t(as.matrix(Q))}
  if (size(Q,2)==1) {Q=t(Q)}
  Q[,2]=Q[,2]/sum(Q[,2])
  
  return(Q)
}
################# Main function ###########################################
fit.scale.mixture=function(y,censoring=y/y,maxitr=50,weight=y/y,con=IQR(y)/20,Q=NULL){
  
  eps=1.1e-14; lik=NULL
  if (is.null(Q))
  {
    weighted_variance=sum(weight*(y-sum(y*weight)/sum(weight))^2)/sum(weight)
    Q=matrix(c(con*5,weighted_variance,0.5,0.5),2,2)
  } 
  
  for (itr in 1:maxitr) {
    lik[itr]=log_lik(y,Q,weight=weight)
    old_Q=Q
    new_phi=unique(find.sigma(y,Q,con=con,weight=weight))
    if (length(new_phi)==0) {
      Q=update.weight(y,Q,weight);Q=Q[Q[,2]>eps,];Q=update.weight(y,Q,weight)
      return(list(Q=Q,lik=lik,lb=con,iter=itr,conv=TRUE))
    }
    if (length(unique(y))<=size(Q,1)+length(new_phi)) {new_phi=new_phi[1]}
    new_Q=cbind(as.matrix(new_phi),matrix(0.001,length(new_phi),1))
    Q=rbind(Q,new_Q) 
    for (i in 1:50){
      Old_Q=Q
      Q=update.weight(y,Q,weight)
      if (i>3 && gradient.max(y,Q,weight)<1.0e-10) break
      if (i>3) {
        Q=Q[Q[,2]>eps,]
        if (size(Q,1)==1) {Q=t(as.matrix(Q))}
        if (size(Q,2)==1) {Q=t(Q)}
      } 
    }
    #Q=Q[Q[,2]>eps,]
    #if (size(Q,1)==1 || size(Q,2)==1) {Q=t(as.matrix(Q))}
    Q=update.weight(y,Q,weight)
    #if (size(Q,1)==1) {Q=t(as.matrix(Q))}
  }
  Q=update.weight(y,Q,weight);Q=Q[Q[,2]>eps,];Q=update.weight(y,Q,weight)
  return(list(Q=Q,lik=lik,lb=con,iter=itr,conv=FALSE))
}
##########################################################################
log_lik=function(y,Q,weight=y/y){ return(sum(weight*log(normal.scale.mix(y,Q))))}
###################################################################
################################################################################
# mixture regression with nonparametric scale mixture error
################################################################
fmix_reg_scalemix=function(formula,data, m,p=NULL,beta=NULL,con=NULL,ini_sigma=NULL,maxitr=100, tRatio = 0.05){

	library(pracma)
	library(quadprog)  
	library(lsei)
	library(mixtools)
	library("RobMixReg")

	formula = as.formula(formula)
	dt = model.frame(formula, data = data)
	y = as.matrix(dt[,1])
	x = as.matrix(dt[,-1])
 
	# step1: initial value

 		 if (is.null(p) | is.null(beta) | is.null(con) ) regmix=TLE(formula, dt, nc = m, tRatio = tRatio, MaxIt = 200)
 		 if (is.null(p)) p=regmix@compcoef[(ncol(x)+3),]
 		 if (is.null(beta)){ 
				if(m != 1){
					beta=regmix@compcoef[1:(ncol(x)+1),]
				}else{					
					beta=matrix(regmix@compcoef[1:(ncol(x)+1),], ncol = 1)
				}
		 }
		 if (is.null(con)) con=regmix@compcoef[(ncol(x)+2),]/20
		 if (is.null(ini_sigma)) {
 			   tmpQ=rep(regmix@compcoef[(ncol(x)+2),],each=2)
 			   tmpQ[2*c(1:m)]=1
		  } else {
		    tmpQ=rep(ini_sigma,each=2)
  		    tmpQ[2*c(1:m)]=1
 		  }
	  beta_ini=beta

 
  n=length(y)
  
	# step2: E-step in the outer EM algorithm and CM-step for pi

  maxQ=50;Q=zeros(maxQ,(2*m))
  Q[1,]=tmpQ
  lik=NULL
  X=model.matrix(~., data=data.frame(x))

  for (i in 1:maxitr){
    for (j in 1:10){
      S=w=zeros(n,m)
      for (h in 1:m){
        tmp_Q=Q[rowSums(Q[,(2*h-1):(2*h)])>0,(2*h-1):(2*h)]
        if (is.vector(tmp_Q)) tmp_Q=t(tmp_Q)
        S[,h]=p[h]*normal.scale.mix(y-X%*%beta[,h],tmp_Q)
      }
      for (h in 1:m){
        w[,h]=S[,h]/rowSums(S)  # posterior
      }
      p=colMeans(w)  # pi

	# step2: E-step in the inner EM algorithm 
      
      z=zeros(n,m)
      for (h in 1:m){
        tmp_Q=Q[rowSums(Q[,(2*h-1):(2*h)])>0,(2*h-1):(2*h)]
        if (is.vector(tmp_Q)) tmp_Q=t(tmp_Q)
        d=dim(tmp_Q)[1]
        for (k in 1:d){
          z[,h]=z[,h]+normal.scale.mix(y-X%*%beta[,h],t(as.matrix(tmp_Q[k,])))/normal.scale.mix(y-X%*%beta[,h],tmp_Q)/tmp_Q[k,1]
        }
      }
      vw=z*w
      
	#step3: CM-step for beta

      for (h in 1:m){
	   	beta[,h]= beta_esti(x, diag(vw[,h]), y)
      }
      
    }


	  # CM-step for Q
    tmp_lik=zeros(n,1)

	 numer <- matrix(0, nrow = nrow(x), ncol = m)

    for (h in 1:m){
      tmp_Q=Q[rowSums(Q[,(2*h-1):(2*h)])>0,(2*h-1):(2*h)]
      if (is.vector(tmp_Q)) tmp_Q=t(tmp_Q)
      result=fit.scale.mixture(y-X%*%beta[,h],weight=w[,h],con=con[h],Q=tmp_Q,maxitr=10)
      dim_Q=dim(result$Q)
      Q[1:dim_Q[1],(2*h-1):(2*h)]=result$Q
      Q[(dim_Q[1]+1):maxQ,(2*h-1):(2*h)]=0

		numer[,h] = normal.scale.mix(y-X%*%beta[,h],result$Q)
      tmp_lik=tmp_lik+p[h]*numer[,h] #normal.scale.mix(y-X%*%beta[,h],result$Q)
    }
    lik[i]=sum(log(tmp_lik))
    if (i>1 && abs(lik[i]-lik[i-1])<1.0e-5) break
  }
  Q=Q[rowSums(Q)!=0,]
  Q = matrix(Q, ncol = 2*m) 

	support = list()
	num_support = NULL

	if(m != 1){
   	for(h in 1 : m){
		   support[[h]] = matrix(Q[,(2*h-1):(2*h)], ncol = 2)
   	   support[[h]] = support[[h]][support[[h]][,2] != 0,]
   	   num_support = append(num_support, nrow(support[[h]]))
		}
	}else{
		support[[1]] = Q
		num_support = nrow(Q)
	}

   if(m != 1){
		deno <- tmp_lik
		G <- matrix(0, nrow = nrow(x), ncol = (m - 1))

		for(h in 1 : (m-1) ){
			G[,h] <- (numer[,h]/deno) - 1
		}

		I <- matrix(0, nrow = nrow(x), ncol = 1)
		for(j in 1 : nrow(x)){
			tryCatch(
				I[j,]<- t(G[j,]) %*% solve(t(G) %*% G) %*% t(t(G[j,])),
				error = function(e){
				library(MASS)
				I[j,]<- t(G[j,]) %*% ginv(t(G) %*% G) %*% t(t(G[j,]))						
				}
			)
		}

		 max_order <- order(I, decreasing = T)
		dim <- (m-1)*(ncol(x)) #* (m-1)

			mix<-list()  #각 experts의 확률밀도값  (분자)
			for(h in 1 : m){
				mix[[h]] <- p[h]*numer[max_order[1:dim],h]
			}

			mix_Sum <- matrix(rep(0,dim),nrow= dim)  #각 experts의 확률밀도값의 합 (분모)
			for(h in 1 : m){
				mix_Sum <- mix_Sum + mix[[h]]
			}

		likep = lik[length(lik)] - sum(log(mix_Sum))
	}else{
		dim <- (ncol(x)) #* (m-1)
		likep = lik[length(lik)]
	}


  #BIC

	df = (m - 1) + m * nrow(beta) + sum(num_support) + (sum(num_support) - m)  
	Bic = -2 * lik[length(lik)] + df * log(n)

  #AIC
	Aic = -2 * lik[length(lik)] + df * 2

  #Icl
  	map = apply(w, 1, which.max)
	penalty = 0
   for(i in 1 : n){
		penalty = penalty + log(w[i,map[i]])
	}
	Icl = Bic - 2 * penalty

  #Final Results 	
  return(list(p=p,beta0=beta_ini,beta=beta,Q=support,con=con,loglik=lik, likep = likep, AIC = Aic, BIC = Bic, ICL = Icl, conv=i!=maxitr, z = w, cluster = map))
}
##############################################################################

