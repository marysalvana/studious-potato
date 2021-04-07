
workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/' 		else 		directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))

sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))

N <- 10
n <- N^2
TT <- 2
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

ind <- 14 
set.seed(ind)
PARAMETER_NONSTAT <- runif(15, -3, 3)

set.seed(3)
PARAMETER_DEFORMATION <- c(runif(2, -8, 8), runif(1, 0, 8), runif(2, -1, 1))

cat('Computing covariances...', '\n')

cov1 <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(1, 0.23, 1, 0.1001, 0.1001, 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT)

set.seed(1)
r1 <- rmvn(500, rep(0, ncol(cov1)), cov1, ncores = 25)

emp_cov <- cov(r1)

NN <- 50
grid_x <- seq(from = -4, to = 5, length.out = NN)
grid_y <- seq(from = -4, to = 5, length.out = NN)
X <- expand.grid(grid_x, grid_y) %>% as.matrix()
nn <- nrow(X)

NEGLOGLIK <- function(p){
	
	w <- c(0.1001, 0.1001)
	#w <- p[1:2]

	#cat(w, '\n')

	Xtarg <- Xtarg_orig <- sim_grid_locations 

	if (TT > 1){
		for (tt in 1:(TT - 1)){
			temp_locs <- cbind(Xtarg_orig[, 1] - tt * w[1], Xtarg_orig[, 2] - tt * w[2])
			Xtarg <- rbind(Xtarg, temp_locs)
		}
	}

	htarg <- as.matrix(dist(rbind(X, Xtarg),diag=TRUE,upper=TRUE))
	htarg <- htarg[1:nn, (nn + 1):(nn + nrow(Xtarg))]
	sigma <- htarg^2 * log(htarg)
	sigma[htarg == 0] <- 0

	beta1 <- matrix(p[1:length(jWarp)], ncol = ncol(sigma), nrow = length(jWarp))

	sub_sigma <- beta1 * sigma[1:length(jWarp), ]
	
	nonparametric_cov_est <- t(sub_sigma) %*% sub_sigma #+ diag(diag(emp_cov), ncol(sigma), ncol(sigma))

	err <- sum((emp_cov - nonparametric_cov_est)^2)

	return(err)
}

jWarp = 1:10
init <- c(rep(0, length(jWarp)))
fit <- optim(par = init, fn = NEGLOGLIK, control = list(trace = 5, maxit = 500))

#10 = 13317.656804

jWarp = 1:15
init <- c(fit$par, rep(0, length(jWarp) - length(fit$par)))
fit <- optim(par = init, fn = NEGLOGLIK, control = list(trace = 5, maxit = 1000))


jWarp = 1:10
init <- c(rep(0, length(jWarp)))
fit <- optim(par = init, fn = NEGLOGLIK, control = list(trace = 5, maxit = 10000))

jWarp = 1:15
init <- c(fit$par, rep(0, length(jWarp) - length(fit$par) + 2))
fit <- optim(par = init, fn = NEGLOGLIK, control = list(trace = 5, maxit = 1000))

reference_locations <- c(1, ceiling(n / 2) + ceiling(N / 2))

pdf(file = paste(root, 'Figures/4_estimation.pdf', sep = ''), width = 12, height = 10)

image.plot(matrix(cov1[55, ], N, N))

dev.off()



NEGLOGLIK <- function(p){
	
	Xtarg <- Xtarg_orig <- sim_grid_locations 

	htarg <- as.matrix(dist(rbind(X, Xtarg),diag=TRUE,upper=TRUE))
	htarg <- htarg[1:nn, (nn + 1):(nn + nrow(Xtarg))]
	sigma <- htarg^2 * log(htarg)
	sigma[htarg == 0] <- 0

	beta1 <- matrix(p[1:length(jWarp)], ncol = ncol(sigma), nrow = length(jWarp))

	sub_sigma <- beta1 * sigma[3 + 1:length(jWarp), ]
	
	nonparametric_cov_est <- t(sub_sigma) %*% sub_sigma #+ diag(diag(emp_cov), ncol(sigma), ncol(sigma))

	err <- sum((emp_cov[1:n, 1:n] - nonparametric_cov_est)^2)

	return(err)
}

jWarp = 1:10
init <- c(rep(0, length(jWarp)))
fit <- optim(par = init, fn = NEGLOGLIK, control = list(trace = 5, maxit = 500))

p <- fit$par

nonparametric_cov_est[1:5, 1:5]
emp_cov[1:5, 1:5]

set.seed(1234)
init <- runif(length(jWarp), -0.1, 0.1)
fit <- optim(par = init, fn = NEGLOGLIK, control = list(trace = 5, maxit = 500))

########  CHANG HSU



source(file = paste(root, "R_codes/cls.r", sep = ''))
library(mapproj)


###################################################
#### Read data
###################################################

loc1 <- sim_grid_locations
dist1<-as.matrix(dist(loc1))

A <- 1:nrow(loc1)

####################################################
#### function for finding alpha_l in stationary process W_l
####################################################
find.alpha<-function(Adist,ASig.em,p){
A.length<-dim(Adist)[1]
alpha.all<-NULL
for(i in 1:A.length){
	for(j in i:A.length){
	
		if (is.na(ASig.em[i,j])==0 & ASig.em[i,j]>0) al.temp<-(-Adist[i,j]/log(ASig.em[i,j]/sqrt(ASig.em[i,i]*ASig.em[j,j])))
    else  al.temp<-NA
		alpha.all<-c(alpha.all,al.temp)
		
	}
}
alpha.all<-alpha.all[is.na(alpha.all)==0]
alpha.set<-as.numeric(quantile(alpha.all,p))
return(alpha.set)
}



#####################################################
#### Gaussian kernal
#####################################################
sd.set<-c(1:4)/8
gg<-NULL
loc.mu<-NULL

for(i in 1:length(sd.set)){
vari<-sd.set[i]^2
loc.mui<-NULL
mu.m<-sd.set[i]
hex.nx<-floor((1+1.5*mu.m)/mu.m)
hex.ny<-floor((1+mu.m*sqrt(3)/2)/(mu.m*sqrt(3)/2))
for(k in -1:(hex.ny)){
	for(j in -1:(hex.nx-k%%2)){
	loc.jk<-c(mu.m*(j+(k%%2)/2),mu.m*k*sqrt(3)/2)
	loc.mui<-rbind(loc.mui,c(sd.set[i],loc.jk))
	temp.g<-diag((loc1[,1:2]-rep(loc.jk,each=n))%*%t(loc1[,1:2]-rep(loc.jk,each=n)))
 	temp.g<-exp(-0.5*temp.g/vari)
	gg<-rbind(gg,temp.g)
	}
	}
	loc.mu<-rbind(loc.mu,loc.mui)
}


####################################################
#### Gaussian kernal matrix for Lasso
####################################################
AM<-gg[,A]
k<-dim(AM)[1]
a<-dim(AM)[2]
si.fix<-NULL
si.temp1<-matrix(as.numeric(outer(c(1:a),c(1:a),"==")),a,a)
si.temp2<-NULL
for(s in 1:a){
	si.temp2<-c(si.temp2,si.temp1[s:a,s])

}
si.fix<-cbind(si.fix,si.temp2)

for(i in 1:k){
	si.temp1<-outer(AM[i,],AM[i,],"*")
	si.temp2<-NULL
for(s in 1:a){
	si.temp2<-c(si.temp2,si.temp1[s:a,s])
}
si.fix<-cbind(si.fix,si.temp2)
}
si.dim<-dim(si.fix)
si.fix<-matrix(as.numeric(si.fix),si.dim[1],si.dim[2])




###################################################
#### finding  alpha_l for W_l
###################################################
Sig.em<-cov1[1:n, 1:n]
p<-c(1:20)/20-1/40
alpha.set<-find.alpha(dist1[A,A],Sig.em,p)

Adist<-dist1[A,A]

###################################################
#### matrix for lasso
###################################################
si<-NULL
for(i in 1:length(alpha.set)){
	exp.temp<-NULL
	for(s in 1:length(A)){
			temp<-exp(-abs(Adist[s:a,s])/alpha.set[i])
			exp.temp<-c(exp.temp,temp)
		}
	si<-cbind(si,exp.temp)
	}
si.dim<-dim(si)
si<-matrix(as.numeric(si),si.dim[1],si.dim[2])
si<-cbind(si.fix,si)



AY.lars<-NULL
loc.y<-NULL
for(i in 1:length(A)){
		AY.lars<-c(AY.lars,Sig.em[(i:a),i])
		loc.y<-rbind(loc.y,cbind(rep(i,(a-i+1)),c(i:a)))
	}
id.train<-c(1:length(AY.lars))[is.na(AY.lars)==0]



###############################################
#### CLS
###############################################
si.re<-si[id.train,]
BY.lars<-cbind(loc.y[id.train,],AY.lars[id.train])
test1<-larspositive(si.re,AY.lars[id.train],type="lasso")
temp1<-cv.lars.positive(si.re,BY.lars)
s.temp<-temp1$fraction[which.min(temp1$cv)]
pre1<-predict.larspositive(test1,s=s.temp,type="coefficient", mode="fraction")$coefficients

#The estimates might be different from those reported in the paper due to the CV method 
#randomly partition the points into K-folds. 
#If you want to get the same estimates as those in our paper, 
#please set s.temp<-0.5050505.


pre1[1]  #parameter estimate:error
M.n<-dim(gg)[1]
est.basis<-cbind(loc.mu,pre1[2:(M.n+1)])
est.basis[est.basis[,4]>0,]   #parameter estimate:basis function
est.stat<-cbind(alpha.set,pre1[-c(1:(M.n+1))])
est.stat[est.stat[,2]>0,] #parameter estimate: stationary component

####est. Sigma by CLS
M<-gg
M.n<-dim(gg)[1]
Sig.g<-(t(M)%*%diag(pre1[(2:(M.n+1))])%*%M)+diag(pre1[1],n)

Sig.exp<-matrix(0,n,n)
for(i in 1:length(alpha.set)){
	Sig.expt<-exp(-dist1/alpha.set[i])*pre1[-c(1:(M.n+1))][i]
	Sig.exp<-Sig.exp+Sig.expt
  }
Sig.CLS<-Sig.g+Sig.exp

