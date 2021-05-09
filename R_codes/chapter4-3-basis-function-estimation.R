
workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/' 		else 		directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
library(future.apply)
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))

sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))
sourceCpp(file = paste(root, "R_codes/Functions/distR.cpp",sep=''))

N <- 20
n <- N^2
TT <- 1
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

########  CHANG HSU

source(file = paste(root, "R_codes/cls.r", sep = ''))


###################################################
#### Read data
###################################################

loc1 <- loc2 <- sim_grid_locations
if(TT > 1){
	for(tt in 1:(TT - 1)){
		loc2 <- rbind(loc2, cbind(sim_grid_locations[, 1] + 0.1001 * tt, sim_grid_locations[, 2] + 0.1001 * tt))
	}
	loc1 <- loc2
}

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
		cat(i / A.length, ', ', j / A.length, '\n')	
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
#sd.set<- seq(0.06, 0.5, length.out = 5)
nu.set <- c(0.5, 1)
gg<-NULL
loc.mu<-NULL
for(i in 1:length(sd.set)){
	cat(i, '\n')
	vari<-sd.set[i]^2
	loc.mui<-NULL
	mu.m<-sd.set[i]
	hex.nx<-floor((1+1.5*mu.m)/mu.m)
	hex.ny<-floor((1+mu.m*sqrt(3)/2)/(mu.m*sqrt(3)/2))
	anchor_locs <- expand.grid(-2:(hex.ny), -2:(hex.nx)) %>% as.matrix()
	loc.jk <- cbind(mu.m * (anchor_locs[, 2] + (anchor_locs[, 1] %% 2) / 2), mu.m * anchor_locs[, 1] * sqrt(3)/2)
	loc.mui <- cbind(rep(sd.set[i], nrow(loc.jk)), loc.jk)
	temp.g <- distR_C(loc.jk, loc1)
	temp.g<-exp(-0.5*temp.g/vari)
	gg<-rbind(gg,temp.g)
	loc.mu<-rbind(loc.mu,loc.mui)
}

shift <- 0
for(nu_val in 1:length(nu.set)){
	shift <- shift - 0.1
	for(i in 1:length(sd.set)){
		cat(i, '\n')
		vari<-sd.set[i]^2
		loc.mui<-NULL
		mu.m<-sd.set[i]
		hex.nx<-floor((1+1.5*mu.m)/mu.m)
		hex.ny<-floor((1+mu.m*sqrt(3)/2)/(mu.m*sqrt(3)/2))
		anchor_locs <- expand.grid(-1:(hex.ny), -1:(hex.nx)) %>% as.matrix()
		loc.jk <- cbind(mu.m * (anchor_locs[, 2] + (anchor_locs[, 1] %% 2) / 2) + shift, mu.m * anchor_locs[, 1] * sqrt(3)/2 + shift)
		loc.mui <- cbind(rep(sd.set[i], nrow(loc.jk)), loc.jk)
		temp.g <- distR_C(loc.jk, loc1)
		temp.g<- Matern(temp.g, range = vari, nu = nu.set[nu_val])
		gg<-rbind(gg,temp.g)
		loc.mu<-rbind(loc.mu,loc.mui)
	}
}

ind <- which(loc.mu[, 2] > 0.5 | loc.mu[, 3] > 0.5)
gg <- gg[-ind, ]

####################################################
#### Gaussian kernal matrix for Lasso
####################################################

AM<-gg[,A]
k<-dim(AM)[1]
a<-dim(AM)[2]
si.fix<-NULL
si.temp1<-matrix(as.numeric(outer(c(1:a),c(1:a),"==")),a,a)
si.temp2<-NULL
si.temp3 <- matrix(1:(a^2), a, a)
si.temp4<-NULL
for(s in 1:a){
	cat('s = ', s/a, '\n')	
	si.temp2<-c(si.temp2,si.temp1[s:a,s])
	si.temp4<-c(si.temp4,si.temp3[s:a,s])

}
si.fix<-cbind(si.fix,si.temp2)

si.fix_new <- future_apply(AM, 1, function(x) x %o% x )
si.fix <- cbind(si.fix, si.fix_new[si.temp4, ])
si.dim<-dim(si.fix)
si <- si.fix<-matrix(as.numeric(si.fix),si.dim[1],si.dim[2])

Sig.em<-cov1

AY.lars<-NULL
loc.y<-NULL
for(i in 1:length(A)){
		cat('i = ', i/length(A), '\n')	
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


####est. Sigma by CLS
M<-gg
M.n<-dim(gg)[1]
Sig.CLS <- Sig.g<-(t(M)%*%diag(pre1[(2:(M.n+1))])%*%M) +diag(pre1[1],n)

###################################################################################################################

reference_locations <- c(1, ceiling(n / 2) + ceiling(N / 2))

pdf(file = paste(root, 'Figures/4-basis-function-estimation.pdf', sep = ''), width = 15, height = 15)

par(mfrow = c(2, 2))

par(mai=c(1.5, 1.5, 1.5, 1.5))
par(pty="s") 
image.plot(matrix(cov1[reference_locations[2], 1:n], N, N), zlim = c(0, max(Sig.CLS)))
mtext('Reference Location 1', side = 2, line = 3, cex = 1.5, font = 2, col="#000080")
mtext('Theoretical Covariance at t = 1', side = 3, line = 2, cex = 1.5, font = 2, col="#000080")

par(mai=c(1.5, 1.5, 1.5, 1.5))
par(pty="s") 
image.plot(matrix(Sig.CLS[reference_locations[2], 1:n], N, N), zlim = c(0, max(Sig.CLS)))
mtext('Estimated Covariance at t = 1', side = 3, line = 2, cex = 1.5, font = 2, col="#000080")

par(mai=c(1.5, 1.5, 1.5, 1.5))
par(pty="s") 
image.plot(matrix(cov1[reference_locations[1], 1:n], N, N), zlim = c(0, max(Sig.CLS)))
mtext('Reference Location 2', side = 2, line = 3, cex = 1.5, font = 2, col="#000080")

par(mai=c(1.5, 1.5, 1.5, 1.5))
par(pty="s") 
image.plot(matrix(Sig.CLS[reference_locations[1], 1:n], N, N), zlim = c(0, max(Sig.CLS)))

dev.off()


pdf(file = paste(root, 'Figures/4-basis-function-location.pdf', sep = ''), width = 15, height = 15)
plot(loc.mu[, 2:3])
dev.off()

NEGLOGLIK <- function(p){
	
	Sig.g<-(t(M)%*%diag(p[(2:(M.n+1))])%*%M) +diag(p[1],n)

	Sig.CLS <- theo_cov<-Sig.g+Sig.exp
	
	err <- sum((theo_cov - cov1)^2/(1.001 - cov1))

	return(err)
}

M<-gg
M.n<-dim(gg)[1]

#set.seed(1234)
#init <- runif(M.n+1)
init <- pre1[(1:(M.n+1))]
fit <- optim(par = init, fn = NEGLOGLIK, control = list(trace = 5, maxit = 1500))

p <- fit$par



###################################################################################


NN <- 30
grid_x <- seq(from = 0, to = 1, length.out = NN)
grid_y <- seq(from = 0, to = 1, length.out = NN)
X <- expand.grid(grid_x, grid_y) %>% as.matrix()

nn <- nrow(X)
## calculating bending energy matrix B;
# See the texts by Bookstein (1991) or Dryden and Mardia (1998)

h <- as.matrix(dist(X, diag=TRUE,upper=TRUE))
K <- h^2 * log(h)
diag(K) <- 0
one <- rep(1, nn)
Gamma <- cbind(K,one,X)
Gamma <- rbind(Gamma,c(one,0,0,0))
Gamma <- rbind(Gamma,cbind(t(X),matrix(0,2,3)))
Ginv <- solve(Gamma)
Ginv <- (Ginv + t(Ginv))/2  # make exactly symmetric prior to eigen

B <- Ginv[1:nn, 1:nn]
Beig <- eigen(B)
g <- Beig$vectors
l <- Beig$values
g <- g[,order(l)]
l <- l[order(l)]

NEGLOGLIK <- function(p){

	#w <- p[length(jWarp) + 1:2]
	Xtarg <- Xtarg_orig <- sim_grid_locations 

	#if (TT > 1){
	#	for (tt in 1:(TT - 1)){
	#		temp_locs <- cbind(Xtarg_orig[, 1] - tt * w[1], Xtarg_orig[, 2] - tt * w[2])
	#		Xtarg <- rbind(Xtarg, temp_locs)
	#	}
	#}

	htarg <- as.matrix(dist(rbind(X, Xtarg),diag=TRUE,upper=TRUE))
	htarg <- htarg[1:nn, (nn + 1):(nn + nrow(Xtarg))]
	sigma <- htarg^2 * log(htarg)
	sigma[htarg == 0] <- 0

	beta1 <- p[1:length(jWarp)]
  
	parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T))

	spline1 <- c(parWarpsSum %*% sigma)
	nonparametric_cov_est <- outer(spline1, spline1, FUN = "*") #+ Sig.orig
	
	err <- sum((nonparametric_cov_est[reference_locations, ] - cov1[reference_locations, ])^2)
	return(err)
}

jWarp = 1:3
init2 <- rep(0, length(jWarp))
#set.seed(1234)
#init2 <- runif(length(jWarp), -1, 1)
fit2 <- optim(par = init2, fn = NEGLOGLIK, control = list(trace = 5, maxit = 2000))


jWarp = 1:4
fit2 <- optim(par = c(fit2$par, 0), fn = NEGLOGLIK, control = list(trace = 5, maxit = 2000))

for(tt in 1:1000){
	fit2 <- optim(par = fit2$par, fn = NEGLOGLIK, control = list(trace = 5, maxit = 500)) 
}

p <- fit2$par

Xtarg <- Xtarg_orig <- sim_grid_locations 

htarg <- as.matrix(dist(rbind(X, Xtarg),diag=TRUE,upper=TRUE))
htarg <- htarg[1:nn, (nn + 1):(nn + nrow(Xtarg))]
sigma <- htarg^2 * log(htarg)
sigma[htarg == 0] <- 0

beta1 <- p[1:length(jWarp)]

parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T))

spline1 <- c(parWarpsSum %*% sigma)
Sig.CLS <- nonparametric_cov_est <- outer(spline1, spline1, FUN = "*")# + Sig.orig

Sig.CLS[1:5, 1:5]

reference_locations <- c(1, ceiling(n / 2) + ceiling(N / 2))

pdf(file = paste(root, 'Figures/4-basis-function-estimation.pdf', sep = ''), width = 15, height = 15)

par(mfrow = c(2, 2))

par(mai=c(1.5, 1.5, 1.5, 1.5))
par(pty="s") 
image.plot(matrix(cov1[reference_locations[2], 1:n], N, N), zlim = c(0, max(Sig.CLS)))
mtext('Reference Location 1', side = 2, line = 3, cex = 1.5, font = 2, col="#000080")
mtext('Theoretical Covariance at t = 1', side = 3, line = 2, cex = 1.5, font = 2, col="#000080")

par(mai=c(1.5, 1.5, 1.5, 1.5))
par(pty="s") 
image.plot(matrix(Sig.CLS[reference_locations[2], 1:n], N, N), zlim = c(0, max(Sig.CLS)))
mtext('Estimated Covariance at t = 1', side = 3, line = 2, cex = 1.5, font = 2, col="#000080")

par(mai=c(1.5, 1.5, 1.5, 1.5))
par(pty="s") 
image.plot(matrix(cov1[reference_locations[1], 1:n], N, N), zlim = c(0, max(Sig.CLS)))
mtext('Reference Location 2', side = 2, line = 3, cex = 1.5, font = 2, col="#000080")

par(mai=c(1.5, 1.5, 1.5, 1.5))
par(pty="s") 
image.plot(matrix(Sig.CLS[reference_locations[1], 1:n], N, N), zlim = c(0, max(Sig.CLS)))

dev.off()


######################################################################################################

NEGLOGLIK <- function(p){

	#w <- p[length(jWarp) + 1:2]
	Xtarg <- Xtarg_orig <- sim_grid_locations 

	#if (TT > 1){
	#	for (tt in 1:(TT - 1)){
	#		temp_locs <- cbind(Xtarg_orig[, 1] - tt * w[1], Xtarg_orig[, 2] - tt * w[2])
	#		Xtarg <- rbind(Xtarg, temp_locs)
	#	}
	#}

	htarg <- as.matrix(dist(rbind(X, Xtarg),diag=TRUE,upper=TRUE))
	htarg <- htarg[1:nn, (nn + 1):(nn + nrow(Xtarg))]
	sigma <- htarg^2 * log(htarg)
	sigma[htarg == 0] <- 0

	beta1 <- p[1 + 1:length(jWarp)]
	beta2 <- p[1 + length(jWarp) + 1:length(jWarp)]
	beta3 <- p[1 + 2 * length(jWarp) + 1:length(jWarp)]
	beta4 <- p[1 + 3 * length(jWarp) + 1:length(jWarp)]
	beta5 <- p[1 + 4 * length(jWarp) + 1:length(jWarp)]
	beta6 <- p[1 + 5 * length(jWarp) + 1:length(jWarp)]
	beta7 <- p[1 + 6 * length(jWarp) + 1:length(jWarp)]
  
	parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T))
	spline1 <- c(parWarpsSum %*% sigma)
	parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T))
	spline2 <- c(parWarpsSum %*% sigma)
	parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta3, ncol=length(beta3), nrow=nrow(X), byrow=T))
	spline3 <- c(parWarpsSum %*% sigma)
	parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta4, ncol=length(beta4), nrow=nrow(X), byrow=T))
	spline4 <- c(parWarpsSum %*% sigma)
	parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta5, ncol=length(beta5), nrow=nrow(X), byrow=T))
	spline5 <- c(parWarpsSum %*% sigma)
	parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta6, ncol=length(beta6), nrow=nrow(X), byrow=T))
	spline6 <- c(parWarpsSum %*% sigma)
	parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta7, ncol=length(beta7), nrow=nrow(X), byrow=T))
	spline7 <- c(parWarpsSum %*% sigma)

	nonparametric_cov_est <- spline1 %o% spline1 + spline2 %o% spline2 + spline3 %o% spline3 + spline4 %o% spline4 + spline5 %o% spline5 + spline6 %o% spline6 + spline7 %o% spline7 + p[1] * diag(n)
	
	err <- sum((nonparametric_cov_est - cov1)^2)
	return(err)
}

jWarp = 1:5
init2 <- c(0.5, rep(0, 7 * length(jWarp)))
fit2 <- optim(par = init2, fn = NEGLOGLIK, control = list(trace = 5, maxit = 10000))
for(tt in 1:1000){
	fit2 <- optim(par = fit2$par, fn = NEGLOGLIK, control = list(trace = 5, maxit = 1000)) 
}


p <- fit2$par
Xtarg <- Xtarg_orig <- sim_grid_locations 

htarg <- as.matrix(dist(rbind(X, Xtarg),diag=TRUE,upper=TRUE))
htarg <- htarg[1:nn, (nn + 1):(nn + nrow(Xtarg))]
sigma <- htarg^2 * log(htarg)
sigma[htarg == 0] <- 0

beta2 <- p[1 + length(jWarp) + 1:length(jWarp)]
beta3 <- p[1 + 2 * length(jWarp) + 1:length(jWarp)]
beta4 <- p[1 + 3 * length(jWarp) + 1:length(jWarp)]
beta5 <- p[1 + 4 * length(jWarp) + 1:length(jWarp)]

parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T))
spline1 <- c(parWarpsSum %*% sigma)
parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T))
spline2 <- c(parWarpsSum %*% sigma)
parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta3, ncol=length(beta3), nrow=nrow(X), byrow=T))
spline3 <- c(parWarpsSum %*% sigma)
parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta4, ncol=length(beta4), nrow=nrow(X), byrow=T))
spline4 <- c(parWarpsSum %*% sigma)
parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta5, ncol=length(beta5), nrow=nrow(X), byrow=T))
spline5 <- c(parWarpsSum %*% sigma)

Sig.CLS <- nonparametric_cov_est <- spline1 %o% spline1 + spline2 %o% spline2 + spline3 %o% spline3 + spline4 %o% spline4 + spline5 %o% spline5 + p[1] * diag(n)

Sig.CLS[1:5, 1:5]

pdf(file = paste(root, 'Figures/4-basis-function-estimation.pdf', sep = ''), width = 15, height = 15)

par(mfrow = c(2, 2))

image.plot(matrix(cov1[reference_locations[2], 1:n], N, N))

image.plot(matrix(Sig.CLS[reference_locations[2], 1:n], N, N))

dev.off()

#############################################################################################################

NEGLOGLIK <- function(p){

	Xtarg <- Xtarg_orig <- sim_grid_locations 

	htarg <- as.matrix(dist(rbind(X, Xtarg),diag=TRUE,upper=TRUE))
	htarg <- htarg[1:nn, (nn + 1):(nn + nrow(Xtarg))]
	sigma <- htarg^2 * log(htarg)
	sigma[htarg == 0] <- 0

	beta1 <- p[1:length(jWarp)]
	beta2 <- p[length(jWarp) + 1:length(jWarp)]
	beta3 <- p[2 * length(jWarp) + 1:length(jWarp)]
	beta4 <- p[3 * length(jWarp) + 1:length(jWarp)]
	beta5 <- p[4 * length(jWarp) + 1:length(jWarp)]
	beta6 <- p[5 * length(jWarp) + 1:length(jWarp)]
	beta7 <- p[6 * length(jWarp) + 1:length(jWarp)]
  
	parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T))
	spline1 <- c(parWarpsSum %*% sigma)
	parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T))
	spline2 <- c(parWarpsSum %*% sigma)
	parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta3, ncol=length(beta3), nrow=nrow(X), byrow=T))
	spline3 <- c(parWarpsSum %*% sigma)
	parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta4, ncol=length(beta4), nrow=nrow(X), byrow=T))
	spline4 <- c(parWarpsSum %*% sigma)
	parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta5, ncol=length(beta5), nrow=nrow(X), byrow=T))
	spline5 <- c(parWarpsSum %*% sigma)
	parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta6, ncol=length(beta6), nrow=nrow(X), byrow=T))
	spline6 <- c(parWarpsSum %*% sigma)
	parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta7, ncol=length(beta7), nrow=nrow(X), byrow=T))
	spline7 <- c(parWarpsSum %*% sigma)

	nonparametric_cov_est <- spline1 %o% spline1 + spline2 %o% spline2 + spline3 %o% spline3 + spline4 %o% spline4 + spline5 %o% spline5 + spline6 %o% spline6 + spline7 %o% spline7 + Sig.orig
	
	err <- sum((nonparametric_cov_est - cov1)^2)
	return(err)
}

jWarp = 1:5
init2 <- rep(0, 7 * length(jWarp))
fit2 <- optim(par = init2, fn = NEGLOGLIK, control = list(trace = 5, maxit = 10000))



