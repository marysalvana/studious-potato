
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

########  CHANG HSU

source(file = paste(root, "R_codes/cls.r", sep = ''))

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

#sd.set<-c(1:4)/8
sd.set<- seq(0.08, 0.5, length.out = 5)
nu.set <- c(0.5, 1, 1.5)
gg<-NULL
loc.mu<-NULL
shift <- 0
for(nu_val in 1:length(nu.set)){
	shift <- shift + 0.2
	for(i in 1:length(sd.set)){
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
Sig.CLS <- (t(M)%*%diag(pre1[(2:(M.n+1))])%*%M) +diag(pre1[1],n)

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

###################################################################################

NEGLOGLIK <- function(p){
	
	loc2 <- sim_grid_locations
	for(tt in 1:(TT - 1)){
		loc2 <- rbind(loc2, cbind(sim_grid_locations[, 1] - p[1] * tt, sim_grid_locations[, 2] - p[2] * tt))
	}

	dist1<-as.matrix(dist(loc2))

	A <- 1:nrow(loc2)

	sd.set<- seq(0.08, 0.5, length.out = 5)
	nu.set <- c(0.5, 1, 1.5)
	gg<-NULL
	loc.mu<-NULL
	shift <- 0
	for(nu_val in 1:length(nu.set)){
		shift <- shift + 0.2
		for(i in 1:length(sd.set)){
			vari<-sd.set[i]^2
			loc.mui<-NULL
			mu.m<-sd.set[i]
			hex.nx<-floor((1+1.5*mu.m)/mu.m)
			hex.ny<-floor((1+mu.m*sqrt(3)/2)/(mu.m*sqrt(3)/2))
			anchor_locs <- expand.grid(-1:(hex.ny), -1:(hex.nx)) %>% as.matrix()
			loc.jk <- cbind(mu.m * (anchor_locs[, 2] + (anchor_locs[, 1] %% 2) / 2) + shift, mu.m * anchor_locs[, 1] * sqrt(3)/2 + shift)
			loc.mui <- cbind(rep(sd.set[i], nrow(loc.jk)), loc.jk)
			temp.g <- distR_C(loc.jk, loc2)
			temp.g<- Matern(temp.g, range = vari, nu = nu.set[nu_val])
			gg<-rbind(gg,temp.g)
			loc.mu<-rbind(loc.mu,loc.mui)
		}
	}

	M<-gg
	M.n<-dim(gg)[1]

	Sig.CLS <- (t(M)%*%diag(pre1[(2:(M.n+1))])%*%M) +diag(pre1[1], n * TT)

	err <- sum((Sig.CLS - cov1)^2)

	return(err)
}

init <- c(0, 0)
fit <- optim(par = init, fn = NEGLOGLIK, control = list(trace = 5, maxit = 5000))

p <- fit$par

loc2 <- sim_grid_locations
for(tt in 1:(TT - 1)){
	loc2 <- rbind(loc2, cbind(sim_grid_locations[, 1] - p[1] * tt, sim_grid_locations[, 2] - p[2] * tt))
}

dist1<-as.matrix(dist(loc2))

A <- 1:nrow(loc2)

sd.set<- seq(0.08, 0.5, length.out = 5)
nu.set <- c(0.5, 1, 1.5)
gg<-NULL
loc.mu<-NULL
shift <- 0
for(nu_val in 1:length(nu.set)){
	shift <- shift + 0.2
	for(i in 1:length(sd.set)){
		vari<-sd.set[i]^2
		loc.mui<-NULL
		mu.m<-sd.set[i]
		hex.nx<-floor((1+1.5*mu.m)/mu.m)
		hex.ny<-floor((1+mu.m*sqrt(3)/2)/(mu.m*sqrt(3)/2))
		anchor_locs <- expand.grid(-1:(hex.ny), -1:(hex.nx)) %>% as.matrix()
		loc.jk <- cbind(mu.m * (anchor_locs[, 2] + (anchor_locs[, 1] %% 2) / 2) + shift, mu.m * anchor_locs[, 1] * sqrt(3)/2 + shift)
		loc.mui <- cbind(rep(sd.set[i], nrow(loc.jk)), loc.jk)
		temp.g <- distR_C(loc.jk, loc2)
		temp.g<- Matern(temp.g, range = vari, nu = nu.set[nu_val])
		gg<-rbind(gg,temp.g)
		loc.mu<-rbind(loc.mu,loc.mui)
	}
}


AM<-gg[,A]
k<-dim(AM)[1]
a<-dim(AM)[2]
si.fix<-NULL
si.temp1<-matrix(as.numeric(outer(c(1:a),c(1:a),"==")),a,a)
si.temp2<-NULL
si.temp3 <- matrix(1:(a^2), a, a)
si.temp4<-NULL
for(s in 1:a){
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
Sig.CLS <-(t(M)%*%diag(pre1[(2:(M.n+1))])%*%M) +diag(pre1[1],n * TT)

fit <- optim(par = fit$par, fn = NEGLOGLIK, control = list(trace = 5, maxit = 5000))

p <- fit$par
TT <- 5
loc2 <- sim_grid_locations
for(tt in 1:(TT - 1)){
	loc2 <- rbind(loc2, cbind(sim_grid_locations[, 1] - p[1] * tt, sim_grid_locations[, 2] - p[2] * tt))
}

dist1<-as.matrix(dist(loc2))

A <- 1:nrow(loc2)

sd.set<- seq(0.08, 0.5, length.out = 5)
nu.set <- c(0.5, 1, 1.5)
gg<-NULL
loc.mu<-NULL
shift <- 0
for(nu_val in 1:length(nu.set)){
	shift <- shift + 0.2
	for(i in 1:length(sd.set)){
		vari<-sd.set[i]^2
		loc.mui<-NULL
		mu.m<-sd.set[i]
		hex.nx<-floor((1+1.5*mu.m)/mu.m)
		hex.ny<-floor((1+mu.m*sqrt(3)/2)/(mu.m*sqrt(3)/2))
		anchor_locs <- expand.grid(-1:(hex.ny), -1:(hex.nx)) %>% as.matrix()
		loc.jk <- cbind(mu.m * (anchor_locs[, 2] + (anchor_locs[, 1] %% 2) / 2) + shift, mu.m * anchor_locs[, 1] * sqrt(3)/2 + shift)
		loc.mui <- cbind(rep(sd.set[i], nrow(loc.jk)), loc.jk)
		temp.g <- distR_C(loc.jk, loc2)
		temp.g<- Matern(temp.g, range = vari, nu = nu.set[nu_val])
		gg<-rbind(gg,temp.g)
		loc.mu<-rbind(loc.mu,loc.mui)
	}
}

####est. Sigma by CLS
M<-gg
M.n<-dim(gg)[1]
Sig.CLS <-(t(M)%*%diag(pre1[(2:(M.n+1))])%*%M) +diag(pre1[1],n * TT)

#########################   TESTING   ########################

cov1 <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(1, 0.23, 1, 0.1001, 0.1001, 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT)
cov2 <- Sig.CLS

set.seed(1234)
r1 <- mvrnorm(500, rep(0, ncol(cov1)), cov1)
set.seed(1234)
r2 <- mvrnorm(500, rep(0, ncol(Sig.CLS)), Sig.CLS)
set.seed(5678)
r3 <- mvrnorm(500, rep(0, ncol(Sig.CLS)), Sig.CLS)

EMPCOV <- array(, dim = c(dim(cov1), 2))

EMPCOV[, , 1] <- cov(r1)
EMPCOV[, , 2] <- cov(r2)
EMPCOV[, , 3] <- cov(r3)

adj_mu <- 0.1001
adj_sig <- 0.001


locs_sub_index <- which(sim_grid_locations[, 1] >= 0.25 & sim_grid_locations[, 1] <= 0.75 & sim_grid_locations[, 2] >= 0.25 & sim_grid_locations[, 2] <= 0.75)
locs_sub_length <- length(locs_sub_index)

DIFF_ARRAY_EMP <- array(, dim = c(locs_sub_length, TT - 1, length(adj_mu), dim(EMPCOV)[3]))

for(model in 1:dim(EMPCOV)[3]){

	empcov <- EMPCOV[, , model]

	for(m in 1:length(adj_mu)){

		set.seed(1234)
		WIND_SIMULATED <- matrix(mvrnorm(n_sim, mu = wind_mu + adj_mu[m], Sigma = matrix(wind_sigma * adj_sig[m], ncol = 2, nrow = 2)), ncol = 2, byrow = T)

		count <- 1
		diff_cov_emp <- diff_cov_theo <- matrix(, ncol = 4, nrow = locs_sub_length)
		for(l2 in locs_sub_index){
			diff_cov_emp_temp <- diff_cov_theo_temp <- NULL
			for(t2 in 1:(TT - 1)){
				cov_purely_time_emp <- empcov[l2, l2 + n * t2]
				cov_purely_space_emp_temp <- cov_purely_space_theo_temp <- 0
				for(k in 1:n_sim){
					wind <- WIND_SIMULATED[k, ]
					new_loc <- matrix(c(sim_grid_locations[l2, 1] - wind[1] * t2, sim_grid_locations[l2, 2] - wind[2] * t2), ncol = 2)
					find_new_loc_index <- which.min(distR_C(sim_grid_locations, new_loc))[1]

					cov_purely_space_emp_temp <- cov_purely_space_emp_temp + empcov[l2, find_new_loc_index]
				}
				cov_purely_space_emp <- cov_purely_space_emp_temp / n_sim
				diff_cov_emp_temp <- c(diff_cov_emp_temp, cov_purely_time_emp - cov_purely_space_emp)
			}
			diff_cov_emp[count, ] <- diff_cov_emp_temp
			count <- count + 1
		}

		DIFF_ARRAY_EMP[, , m, model] <- diff_cov_emp
	}
}

dataX <- t(DIFF_ARRAY_EMP[, , 1, 1])
dataRef <- t(DIFF_ARRAY_EMP[, , 1, 2])
g <- t(DIFF_ARRAY_EMP[, , 1, 3])

rank.value=rankTest(dataX, g, dataRef)
rank.value

