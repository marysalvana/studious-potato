
workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/' 		else 		directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))

sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))

N <- 10
n <- N^2
TT <- 10
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
r1 <- rmvn(100, rep(0, ncol(cov1)), cov1, ncores = 25)

emp_cov <- cov(r1)

NN <- 50
grid_x <- seq(from = -3, to = 3, length.out = NN)
grid_y <- seq(from = -3, to = 3, length.out = NN)
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

	w <- p[length(jWarp) + 1:2]
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

	beta1 <- p[1:length(jWarp)]
  
	parWarpsSum <- rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T))

	nonparametric_cov_est <- c(parWarpsSum %*% sigma)

	err <- sum((nonstat_est1 - nonstat_est1_targ)^2 + (nonstat_est2 - nonstat_est2_targ)^2 + (nonstat_est3 - nonstat_est3_targ)^2)

	return(err)
}

jWarp = 1:50
init2 <- rep(0, 3 * length(jWarp))
fit2 <- optim(par = init2, fn = NEGLOGLIK2_SCRATCH, control = list(trace = 5, maxit = 1000))
for(tt in 1:1000){
	fit2 <- optim(par = fit2$par, fn = NEGLOGLIK2_SCRATCH, control = list(trace = 5, maxit = 500)) 
	plot_spatially_varying_scratch(fit2$par)
}
