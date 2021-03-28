
workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/' 		else 		directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))

sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))

N <- 10
n <- N^2
TT <- 5
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
grid_x <- seq(from = -3, to = 3, length.out = NN)
grid_y <- seq(from = -3, to = 3, length.out = NN)
X <- expand.grid(grid_x, grid_y) %>% as.matrix()

NEGLOGLIK <- function(p){
	
	w <- p[1:2]

	cat(w, '\n')

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

	beta1 <- matrix(p[2 + 1:length(jWarp)], ncol = ncol(sigma), nrow = length(jWarp))

	sub_sigma <- beta1 * sigma[1:length(jWarp), ]
	
	nonparametric_cov_est <- t(sub_sigma) %*% sub_sigma + diag(diag(emp_cov), ncol(sigma), ncol(sigma))

	err <- sum((emp_cov - nonparametric_cov_est)^2)

	return(err)
}

jWarp = 1:10
init <- c(0.1, 0.1, rep(0, length(jWarp)))
fit <- optim(par = init, fn = NEGLOGLIK, control = list(trace = 5, maxit = 10000))

jWarp = 1:15
init <- c(fit$par, rep(0, length(jWarp) - length(fit$par) + 2))
fit <- optim(par = init, fn = NEGLOGLIK, control = list(trace = 5, maxit = 1000))

