workstation = F

if(workstation){
	directory <- '/home/salvanmo/Desktop/'
	root <- paste(directory, 'studious-potato/', sep = '')
	source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
	sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))
}else{
	directory <- '/ibex/scratch/salvanmo/'
	root <- paste(directory, 'studious-potato/', sep = '')
	source(file = paste(root, "R_codes/Functions/load_packages-IBEX.R", sep = ''))
	sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2-IBEX.cpp",sep=''))
}

source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))

start_time <- Sys.time()



distributed = T



AREA <- 'SAUDI'

DATA <- data_format(aggregate = 1, area = AREA)
locs <- as.matrix(DATA[[2]])
locs <- cbind((locs[, 1] - mean(locs[, 1])) / sd(locs[, 1]), (locs[, 2] - mean(locs[, 2])) / sd(locs[, 2]))



#######################################################################################


velocity_mu_config = 2
velocity_var_config = 2

mu_k <- c(0, 0.2001)
var_k <- c(0.0001, 0.01, 1)

WIND <- WIND_MU <- rep(mu_k[velocity_mu_config], 2)
WIND_VAR <- matrix(var_k[velocity_var_config] * diag(2), 2, 2)



N <- 50
n <- N^2
TT <- 5
grid_x <- seq(from = min(locs[, 1]), to = max(locs[, 1]), length.out = N)
grid_y <- seq(from = min(locs[, 2]), to = max(locs[, 2]), length.out = N)
sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

reference_locations <- c(2 * N + 5, ceiling(n / 2) + ceiling(N / 2))


NN <- 50
grid_x <- seq(from = -3, to = 3, length.out = NN)
grid_y <- seq(from = -3, to = 3, length.out = NN)
X <- expand.grid(grid_x, grid_y) %>% as.matrix()
nn <- nrow(X)

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

Xtarg <- sim_grid_locations 

if (TT > 1){
	for (tt in 1:(TT - 1)){
		temp_locs <- cbind(sim_grid_locations[, 1] - tt * WIND[1], sim_grid_locations[, 2] - tt * WIND[2])
		Xtarg <- rbind(Xtarg, temp_locs)
	}
}

htarg <- as.matrix(dist(rbind(X, Xtarg),diag=TRUE,upper=TRUE))
htarg <- htarg[1:nn, (nn + 1):(nn + nrow(Xtarg))]
sigma <- htarg^2 * log(htarg)
sigma[htarg == 0] <- 0


#p <- fit1$par
p <- c(0.0093092093, 0.0149825632, 0.0182713056, -0.0096880506, 0.0021445396, -0.0160542132, -0.0061874637, 0.0043052011, 0.0268706150, -0.0699589454, 0.0101119343, 0.0264518611, 0.0245503625, 0.0076773012, -0.0016969535, 0.0099870058, -0.0210626524, -0.0009397276, -0.0097587524, 0.0122824133, 0.0105428080, -0.0627256549, -0.0065139657, 0.0124965931, 0.0097002523, -0.0039516145, 0.0145425766, 0.0372800506, 0.0221969349, -0.0242916859, 0.0056526464, 0.0312454629, -0.0022965453)

jWarp = 1:10
theta <- exp(p[1:3])

beta1 <- p[3 + 1:length(jWarp)]
beta2 <- p[3 + length(jWarp) + 1:length(jWarp)]
beta3 <- p[3 + 2 * length(jWarp) + 1:length(jWarp)]

parWarpsSum <- cbind(rowSums( g[,3+jWarp] * matrix(beta3, ncol=length(beta3), nrow=nrow(X), byrow=T)),
	       rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T)),
		rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T)))

PARAMETER_NONSTAT <- t(sigma) %*% parWarpsSum

#ind <- 72 #40, 54 #3, 4, 6, 8, 9, 14, 16
#set.seed(ind)
#PARAMETER_NONSTAT <- runif(12, -1, 1)
#PARAMETER_NONSTAT <- c(0, PARAMETER_NONSTAT[1:4], 0, PARAMETER_NONSTAT[5:8], 0, PARAMETER_NONSTAT[9:12])
#PARAMETER_NONSTAT <- runif(15, -1, 1)



#######################################################################################



cat('Computing covariances...', '\n')

if(!distributed){

	cov1 <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(1, 0.23, 1, WIND, 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT, FITTING = T)

	cat('Generating realizations...', '\n')

	set.seed(1)
	r1 <- rmvn(100, rep(0, n * TT), cov1[["covariance"]], ncores = number_of_cores_to_use)

	cat('Saving the values...', '\n')

	write.table(cov1[["covariance"]][reference_locations, ], file = paste(root, "Data/univariate-nonstationary/cov-example-1-velocity_mu_config", velocity_mu_config, "_velocity_var_config_", velocity_var_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

}else{

	cores=detectCores()
	number_of_cores_to_use = 39
	#number_of_cores_to_use = cores[1]-1 # not to overload the computer

	cat('Registering', number_of_cores_to_use, 'cores...', '\n')

	cl <- makeCluster(number_of_cores_to_use) 
	registerDoParallel(cl)

	number_of_chunks = number_of_cores_to_use

	clusterExport(cl, "root")
	clusterEvalQ(cl, source(paste(root, "R_codes/Functions/cov_func.R", sep = '')))
	
	if(workstation){
		clusterEvalQ(cl, source(paste(root, "R_codes/Functions/load_packages.R", sep = '')))
		clusterEvalQ(cl, sourceCpp(paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp", sep = '')))
	}else{
		clusterEvalQ(cl, source(paste(root, "R_codes/Functions/load_packages-IBEX.R", sep = '')))
		clusterEvalQ(cl, sourceCpp(paste(root, "R_codes/Functions/spatially_varying_parameters2-IBEX.cpp", sep = '')))
	}

	clusterExport(cl, c("PARAMETER_NONSTAT", "TT"), envir = environment())

	cat('Simulating wind values...', '\n')

	set.seed(1234)
	wind_vals <- mvrnorm(100, WIND_MU, WIND_VAR)

	cat('Distributing computations over', number_of_cores_to_use, 'cores...', '\n')

	output <- foreach(i=1:nrow(wind_vals), .combine=cbind, .packages = "Rcpp", .noexport = "SPATIALLY_VARYING_PARAMETERS_FOR_FITTING_PARALLEL") %dopar% {
		
		COVARIANCE <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(1, 0.23, 1, wind_vals[i, ], 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT, FITTING = T, PARALLEL = T)

		return(c(COVARIANCE))
	}

	stopCluster(cl)



	cov1 <- matrix(rowSums(output), n * TT, n * TT) / nrow(wind_vals) 

	cat('Generating realizations...', '\n')

	set.seed(1)
	r1 <- rmvn(100, rep(0, n * TT), cov1, ncores = number_of_cores_to_use)

	cat('Saving the values...', '\n')

	write.table(cov1[reference_locations, ], file = paste(root, "Data/univariate-nonstationary/cov-example-1-velocity_mu_config", velocity_mu_config, "_velocity_var_config_", velocity_var_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

}

write.table(r1[1:10, ], file = paste(root, "Data/univariate-nonstationary/realizations-example-1-velocity_mu_config", velocity_mu_config, "_velocity_var_config_", velocity_var_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

end_time <- Sys.time()

cat('DONE.', '\n')
cat("Textfiles are saved in ", paste(root, 'Data/univariate-nonstationary/', sep = ''), '\n')
cat(end_time - start_time, 'seconds', '\n')



#######################################################################################

set.seed(3)
PARAMETER_DEFORMATION <- c(runif(2, -8, 8), runif(1, 0, 8), runif(2, -1, 1))

#cov2 <- MATERN_UNI_DEFORMATION(PARAMETER = c(1, 0.23, 1, 0.1001, 0.1001, 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_DEFORMATION = PARAMETER_DEFORMATION)

