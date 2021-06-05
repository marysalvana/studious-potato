workstation = T

if(workstation){
	directory <- '/home/salvanmo/Desktop/'
	root <- paste(directory, 'studious-potato/', sep = '')
	source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
	sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))
	number_of_cores_to_use = 24
}else{
	directory <- '/ibex/scratch/salvanmo/'
	root <- paste(directory, 'studious-potato/', sep = '')
	source(file = paste(root, "R_codes/Functions/load_packages-IBEX.R", sep = ''))
	sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2-IBEX.cpp",sep=''))
	number_of_cores_to_use = 39
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

model = 2
velocity_mu_config = 2
velocity_var_config = 1

rho_config = 3


mu_k <- c(0, 0.3001)
var_k <- c(0.001, 0.1, 1)

WIND <- WIND_MU <- rep(mu_k[velocity_mu_config], 2)
WIND_VAR <- matrix(var_k[velocity_var_config] * diag(2), 2, 2)

rho_k <- c(-0.5, 0, 0.5)
VARIABLE_RHO <- rho_k[rho_config]


N <- 20
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


if(model == 1){

	#p <- fit1$par
	#p <- c(0.0093092093, 0.0149825632, 0.0182713056, -0.0096880506, 0.0021445396, -0.0160542132, -0.0061874637, 0.0043052011, 0.0268706150, -0.0699589454, 0.0101119343, 0.0264518611, 0.0245503625, 0.0076773012, -0.0016969535, 0.0099870058, -0.0210626524, -0.0009397276, -0.0097587524, 0.0122824133, 0.0105428080, -0.0627256549, -0.0065139657, 0.0124965931, 0.0097002523, -0.0039516145, 0.0145425766, 0.0372800506, 0.0221969349, -0.0242916859, 0.0056526464, 0.0312454629, -0.0022965453)
	
	p <- c(0.0071644151, 0.0094228436, -0.0094357164, 0.0011911244, -0.0001946439, 0.0044519415, -0.0029900530, 0.0024820814, -0.0077504082, 0.0037965748, -0.0691636929, 0.0016734147, 0.0056327376, 0.0014253809, -0.0057678577, -0.0030793794, 0.0087365585, 0.0058047295, 0.0040945055, -0.0041464973, 0.0366414100, 0.0088822294, 0.0167993236, 0.0064401123, 0.0090248538, 0.0084194860, 0.0063182756, 0.0122171447, 0.0069259888, 0.0068051953, -0.0066732051, 0.0201189423, 0.0066393622)

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

		write.table(cov1[["covariance"]][reference_locations, ], file = paste(root, "Data/univariate-nonstationary/cov-example-1-velocity_mu_config_", velocity_mu_config, "_velocity_var_config_", velocity_var_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

	}else{

		cores=detectCores()

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

		write.table(cov1[reference_locations, ], file = paste(root, "Data/univariate-nonstationary/cov-example-1-velocity_mu_config_", velocity_mu_config, "_velocity_var_config_", velocity_var_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

	}

	write.table(r1[1:10, ], file = paste(root, "Data/univariate-nonstationary/realizations-example-1-velocity_mu_config_", velocity_mu_config, "_velocity_var_config_", velocity_var_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

}else if(model == 2){

	#set.seed(3)
	#PARAMETER_DEFORMATION <- c(runif(2, -8, 8), runif(1, 0, 8), runif(2, -1, 1))

	#p <- fit1$par

	#p <- c(0.0093092093, 0.0149825632, 0.0182713056, -0.0096880506, 0.0021445396, -0.0160542132, -0.0061874637, 0.0043052011, 0.0268706150, -0.0699589454, 0.0101119343, 0.0264518611, 0.0245503625, 0.0076773012, -0.0016969535, 0.0099870058, -0.0210626524, -0.0009397276, -0.0097587524, 0.0122824133, 0.0105428080, -0.0627256549, -0.0065139657, 0.0124965931, 0.0097002523, -0.0039516145, 0.0145425766, 0.0372800506,  0.0221969349, -0.0242916859, 0.0056526464, 0.0312454629, -0.0022965453)

	p <- c(0.067172812, -0.141482629, -0.123072246, -0.020615311, 0.026682032, -0.064454540, 0.006658298, -0.007649700, 0.047388007, -0.004046326, 0.035169697, -0.032544016, -0.009401134, 0.036043187, 0.004287410, -0.010179291, 0.016116570, 0.022926999, -0.105069368, 0.117549831, -0.017031510, 0.094188039, 0.057326933)

        jWarp = 1:10
	theta <- exp(p[1:3])

	beta1 <- p[3 + 1:length(jWarp)]
	beta2 <- p[3 + length(jWarp) + 1:length(jWarp)]

	parWarpsSum <- cbind(rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T)),
		       rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T)))

	PARAMETER_DEFORMATION <- t(sigma) %*% parWarpsSum



	#######################################################################################



	cat('Computing covariances...', '\n')


	if(!distributed){
		cov2 <- MATERN_UNI_DEFORMATION(PARAMETER = c(1, 0.23, 1, WIND, 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_DEFORMATION = PARAMETER_DEFORMATION, FITTING = T)

		cat('Generating realizations...', '\n')

		set.seed(1)
		r2 <- rmvn(100, rep(0, n * TT), cov2[["covariance"]], ncores = number_of_cores_to_use)

		cat('Saving the values...', '\n')

		write.table(cov2[["covariance"]][reference_locations, ], file = paste(root, "Data/univariate-nonstationary/cov-example-2-velocity_mu_config_", velocity_mu_config, "_velocity_var_config_", velocity_var_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

	}else{

		cores=detectCores()

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

		clusterExport(cl, c("PARAMETER_DEFORMATION", "TT"), envir = environment())

		cat('Simulating wind values...', '\n')

		set.seed(1234)
		wind_vals <- mvrnorm(100, WIND_MU, WIND_VAR)

		cat('Distributing computations over', number_of_cores_to_use, 'cores...', '\n')

		output <- foreach(i=1:nrow(wind_vals), .combine=cbind, .packages = "Rcpp", .noexport = "DEFORMATION_FOR_FITTING_PARALLEL") %dopar% {
			
			COVARIANCE <- MATERN_UNI_DEFORMATION(PARAMETER = c(1, 0.23, 1, wind_vals[i, ], 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_DEFORMATION = PARAMETER_DEFORMATION, FITTING = T, PARALLEL = T)

			return(c(COVARIANCE))
		}

		stopCluster(cl)



		cov2 <- matrix(rowSums(output), n * TT, n * TT) / nrow(wind_vals) 

		print(range(cov2))
		cat('Generating realizations...', '\n')

		set.seed(1)
		r2 <- rmvn(100, rep(0, n * TT), cov2, ncores = number_of_cores_to_use)

		cat('Saving the values...', '\n')

		write.table(cov2[reference_locations, ], file = paste(root, "Data/univariate-nonstationary/cov-example-2-velocity_mu_config_", velocity_mu_config, "_velocity_var_config_", velocity_var_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)
	}

	write.table(r2[1:10, ], file = paste(root, "Data/univariate-nonstationary/realizations-example-2-velocity_mu_config_", velocity_mu_config, "_velocity_var_config_", velocity_var_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

}else if(model == 3){


	cat('Computing covariances...', '\n')

	PARAMETER_NONSTAT2 <- matrix(0, ncol = ncol(PARAMETER_NONSTAT), nrow = nrow(PARAMETER_NONSTAT))


        if(!distributed){

		#cov3_uni11 <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(1, 0.23, 0.5, WIND, 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT, FITTING = T)
		#cov3_uni22 <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(1, 0.23, 1, WIND, 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT, FITTING = T)
		#cov3_uni12 <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(VARIABLE_RHO, 0.23, 0.5 * (0.5 + 1), WIND, 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT, FITTING = T)

		#set.seed(1)
		#cov3_theo <- rbind(cbind(cov3_uni11[['covariance']], cov3_uni12[['covariance']]), cbind(t(cov3_uni12[['covariance']]), cov3_uni22[['covariance']]))
		#r3 <- rmvn(100, rep(0, n * TT * 2), cov3_theo, ncores = number_of_cores_to_use)

		cov3 <- MULTIVARIATE_MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(1, 1, 0.23, 0.5, 1, VARIABLE_RHO, WIND, 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT, PARAMETER_NONSTAT2 = PARAMETER_NONSTAT2, FITTING = T, PARALLEL = T)

		set.seed(1)
		r3 <- rmvn(100, rep(0, n * TT * 2), cov3, ncores = number_of_cores_to_use)
	}else{
		cores=detectCores()

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

		output22 <- foreach(i=1:nrow(wind_vals), .combine=cbind, .packages = "Rcpp", .noexport = "SPATIALLY_VARYING_PARAMETERS_FOR_FITTING_PARALLEL") %dopar% {
			
			COVARIANCE <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(1, 0.23, 0.5, wind_vals[i, ], 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT, FITTING = T, PARALLEL = T)

			return(c(COVARIANCE))
		}

		output12 <- foreach(i=1:nrow(wind_vals), .combine=cbind, .packages = "Rcpp", .noexport = "SPATIALLY_VARYING_PARAMETERS_FOR_FITTING_PARALLEL") %dopar% {
			
			COVARIANCE <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(VARIABLE_RHO, 0.23, 0.5 * (0.5 + 1), wind_vals[i, ], 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT, FITTING = T, PARALLEL = T)

			return(c(COVARIANCE))
		}

		stopCluster(cl)



		cov11 <- matrix(rowSums(output), n * TT, n * TT) / nrow(wind_vals) 
		cov22 <- matrix(rowSums(output22), n * TT, n * TT) / nrow(wind_vals) 
		cov12 <- matrix(rowSums(output12), n * TT, n * TT) / nrow(wind_vals) 

		cov3 <- rbind(cbind(cov11, cov12), cbind(t(cov12), cov22))

		cat('Generating realizations...', '\n')

		set.seed(1)
		r3 <- rmvn(100, rep(0, n * TT * 2), cov3, ncores = number_of_cores_to_use)

		cat('Saving the values...', '\n')
	
		write.table(cov3[reference_locations, ], file = paste(root, "Data/univariate-nonstationary/cov-example-3-velocity_mu_config_", velocity_mu_config, "_velocity_var_config_", velocity_var_config, "_rho_config_", rho_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)
	}
	write.table(r3[1:10, ], file = paste(root, "Data/univariate-nonstationary/realizations-example-3-velocity_mu_config_", velocity_mu_config, "_velocity_var_config_", velocity_var_config, "_rho_config_", rho_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

}else if(model == 4){

}

end_time <- Sys.time()

cat('DONE.', '\n')
cat("Textfiles are saved in ", paste(root, 'Data/univariate-nonstationary/', sep = ''), '\n')
cat(end_time - start_time, 'seconds', '\n')

