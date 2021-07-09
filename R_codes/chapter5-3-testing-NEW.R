


source("./pkg-config.R")



DISTRIBUTED = T
MODEL = 1



WIND <- WIND_MU <- rep(0.3001, 2)
WIND_VAR <- matrix(0.1 * diag(2), 2, 2)



N <- 20
n <- N^2
TT <- 5
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
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




if(MODEL == 1){



	p <- c(-0.0445389018, 0.0133324956, 0.1284334294, 0.0122780156, -0.0023861828, -0.0032722984, -0.0060019839, 0.0086527895, -0.0111882022, -0.0349673412, 0.0133417022, -0.0887673908, -0.0275147364, -0.0146839115, -0.0024901228, 0.0041304654, 0.0018204541, -0.0022792563, -0.0189851402, -0.0039365573, 0.0036359306, -0.0121414575, 0.0749431348, -0.0044259752, 0.0033292430, 0.0084257633, -0.0010380224, 0.0100569830, -0.0001364005, 0.0464309935, 0.0040019716, 0.0194226880, 0.0189757404)

	jWarp = 1:10
	theta <- exp(p[1:3])

	beta1 <- p[3 + 1:length(jWarp)]
	beta2 <- p[3 + length(jWarp) + 1:length(jWarp)]
	beta3 <- p[3 + 2 * length(jWarp) + 1:length(jWarp)]

	parWarpsSum <- cbind(rowSums( g[,3+jWarp] * matrix(beta3, ncol=length(beta3), nrow=nrow(X), byrow=T)),
		       rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T)),
			rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T)))

	PARAMETER_NONSTAT <- t(sigma) %*% parWarpsSum



	cat('Computing covariances...', '\n')

	if(!DISTRIBUTED){

		cov1 <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(1, 0.23, 1, WIND, 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT, FITTING = T)

		cat('Generating realizations...', '\n')

		set.seed(1)
		r1 <- rmvn(10, rep(0, n * TT), cov1[["covariance"]], ncores = number_of_cores_to_use)

		cat('Saving the values...', '\n')

		write.table(cov1[["covariance"]][reference_locations, ], file = paste(root, "Data/nonstationary-taylor-hypothesis/cov-example-1-velocity_mu_config_", velocity_mu_config, "_velocity_var_config_", velocity_var_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

	}else{

		cores=detectCores()

		number_of_cores_to_use = cores[1]-1 # not to overload the computer
		cat('Registering', number_of_cores_to_use, 'cores...', '\n')

		cl <- makeCluster(number_of_cores_to_use) 
		registerDoParallel(cl)

		clusterEvalQ(cl, source("./pkg-config.R"))
		clusterExport(cl, c("PARAMETER_NONSTAT", "TT"), envir = environment())

		cat('Simulating wind values...', '\n')

		set.seed(1234)
		wind_vals <- mvrnorm(100, WIND_MU, WIND_VAR)

		cat('Distributing computations over', number_of_cores_to_use, 'cores...', '\n')

		output <- foreach(i=1:nrow(wind_vals), .combine='+', .packages = "Rcpp", .noexport = "SPATIALLY_VARYING_PARAMETERS_FOR_FITTING_PARALLEL") %dopar% {
			
			COVARIANCE <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(1, 0.23, 1, wind_vals[i, ], 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT, FITTING = T, PARALLEL = T)

			return(c(COVARIANCE))
		}

		stopCluster(cl)



		cov1 <- matrix(output, n * TT, n * TT) / nrow(wind_vals) 

		cat('Generating realizations...', '\n')

		set.seed(1)
		r1 <- rmvn(300, rep(0, n * TT), cov1, ncores = number_of_cores_to_use)

		cat('Computing empirical covariance...', '\n')

		empcov <- cov(r1)

	}


}else if(MODEL == 2){

	p <- c(0.067172812, -0.141482629, -0.123072246, -0.020615311, 0.026682032, -0.064454540, 0.006658298, -0.007649700, 0.047388007, -0.004046326, 0.035169697, -0.032544016, -0.009401134, 0.036043187, 0.004287410, -0.010179291, 0.016116570, 0.022926999, -0.105069368, 0.117549831, -0.017031510, 0.094188039, 0.057326933)

	jWarp = 1:10
	theta <- exp(p[1:3])

	beta1 <- p[3 + 1:length(jWarp)]
	beta2 <- p[3 + length(jWarp) + 1:length(jWarp)]

	parWarpsSum <- cbind(rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T)),
		       rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T)))

	PARAMETER_DEFORMATION <- t(sigma) %*% parWarpsSum



	cat('Computing covariances...', '\n')


	if(!DISTRIBUTED){
		cov2 <- MATERN_UNI_DEFORMATION(PARAMETER = c(1, 0.23, 1, WIND, 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_DEFORMATION = PARAMETER_DEFORMATION, FITTING = T)

		cat('Generating realizations...', '\n')

		set.seed(1)
		r2 <- rmvn(10, rep(0, n * TT), cov2[["covariance"]], ncores = number_of_cores_to_use)

		cat('Saving the values...', '\n')

		write.table(cov2[["covariance"]][reference_locations, ], file = paste(root, "Data/nonstationary-taylor-hypothesis/cov-example-2-velocity_mu_config_", velocity_mu_config, "_velocity_var_config_", velocity_var_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

	}else{

		cores=detectCores()

		#number_of_cores_to_use = cores[1]-1 # not to overload the computer
		cat('Registering', number_of_cores_to_use, 'cores...', '\n')

		cl <- makeCluster(number_of_cores_to_use) 
		registerDoParallel(cl)

		clusterEvalQ(cl, source("./pkg-config.R"))
		clusterExport(cl, c("PARAMETER_DEFORMATION", "TT"), envir = environment())

		cat('Simulating wind values...', '\n')

		set.seed(1234)
		wind_vals <- mvrnorm(100, WIND_MU, WIND_VAR)

		cat('Distributing computations over', number_of_cores_to_use, 'cores...', '\n')

		output <- foreach(i=1:nrow(wind_vals), .combine='+', .packages = "Rcpp", .noexport = "DEFORMATION_FOR_FITTING_PARALLEL") %dopar% {
			
			COVARIANCE <- MATERN_UNI_DEFORMATION(PARAMETER = c(1, 0.23, 1, wind_vals[i, ], 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_DEFORMATION = PARAMETER_DEFORMATION, FITTING = T, PARALLEL = T)

			return(c(COVARIANCE))
		}

		stopCluster(cl)



		cov2 <- matrix(output, n * TT, n * TT) / nrow(wind_vals) 

		cat('Generating realizations...', '\n')

		set.seed(1)
		r2 <- rmvn(10, rep(0, n * TT), cov2, ncores = number_of_cores_to_use)

		cat('Saving the values...', '\n')

	}


}

