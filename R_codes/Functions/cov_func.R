MATERN_UNI_STATIONARY <- function(PARAMETER, LOCATION){

        nu <- PARAMETER[3]
        beta <- PARAMETER[2]
        sigma2 <- PARAMETER[1]

	dist0 <- dist(x = LOCATION, diag = TRUE, upper = TRUE) %>% as.matrix()

	S <- ifelse(dist0 != 0, (dist0 / beta)^nu * besselK(dist0 / beta, nu) / (2^(nu - 1) * gamma(nu)) * sigma2, sigma2)

  	return(S)
}

MATERN_UNI_SPATIALLY_VARYING_PARAMETERS <- function(PARAMETER, PARAMETER_NONSTAT, LOCATION, TIME, N_SIM = NULL, FITTING = F, PARALLEL = F) {

	n <- nrow(LOCATION)

	LOCATION_NEW <- NULL

        for(tt in 0:(TIME - 1)){
                LOCATION_NEW <- rbind(LOCATION_NEW, cbind(LOCATION, rep(tt, n)))
        }

	if(!is.null(N_SIM)){
		set.seed(1234)
		WIND_SIMULATED <- matrix(mvrnorm(n_sim, mu = PARAMETER[4:5], Sigma = matrix(PARAMETER[6:9], ncol = 2, nrow = 2)), ncol = 2, byrow = T)
	}else{
		WIND_SIMULATED <- matrix(PARAMETER[4:5], ncol = 2, byrow = T)	
	}

	if(!FITTING & !PARALLEL){
		SIGMA <- SPATIALLY_VARYING_PARAMETERS(Loc = LOCATION_NEW, param = PARAMETER[1:3], wind = WIND_SIMULATED, param_nonstat = PARAMETER_NONSTAT, time = TIME)

		FINAL_DATA <- list("covariance" = SIGMA[[1]], "parameters" = SIGMA[[2]])
	}else if(!FITTING & PARALLEL){
                SIGMA <- SPATIALLY_VARYING_PARAMETERS_PARALLEL(Loc = LOCATION_NEW, param = PARAMETER[1:3], wind = WIND_SIMULATED, param_nonstat = PARAMETER_NONSTAT, time = TIME)

                FINAL_DATA <- SIGMA
        }else if(FITTING & !PARALLEL){
		SIGMA <- SPATIALLY_VARYING_PARAMETERS_FOR_FITTING(Loc = LOCATION_NEW, param = PARAMETER[1:3], wind = WIND_SIMULATED, param_nonstat = PARAMETER_NONSTAT, time = TIME)

		FINAL_DATA <- list("covariance" = SIGMA[[1]])
	}else if(FITTING & PARALLEL){
                SIGMA <- SPATIALLY_VARYING_PARAMETERS_FOR_FITTING_PARALLEL(Loc = LOCATION_NEW, param = PARAMETER[1:3], wind = WIND_SIMULATED, param_nonstat = PARAMETER_NONSTAT, time = TIME)

                FINAL_DATA <- SIGMA
        }
        return(FINAL_DATA) 
}

MATERN_UNI_DEFORMATION <- function(PARAMETER, PARAMETER_DEFORMATION, LOCATION, TIME, N_SIM = NULL) {

	n <- nrow(LOCATION)

	LOCATION_NEW <- NULL

        for(tt in 0:(TIME - 1)){
                LOCATION_NEW <- rbind(LOCATION_NEW, cbind(LOCATION, rep(tt, n)))
        }

	if(!is.null(N_SIM)){
		set.seed(1234)
		WIND_SIMULATED <- matrix(mvrnorm(n_sim, mu = PARAMETER[4:5], Sigma = matrix(PARAMETER[6:9], ncol = 2, nrow = 2)), ncol = 2, byrow = T)
	}else{
		WIND_SIMULATED <- matrix(PARAMETER[4:5], ncol = 2, byrow = T)	
	}

	SIGMA <- DEFORMATION(Loc = LOCATION_NEW, param = PARAMETER[1:3], wind = WIND_SIMULATED, param_nonstat = PARAMETER_DEFORMATION)
 	
        return(SIGMA) 
}

MATERN_UNI_POINT_SOURCE_DEFORMATION <- function(PARAMETER, PARAMETER_DEFORMATION, LOCATION1, LOCATION2) {

	SIGMA <- POINT_SOURCE_DEFORMATION(Loc1 = LOCATION1, Loc2 = LOCATION2, param = PARAMETER[1:3], param_nonstat = PARAMETER_DEFORMATION)
 	
        return(SIGMA) 
}


#---------STATIONARY------------#

#-------UNIVARIATE----------#

schlather_stationary <- function(theta, wind_mu, wind_var, LOCS){

	w <- wind_mu
        d <- length(w) 
        Sigma <- wind_var

	n <- nrow(LOCS)

        SS <- list()

        dist0 <- parDist(x = LOCS, method = "euclidean") %>% as.matrix()

	tempcov <- exp(-dist0^2 / theta[2])

	SS[[1]] <- tempcov 

	for(tt in 1:(TT - 1)){
		denom <- (diag(d) + Sigma * tt^2 )
		denom_inv <- solve(denom)
		tempcov <- apply(LOCS, 1, quadratic_func_FINAL, w = w, SIGS_INV = denom_inv, tt, scale = theta[2], LOCS = LOCS) / sqrt(det(denom)) 

		SS[[tt + 1]] <- tempcov
	}
        S <- theta[1] * toeplitz_mat(SS)
	return(S)
}

gaussian_spatially_varying_covariates_uni_cov_ST_for_fitting <- function(theta, LOCS, TIME_PARAM, NONSTATPARAMS){

	beta <- theta[2]
        sigma2 <- theta[1]

        locs <- coords <- LOCS

  	S <- spatially_varying_parameters3_for_fitting(Loc = locs, param = c(sigma2, beta, TIME_PARAM), NONSTATPARAMS)

  	return(S)
}

gaussian_st_uni_cov_NS_for_fitting <- function(theta, max_time_lag = 0, LOCS, NONSTATPARAMS){

        sigma_t <- theta[3]
        sigma2 <- theta[1]

	SS <- list()
	S <- SS[[1]] <- gaussian_spatially_varying_covariates_uni_cov_ST_for_fitting(theta, LOCS, TIME_PARAM = 1, NONSTATPARAMS)
	
	if(max_time_lag > 0){

		for(tt in 1:max_time_lag){
			
			SS[[tt + 1]] <- gaussian_spatially_varying_covariates_uni_cov_ST_for_fitting(theta, LOCS, TIME_PARAM = (1 + sigma_t * tt^2), NONSTATPARAMS) / (1 + sigma_t * tt^2) 
		}

		S <- sigma2 * toeplitz_mat(SS)
	}

  	return(S)
}

gaussian_spatially_varying_covariates_uni_cov_ST <- function(theta, LOCS, TIME_PARAM){

	beta <- theta[2]
        sigma2 <- theta[1]

        locs <- coords <- LOCS

  	S <- spatially_varying_parameters3(Loc = locs, param = c(sigma2, beta, TIME_PARAM))

  	return(S)
}

gaussian_st_uni_cov_NS <- function(theta, max_time_lag = 0, LOCS){

        sigma_t <- theta[3]
        sigma2 <- theta[1]

	SS <- list()
	SS[[1]] <- gaussian_spatially_varying_covariates_uni_cov_ST(theta, LOCS, TIME_PARAM = 1)

	for(tt in 1:max_time_lag){
		
		SS[[tt + 1]] <- gaussian_spatially_varying_covariates_uni_cov_ST(theta, LOCS, TIME_PARAM = (1 + sigma_t * tt^2)) / (1 + sigma_t * tt^2) 
	}

        S <- sigma2 * toeplitz_mat(SS)

  	return(S)
}

matern_spatially_varying_covariates_uni_cov_ST <- function(theta, LOCS, TIME_PARAM){

	nu <- theta[2]
        sigma2 <- theta[1]

        locs <- coords <- LOCS

  	S <- spatially_varying_parameters3(Loc = locs, param = c(sigma2, nu, TIME_PARAM))

  	return(S)
}

matern_st_uni_cov_NS <- function(theta, max_time_lag = 0, LOCS){

        sigma_t <- theta[3]
        sigma2 <- theta[1]

	SS <- list()
	#SS[[1]] <- Matern(dist1, range = 1, nu = nu)
	SS[[1]] <- matern_spatially_varying_covariates_uni_cov_ST(theta, LOCS, TIME_PARAM = 1)

	for(tt in 1:max_time_lag){
		
		#SS[[tt + 1]] <-  Matern(dist2, range = 1, nu = nu) / (1 + sigma_t * tt^2)
		SS[[tt + 1]] <- matern_spatially_varying_covariates_uni_cov_ST(theta, LOCS, TIME_PARAM = (1 + sigma_t * tt^2)) / (1 + sigma_t * tt^2) 
	}
        S <- sigma2 * toeplitz_mat(SS)

  	return(S)
}

gaussian_uni_cov <- function(theta, LOCS){

	nug <- theta[3]
        beta <- theta[2]
        sigma2 <- theta[1]

        dist0 <- parDist(x = LOCS, method = "euclidean") %>% as.matrix()

	S <- ifelse(dist0 != 0, sigma2 * exp( -dist0^2 / beta), sigma2 + nug)

  	return(S)
}

quadratic_func <- function(H_MAT, MU, SIGMA, U, scale){
	
	d <- length(MU)	
	DENOM <- diag(d) + SIGMA * U^2
	DENOM_INV <- solve(DENOM)

	H <- matrix(H_MAT - MU * U, nrow = 1)
	new_dist <- H %*% DENOM_INV %*% t(H)
	cov_val <- ifelse(H_MAT[1] == 0 & H_MAT[2] == 0 & U == 0, 1, exp( -new_dist / scale) / sqrt(det(DENOM)))
	return(cov_val)

}
quadratic_func_FINAL <- function(loc1, w, SIGS_INV, tt, scale, LOCS){
	
	H <- cbind(loc1[1] - LOCS[, 1] - w[1] * tt, loc1[2] - LOCS[, 2] - w[2] * tt)

	test <- diag(H %*% SIGS_INV %*% t(H))	
	cov_val <- exp( - test / scale) 

	return(cov_val)

}
quadratic_func_slow <- function(H_MAT, MU, DENOM_INV, U, scale){
	
	H <- matrix(H_MAT - MU * U, nrow = 1)
	new_dist <- H %*% DENOM_INV %*% t(H)
	cov_val <- exp( -new_dist / scale) 
	return(cov_val)

}

frozen_matern_uni_cov <- function(theta, wind, max_time_lag = 0, LOCS){

        nu <- theta[3]
        beta <- theta[2]
        sigma2 <- theta[1]

	w <- wind

        locs <- coords <- LOCS

        if (max_time_lag > 0){
                for (tt in 1:max_time_lag){
                        temp_locs <- cbind(coords[, 1] - tt * w[1], coords[, 2] - tt * w[2])
                        locs <- rbind(locs, temp_locs)
                }
        }
   	
        dist0 <- parDist(x = locs, method = "euclidean") %>% as.matrix()

	S <- ifelse(dist0 != 0, (dist0 / beta)^nu * besselK(dist0 / beta, nu) / (2^(nu - 1) * gamma(nu)), sigma2)

  	return(S)
}


gneiting_st_uni_cov <- function(theta, max_time_lag = 0, LOCS){

        sigma_t <- theta[3]
        #nu <- theta[3]
        beta <- theta[2]
        sigma2 <- theta[1]

        locs <- coords <- LOCS

        #dist0 <- rdist.earth(locs)
        dist0 <- parDist(x = locs, method = "euclidean") %>% as.matrix()
	dist1 <- (dist0^2 / beta)
	SS <- list()
	#SS[[1]] <- Matern(dist1, range = 1, nu = nu)
	SS[[1]] <- exp(-dist1)

	for(tt in 1:max_time_lag){
		
		dist2 <- dist1 / (1 + sigma_t * tt^2)
		#SS[[tt + 1]] <-  Matern(dist2, range = 1, nu = nu) / (1 + sigma_t * tt^2)
		SS[[tt + 1]] <-  exp(-dist2) / (1 + sigma_t * tt^2)
	}
        S <- sigma2 * toeplitz_mat(SS)

  	return(S)
}

nonfrozen_schlather_uni_cov_ORIG <- function(theta, wind, wind_var, max_time_lag, LOCS, h_mat){

	w <- wind
        d <- length(w) 
        Sigma <- matrix(wind_var, ncol = d)

        loc <- coords <- LOCS

	n <- nrow(loc)

        SS <- list()

        dist0 <- parDist(x = loc, method = "euclidean") %>% as.matrix()

	tempcov <- exp(-dist0^2 / theta[2])

	SS[[1]] <- tempcov 

	for(tt in 1:max_time_lag){
		denom <- (diag(d) + Sigma * tt^2 )
		denom_inv <- solve(denom)

		tempcov <- apply(h_mat, 1, quadratic_func_slow, MU = w, DENOM_INV = denom_inv, U = tt, scale = theta[2]) 
		tempcov <- matrix(tempcov, ncol = n, nrow = n) / sqrt(det(denom))

		SS[[tt + 1]] <- tempcov
	}
        S <- theta[1] * toeplitz_mat(SS)
        return(S)
}

array_split <- function(data, number_of_chunks) {
  	# [Partition matrix into N equally-sized chunks with R](https://stackoverflow.com/a/45198299/395857)
  	rowIdx <- seq_len(nrow(data))
  	lapply(split(rowIdx, cut(rowIdx, pretty(rowIdx, number_of_chunks))), function(x) data[x, ])
}


nonfrozen_schlather_uni_cov_distributed <- function(theta, locs, TT) {

  	# setup parallel backend to use many processors
  	cores=detectCores()
	number_of_cores_to_use = cores[1]-1 # not to overload the computer
	#cat(paste('number_of_cores_to_use:',number_of_cores_to_use))
	cl <- makeCluster(number_of_cores_to_use) 
	clusterExport(cl=cl, varlist=c('weights'))
	registerDoParallel(cl)

	#cat('...Starting array split')
	number_of_chunks = number_of_cores_to_use

	n <- nrow(locs)

	LOCS <- NULL

        for(tt in 0:(TT - 1)){
                LOCS <- rbind(LOCS, cbind(locs, rep(tt, n)))
        }

	n <- nrow(LOCS)

	data_chunks = array_split(data = LOCS, number_of_chunks=number_of_chunks)

	clusterEvalQ(cl, source("/home/salvanmo/Desktop/univariate-nonstationary-lagrangian/R_codes/Functions/cov_func.R"))

	#nug <- theta[3]
        scale <- theta[2]
        sigma2 <- theta[1]

	MU <- w <- theta[3:4]
        d <- length(w) 
        SIGMA <- matrix(c(theta[5], theta[7], theta[7], theta[6]), ncol = d)

	clusterExport(cl, c("scale", "MU", "SIGMA"), envir = environment())

  	#cat('...Starting foreach initialization')
  	output <- foreach(i=1:length(data_chunks), .combine=rbind) %dopar% {

		data_temporary = data_chunks[[i]]
		data_full <- LOCS		

		output_temporary = matrix(0, nrow = nrow(data_temporary) * nrow(data_full), ncol = 1)
		for(i in 1:length(data_temporary[, 1])) {
			for(j in 1:length(data_full[, 1])) {
				h_temp <- data_temporary[i, -ncol(data_temporary)] - data_full[j, -ncol(data_temporary)]
				u_temp <- data_temporary[i, ncol(data_temporary)] - data_full[j, ncol(data_temporary)]
				cov_val <- quadratic_func(h_temp, MU = MU, SIGMA = SIGMA, U = u_temp, scale = scale)
				output_temporary[(i - 1) * length(data_full[, 1]) + j, ] = cov_val
			}
		}
		return(output_temporary)
  	}

  	# stop cluster
  	#cat('...Stop cluster')

  	stopCluster(cl)
	
	S <- sigma2 * matrix(output, ncol = n, nrow = n, byrow = T) 
 	#S <- S + nug * diag(n)
        return(S) 
}

nonfrozen_matern_uni_cov_distributed_montecarlo <- function(theta, locs, TT) {

  	# setup parallel backend to use many processors
  	cores=detectCores()
	number_of_cores_to_use = cores[1]-1 # not to overload the computer
	#cat(paste('number_of_cores_to_use:',number_of_cores_to_use))
	cl <- makeCluster(number_of_cores_to_use) 
	registerDoParallel(cl)

	n <- nrow(locs) 

	clusterExport(cl, c("theta", "root"), envir = environment())
	clusterEvalQ(cl, source(paste(root, "R_codes/Functions/load_packages.R", sep = '')))
	clusterEvalQ(cl, sourceCpp(paste(root, "R_codes/Functions/nonfrozen_integration.cpp", sep = '')))

  	#cat('...Starting foreach initialization')
  	system.time(output <- foreach(i=1:40, .combine = '+', .packages = "MASS", .noexport = "nonfrozen_matern_uni_cpp_distributed") %dopar% {

		nu <- 1
		beta <- 0.23
		sigma2 <- 1

		S <- matrix(0, ncol = n * TT, nrow = n * TT)

		for(round in 1:100){

			LOCS <- coords <- locs

			for (tt in 1:(TT - 1)){
				wind <- mvrnorm(n, mu = rep(0.1, 2), Sigma = 0.001 * diag(2))
				temp_locs <- cbind(coords[, 1] - tt * wind[, 1], coords[, 2] - tt * wind[, 2])
				LOCS <- rbind(LOCS, temp_locs)
			}
			
			dist0 <- dist(x = LOCS, diag = TRUE, upper = TRUE) %>% as.matrix()

			S <- S + ifelse(dist0 != 0, (dist0 / beta)^nu * besselK(dist0 / beta, nu) / (2^(nu - 1) * gamma(nu)), sigma2)
		}

		return(S)
  	})

  	# stop cluster
  	#cat('...Stop cluster')

  	stopCluster(cl)
	
	#S <- matrix(output, ncol = nrow(LOCS), nrow = nrow(LOCS), byrow = T) 
 	
        return(output) 
}

nonfrozen_matern_uni_cov_distributed_nonstationary <- function(theta, locs, TT) {

  	# setup parallel backend to use many processors
  	cores=detectCores()
	number_of_cores_to_use = cores[1]-1 # not to overload the computer
	#cat(paste('number_of_cores_to_use:',number_of_cores_to_use))
	cl <- makeCluster(number_of_cores_to_use) 
	registerDoParallel(cl)

	#cat('...Starting array split')
	number_of_chunks = number_of_cores_to_use
	n <- nrow(locs)

	LOCS <- NULL

        for(tt in 0:(TT - 1)){
                LOCS <- rbind(LOCS, cbind(locs, rep(tt, n)))
        }

	data_chunks = array_split(data = LOCS, number_of_chunks=number_of_chunks)

	clusterExport(cl, c("theta", "root"), envir = environment())
	clusterEvalQ(cl, source(paste(root, "R_codes/Functions/load_packages.R", sep = '')))
	clusterEvalQ(cl, sourceCpp(paste(root, "R_codes/Functions/nonfrozen_integration.cpp", sep = '')))

  	#cat('...Starting foreach initialization')
  	system.time(output <- foreach(i=1:length(data_chunks), .combine=rbind, .packages = "Rcpp", .noexport = "nonfrozen_matern_uni_cpp_distributed") %dopar% {
	
		data_temporary = data_chunks[[i]]
		data_full <- LOCS		
		output_temporary <- nonfrozen_matern_uni_cpp_distributed_nonstat(data_temporary, data_full, theta)

		return(output_temporary)
  	})

  	# stop cluster
  	#cat('...Stop cluster')

  	stopCluster(cl)
	
	#S <- matrix(output, ncol = nrow(LOCS), nrow = nrow(LOCS), byrow = T) 
 	
        return(output) 
}

nonfrozen_matern_uni_cov_distributed_deform <- function(theta, locs, TT) {

  	# setup parallel backend to use many processors
  	cores=detectCores()
	number_of_cores_to_use = cores[1]-1 # not to overload the computer
	#cat(paste('number_of_cores_to_use:',number_of_cores_to_use))
	cl <- makeCluster(number_of_cores_to_use) 
	registerDoParallel(cl)

	#cat('...Starting array split')
	number_of_chunks = number_of_cores_to_use
	n <- nrow(locs)

	LOCS <- NULL

        for(tt in 0:(TT - 1)){
                LOCS <- rbind(LOCS, cbind(locs, rep(tt, n)))
        }

	data_chunks = array_split(data = LOCS, number_of_chunks=number_of_chunks)

	clusterExport(cl, c("theta", "root"), envir = environment())
	clusterEvalQ(cl, source(paste(root, "R_codes/Functions/load_packages.R", sep = '')))
	clusterEvalQ(cl, sourceCpp(paste(root, "R_codes/Functions/nonfrozen_integration.cpp", sep = '')))

  	#cat('...Starting foreach initialization')
  	system.time(output <- foreach(i=1:length(data_chunks), .combine=rbind, .packages = "Rcpp", .noexport = "nonfrozen_matern_uni_cpp_distributed") %dopar% {
	
		data_temporary = data_chunks[[i]]
		data_full <- LOCS		
		output_temporary <- nonfrozen_matern_uni_cpp_distributed_deform(data_temporary, data_full, theta)

		return(output_temporary)
  	})

  	# stop cluster
  	#cat('...Stop cluster')

  	stopCluster(cl)
	
	#S <- matrix(output, ncol = nrow(LOCS), nrow = nrow(LOCS), byrow = T) 
 	
        return(output) 
}

nonfrozen_matern_uni_cov_distributed <- function(theta, locs, TT) {

  	# setup parallel backend to use many processors
  	cores=detectCores()
	number_of_cores_to_use = cores[1]-1 # not to overload the computer
	#cat(paste('number_of_cores_to_use:',number_of_cores_to_use))
	cl <- makeCluster(number_of_cores_to_use) 
	registerDoParallel(cl)

	#cat('...Starting array split')
	number_of_chunks = number_of_cores_to_use
	n <- nrow(locs)

	LOCS <- NULL

        for(tt in 0:(TT - 1)){
                LOCS <- rbind(LOCS, cbind(locs, rep(tt, n)))
        }

	data_chunks = array_split(data = LOCS, number_of_chunks=number_of_chunks)

	clusterExport(cl, c("theta", "root"), envir = environment())
	clusterEvalQ(cl, source(paste(root, "R_codes/Functions/load_packages.R", sep = '')))
	clusterEvalQ(cl, sourceCpp(paste(root, "R_codes/Functions/nonfrozen_integration.cpp", sep = '')))

  	#cat('...Starting foreach initialization')
  	system.time(output <- foreach(i=1:length(data_chunks), .combine=rbind, .packages = "Rcpp", .noexport = "nonfrozen_matern_uni_cpp_distributed") %dopar% {
		
		data_temporary = data_chunks[[i]]
		data_full <- LOCS		
		
		output_temporary <- nonfrozen_matern_uni_cpp_distributed(data_temporary, data_full, theta)

		return(output_temporary)
  	})

  	# stop cluster
  	#cat('...Stop cluster')

  	stopCluster(cl)
	
	#S <- matrix(output, ncol = nrow(LOCS), nrow = nrow(LOCS), byrow = T) 
 	
        return(output) 
}

nonfrozen_matern_uni_cov_distributed_SLOW <- function(theta, locs, TT) {

  	# setup parallel backend to use many processors
  	cores=detectCores()
	number_of_cores_to_use = cores[1]-1 # not to overload the computer
	#cat(paste('number_of_cores_to_use:',number_of_cores_to_use))
	cl <- makeCluster(number_of_cores_to_use) 
	registerDoParallel(cl)

	#cat('...Starting array split')
	number_of_chunks = number_of_cores_to_use
	n <- nrow(locs)

	LOCS <- NULL

        for(tt in 0:(TT - 1)){
                LOCS <- rbind(LOCS, cbind(locs, rep(tt, n)))
        }

	data_chunks = array_split(data = LOCS, number_of_chunks=number_of_chunks)

	clusterExport(cl, c("theta", "root"), envir = environment())
	clusterEvalQ(cl, source(paste(root, "R_codes/Functions/load_packages.R", sep = '')))
	clusterEvalQ(cl, sourceCpp(paste(root, "R_codes/Functions/nonfrozen_integration.cpp", sep = '')))

  	#cat('...Starting foreach initialization')
  	output <- foreach(i=1:length(data_chunks), .combine=rbind, .packages = "Rcpp", .noexport = "nonfrozen_matern_uni_cpp_distributed") %dopar% {
		
		data_temporary = data_chunks[[i]]
		data_full <- LOCS		

		output_temporary = matrix(0, nrow = nrow(data_temporary) * nrow(data_full), ncol = 1)
		for(i in 1:length(data_temporary[, 1])) {
			cov_val <- nonfrozen_matern_uni_cpp_distributed(matrix(data_temporary[i, ], ncol = 3), data_full, theta)
			output_temporary[(i - 1) * length(data_full[, 1]) + 1:nrow(data_full), ] = cov_val
	
		}
		return(output_temporary)
  	}

  	# stop cluster
  	#cat('...Stop cluster')

  	stopCluster(cl)
	
	S <- matrix(output, ncol = nrow(LOCS), nrow = nrow(LOCS), byrow = T) 
 	
        return(S) 
}

nonfrozen_matern_uni_cov_distributed_FAST <- function(theta, locs, TT) {

	Sys.getenv("SLURM_NODELIST") #get names of nodes

	nodelist <- unlist(strsplit(Sys.getenv("SLURM_NODELIST"), ","))

	#cl <- makeCluster(rep(nodelist , each = 10))

  	# setup parallel backend to use many processors
  	cores=detectCores()
	number_of_cores_to_use = cores[1]-1 # not to overload the computer
	#cat(paste('number_of_cores_to_use:',number_of_cores_to_use))
	cl <- makeCluster(number_of_cores_to_use) 
	clusterExport(cl=cl, varlist=c('weights'))
	registerDoParallel(cl)

	#cat('...Starting array split')
	number_of_chunks = number_of_cores_to_use
	n <- nrow(locs)

	LOCS <- NULL

        for(tt in 0:(TT - 1)){
                LOCS <- rbind(LOCS, cbind(locs, rep(tt, n)))
        }

	data_chunks = array_split(data = LOCS, number_of_chunks=number_of_chunks)
	#data_chunks = array_split(data = LOCS, number_of_chunks=15)

	chunks_pairs <- list()

	for(aa in 1:length(data_chunks)){
		for(bb in 1:length(data_chunks)){
			chunks_pairs_temp <- list()
			chunks_pairs_temp[[1]] <- data_chunks[[aa]]
			chunks_pairs_temp[[2]] <- data_chunks[[bb]]
			chunks_pairs[[ (aa - 1) * length(data_chunks) + bb]] <- chunks_pairs_temp
		}
	}

	clusterExport(cl, c("theta", "root"), envir = environment())
	clusterEvalQ(cl, source(paste(root, "R_codes/Functions/load_packages.R", sep = '')))
	clusterEvalQ(cl, sourceCpp(paste(root, "R_codes/Functions/nonfrozen_integration.cpp", sep = '')))
	clusterEvalQ(cl, source(paste(root, "R_codes/Functions/cov_func.R", sep = '')))

	
  	#cat('...Starting foreach initialization')
  	system.time(output <- foreach(i=1:length(chunks_pairs), .combine=rbind, .packages = "Rcpp", .noexport = "nonfrozen_matern_uni_cpp_distributed") %dopar% {
		
		data_temporary = chunks_pairs[[i]][[1]]
		data_full <- chunks_pairs[[i]][[2]]	

		output_temporary <- nonfrozen_matern_uni_cpp_distributed(data_temporary, data_full, theta)

		return(output_temporary)
  	})

  	# stop cluster
  	#cat('...Stop cluster')

  	stopCluster(cl)
	
	S <- matrix(output, ncol = nrow(LOCS), nrow = nrow(LOCS), byrow = T) 
 	
        return(S) 
}

nonfrozen_matern_uni_cov_distributed_FAST_NEW <- function(theta, locs, TT) {

  	# setup parallel backend to use many processors
  	cores=detectCores()
	number_of_cores_to_use = cores[1]-1 # not to overload the computer
	#cat(paste('number_of_cores_to_use:',number_of_cores_to_use))
	cl <- makeCluster(number_of_cores_to_use) 
	clusterExport(cl=cl, varlist=c('weights'))
	registerDoParallel(cl)

	#cat('...Starting array split')
	number_of_chunks = number_of_cores_to_use
	n <- nrow(locs)

	LOCS <- NULL

        for(tt in 0:(TT - 1)){
                LOCS <- rbind(LOCS, cbind(locs, rep(tt, n)))
        }

	n <- nrow(LOCS)

	seq_ind <- ceiling(seq(1, number_of_chunks - 1, by = 1) / (number_of_chunks * (number_of_chunks + 1) / 2) * n)
	seq_ind <- c(seq_ind, n - sum(seq_ind))
	seq_ind2 <- c(0, cumsum(seq_ind))

	data_chunks <- list()
	for(aa in 1:length(seq_ind)){
		data_chunks[[aa]] <- LOCS[seq_ind2[aa] + 1:seq_ind[aa], ]
	}

	clusterExport(cl, c("theta", "root"), envir = environment())
	clusterEvalQ(cl, source(paste(root, "R_codes/Functions/load_packages.R", sep = '')))
	clusterEvalQ(cl, sourceCpp(paste(root, "R_codes/Functions/nonfrozen_integration.cpp", sep = '')))

  	output <- foreach(i=1:length(data_chunks), .combine=c, .packages = "Rcpp", .noexport = "nonfrozen_matern_uni_cpp_distributed") %dopar% {
		
		data_temporary = data_chunks[[i]]
		data_full <- LOCS		

		ind <- seq_ind2[i]	

		output_temporary = NULL
		for(k in 1:length(data_temporary[, 1])) {
			cov_val <- nonfrozen_matern_uni_cpp_distributed(matrix(data_temporary[k, ], ncol = 3), matrix(data_full[(ind + k):n, ], ncol = 3), theta)
			cat(length(cov_val), '\n')
			output_temporary <- c(output_temporary, cov_val)
	
		}
		return(output_temporary)
  	}

  	stopCluster(cl)
	
	S <- matrix(0, ncol = n, nrow = n, byrow = T) 

	S[lower.tri(S, diag=T)] <- output
	S <- t(S)
	S[lower.tri(S, diag = T)] <- output 	
	
        return(S) 
}


nonfrozen_matern_uni_cov_distributed_ORIG <- function(theta, wind, wind_var, locs) {

  	# setup parallel backend to use many processors
  	cores=detectCores()
	number_of_cores_to_use = cores[1]-1 # not to overload the computer
	#cat(paste('number_of_cores_to_use:',number_of_cores_to_use))
	cl <- makeCluster(number_of_cores_to_use) 
	clusterExport(cl=cl, varlist=c('weights'))
	registerDoParallel(cl)

	#cat('...Starting array split')
	number_of_chunks = number_of_cores_to_use
	data_chunks = array_split(data = locs, number_of_chunks=number_of_chunks)

	clusterEvalQ(cl, source("/home/salvanmo/Desktop/univariate-nonstationary-lagrangian/R_codes/Functions/cov_func.R"))
	clusterEvalQ(cl, source("/home/salvanmo/Desktop/univariate-nonstationary-lagrangian/R_codes/Functions/num_integ_func.R"))

	n <- nrow(locs)

	nu <- theta[3]
        scale <- theta[2]
        sigma2 <- theta[1]

	MU <- w <- wind
        d <- length(w) 
        SIGMA <- matrix(wind_var, ncol = d)

	clusterExport(cl, c("theta", "MU", "SIGMA"), envir = environment())

  	#cat('...Starting foreach initialization')
  	output <- foreach(i=1:length(data_chunks), .combine=rbind) %dopar% {
		
		data_temporary = data_chunks[[i]]
		data_full <- locs		

		output_temporary = matrix(0, nrow = nrow(data_temporary) * nrow(data_full), ncol = 1)
		for(i in 1:length(data_temporary[, 1])) {
			for(j in 1:length(data_full[, 1])) {
				h_temp <- data_temporary[i, -ncol(data_temporary)] - data_full[j, -ncol(data_temporary)]
				u_temp <- data_temporary[i, ncol(data_temporary)] - data_full[j, ncol(data_temporary)]
				cov_val <- numerical_integ(h_temp[1:2], MU = MU[1:2], SIGMA = SIGMA[1:2, 1:2], U = u_temp, THETA = theta)
				output_temporary[(i - 1) * length(data_full[, 1]) + j, ] = cov_val
			}
		}
		return(output_temporary)
  	}

  	# stop cluster
  	#cat('...Stop cluster')

  	stopCluster(cl)
	
	S <- matrix(output, ncol = n, nrow = n, byrow = T) 
 	
        return(S) 
}

nonfrozen_schlather_uni_cov_distributed_ORIG <- function(theta, wind, wind_var, tt, h_mat) {

  	# setup parallel backend to use many processors
  	cores=detectCores()
	number_of_cores_to_use = cores[1]-1 # not to overload the computer
	#cat(paste('number_of_cores_to_use:',number_of_cores_to_use))
	cl <- makeCluster(number_of_cores_to_use) 
	clusterExport(cl=cl, varlist=c('weights'))
	registerDoParallel(cl)

	#cat('...Starting array split')
	number_of_chunks = number_of_cores_to_use
	data_chunks = array_split(data = h_mat, number_of_chunks=number_of_chunks)

	clusterEvalQ(cl, source("/home/salvanmo/Desktop/univariate-nonstationary-lagrangian/R_codes/Functions/cov_func.R"))

	scale <- theta
	MU <- w <- wind
        d <- length(w) 
        Sigma <- matrix(wind_var, ncol = d)

	U <- tt
	DENOM <- denom <- (diag(d) + Sigma * tt^2 )
	DENOM_INV <- denom_inv <- solve(denom)

	clusterExport(cl, c("scale", "MU", "DENOM_INV", "U"), envir = environment())

  	#cat('...Starting foreach initialization')
  	output <- foreach(i=1:length(data_chunks), .combine=rbind) %dopar% {

		data_temporary = data_chunks[[i]]
		output_temporary = matrix(0, nrow=nrow(data_temporary), ncol = 1)
		for(i in 1:length(data_temporary[,1])) {
			cov_val <- quadratic_func_slow(data_temporary[i, ], MU = MU, DENOM_INV = DENOM_INV, U = U, scale = scale)
			output_temporary[i,] = cov_val
		}
		return(output_temporary)
  	}

  	# stop cluster
  	#cat('...Stop cluster')

  	stopCluster(cl)

  	return(output)
}

nonfrozen_schlather_uni_cov_distributed_FINAL_ORIG <- function(theta, wind, wind_var, max_time_lag, LOCS, h_mat){

	scale <- theta
        loc <- coords <- LOCS

	n <- nrow(loc)

        SS <- list()

        dist0 <- parDist(x = loc, method = "euclidean") %>% as.matrix()

	tempcov <- exp(-dist0^2 / theta)

	SS[[1]] <- tempcov 

	scale <- theta
	MU <- w <- wind
        d <- length(w) 
        Sigma <- matrix(wind_var, ncol = d)

	for(tt in 1:max_time_lag){
		DENOM <- denom <- (diag(d) + Sigma * tt^2 )
		DENOM_INV <- denom_inv <- solve(denom)

		temploc <- nonfrozen_schlather_uni_cov_distributed(theta, wind, wind_var, tt = tt, h_mat)
		temploc <- matrix(temploc, ncol = n, nrow = n, byrow = T)
		tempcov <- ifelse(temploc == 0, 1, temploc / sqrt(det(denom)))

		SS[[tt + 1]] <- tempcov
	}

        S <- toeplitz_mat(SS)

        return(S)
}

nonfrozen_matern_uni_cov_numerical <- function(theta, wind, wind_var, max_time_lag, LOCS, h_mat){

        nu <- theta[1]
        beta <- theta[2]
        var <- theta[3]

	MEAN <- w <- wind
        d <- length(w) 
        SIGMA <- matrix(wind_var, ncol = d)

	clusterExport(cl, c("MEAN", "SIGMA", "nu", "beta", "var"))

	U <- tt <- 1

        loc <- coords <- LOCS

	n <- nrow(loc)

        SS <- list()

        dist0 <- parDist(x = loc, method = "euclidean") %>% as.matrix()

	temploc <- var * Matern(dist0, range = beta, nu = nu) 
	SS[[1]] <- temploc 

	temploc <- adply(h_mat, .margins = 1, .parallel = T, .fun = numerical_integ)
	temploc <- matrix(temploc[, 2], ncol = n, nrow = n, byrow = T)

	SS[[2]] <- temploc

        S <- toeplitz_mat(SS)

        return(S)
}

nonfrozen_matern_uni_cov <- function(theta, wind, wind_var, max_time_lag, LOCS, h_mat){

        nu <- theta[1]
        beta <- theta[2]
        var <- theta[3]

        w <- wind
        Sigma <- matrix(c(wind_var[1:2], wind_var[2:3]), ncol = 2)

	clusterExport(cl, c("nu", "beta", "var", "w", "Sigma"))

        loc <- coords <- LOCS

	n <- nrow(loc)

        SS <- list()

        dist0 <- parDist(x = loc, method = "euclidean") %>% as.matrix()

	temploc <- ifelse(dist0 != 0, var * (dist0 / beta)^nu * besselK(dist0 / beta, nu) / (2^(nu - 1) * gamma(nu)), var)

	SS[[1]] <- temploc 

	temploc <- adply(h_mat, .margins = 1, .parallel = T, .fun = num_integ)
	temploc <- matrix(temploc[, 2], ncol = n, nrow = n, byrow = T)
	tempcov <- ifelse(temploc == 0, var, var * 4 * pi * temploc * beta^2 / gamma(nu))

	SS[[2]] <- tempcov
	
	#stopCluster(cl)
	
        S <- toeplitz_mat(SS)
        return(S)
}

#-------MULTIVARIATE----------#

frozen_matern_cov <- function(theta, wind, max_time_lag = 0, q = 2, LOCS){

	###################												###################
	################### RETURNS a q * nrow(LOCS) * (max_time_lag + 1) x q * nrow(LOCS) * (max_time_lag + 1) matrix 	###################
	###################												###################	

	w <- wind
	
	locs <- coords <- LOCS
  
  	if (max_time_lag > 0){
		for (tt in 1:max_time_lag){
			temp_locs <- cbind(coords[, 1] - tt * w[1], coords[, 2] - tt * w[2])
			locs <- rbind(locs, temp_locs)
		}
  	} 
  
	nu <- theta[4:5]
	beta <- theta[3]
	sigma2 <- theta[1:2]
  
	dist0 <- parDist(x = locs, method = "euclidean") %>% as.matrix()	

	S <- matrix(NA,  q * nrow(dist0), q * nrow(dist0))
  
	for(i in 1:q){
		for(j in 1:i){

			temp <- (i - 1) * nrow(dist0) + 1:nrow(dist0)
	      		temp1 <- (j - 1) * nrow(dist0) + 1:nrow(dist0)
	      
	      		if(i == j){
		
				temp2 <- ifelse(dist0 != 0, sigma2[i] * (dist0 / beta)^nu[i] * besselK(dist0 / beta, nu[i]) / (2^(nu[i] - 1) * gamma(nu[i])), sigma2[i])
				S[temp, temp1] <- temp2
		
	      		}
	      
		      	if(i != j){
			
				nu1 <- nu[i]
				nu2 <- nu[j]
				nu3 <- (nu1 + nu2)/2
			
				rho <- theta[6] * (gamma(nu1 + 3/2) / gamma(nu1))^(1/2) * (gamma(nu2 + 3/2) / gamma(nu2))^(1/2) * gamma(nu3) / (gamma(nu3 + 3/2))
			
				temp3 <- (dist0 / beta)^nu3 * besselK(dist0 / beta, nu3)/(2^(nu3 - 1) * gamma(nu3)) * sqrt(sigma2[i] * sigma2[j]) * rho
				temp3[is.na(temp3)] <- sqrt(sigma2[i] * sigma2[j]) * rho
				S[temp, temp1] <- temp3
				S[temp1, temp] <- t(temp3)
		      }
	    	}
	  }
  
	return(S)
}

frozen_matern_cov_rep_I <- function(theta, wind, max_time_lag = 0, q = 2, LOCS){

	###################												###################
	################### RETURNS a q * nrow(LOCS) * (max_time_lag + 1) x q * nrow(LOCS) * (max_time_lag + 1) matrix 	###################
	###################												###################	

	w <- wind
	
	locs <- coords <- LOCS
  
  	if (max_time_lag > 0){
		for (tt in 1:max_time_lag){
			temp_locs <- cbind(coords[, 1] - tt * w[1], coords[, 2] - tt * w[2])
			locs <- rbind(locs, temp_locs)
		}
  	} 
  
	nu <- theta[4:5]
	beta <- theta[3]
	sigma2 <- theta[1:2]
  
	dist0 <- parDist(x = locs, method = "euclidean") %>% as.matrix()	

	S <- matrix(NA,  q * nrow(dist0), q * nrow(dist0))
  
	for(i in 1:q){
		for(j in 1:i){

			temp <- seq(i, nrow(dist0) * q, by = q)
	      		temp1 <- seq(j, nrow(dist0) * q, by = q)
	      
	      		if(i == j){
		
				temp2 <- ifelse(dist0 != 0, sigma2[i] * (dist0 / beta)^nu[i] * besselK(dist0 / beta, nu[i]) / (2^(nu[i] - 1) * gamma(nu[i])), sigma2[i])
				S[temp, temp1] <- temp2
		
	      		}
	      
		      	if(i != j){
			
				nu1 <- nu[i]
				nu2 <- nu[j]
				nu3 <- (nu1 + nu2)/2
			
				rho <- theta[6] * (gamma(nu1 + 3/2) / gamma(nu1))^(1/2) * (gamma(nu2 + 3/2) / gamma(nu2))^(1/2) * gamma(nu3) / (gamma(nu3 + 3/2))
			
				temp3 <- (dist0 / beta)^nu3 * besselK(dist0 / beta, nu3)/(2^(nu3 - 1) * gamma(nu3)) * sqrt(sigma2[i] * sigma2[j]) * rho
				temp3[is.na(temp3)] <- sqrt(sigma2[i] * sigma2[j]) * rho
				S[temp, temp1] <- temp3
				S[temp1, temp] <- t(temp3)
		      }
	    	}
	  }
  
	return(S)
}

frozen_matern_cov_for_heatmap <- function(theta, wind, q = 2){

	###################															###################
	################### RETURNS a q * nrow(LOCS) * (max_time_lag + 1) x 3 matrix of covariance values with respect to location (0,0) 	###################
	###################															###################	

	N <- 51
        n <- N^2
        TT <- 3
        grid_x <- seq(from = -0.5, to = 0.5, length.out = N)
	sim_grid_locations <- expand.grid(grid_x, grid_x) %>% as.matrix()

	w <- wind
	
	nu <- theta[4:5]
	beta <- theta[3]
	sigma2 <- theta[1:2]

	nu1 <- nu[1]
	nu2 <- nu[2]
	nu3 <- (nu1 + nu2)/2
	rho <- theta[6] * (gamma(nu1 + 3/2) / gamma(nu1))^(1/2) * (gamma(nu2 + 3/2) / gamma(nu2))^(1/2) * gamma(nu3) / (gamma(nu3 + 3/2))

	S <- matrix(NA,  n * TT, 3)

	for(i in 1:3){
		for(tt in 0:(TT - 1)){

			for(l in 1:nrow(sim_grid_locations)){

				temp_locs <- rbind(cbind(0, 0), cbind(sim_grid_locations[l, 1] - tt * w[1], sim_grid_locations[l, 2] - tt * w[2]))
				dist0 <- dist(temp_locs) %>% as.numeric()
				
				if(i < 3) temp2 <- ifelse(dist0 != 0, sigma2[i] * (dist0 / beta)^nu[i] * besselK(dist0 / beta, nu[i]) / (2^(nu[i] - 1) * gamma(nu[i])), sigma2[i])
				else temp2 <- ifelse(dist0 != 0, (dist0 / beta)^nu3 * besselK(dist0 / beta, nu3)/(2^(nu3 - 1) * gamma(nu3)) * sqrt(sigma2[1] * sigma2[2]) * rho, sqrt(sigma2[1] * sigma2[2]) * rho)
				S[tt * n + l, i] <- temp2 
			}
		}
	}

	return(S)
}



#---------NONSTATIONARY------------#

frozen_matern_dimension_expansion_uni_cov <- function(theta, wind, max_time_lag = 0, LOCS){

        nu <- theta[3]
        beta <- theta[2]
        sigma2 <- theta[1]

	w <- wind

        locs <- coords <- LOCS

        if (max_time_lag > 0){
                for (tt in 1:max_time_lag){
                        temp_locs <- cbind(coords[, 1] - tt * w[1], coords[, 2] - tt * w[2], coords[, 3] - tt * w[3])
                        locs <- rbind(locs, temp_locs)
                }
        }
   	
        dist0 <- parDist(x = locs, method = "euclidean") %>% as.matrix()

	S <- ifelse(dist0 != 0, (dist0 / beta)^nu * besselK(dist0 / beta, nu) / (2^(nu - 1) * gamma(nu)), sigma2)

  	return(S)
}

frozen_matern_deform_uni_cov <- function(theta, wind, max_time_lag = 0, LOCS){

        nu <- theta[3]
        beta <- theta[2]
        sigma2 <- theta[1]

	w <- wind

        locs <- coords <- LOCS

        if (max_time_lag > 0){
                for (tt in 1:max_time_lag){
                        temp_locs <- cbind(coords[, 1] - tt * w[1], coords[, 2] - tt * w[2])
                        locs <- rbind(locs, temp_locs)
                }
        }
   	
	new_locations <- c(0.5, 0.5) + (locs - c(0.5, 0.5)) * as.numeric(distR_C(locs, matrix(c(0.5, 0.5), nrow = 1)))

        dist0 <- parDist(x = new_locations, method = "euclidean") %>% as.matrix()

	S <- ifelse(dist0 != 0, (dist0 / beta)^nu * besselK(dist0 / beta, nu) / (2^(nu - 1) * gamma(nu)), sigma2)

  	return(S)
}

frozen_matern_spatially_varying_covariates_uni_cov <- function(theta, wind, max_time_lag = 0, LOCS, NONSTAT_PARAMS){

	nu <- theta[2]
        sigma2 <- theta[1]

        w <- wind

        locs <- coords <- LOCS

        if (max_time_lag > 0){
                for (tt in 1:max_time_lag){
                        temp_locs <- cbind(coords[, 1] - tt * w[1], coords[, 2] - tt * w[2])
                        locs <- rbind(locs, temp_locs)
                }
        }	

  	S <- spatially_varying_parameters2(Loc = locs, param = c(sigma2, nu), Nonstat_params = NONSTAT_PARAMS)

  	return(S)
}

matern_spatially_varying_covariates_uni_cov <- function(theta, LOCS, NONSTAT_PARAMS){

	nu <- theta[3]
	beta <- theta[2]
        sigma2 <- theta[1]

        locs <- coords <- LOCS

  	S <- spatially_varying_parameters2(Loc = locs, param = c(sigma2, beta, nu), Nonstat_params = NONSTAT_PARAMS)

  	return(S)
}

matern_deform_uni_cov <- function(theta, LOCS){

        nu <- theta[3]
        beta <- theta[2]
        sigma2 <- theta[1]

        locs <- coords <- LOCS

	new_locations <- c(0.5, 0.5) + (locs - c(0.5, 0.5)) * as.numeric(distR_C(locs, matrix(c(0.5, 0.5), nrow = 1)))

        dist0 <- parDist(x = new_locations, method = "euclidean") %>% as.matrix()

	S <- ifelse(dist0 != 0, (dist0 / beta)^nu * besselK(dist0 / beta, nu) / (2^(nu - 1) * gamma(nu)), sigma2)

  	return(S)
}


#------------------------------- END ----------------------------#



matern_cov_soph <- function(theta, wind, max_time_lag, q, new_locations = locations, meters = T, nug_eff = F, kap = F){
  
  w <- wind
  
  loc1 <- coords1 <- new_locations
  
  if (max_time_lag == 0){
    loc1 <- loc1
  } else {
    for (tt in 1:max_time_lag){
      temploc <- matrix(, ncol=2, nrow=nrow(coords1))
      for(rr in 1:nrow(coords1)){
        temploc[rr,] <- c(coords1[rr,1] - tt*w[1], coords1[rr,2] - tt*w[2])
      }
      loc1 <- rbind(loc1, temploc)
    }
  }
  
  loc2 <- coords2 <- cbind(new_locations[,1] - kap[1], new_locations[,2] - kap[2])
  
  if (max_time_lag == 0){
    loc2 <- loc2
  } else {
    for (tt in 1:max_time_lag){
      temploc <- matrix(, ncol=2, nrow=nrow(coords2))
      for(rr in 1:nrow(coords2)){
        temploc[rr,] <- c(coords2[rr,1] - tt*w[1], coords2[rr,2] - tt*w[2])
      }
      loc2 <- rbind(loc2, temploc)
    }
  }
  loc <- rbind(loc1, loc1, loc2)
  
  if(meters == T){
    dist0 <- spDists(loc, longlat=F)/1000
  }else{
    dist0 <- spDists(loc, longlat=F)
  }
  
  nu <- theta[1:2]
  beta <- theta[3]
  var <- theta[4:5]
  rho <- theta[6]
  
  if(nug_eff == T){
    nug <- theta[7:8]
  }else{
    nug <- c(0, 0)
  }
  
  S=matrix(NA,  q*dim(dist0)[1], q*dim(dist0)[1])
  
  for(i in 1:q){
    for(j in 1:i){
      
      temp=(i-1)*dim(dist0)[1]+1:dim(dist0)[1]
      temp1=(j-1)*dim(dist0)[1]+1:dim(dist0)[1]
      
      if(i==j){
        
        temp2=ifelse(dist0!=0,var[i]*(dist0/beta)^nu[i] * besselK(dist0/beta, nu[i])/(2^(nu[i]-1)*gamma(nu[i])),var[i]+nug[i])
        S[temp,temp1]=temp2
        
      }
      
      if(i != j){
        
        nu1 <- nu[i]
        nu2 <- nu[j]
        nu3 <- (nu1 + nu2)/2
        
        #rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        
        temp3 <- (dist0/beta)^nu3 * besselK(dist0/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(var[i] * var[j])*rho
        temp3[is.na(temp3)] <- sqrt(var[i] * var[j])*rho
        S[temp,temp1] <- temp3
        S[temp1,temp] <- t(temp3)
      }
    }
  }
  
  S1 <- rbind(cbind(S[1:nrow(loc1), 1:nrow(loc1)], S[1:nrow(loc1), (nrow(loc1)*3 + 1):(nrow(loc1)*4)]),
              cbind(S[(nrow(loc1)*3 + 1):(nrow(loc1)*4), 1:nrow(loc1)], S[(nrow(loc1)*3 + 1):(nrow(loc1)*4), (nrow(loc1)*3 + 1):(nrow(loc1)*4)]))
  return(S1)
}

matern_cov_old <- function(theta, wind, max_time_lag, q, new_locations = locations, meters = T, nug_eff, kap = matrix(c(704400+ 100, 205700 + 100, 100, 100), ncol = 2, nrow = 2, byrow=T)){
  
  w <- wind
  
  loc1 <- coords1 <- cbind(new_locations[,1] + kap[1,1] - kap[2,1], new_locations[,2] + kap[1,2] - kap[2,2])
  
  if (max_time_lag == 0){
    loc1 <- loc1
  } else {
    for (tt in 1:max_time_lag){
      temploc <- matrix(, ncol=2, nrow=nrow(coords1))
      for(rr in 1:nrow(coords1)){
        temploc[rr,] <- c(coords1[rr,1] - tt*w[1], coords1[rr,2] - tt*w[2])
      }
      loc1 <- rbind(loc1, temploc)
    }
  }
  
  loc2 <- coords2 <- cbind(new_locations[,1] - kap[1,1] + kap[2,1], new_locations[,2] - kap[1,2] + kap[2,2])
  
  if (max_time_lag == 0){
    loc2 <- loc2
  } else {
    for (tt in 1:max_time_lag){
      temploc <- matrix(, ncol=2, nrow=nrow(coords2))
      for(rr in 1:nrow(coords2)){
        temploc[rr,] <- c(coords2[rr,1] - tt*w[1], coords2[rr,2] - tt*w[2])
      }
      loc2 <- rbind(loc2, temploc)
    }
  }
  loc <- rbind(loc1, loc2)
  
  if(meters == T){
    dist0 <- spDists(loc, longlat=F)/1000
  }else{
    dist0 <- spDists(loc, longlat=F)
  }
  
  nu <- theta[1:2]
  beta <- theta[3]
  var <- theta[4:5]
  rho <- theta[6]
  
  if(nug_eff == T){
    nug <- theta[7:8]
  }else{
    nug <- c(0, 0)
  }
  
  S=matrix(NA,  q*dim(dist0)[1], q*dim(dist0)[1])
  
  for(i in 1:q){
    for(j in 1:i){
      
      temp=(i-1)*dim(dist0)[1]+1:dim(dist0)[1]
      temp1=(j-1)*dim(dist0)[1]+1:dim(dist0)[1]
      
      if(i==j){
        
        temp2=ifelse(dist0!=0,var[i]*(dist0/beta)^nu[i] * besselK(dist0/beta, nu[i])/(2^(nu[i]-1)*gamma(nu[i])),var[i]+nug[i])
        S[temp,temp1]=temp2
        
      }
      
      if(i != j){
        
        nu1 <- nu[i]
        nu2 <- nu[j]
        nu3 <- (nu1 + nu2)/2
        
        #rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        
        temp3 <- (dist0/beta)^nu3 * besselK(dist0/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(var[i] * var[j])*rho
        temp3[is.na(temp3)] <- sqrt(var[i] * var[j])*rho
        S[temp,temp1] <- temp3
        S[temp1,temp] <- t(temp3)
      }
    }
  }
  
  S1 <- rbind(cbind(S[1:nrow(loc1), 1:nrow(loc1)], S[1:nrow(loc1), (nrow(loc1)*3 + 1):(nrow(loc1)*4)]),
              cbind(S[(nrow(loc1)*3 + 1):(nrow(loc1)*4), 1:nrow(loc1)], S[(nrow(loc1)*3 + 1):(nrow(loc1)*4), (nrow(loc1)*3 + 1):(nrow(loc1)*4)]))
  return(S1)
}

matern_random_cov <- function(theta, wind, wind_var, max_time_lag, q, new_locations, meters = T, nug_eff, kap){
  
  nu <- theta[1:2]
  beta <- theta[3]
  rho <- theta[6]
  var <- theta[4:5]
  if(nug_eff == T){
    nug <- theta[7:8]
  }else{
    nug <- c(0, 0)
  }
  
  if(meters == T){
    w <- wind/1000
    Sigma <- matrix(c(wind_var[1:2], wind_var[2:3]), ncol = 2)/1000
    loc <- coords <- new_locations/1000
    loc2 <- coords2 <- cbind(new_locations[,1] - kap[1], new_locations[,2] - kap[2])/1000
  }else{
    w <- wind
    Sigma <- matrix(c(wind_var[1:2], wind_var[2:3]), ncol=2)
    loc <- coords <- new_locations
    loc2 <- coords2 <- cbind(new_locations[,1] - kap[1], new_locations[,2] - kap[2])/1000
  }
  
  SS <- list()
  
  S <- matrix(NA,  q*nrow(coords)*(max_time_lag + 1), q*nrow(coords)*(max_time_lag + 1))
  
  for(i in 1:q){
    for(j in 1:i){
      
      temp2 <- (i-1)*(nrow(coords)*(max_time_lag + 1)) + 1:(nrow(coords)*(max_time_lag + 1))
      temp1 <- (j-1)*(nrow(coords)*(max_time_lag + 1)) + 1:(nrow(coords)*(max_time_lag + 1))
      
      if(i == j){
        
        for(tt in 0:max_time_lag){
          temploc <- matrix(, ncol=nrow(coords), nrow=nrow(coords))
          for(rr in 1:nrow(coords)){
            for(ss in 1:nrow(coords)){
              cat(tt,rr,ss,'\n')
              h <- c(coords[rr,1]-coords[ss,1],coords[rr,2]-coords[ss,2])
              emp_cov1 <- c(h[1], h[2], tt)
              Int.func <- function(c, hvec){   
                y.fun  <- function(y) y^(nu[i])*exp(-y)*dmvn(X = hvec[1:2], mu = hvec[3]*w, sigma = (hvec[3]^2*Sigma + beta^2*2*y*diag(2)))
                sapply(c, y.fun)
              }
              lai <- function(xxxx) integrate(Int.func, lower = 0, upper = Inf, hvec = xxxx, abs.tol = 1e-18, rel.tol = 1e-18)$val
              temp <- lai(emp_cov1) 
              
              temploc[rr,ss] <- ifelse(tt == 0 & h[1] == 0 & h[2] == 0, var[i], var[i]*4*pi*temp*beta^2/gamma(nu[i]))
            }
          }
          SS[[tt + 1]] <- temploc
        }
        S2 <- toeplitz_mat(SS)
        S[temp2,temp1] <- S2 
        
      }
      
      if(i != j){
        
        nu1 <- nu[i]
        nu2 <- nu[j]
        nu3 <- (nu1 + nu2)/2
        
        #rho=rot*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        
        for(tt in 0:max_time_lag){
          temploc <- matrix(, ncol = nrow(coords), nrow = nrow(coords))
          for(rr in 1:nrow(coords)){
            for(ss in 1:nrow(coords)){
              cat(tt, rr, ss,'\n')
              
              h <- c(coords2[rr,1] - coords2[ss,1], coords2[rr,2] - coords2[ss,2])
              emp_cov1 <- c(h[1], h[2], tt)
              Int.func <- function(c, hvec){   
                y.fun  <- function(y) y^(nu3)*exp(-y)*dmvn(X = hvec[1:2], mu = hvec[3]*w, sigma = (hvec[3]^2*Sigma + beta^2*2*y*diag(2)))
                sapply(c, y.fun)
              }
              lai <- function(xxxx) integrate(Int.func, lower = 0, upper = Inf, hvec = xxxx, abs.tol = 1e-18, rel.tol = 1e-18)$val
              temp <- lai(emp_cov1) 
              
              temploc[rr, ss] <- ifelse(tt == 0 & h[1] == 0 & h[2] == 0, sqrt(var[i] * var[j])*rho, sqrt(var[i] * var[j])*rho*4*pi*temp*beta^2/gamma(nu3))
            }
          }
          SS[[tt + 1]] <- temploc
        }
        S2 <- toeplitz_mat(SS)
        S[temp2, temp1] <- S2
        S[temp1, temp2] <- t(S2)
      }
    }
  }
  return(S)
}

lmc_cov <- function(theta, wind, max_time_lag, q, new_locations = locations, meters = T, nug_eff){
  
  nu <- theta[1:2]
  beta <- theta[3:4]
  var <- theta[5:6]
  
  if(nug_eff == T){
    nug <- theta[7:8]
  }else{
    nug <- c(0, 0)
  }
  
  alpha <- matrix(c(theta[7], theta[8], theta[9], theta[10]), ncol=2, byrow=T)
  
  S <- list()
  for(i in 1:q){
    
    if(meters == T){
      w <- matrix(wind, ncol = 2, byrow = T)/1000
      loc <- coords <- locations/1000
    }else{
      w <- matrix(wind, ncol = 2, byrow = T)
      loc <- coords <- locations
    }
    
    if (max_time_lag == 0){
      loc <- coords
    } else {
      for (tt in 1:max_time_lag){
        temploc <- matrix(, ncol=2, nrow=nrow(coords))
        for(rr in 1:nrow(coords)){
          temploc[rr,] <- c(coords[rr,1] - tt*w[i,1], coords[rr,2] - tt*w[i,2])
        }
        loc <- rbind(loc, temploc)
      }
    }
    dist0 <- spDists(loc, longlat = F)
  
    SS <- ifelse(dist0 != 0, var[i]*(dist0/beta[i])^nu[i] * besselK(dist0/beta[i], nu[i])/(2^(nu[i] - 1)*gamma(nu[i])), var[i] + nug[i])
    
    S[[i]] <- SS
  }
  S1 <- rbind(cbind(alpha[1,1]^2*S[[1]] + alpha[1,2]^2*S[[2]], alpha[1,1]*alpha[2,1]*S[[1]] + alpha[1,2]*alpha[2,2]*S[[2]]),
              cbind(alpha[1,1]*alpha[2,1]*S[[1]] + alpha[1,2]*alpha[2,2]*S[[2]], alpha[2,1]^2*S[[1]] + alpha[2,2]^2*S[[2]]))
  
  return(S1)
}

lmc_random_cov <- function(theta, wind, wind_var, max_time_lag, q, new_locations, meters = T, nug_eff){
  
  nu <- theta[1:2]
  beta <- theta[3:4]
  var <- theta[5:6]
  
  if(nug_eff == T){
    nug <- theta[7:8]
  }else{
    nug <- c(0, 0)
  }
  
  alpha <- matrix(c(theta[7], theta[8], theta[9], theta[10]), ncol=2, byrow=T)
  
  if(meters == T){
    w <- matrix(wind, ncol = 2, byrow = T)/1000
    loc <- coords <- locations/1000
    sigma <- wind_var
    sigma[[1]] <- sigma[[1]]/1000
    sigma[[2]] <- sigma[[2]]/1000
  }else{
    w <- matrix(wind, ncol = 2, byrow = T)
    loc <- coords <- locations
    sigma <- wind_var
  }
  
  S <- list()
  temploc <- denom <- list()
  
  for(i in 1:q){
    
    temploc[[1]] <- loc <- spDists(coords, longlat = F)
    denom[[1]] <- matrix(1, ncol = ncol(temploc[[1]]), nrow = nrow(temploc[[1]]))
    
    if (max_time_lag == 0){
      loc <- spDists(coords, longlat = F)
    } else {
      for (tt in 1:max_time_lag){
        
        temploc.temp <- matrix(, ncol=nrow(coords), nrow=nrow(coords))
        for(rr in 1:nrow(coords)){
          for(ss in 1:nrow(coords)){
            temploc.temp[rr,ss] <- sqrt((coords[rr,] - coords[ss,] - tt*w[i,])%*%solve(diag(2) + sigma[[i]])%*%matrix((coords[rr,] - coords[ss,] - tt*w[i,]), ncol=1))
          }
        }
        
        temploc[[tt + 1]] <- temploc.temp
        denom[[tt + 1]] <- sqrt(det(diag(2) + tt^2*sigma[[i]]))
      }
    }
    dist0 <- toeplitz_mat(temploc)
    denom.fin <- toeplitz_mat(denom)
    
    SS <- ifelse(dist0 != 0, var[i]*(dist0/beta[i])^nu[i] * besselK(dist0/beta[i], nu[i])/(2^(nu[i]-1)*gamma(nu[i])), var[i] + nug[i])
    
    S[[i]] <- SS/denom.fin
  }
  S1 <- rbind(cbind(alpha[1,1]^2*S[[1]] + alpha[1,2]^2*S[[2]], alpha[1,1]*alpha[2,1]*S[[1]] + alpha[1,2]*alpha[2,2]*S[[2]]),
              cbind(alpha[1,1]*alpha[2,1]*S[[1]] + alpha[1,2]*alpha[2,2]*S[[2]], alpha[2,1]^2*S[[1]] + alpha[2,2]^2*S[[2]]))
  
  return(S1)
}

matern_allard <- function(theta, max_time_lag, q, new_locations = locations, meters = T, nug_eff){
  
   loc <- coords <- new_locations
   tloc <- list()
  
  if (max_time_lag == 0){
    loc <- loc
    tloc.temp <- matrix(0, ncol = nrow(coords), nrow = nrow(coords))
    tloc[[1]] <- tloc.temp
  } else {
    tloc.temp <- matrix(0, ncol = nrow(coords), nrow = nrow(coords))
    tloc[[1]] <- tloc.temp
    for (tt in 1:max_time_lag){
      loc <- rbind(loc, coords)
      tloc.temp <- matrix(tt, ncol = nrow(coords), nrow = nrow(coords))
      tloc[[tt + 1]] <- tloc.temp
    }
  }
  
  if(meters == T){
    dist1 <- spDists(loc, longlat = F)/1000
  }else{
    dist1 <- spDists(loc, longlat = F)
  }
  
  nu <- theta[1:2]
  beta <- theta[3]
  var <- theta[4:5]
  rho <- theta[6]
  
  if(nug_eff == T){
    nug <- theta[7:8]
  }else{
    nug <- c(0, 0)
    alpha <- theta[7]
    a <- 1
    b <- theta[8]
  }
  
  dist2 <- toeplitz_mat(tloc)
  dist0 <- dist1/(alpha*dist2^(2*a) + 1)^(b/2)
  
  S=matrix(NA,  q*dim(dist0)[1], q*dim(dist0)[1])
  
  for(i in 1:q){
    for(j in 1:i){
      
      temp=(i-1)*dim(dist0)[1]+1:dim(dist0)[1]
      temp1=(j-1)*dim(dist0)[1]+1:dim(dist0)[1]
      
      if(i == j){
        
        temp2=ifelse(dist0 != 0, var[i]*(dist0/beta)^nu[i] * besselK(dist0/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i]))/((alpha*(dist2)^(2*a)+1)),(var[i]+nug[i])/(alpha*(dist2)^(2*a)+1))
        #temp2=ifelse(dist0!=0,var[i]*(dist0/beta[i])^nu[i] * besselK(dist0/beta[i],nu[i])/(2^(nu[i]-1)*gamma(nu[i])),var[i]+nug[i])
        
        S[temp,temp1]=temp2
        
      }
      
      if(i != j){
        
        nu1 <- nu[i]
        nu2 <- nu[j]
        nu3 <- (nu1+nu2)/2
        
        #rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        
        lai=ifelse(dist0 != 0, (dist0/beta)^nu3 * besselK(dist0/beta, nu3)/(2^(nu3 - 1)*gamma(nu3))*sqrt(var[i] * var[j])*rho/((alpha*(dist2)^(2*a) + 1)), sqrt(var[i] * var[j])*rho/((alpha*(dist2)^(2*a)+1)))
        S[temp,temp1] <- lai
        S[temp1,temp] <- t(lai)
        
      }
    }
  }
  return(S)
}

cov_lagrangian <- function(theta, wind, max_time_lag, q, new_locations = locations, meters = T, nug_eff){
  
  w <- wind
  
  loc <- coords <- new_locations
  
  if (max_time_lag == 0){
    loc <- loc
  } else {
    for (tt in 1:max_time_lag){
      temploc <- matrix(, ncol=2, nrow=nrow(coords))
      for(rr in 1:nrow(coords1)){
        temploc[rr,] <- c(coords[rr,1] - tt*w[1], coords[rr,2] - tt*w[2])
      }
      loc <- rbind(loc, temploc)
    }
  }
  
  if(meters == T){
    dist0 <- spDists(loc, longlat = F)/1000
  }else{
    dist0 <- spDists(loc, longlat = F)
  }
  
  mu <- theta[1:2]
  nu <- theta[3]
  scale <- theta[4]
  beta <- c(1,1)
  
  S <- matrix(NA,  q*dim(dist0)[1], q*dim(dist0)[1])
  
  for(i in 1:q){
    for(j in 1:i){
      
      temp <- (i - 1)*dim(dist0)[1] + 1:dim(dist0)[1]
      temp1 <- (j - 1)*dim(dist0)[1] + 1:dim(dist0)[1]
      
      if(i == j){
        
        temp2 <- pmax((1 - dist0/scale), 0)^(nu + mu[i])
        S[temp, temp1] <- temp2
        
      }
      
      if(i != j){
        
        mu1 <- mu[i]
        mu2 <- mu[j]
        mu3 <- (mu1 + mu2)/2
        
        beta3 <- (gamma(1 + mu3)/gamma(1 + nu + mu3))*sqrt((gamma(1 + nu + mu1)*gamma(1 + nu + mu2))/(gamma(1 + mu1)*gamma(1 + mu2)))
        
        lai <- beta3*pmax((1 - dist0/scale), 0)^(nu + mu3)
        S[temp, temp1] <- lai
        S[temp1, temp] <- t(lai)
        
      }
    }
  }
  return(S)
}


#---------NONSTATIONARY---------#

matern_cov_regular_grid <-function(theta,wind,time){
  
  w <- wind
  
  t=time
  q=2
  
  # create a spatial autocorrelation signature
  # coordinate list
  
  loc <- coords <- sim_grid_locations
  
  n <- nrow(sim_grid_locations)
  
  if (t==1){
    loc <- loc
  } else {
    for (tt in 1:(t-1)){
      temploc <- matrix(,ncol=2,nrow=nrow(coords))
      for(rr in 1:nrow(coords)){
        temploc[rr,] <- c(coords[rr,1]-tt*w[1],coords[rr,2]-tt*w[2])
      }
      loc <- rbind(loc, temploc)
    }
  }
  
  locations <- loc
  
  theta2 <- function (n,beta0,beta1,beta2,beta3,beta4) {
    theta3 <- beta0 + beta1*(locations[,1] - .5) + beta2*(locations[,2]-.5) + 
      beta3*(locations[,1] - .5)^2 + beta4*(locations[,2] - .5)^2
    theta3 <- matrix(theta3,nrow=nrow(locations),ncol=1)
    return(theta3)
  }
  
  #log.lam1.1<-theta2(n,-3,1,1,-6,-7)
  #log.lam1.2<-theta2(n,-5,1,1,6,-4)
  #logit.phi.1<-theta2(n,0,1,-2,0,1)
  
  #log.lam2.1<-theta2(n,-1.65,0.5,0.5,0,0)
  #log.lam2.2<-theta2(n,-2.8,-1,2,0,-7)
  #logit.phi.2<-theta2(n,-3,-1,2,0,-1)
  
  #log.lam2.1<-theta2(n,-3,-1,-1,-6,-7)
  #log.lam2.2<-theta2(n,-5,-1,-1,6,-4)
  #logit.phi.2<-theta2(n,0,-1,-2,0,1)
  
  log.lam1.1<-theta2(n,-3,1,1,-6,-7)
  log.lam1.2<-theta2(n,-5,1,1,2,-12)
  logit.phi.1<-theta2(n,0,1,-2,0,1)
  
  log.lam2.1<-theta2(n,-3,-1,-1,-6,-7)
  log.lam2.2<-theta2(n,-5,-1,-1,2,-12)
  logit.phi.2<-theta2(n,0,-1,-2,0,1)
  
  KERNEL_LIST <- list()
  
  kernel.local <- array(0, dim = c(2, 2, nrow(locations)))
  for(i in 1:nrow(locations)){
    lam1 <- exp(log.lam1.1[i,])
    lam2 <- exp(log.lam1.2[i,])
    phi <- (pi/2)*exp(logit.phi.1[i,])/(1+exp(logit.phi.1[i,]))
    Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
    Dmat <- diag(c(lam1, lam2))
    Sigma <- Pmat %*% Dmat %*% t(Pmat)
    kernel.local[, ,i] <-  Sigma
  }
  
  KERNEL_LIST[[1]] <- kernel.local
  
  for(i in 1:nrow(locations)){
    lam1 <- exp(log.lam2.1[i,])
    lam2 <- exp(log.lam2.2[i,])
    phi <- (pi/2)*exp(logit.phi.2[i,])/(1+exp(logit.phi.2[i,]))
    Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
    Dmat <- diag(c(lam1, lam2))
    Sigma <- Pmat %*% Dmat %*% t(Pmat)
    kernel.local[, ,i] <-  Sigma
  }
  
  KERNEL_LIST[[2]] <- kernel.local
  
  ##Calculate Matern form Nonstationary Covariance function 
  FIN_Sigma.mat <- list()
  dist0 <- list()
  
  for(KK in 1:2){
    Sigma.mat <- matrix(rep(NA, (n*t)^2), nrow = n*t)
    Q.mat <- matrix(rep(NA, (n*t)^2), nrow = n*t)
    Inv_ij <- matrix(rep(NA,4),2,2)
    
    for (i in 1:nrow(locations)) {
      #Sigma.mat[i, i] <- 1
      #Q.mat[i, i] <- 0
      Kernel_i <- KERNEL_LIST[[KK]][, , i]
      det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
      for (j in 1:nrow(locations)) {
        Kernel_j <- KERNEL_LIST[[KK]][, , j]
        det_j <- Kernel_j[1,1] * Kernel_j[2,2] - Kernel_j[1,2] * Kernel_j[2,1]
        Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
        Inv_ij[1,1] <- Kernel_ij[2,2] 
        Inv_ij[2,2] <- Kernel_ij[1,1] 
        Inv_ij[2,1] <- - Kernel_ij[2,1] 
        Inv_ij[1,2] <- - Kernel_ij[1,2] 
        det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - Kernel_ij[1,2] * Kernel_ij[2,1]
        x <- c(locations[i,1] - locations[j,1], locations[i,2] - locations[j,2])
        Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij)
        Q.mat[i, j] <- sqrt(t(x) %*% Inv_ij %*% x/det_ij)
        #Sigma.mat[j, i] <- Sigma.mat[i, j]
        #Q.mat[j, i] <- Q.mat[i, j]
      }
    }
    FIN_Sigma.mat[[KK]] <- Sigma.mat
    dist0[[KK]] <- Q.mat
  }
  
  for (i in 1:nrow(locations)) {
    #Sigma.mat[i, i] <- 1
    #Q.mat[i, i] <- 0
    Kernel_i <- KERNEL_LIST[[1]][, , i]
    det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
    for (j in 1:nrow(locations)) {
      Kernel_j <- KERNEL_LIST[[2]][, , j]
      det_j <- Kernel_j[1,1] * Kernel_j[2,2] - Kernel_j[1,2] * Kernel_j[2,1]
      Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
      Inv_ij[1,1] <- Kernel_ij[2,2] 
      Inv_ij[2,2] <- Kernel_ij[1,1] 
      Inv_ij[2,1] <- - Kernel_ij[2,1] 
      Inv_ij[1,2] <- - Kernel_ij[1,2] 
      det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - Kernel_ij[1,2] * Kernel_ij[2,1]
      x <- c(locations[i,1] - locations[j,1], locations[i,2] - locations[j,2])
      Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij)
      Q.mat[i, j] <- sqrt(t(x) %*% Inv_ij %*% x/det_ij)
      #Sigma.mat[j, i] <- Sigma.mat[i, j]
      #Q.mat[j, i] <- Q.mat[i, j]
    }
  }
  FIN_Sigma.mat[[3]] <- Sigma.mat
  dist0[[3]] <- Q.mat
  
  for (i in 1:nrow(locations)) {
    #Sigma.mat[i, i] <- 1
    #Q.mat[i, i] <- 0
    Kernel_i <- KERNEL_LIST[[2]][, , i]
    det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
    for (j in 1:nrow(locations)) {
      Kernel_j <- KERNEL_LIST[[1]][, , j]
      det_j <- Kernel_j[1,1] * Kernel_j[2,2] - Kernel_j[1,2] * Kernel_j[2,1]
      Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
      Inv_ij[1,1] <- Kernel_ij[2,2] 
      Inv_ij[2,2] <- Kernel_ij[1,1] 
      Inv_ij[2,1] <- - Kernel_ij[2,1] 
      Inv_ij[1,2] <- - Kernel_ij[1,2] 
      det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - Kernel_ij[1,2] * Kernel_ij[2,1]
      x <- c(locations[i,1] - locations[j,1], locations[i,2] - locations[j,2])
      Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij)
      Q.mat[i, j] <- sqrt(t(x) %*% Inv_ij %*% x/det_ij)
      #Sigma.mat[j, i] <- Sigma.mat[i, j]
      #Q.mat[j, i] <- Q.mat[i, j]
    }
  }
  FIN_Sigma.mat[[4]] <- Sigma.mat
  dist0[[4]] <- Q.mat
  
  nu=theta[1:2]
  beta=theta[3]
  rot=theta[4]
  Beta=matrix(0,q,q)
  diag(Beta)=1
  Beta[2,1]=rot
  Beta[1,2]=rot
  
  S=matrix(NA,  q*dim(dist0[[1]])[1], q*dim(dist0[[1]])[1])
  
  for(i in 1:q){
    for(j in 1:q){
      
      temp=(i-1)*dim(dist0[[1]])[1]+1:dim(dist0[[1]])[1]
      temp1=(j-1)*dim(dist0[[1]])[1]+1:dim(dist0[[1]])[1]
      
      if(i==j){
        
        temp2=ifelse(dist0[[i]]!=0,FIN_Sigma.mat[[i]]*(dist0[[i]]/beta)^nu[i] * besselK(dist0[[i]]/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i])),FIN_Sigma.mat[[i]])
        S[temp,temp1]=temp2
        
      }
      
      if(i!=j & i<j){
        
        nu1=nu[i]
        nu2=nu[j]
        nu3=(nu[i]+nu[j])/2
        
        rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        #rho <- Beta[i,j]*gamma(nu3)/sqrt(gamma(nu1)*gamma(nu2))
        
        lai=ifelse(dist0[[3]]!=0 ,(dist0[[3]]/beta)^nu3 * besselK(dist0[[3]]/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho,sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho)
        S[temp,temp1]=lai
        #S[temp1,temp]=t(lai)
      }
      
      if(i!=j & i>j){
        
        nu1=nu[i]
        nu2=nu[j]
        nu3=(nu[i]+nu[j])/2
        
        rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        #rho <- Beta[i,j]*gamma(nu3)/sqrt(gamma(nu1)*gamma(nu2))
        
        lai=ifelse(dist0[[4]]!=0 ,(dist0[[4]]/beta)^nu3 * besselK(dist0[[4]]/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho,sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho)
        S[temp,temp1]=lai
        #S[temp1,temp]=t(lai)
      }
    }
  }
  
  return(S)
}

matern_cov_regular_grid_v4 <-function(theta,wind,time){
  
  w <- wind
  
  t=time
  q=2
  
  # create a spatial autocorrelation signature
  # coordinate list
  
  loc <- coords <- sim_grid_locations
  
  n <- nrow(sim_grid_locations)
  
  if (t==1){
    loc <- loc
  } else {
    for (tt in 1:(t-1)){
      temploc <- matrix(,ncol=2,nrow=nrow(coords))
      for(rr in 1:nrow(coords)){
        temploc[rr,] <- c(coords[rr,1]-tt*w[1],coords[rr,2]-tt*w[2])
      }
      loc <- rbind(loc, temploc)
    }
  }
  
  locations <- loc
  
  theta2 <- function (n,beta0,beta1,beta2,beta3,beta4) {
    theta3 <- beta0 + beta1*(locations[,1] - .5) + beta2*(locations[,2]-.5) + 
      beta3*(locations[,1] - .5)^2 + beta4*(locations[,2] - .5)^2
    theta3 <- matrix(theta3,nrow=nrow(locations),ncol=1)
    return(theta3)
  }
  
  #log.lam1.1<-theta2(n,-3,1,1,-6,-7)
  #log.lam1.2<-theta2(n,-5,1,1,6,-4)
  #logit.phi.1<-theta2(n,0,1,-2,0,1)
  
  #log.lam2.1<-theta2(n,-1.65,0.5,0.5,0,0)
  #log.lam2.2<-theta2(n,-2.8,-1,2,0,-7)
  #logit.phi.2<-theta2(n,-3,-1,2,0,-1)
  
  #log.lam2.1<-theta2(n,-3,-1,-1,-6,-7)
  #log.lam2.2<-theta2(n,-5,-1,-1,6,-4)
  #logit.phi.2<-theta2(n,0,-1,-2,0,1)
  
  log.lam1.1<-theta2(n,-3,1,1,-6,-7)
  log.lam1.2<-theta2(n,-5,1,1,2,-12)
  logit.phi.1<-theta2(n,0,1,-2,0,1)
  
  log.lam2.1<-theta2(n,-3,-1,-1,-6,-7)
  log.lam2.2<-theta2(n,-5,-1,-1,2,-12)
  logit.phi.2<-theta2(n,0,-1,-2,0,1)
  
  KERNEL_LIST <- list()
  
  kernel.local <- array(0, dim = c(2, 2, nrow(locations)))
  for(i in 1:nrow(locations)){
    lam1 <- exp(log.lam1.1[i,])
    lam2 <- exp(log.lam1.2[i,])
    phi <- (pi/2)*exp(logit.phi.1[i,])/(1+exp(logit.phi.1[i,]))
    Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
    Dmat <- diag(c(lam1, lam2))
    Sigma <- Pmat %*% Dmat %*% t(Pmat)
    kernel.local[, ,i] <-  Sigma
  }
  
  KERNEL_LIST[[1]] <- kernel.local
  
  for(i in 1:nrow(locations)){
    lam1 <- exp(log.lam2.1[i,])
    lam2 <- exp(log.lam2.2[i,])
    phi <- (pi/2)*exp(logit.phi.2[i,])/(1+exp(logit.phi.2[i,]))
    Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
    Dmat <- diag(c(lam1, lam2))
    Sigma <- Pmat %*% Dmat %*% t(Pmat)
    kernel.local[, ,i] <-  Sigma
  }
  
  KERNEL_LIST[[2]] <- kernel.local
  
  ##Calculate Matern form Nonstationary Covariance function 
  FIN_Sigma.mat <- list()
  dist0 <- list()
  
  for(KK in 1:2){
    Sigma.mat <- matrix(rep(NA, (n*t)^2), nrow = n*t)
    Q.mat <- matrix(rep(NA, (n*t)^2), nrow = n*t)
    Inv_ij <- matrix(rep(NA,4),2,2)
    
    for (i in 1:nrow(locations)) {
      #Sigma.mat[i, i] <- 1
      #Q.mat[i, i] <- 0
      Kernel_i <- KERNEL_LIST[[KK]][, , i]
      det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
      for (j in 1:nrow(locations)) {
        Kernel_j <- KERNEL_LIST[[KK]][, , j]
        det_j <- Kernel_j[1,1] * Kernel_j[2,2] - Kernel_j[1,2] * Kernel_j[2,1]
        Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
        Inv_ij[1,1] <- Kernel_ij[2,2] 
        Inv_ij[2,2] <- Kernel_ij[1,1] 
        Inv_ij[2,1] <- - Kernel_ij[2,1] 
        Inv_ij[1,2] <- - Kernel_ij[1,2] 
        det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - Kernel_ij[1,2] * Kernel_ij[2,1]
        x <- c(locations[i,1] - locations[j,1], locations[i,2] - locations[j,2])
        Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij)
        Q.mat[i, j] <- sqrt(t(x) %*% Inv_ij %*% x/det_ij)
        #Sigma.mat[j, i] <- Sigma.mat[i, j]
        #Q.mat[j, i] <- Q.mat[i, j]
      }
    }
    FIN_Sigma.mat[[KK]] <- Sigma.mat
    dist0[[KK]] <- Q.mat
  }
  
  for (i in 1:nrow(locations)) {
    #Sigma.mat[i, i] <- 1
    #Q.mat[i, i] <- 0
    Kernel_i <- KERNEL_LIST[[1]][, , i]
    det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
    for (j in 1:nrow(locations)) {
      Kernel_j <- KERNEL_LIST[[2]][, , j]
      det_j <- Kernel_j[1,1] * Kernel_j[2,2] - Kernel_j[1,2] * Kernel_j[2,1]
      Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
      Inv_ij[1,1] <- Kernel_ij[2,2] 
      Inv_ij[2,2] <- Kernel_ij[1,1] 
      Inv_ij[2,1] <- - Kernel_ij[2,1] 
      Inv_ij[1,2] <- - Kernel_ij[1,2] 
      det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - Kernel_ij[1,2] * Kernel_ij[2,1]
      x <- c(locations[i,1] - locations[j,1], locations[i,2] - locations[j,2])
      Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij)
      Q.mat[i, j] <- sqrt(t(x) %*% Inv_ij %*% x/det_ij)
      #Sigma.mat[j, i] <- Sigma.mat[i, j]
      #Q.mat[j, i] <- Q.mat[i, j]
    }
  }
  FIN_Sigma.mat[[3]] <- Sigma.mat
  dist0[[3]] <- Q.mat
  
  for (i in 1:nrow(locations)) {
    #Sigma.mat[i, i] <- 1
    #Q.mat[i, i] <- 0
    Kernel_i <- KERNEL_LIST[[2]][, , i]
    det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
    for (j in 1:nrow(locations)) {
      Kernel_j <- KERNEL_LIST[[1]][, , j]
      det_j <- Kernel_j[1,1] * Kernel_j[2,2] - Kernel_j[1,2] * Kernel_j[2,1]
      Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
      Inv_ij[1,1] <- Kernel_ij[2,2] 
      Inv_ij[2,2] <- Kernel_ij[1,1] 
      Inv_ij[2,1] <- - Kernel_ij[2,1] 
      Inv_ij[1,2] <- - Kernel_ij[1,2] 
      det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - Kernel_ij[1,2] * Kernel_ij[2,1]
      x <- c(locations[i,1] - locations[j,1], locations[i,2] - locations[j,2])
      Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij)
      Q.mat[i, j] <- sqrt(t(x) %*% Inv_ij %*% x/det_ij)
      #Sigma.mat[j, i] <- Sigma.mat[i, j]
      #Q.mat[j, i] <- Q.mat[i, j]
    }
  }
  FIN_Sigma.mat[[4]] <- Sigma.mat
  dist0[[4]] <- Q.mat
  
  nu=theta[1:2]
  beta=theta[3]
  rot=theta[4]
  Beta=matrix(0,q,q)
  diag(Beta)=1
  Beta[2,1]=rot
  Beta[1,2]=rot
  
  S=matrix(NA,  q*dim(dist0[[1]])[1], q*dim(dist0[[1]])[1])
  
  for(i in 1:q){
    for(j in 1:q){
      
      temp=(i-1)*dim(dist0[[1]])[1]+1:dim(dist0[[1]])[1]
      temp1=(j-1)*dim(dist0[[1]])[1]+1:dim(dist0[[1]])[1]
      
      if(i==j){
        
        temp2=ifelse(dist0[[i]]!=0,FIN_Sigma.mat[[i]]*(dist0[[i]]/beta)^nu[i] * besselK(dist0[[i]]/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i])),FIN_Sigma.mat[[i]])
        S[temp,temp1]=temp2
        
      }
      
      if(i!=j & i<j){
        
        nu1=nu[i]
        nu2=nu[j]
        nu3=(nu[i]+nu[j])/2
        
        rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        #rho <- Beta[i,j]*gamma(nu3)/sqrt(gamma(nu1)*gamma(nu2))
        
        lai=ifelse(dist0[[3]]!=0 ,(dist0[[3]]/beta)^nu3 * besselK(dist0[[3]]/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho,sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho)
        S[temp,temp1]=lai
        #S[temp1,temp]=t(lai)
      }
      
      if(i!=j & i>j){
        
        nu1=nu[i]
        nu2=nu[j]
        nu3=(nu[i]+nu[j])/2
        
        rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        #rho <- Beta[i,j]*gamma(nu3)/sqrt(gamma(nu1)*gamma(nu2))
        
        lai=ifelse(dist0[[4]]!=0 ,(dist0[[4]]/beta)^nu3 * besselK(dist0[[4]]/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho,sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho)
        S[temp,temp1]=lai
        #S[temp1,temp]=t(lai)
      }
    }
  }
  
  return(S)
}
#--------------------------------------------------------------------------#

matern_cov_regular_grid_v2_for_estimation_sim_step1 <-function(theta,Q.mat1,Q.mat2,Q.mat3){
  q=2
  dist0 <- list()
  
  dist0[[1]] <- Q.mat1
  dist0[[2]] <- Q.mat2
  dist0[[3]] <- Q.mat3
  
  nu=theta[1:2]
  beta=theta[3]
  rot=theta[4]
  Beta=matrix(0,q,q)
  diag(Beta)=1
  Beta[2,1]=rot
  Beta[1,2]=rot
  
  var=theta[5:6]
  
  S=matrix(NA,  q*dim(dist0[[1]])[1], q*dim(dist0[[1]])[1])
  
  for(i in 1:q){
    for(j in 1:i){
      
      temp=(i-1)*dim(dist0[[1]])[1]+1:dim(dist0[[1]])[1]
      temp1=(j-1)*dim(dist0[[1]])[1]+1:dim(dist0[[1]])[1]
      
      if(i==j){
        
        temp2=ifelse(dist0[[i]]!=0,(dist0[[i]]/beta)^nu[i] * besselK(dist0[[i]]/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i])),var[i])
        #diag(temp2)=var[i]+nug[i]
        S[temp,temp1]=temp2
        
      }
      
      if(i !=j){
        
        nu1=nu[i]
        nu2=nu[j]
        nu3=(nu[i]+nu[j])/2
        
        rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        
        lai=ifelse(dist0[[3]]!=0 ,(dist0[[3]]/beta)^nu3 * besselK(dist0[[3]]/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(var[i]*var[j])*rho,sqrt(var[i]*var[j])*rho)
        S[temp,temp1]=lai
        S[temp1,temp]=t(lai)
      }
    }
  }
  
  return(S)
}

matern_cov_regular_grid_for_estimation_sim_v2 <-function(theta,Q.mat){
  q=2
  dist0 <- Q.mat
  
  nu=theta[1:2]
  beta=theta[3]
  rot=theta[4]
  Beta=matrix(0,q,q)
  diag(Beta)=1
  Beta[2,1]=rot
  Beta[1,2]=rot
  
  var=theta[5:6]
  
  S=matrix(NA,  q*dim(dist0)[1], q*dim(dist0)[1])
  
  for(i in 1:q){
    for(j in 1:i){
      
      temp=(i-1)*dim(dist0)[1]+1:dim(dist0)[1]
      temp1=(j-1)*dim(dist0)[1]+1:dim(dist0)[1]
      
      if(i==j){
        
        temp2=ifelse(dist0!=0,(dist0/beta)^nu[i] * besselK(dist0/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i])),var[i])
        #diag(temp2)=var[i]+nug[i]
        S[temp,temp1]=temp2
        
      }
      
      if(i !=j){
        
        nu1=nu[i]
        nu2=nu[j]
        nu3=(nu[i]+nu[j])/2
        
        rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        
        lai=ifelse(dist0!=0 ,(dist0/beta)^nu3 * besselK(dist0/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(var[i]*var[j])*rho,sqrt(var[i]*var[j])*rho)
        S[temp,temp1]=lai
        S[temp1,temp]=t(lai)
      }
    }
  }
  
  return(S)
}
