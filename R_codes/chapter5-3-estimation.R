
source("./pkg-config.R")



DISTRIBUTED = T
MODEL = 1
ESTIMATION = T
PLOT = F
PLOT_RESULTS = F

args = commandArgs(trailingOnly = TRUE)
REP = as.numeric(args[1])
velocity_mu_config = as.numeric(args[2])
velocity_var_config = as.numeric(args[3])

mu_k <- c(0, 0.0501)
var_k <- c(0.01, 0.1, 1)

WIND <- WIND_MU <- rep(mu_k[velocity_mu_config], 2)
WIND_VAR <- matrix(var_k[velocity_var_config] * diag(2), 2, 2)


number_of_cores_to_use <- 39


N <- 20
TT <- 5
grid_x <- seq(from = 0, to = 1, by = 0.05)
grid_y <- seq(from = 0, to = 1, by = 0.05)
sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()
n <- nrow(sim_grid_locations)

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

n_sim <- 100

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

		#number_of_cores_to_use = cores[1]-1 # not to overload the computer
		cat('Registering', number_of_cores_to_use, 'cores...', '\n')

		cl <- makeCluster(number_of_cores_to_use) 
		registerDoParallel(cl)

		clusterEvalQ(cl, source("./pkg-config.R"))
		clusterExport(cl, c("PARAMETER_NONSTAT", "TT"), envir = environment())

		cat('Simulating wind values...', '\n')

		set.seed(1234)
		wind_vals <- mvrnorm(n_sim, WIND_MU, WIND_VAR)

		cat('Distributing computations over', number_of_cores_to_use, 'cores...', '\n')

		output <- foreach(i=1:nrow(wind_vals), .combine='+', .packages = "Rcpp", .noexport = "SPATIALLY_VARYING_PARAMETERS_FOR_FITTING_PARALLEL") %dopar% {
			
			COVARIANCE <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(1, 0.23, 1, wind_vals[i, ], 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT, FITTING = T, PARALLEL = T)

			return(c(COVARIANCE))
		}

		stopCluster(cl)
	

		cov1 <- matrix(output, n * TT, n * TT) / nrow(wind_vals) 


	}


	cat('Generating realizations...', '\n')

	set.seed(1)
	r1 <- rmvn(1000, rep(0, n * TT), cov1, ncores = number_of_cores_to_use)


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



if(ESTIMATION){


	locs_sub_index <- which(sim_grid_locations[, 1] >= 0.25 & sim_grid_locations[, 1] <= 0.75 & sim_grid_locations[, 2] >= 0.25 & sim_grid_locations[, 2] <= 0.75)
	locs_sub_length <- length(locs_sub_index)


	NEGLOGLIK_NONPARAMETRIC <- function(p){

		wind_mu <- p[1:2]

		wind_var_chol <- matrix(c(p[3], p[4], 0, p[5]), ncol = 2, byrow = T)
        	wind_var <- t(wind_var_chol) %*% wind_var_chol


		set.seed(1234)
		est_wind_vals <- matrix(mvrnorm(n_sim, mu = wind_mu, Sigma = wind_var), ncol = 2, byrow = T)


		diff_cov_emp <- 0

		for(l2 in locs_sub_index){
			cov_purely_space_emp <- cov_purely_time_emp <- matrix(, TT, TT)
			for(t1 in 1:TT){
				for(t2 in t1:TT){
					cov_purely_time_emp[t1, t2] <- empcov[l2 + n * (t1 - 1), l2 + n * (t2 - 1)]
					cov_purely_space_emp_temp <- 0
					for(k in 1:n_sim){
						wind <- est_wind_vals[k, ]
						new_loc <- matrix(c(sim_grid_locations[l2, 1] - wind[1] * (t2 - 1), sim_grid_locations[l2, 2] - wind[2] * (t2 - 1)), ncol = 2)
						find_new_loc_index <- which.min(distR_C(cbind(sim_grid_locations[, 1] - wind[1] * (t1 - 1), sim_grid_locations[, 2] - wind[2] * (t1 - 1)), new_loc))[1]

						cov_purely_space_emp_temp <- cov_purely_space_emp_temp + empcov[l2 + n * (t1 - 1), find_new_loc_index + n * (t1 - 1)]
					}
					cov_purely_space_emp[t1, t2] <- cov_purely_space_emp_temp / n_sim
					if(t2 != t1){
						cov_purely_space_emp[t2, t1] <- cov_purely_space_emp[t1, t2]
						cov_purely_time_emp[t2, t1] <- cov_purely_time_emp[t1, t2]
					}
				}
			}

			diff_cov_emp <- diff_cov_emp + sum((cov_purely_time_emp - cov_purely_space_emp)^2)


		}


		return(diff_cov_emp)
	}


	init <- c(0.2, 0.2, 0.1, 0, 0.1)



	set.seed(REP)
	r1 <- rmvn(1000, rep(0, n * TT), cov1, ncores = number_of_cores_to_use)
	empcov <- cov(r1)

	fit1 <- optim(par = init, fn = NEGLOGLIK_NONPARAMETRIC, control = list(trace = 5, maxit = 3000)) #

	for(aa in 1:5){

		cat("SIMULATION REP: ", REP, " --- ESTIMATION REP: ", aa, '\n')
		fit1 <- optim(par = fit1$par, fn = NEGLOGLIK_NONPARAMETRIC, control = list(trace = 5, maxit = 3000)) #
	}

	p <- fit1$par
	wind_mu <- p[1:2]

	wind_var_chol <- matrix(c(p[3], p[4], 0, p[5]), ncol = 2, byrow = T)
	wind_var <- t(wind_var_chol) %*% wind_var_chol


	est_params <- matrix(c(REP, velocity_mu_config, velocity_var_config, wind_mu, wind_var[1, 1], wind_var[2, 2], wind_var[1, 2]), nrow = 1)



	write.table(est_params, file = paste(root, "Results/5-boxplots-estimated-advection-parameters", sep = ""), sep = " ", append = TRUE, row.names = FALSE, col.names = FALSE)

}


if(PLOT){

	cat('PLOTTING COVARIANCE . . . ', '\n')

	Ref_loc <- c(2 * N + 5, ceiling(n / 2) + ceiling(N / 2))



	if(MODEL == 1){



		###########   SPATIALLY VARYING PARAMETERS MODEL   ###########



		cov_example <- read.table(paste(root, 'Data/nonstationary-taylor-hypothesis/cov-example-1-velocity_mu_config_2_velocity_var_config_1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

		plot_univariate_nonstationary_covariance_heatmap(covariance = cov_example, covariance_est = Sig.CLS[reference_locations, ], locations = sim_grid_locations, reference_locations = Ref_loc, '5-nonstationary-cov1-heatmap-nonparametric-estimation.pdf')



	}else if(MODEL == 2){



	###########    DEFORMATION MODEL   ###########



		cov_example <- read.table(paste(root, 'Data/nonstationary-taylor-hypothesis/cov-example-2-velocity_mu_config_2_velocity_var_config_1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

		plot_univariate_nonstationary_covariance_heatmap(covariance = cov_example, covariance_est = Sig.CLS[reference_locations, ], locations = sim_grid_locations, reference_locations = Ref_loc, '5-nonstationary-cov2-heatmap-nonparametric-estimation.pdf')



	}

}


if(PLOT_RESULTS){

	params_temp <- read.table(paste(root, "Results/5-boxplots-estimated-advection-parameters", sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	params <- params_temp[, -(1:3)]

	
	# width of each boxplot is 0.8
	x0s <- 1:ncol(params) - 0.4
	x1s <- 1:ncol(params) + 0.4
	# these are the y-coordinates for the horizontal lines
	# that you need to set to the desired values.
	y0s <- matrix(c(WIND_MU, WIND_VAR[1, 1], WIND_VAR[2, 2], WIND_VAR[1, 2]), nrow = 1, ncol = ncol(params), byrow = T)


	pdf(file = paste(root, 'Figures/5-boxplots-estimated-advection-parameters.pdf', sep = ''), width = 25, height = 10)

	split.screen( rbind(c(0.08,0.98,0.1,0.95), c(0.98,0.99,0.1,0.95)))
	split.screen( figs = c( 2, 3 ), screen = 1 )

	hr_label <- c('Weak', 'Moderate', 'Strong')
	mod_label <- c('A', 'B')

	for(model in 1:2){

		for(m in 1:3){
		
			screen((model - 1) * 3 + 2 + m)
			par(mai=c(0.2,0.2,0.2,0.2))
			
			boxplot(params[, 1:5], xaxt = 'n', yaxt = 'n')
			segments(x0 = x0s, x1 = x1s, y0 = y0s[1:5], col = "red", lwd = 2)

			if(m == 1){
				#mtext(expression(hat(f)), side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
				text(-0.3, 4, mod_label[model], col = 'blue', xpd = NA, cex = 4, font = 2)
				#mtext(mod_label[model], side = 2, line = 8, adj = 0.5, cex = 4, font = 2, col = 'blue')
				axis(2, cex.axis = 2)
			}

			if(model == 1){
				mtext(hr_label[m], side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
			}else{
				axis(1, at = seq(1, 5, by = 1), cex.axis = 2, mgp = c(1, 1.5, 0))
			}
		}
	}				

	close.screen( all=TRUE)

	dev.off()



}


