


source("./pkg-config.R")



DISTRIBUTED = T
TESTING_REFERENCE = F
PLOT = F
PLOT_MANUSCRIPT = T
SIMULATION1 = T
SIMULATION2 = F



WIND <- WIND_MU <- rep(0.0501, 2)
WIND_VAR <- matrix(0.0001 * diag(2), 2, 2)



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


n_sim <- 500

adj_mu <- c(0, 0, 0, 0)
adj_sig <- c(1, 10, 100, 1000)
adj_alpha <- c(0, 0.1, 0.5, 1)


locs_sub_index <- which(sim_grid_locations[, 1] >= 0.25 & sim_grid_locations[, 1] <= 0.75 & sim_grid_locations[, 2] >= 0.25 & sim_grid_locations[, 2] <= 0.75)
locs_sub_length <- length(locs_sub_index)

DIFF_ARRAY_EMP <- DIFF_ARRAY_THEO <- array(, dim = c(locs_sub_length, TT - 1, length(adj_mu), 2))

for(MODEL in 1:2){

	for(m in 1:length(adj_mu)){

		cat("MODEL: ", MODEL, ";  m: ", m, "\n")

		cat('Simulating wind values...', '\n')

		set.seed(1234)

		if(SIMULATION1){
			WIND_SIMULATED <- matrix(mvrnorm(n_sim, mu = WIND_MU + adj_mu[m], Sigma = matrix(WIND_VAR * adj_sig[m], ncol = 2, nrow = 2)), ncol = 2, byrow = T)
		}else if(SIMULATION2){
			WIND_SIMULATED <- rmsn(n = n_sim, xi = WIND_MU, WIND_VAR, alpha = c(0, 0) + adj_alpha[m])
		}



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
				clusterExport(cl, c("PARAMETER_NONSTAT", "TT", "WIND_SIMULATED"), envir = environment())


				cat('Distributing computations over', number_of_cores_to_use, 'cores...', '\n')

				output <- foreach(i=1:nrow(WIND_SIMULATED), .combine='+', .packages = "Rcpp", .noexport = "SPATIALLY_VARYING_PARAMETERS_FOR_FITTING_PARALLEL") %dopar% {
					
					COVARIANCE <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(1, 0.23, 1, WIND_SIMULATED[i, ], 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT, FITTING = T, PARALLEL = T)

					return(c(COVARIANCE))
				}

				stopCluster(cl)

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

				number_of_cores_to_use = cores[1]-1 # not to overload the computer
				cat('Registering', number_of_cores_to_use, 'cores...', '\n')

				cl <- makeCluster(number_of_cores_to_use) 
				registerDoParallel(cl)

				clusterEvalQ(cl, source("./pkg-config.R"))
				clusterExport(cl, c("PARAMETER_DEFORMATION", "TT", "WIND_SIMULATED"), envir = environment())


				cat('Distributing computations over', number_of_cores_to_use, 'cores...', '\n')

				output <- foreach(i=1:nrow(WIND_SIMULATED), .combine='+', .packages = "Rcpp", .noexport = "DEFORMATION_FOR_FITTING_PARALLEL") %dopar% {
					
					COVARIANCE <- MATERN_UNI_DEFORMATION(PARAMETER = c(1, 0.23, 1, WIND_SIMULATED[i, ], 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_DEFORMATION = PARAMETER_DEFORMATION, FITTING = T, PARALLEL = T)

					return(c(COVARIANCE))
				}

				stopCluster(cl)


			}


		}

		theocov <- matrix(output, n * TT, n * TT) / nrow(WIND_SIMULATED) 

		cat('Generating realizations...', '\n')

		set.seed(1234)
		r <- rmvn(1000, rep(0, n * TT), theocov, ncores = number_of_cores_to_use)


		empcov <- cov(r)


		############################################################################################

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

		fit1 <- optim(par = init, fn = NEGLOGLIK_NONPARAMETRIC, control = list(trace = 5, maxit = 3000)) #

		p <- fit1$par
		EST_WIND_MU <- p[1:2]

		wind_var_chol <- matrix(c(p[3], p[4], 0, p[5]), ncol = 2, byrow = T)
		EST_WIND_VAR <- t(wind_var_chol) %*% wind_var_chol


		############################################################################################

		set.seed(1234)
		est_wind_vals <- mvrnorm(n_sim, EST_WIND_MU, EST_WIND_VAR)


		count <- 1
		diff_cov_emp <- diff_cov_theo <- matrix(, ncol = 4, nrow = locs_sub_length)
		for(l2 in locs_sub_index){
			diff_cov_emp_temp <- diff_cov_theo_temp <- NULL
			for(t2 in 1:(TT - 1)){
				cov_purely_time_emp <- empcov[l2, l2 + n * t2]
				cov_purely_time_theo <- theocov[l2, l2 + n * t2]
				cov_purely_space_emp_temp <- cov_purely_space_theo_temp <- 0
				for(k in 1:n_sim){
					wind <- est_wind_vals[k, ]
					new_loc <- matrix(c(sim_grid_locations[l2, 1] - wind[1] * t2, sim_grid_locations[l2, 2] - wind[2] * t2), ncol = 2)
					find_new_loc_index <- which.min(distR_C(sim_grid_locations, new_loc))[1]

					cov_purely_space_emp_temp <- cov_purely_space_emp_temp + empcov[l2, find_new_loc_index]
					cov_purely_space_theo_temp <- cov_purely_space_theo_temp + theocov[l2, find_new_loc_index]
				}
				cov_purely_space_emp <- cov_purely_space_emp_temp / n_sim
				cov_purely_space_theo <- cov_purely_space_theo_temp / n_sim
				diff_cov_emp_temp <- c(diff_cov_emp_temp, cov_purely_time_emp - cov_purely_space_emp)
				diff_cov_theo_temp <- c(diff_cov_theo_temp, cov_purely_time_theo - cov_purely_space_theo)
			}
			diff_cov_emp[count, ] <- diff_cov_emp_temp
			diff_cov_theo[count, ] <- diff_cov_theo_temp
			count <- count + 1
		}

		DIFF_ARRAY_EMP[, , m, MODEL] <- diff_cov_emp
		DIFF_ARRAY_THEO[, , m, MODEL] <- diff_cov_theo
	}
}



if(TESTING_REFERENCE){


	cat('Computing spatio-temporal empirical covariance...', '\n')

	empcov <- cov(r)



	cat('Simulating velocity from estimated distribution of advection velocity...', '\n')

	EST_WIND_MU <- WIND_MU #+ 0.01
	EST_WIND_VAR <- WIND_VAR #+ 0.01

	W_vals <- c()

	for(rep in 1:500){
		set.seed(rep)
		est_wind_vals <- mvrnorm(n_sim, EST_WIND_MU, EST_WIND_VAR)



		locs_sub_index <- which(sim_grid_locations[, 1] >= 0.25 & sim_grid_locations[, 1] <= 0.75 & sim_grid_locations[, 2] >= 0.25 & sim_grid_locations[, 2] <= 0.75)
		locs_sub_length <- length(locs_sub_index)

		lag_max <- TT - 1
		f_theo <- f_ref <- f_emp <- matrix(, ncol = length(locs_sub_index), nrow = lag_max)
		#f <- matrix(, ncol = length(locs_sub_index) * (TT - lag_max))


		ct <- 0
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

			diff_cov_emp <- cov_purely_time_emp - cov_purely_space_emp

			#D <- eigen(cov_purely_space_emp)$val
			#D[which(D < 0)] <- 0
			#U <- eigen(cov_purely_space_emp)$vectors
			#cov_purely_time_emp <- U %*% diag(D) %*% t(U)

			set.seed((rep - 1) * 1000 + l2)
			r1_time <- mvrnorm(200000, rep(0, TT), cov_purely_space_emp)


			#cat('Computing purely temporal empirical covariance...', '\n')

			empcov_time <- cov(r1_time[1:100000, ])


			#cat('Computing difference between purely spatial and purely temporal empirical covariance...', '\n')

			diff_cov_theo <- empcov_time - cov_purely_space_emp

			empcov_time <- cov(r1_time[100000 + 1:100000, ])

			diff_cov_ref <- empcov_time - cov_purely_space_emp

			ct <- ct + 1
			f_theo[, ct] <- diff_cov_theo[1, 2:TT]					
			f_ref[, ct] <- diff_cov_ref[1, 2:TT]					
			f_emp[, ct] <- diff_cov_emp[1, 2:TT]					

		}

		ff <- fbplot(f_emp, plot=F);
		nonOut <- setdiff(1:ct, ff$outpoint);
		f_emp <- f_emp[, nonOut];

		ff <- fbplot(f_theo, plot=F);
		nonOut <- setdiff(1:ct, ff$outpoint);
		f_theo <- f_theo[, nonOut];

		ff <- fbplot(f_ref, plot=F);
		nonOut <- setdiff(1:ct, ff$outpoint);
		f_ref <- f_ref[, nonOut];


		cat('Computing ranks...', '\n')



		n1=dim(f_theo)[2];
		m=dim(f_emp)[2];
		r=dim(f_ref)[2];
		order=integer(n1 + m);
		for(i in 1:m){
			sample <- cbind(f_ref, f_emp[, i]);
			result <- fbplot(sample, plot = F, method = "Both");
			order[i] <- sum(result$depth[1:r] <= result$depth[r + 1])
		}

		for(i in 1:n1){
			sample <- cbind(f_ref, f_theo[, i]);
			result <- fbplot(sample, plot = F, method = "Both");
			order[i + m] <- sum(result$depth[1:r] <= result$depth[r + 1])
		}
		rk <- (rank(order) - 1) / (n1 + m - 1);
		W <- mean(rk[1:m]);

		cat("rep: ", rep, '   Rank value = ', W, '\n')
		W_vals[rep] <- W
	}

	write.table(W_vals, file = paste(root, "Data/nonstationary-taylor-hypothesis/5-bootstrap-W-values-model-", MODEL, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)



}

if(PLOT_MANUSCRIPT){

	if(SIMULATION1){
		pdf(file = paste(root, 'Figures/5-test-functions-simulation1.pdf', sep = ''), width = 25, height = 10)
	}else if(SIMULATION2){
		pdf(file = paste(root, 'Figures/5-test-functions-simulation2.pdf', sep = ''), width = 25, height = 10)
	}

	split.screen( rbind(c(0.08,0.98,0.1,0.95), c(0.98,0.99,0.1,0.95)))
	split.screen( figs = c( 2, 4 ), screen = 1 )

	hr_label <- c('i', 'ii', 'iii', 'iv')
	mod_label <- c('A', 'B')

	for(model in 1:2){

		for(m in 1:length(adj_mu)){
		
			screen((model - 1) * 4 + 2 + m)
			par(mai=c(0.2,0.2,0.2,0.2))
			
			fbplot(t(DIFF_ARRAY_EMP[, , m, model]), method='MBD', ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', ylim = range(DIFF_ARRAY_EMP[, , , model]))
			abline(h = 0, col = 3, lty = 2, lwd = 5)

			if(m == 1){
				mtext(expression(hat(f)), side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
				text(-0.275, 0, mod_label[model], col = 'blue', xpd = NA, cex = 4, font = 2)
				axis(2, cex.axis = 2)
			}

			if(model == 1){
				mtext(hr_label[m], side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
			}else{
				mtext(expression(t^'*'), side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
				axis(1, at = seq(1, 4, by = 1), cex.axis = 2, mgp = c(1, 1.5, 0))
			}
		}
	}				

	close.screen( all=TRUE)

	dev.off()

}


if(PLOT){



	pdf(file = '/home/salvanmo/Desktop/studious-potato/Figures/5-bootstrap-histogram.pdf', width = 20, height = 15)
	hist(W_vals)
	dev.off()

	pdf(file = '/home/salvanmo/Desktop/studious-potato/Figures/5-skew-normal-advection.pdf', width = 20, height = 15)
	plot(WIND_SIMULATED)
	dev.off()

	pdf(file = '/home/salvanmo/Desktop/studious-potato/Figures/5-3-test-functions.pdf', width = 20, height = 15)

	par(mfrow = c(3, 1))
	fbplot(f_emp, method='MBD', ylab = '', xlab = '')
	fbplot(f_ref, method='MBD', ylab = '', xlab = '')
	fbplot(f_theo, method='MBD', ylab = '', xlab = '')
		
	dev.off()


}
