workstation = T

if(workstation){
	directory <- '/home/salvanmo/Desktop/'
	root <- paste(directory, 'studious-potato/', sep = '')
	source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
	number_of_cores_to_use = 24

	model = 1
	velocity_mu_config = 2
	velocity_var_config = 1

}else{
	directory <- '/ibex/scratch/salvanmo/'
	root <- paste(directory, 'studious-potato/', sep = '')
	source(file = paste(root, "R_codes/Functions/load_packages-IBEX.R", sep = ''))
	number_of_cores_to_use = 39

	args <- commandArgs(trailingOnly = TRUE)

	model = as.numeric(args[1])
	velocity_mu_config = as.numeric(args[2])
	velocity_var_config = as.numeric(args[3])
}

source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))
sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2-IBEX.cpp",sep=''))


source(file = paste(root, "R_codes/Functions/cls.r", sep = ''))
sourceCpp(file = paste(root, "R_codes/Functions/distR.cpp",sep=''))



start_time <- Sys.time()



SIMULATE = F
PLOT = F
NONPARAMETRIC_ESTIMATION = T



AREA <- 'SAUDI'

DATA <- data_format(aggregate = 1, area = AREA)
locs <- as.matrix(DATA[[2]])
locs <- cbind((locs[, 1] - mean(locs[, 1])) / sd(locs[, 1]), (locs[, 2] - mean(locs[, 2])) / sd(locs[, 2]))



#######################################################################################

cat("model:", model, "velocity_mu_config", velocity_mu_config, "velocity_var_config", velocity_var_config, '\n')



mu_k <- c(0, 0.3001)
var_k <- c(0.001, 0.1, 1)

WIND <- WIND_MU <- rep(mu_k[velocity_mu_config], 2)
WIND_VAR <- matrix(var_k[velocity_var_config] * diag(2), 2, 2)


N <- 20
n <- N^2
TT <- 5
grid_x <- seq(from = min(locs[, 1]), to = max(locs[, 1]), length.out = N)
grid_y <- seq(from = min(locs[, 2]), to = max(locs[, 2]), length.out = N)
sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

reference_locations <- c(2 * N + 5, ceiling(n / 2) + ceiling(N / 2))


if(SIMULATE){

	distributed = T

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



		#######################################################################################



		cat('Computing covariances...', '\n')

		if(!distributed){

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

			output <- foreach(i=1:nrow(wind_vals), .combine='+', .packages = "Rcpp", .noexport = "SPATIALLY_VARYING_PARAMETERS_FOR_FITTING_PARALLEL") %dopar% {
				
				COVARIANCE <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(1, 0.23, 1, wind_vals[i, ], 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT, FITTING = T, PARALLEL = T)

				return(c(COVARIANCE))
			}

			stopCluster(cl)



			cov1 <- matrix(output, n * TT, n * TT) / nrow(wind_vals) 

			cat('Generating realizations...', '\n')

			set.seed(1)
			r1 <- rmvn(10, rep(0, n * TT), cov1, ncores = number_of_cores_to_use)

			cat('Saving the values...', '\n')

			write.table(cov1[reference_locations, ], file = paste(root, "Data/nonstationary-taylor-hypothesis/cov-example-1-velocity_mu_config_", velocity_mu_config, "_velocity_var_config_", velocity_var_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

		}

		write.table(r1[1:10, ], file = paste(root, "Data/nonstationary-taylor-hypothesis/realizations-example-1-velocity_mu_config_", velocity_mu_config, "_velocity_var_config_", velocity_var_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

	}else if(model == 2){

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
			r2 <- rmvn(10, rep(0, n * TT), cov2[["covariance"]], ncores = number_of_cores_to_use)

			cat('Saving the values...', '\n')

			write.table(cov2[["covariance"]][reference_locations, ], file = paste(root, "Data/nonstationary-taylor-hypothesis/cov-example-2-velocity_mu_config_", velocity_mu_config, "_velocity_var_config_", velocity_var_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

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

			write.table(cov2[reference_locations, ], file = paste(root, "Data/nonstationary-taylor-hypothesis/cov-example-2-velocity_mu_config_", velocity_mu_config, "_velocity_var_config_", velocity_var_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)
		}

		write.table(r2[1:10, ], file = paste(root, "Data/nonstationary-taylor-hypothesis/realizations-example-2-velocity_mu_config_", velocity_mu_config, "_velocity_var_config_", velocity_var_config, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

	}

	end_time <- Sys.time()

	cat('DONE.', '\n')
	cat("Textfiles are saved in ", paste(root, 'Data/nonstationary-taylor-hypothesis/', sep = ''), '\n')
	cat(end_time - start_time, 'seconds', '\n')

}

if(NONPARAMETRIC_ESTIMATION){

	locs1 <- sim_grid_locations

	dist1<-as.matrix(dist(locs1))
	A <- 1:nrow(locs1)

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
			temp.g <- distR_C(loc.jk, locs1)
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

	Sig.em<-cov1[1:n, 1:n]

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

}


if(PLOT){

	cat('PLOTTING COVARIANCE . . . ', '\n')

	Ref_loc <- c(2 * N + 5, ceiling(n / 2) + ceiling(N / 2))



	if(model == 1){



		###########   SPATIALLY VARYING PARAMETERS MODEL   ###########



		cov_example <- read.table(paste(root, 'Data/nonstationary-taylor-hypothesis/cov-example-1-velocity_mu_config_2_velocity_var_config_1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

		plot_univariate_nonstationary_covariance_heatmap(covariance = cov_example, covariance_est = Sig.CLS[reference_locations, ], locations = sim_grid_locations, reference_locations = Ref_loc, '5-nonstationary-cov1-heatmap-nonparametric-estimation.pdf')



	}else if(model == 2){



	###########    DEFORMATION MODEL   ###########



		cov_example <- read.table(paste(root, 'Data/nonstationary-taylor-hypothesis/cov-example-2-velocity_mu_config_2_velocity_var_config_1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

		plot_univariate_nonstationary_covariance_heatmap(covariance = cov_example, covariance_est = Sig.CLS[reference_locations, ], locations = sim_grid_locations, reference_locations = Ref_loc, '5-nonstationary-cov2-heatmap-nonparametric-estimation.pdf')



	}



}

