
workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/' 		else 		directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))

sourceCpp(file=paste(root,"R_codes/Functions/distR.cpp",sep=''))

N <- 10
n <- N^2
TT <- 10
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

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

NEGLOGLIK2_SCRATCH <- function(p){

	w <- p[4:5]
	Xtarg <- Xtarg_orig <- locs 

	if (TT > 1){
		for (tt in 1:(TT - 1)){
			temp_locs <- cbind(Xtarg_orig[, 1] - tt * w[1], Xtarg_orig[, 2] - tt * w[2])
			Xtarg <- rbind(Xtarg, temp_locs)
		}
	}

	Ytarg <- c(0.5, 0.5) + (Xtarg - c(0.5, 0.5)) * as.numeric(distR_C(Xtarg, matrix(c(0.5, 0.5), nrow = 1)))
	#Ytarg[, 1] <- (Ytarg[, 1] - mean(Ytarg[, 1])) / sd(Ytarg[, 1])
	#Ytarg[, 2] <- (Ytarg[, 2] - mean(Ytarg[, 2])) / sd(Ytarg[, 2])

	beta1 <- beta_est[1:length(jWarp)]
	beta2 <- beta_est[length(jWarp) + 1:length(jWarp)]
  
	parWarpsSum <- cbind(rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T)),
		       rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T)))

	htarg <- as.matrix(dist(rbind(X, Xtarg),diag=TRUE,upper=TRUE))
	htarg <- htarg[1:nn, (nn + 1):(nn + nrow(Xtarg))]
	sigma <- htarg^2 * log(htarg)
	sigma[htarg == 0] <- 0

	Y <- Xtarg
	Y[, 1] <- Y[, 1] + c(parWarpsSum[, 1] %*% sigma)
	Y[, 2] <- Y[, 2] + c(parWarpsSum[, 2] %*% sigma)

	beta1 <- p[1:length(jWarp)]
	beta2 <- p[length(jWarp) + 1:length(jWarp)]
	beta3 <- p[2 * length(jWarp) + 1:length(jWarp)]
  
	parWarpsSum <- cbind(rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T)),
		       rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T)), 
			rowSums( g[,3+jWarp] * matrix(beta3, ncol=length(beta3), nrow=nrow(X), byrow=T)))

	nonstat_est1 <- c(parWarpsSum[, 1] %*% sigma)
	nonstat_est2 <- c(parWarpsSum[, 2] %*% sigma)
	nonstat_est3 <- c(parWarpsSum[, 3] %*% sigma)

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
