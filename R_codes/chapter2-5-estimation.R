


source("./pkg-config.R")



## INPUT: N x T matrix of log PM2.5 concentrations, where N is the number of spatial locations and T is the number of temporal locations 
## INPUT: N x 2 matrix of locations containing longitude and latitude

## OUTPUT: textfile of large training and testing datasets (measurements, spatial locations, temporal locations)



AREA <- 'SAUDI'

DATA <- data_format(aggregate = 1, area = AREA)

locs <- as.matrix(DATA[[2]])
locs <- cbind((locs[, 1] - mean(locs[, 1])) / sd(locs[, 1]), (locs[, 2] - mean(locs[, 2])) / sd(locs[, 2]))


####################################################################################################


TT <- 2

WIND <- rep(0.1, 2)

NN <- 50
grid_x <- seq(from = -3, to = 3, length.out = NN)
grid_y <- seq(from = -3, to = 3, length.out = NN)
X <- expand.grid(grid_x, grid_y) %>% as.matrix()

#X <- rbind(locs, X)
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

Xtarg <- locs

if (TT > 1){
	for (tt in 1:(TT - 1)){
		temp_locs <- cbind(locs[, 1] - tt * WIND[1], locs[, 2] - tt * WIND[2])
		Xtarg <- rbind(Xtarg, temp_locs)
	}
}

htarg <- as.matrix(dist(rbind(X, Xtarg),diag=TRUE,upper=TRUE))
htarg <- htarg[1:nn, (nn + 1):(nn + nrow(Xtarg))]
sigma <- htarg^2 * log(htarg)
sigma[htarg == 0] <- 0

Z_rand_sample <- matrix(c(DATA[[1]][1, ], DATA[[1]][2, ]), nrow = 1)
#Z_rand_sample <- matrix(DATA[[1]][seq(1, 35040, by = 8760), ], nrow = 1)



####################################################################################################



NEGLOGLIK_DEFORM <- function(p){

	theta <- exp(p[1:3])

	beta1 <- p[3 + 1:length(jWarp)]
	beta2 <- p[3 + length(jWarp) + 1:length(jWarp)]
  
	parWarpsSum <- cbind(rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T)),
		       rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T)))

	Y <- Xtarg
	Y[, 1] <- Y[, 1] + c(parWarpsSum[, 1] %*% sigma)
	Y[, 2] <- Y[, 2] + c(parWarpsSum[, 2] %*% sigma)

	dist0 <- as.matrix(dist(Y, diag = TRUE, upper = TRUE))	

	Sigma <- theta[1] * Matern(dist0, range = theta[2], nu = theta[3])
	cholmat <- t(cholesky(Sigma, parallel = TRUE))
	z <- forwardsolve(cholmat, t(Z_rand_sample))
	logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
	out  <- 1/2 * logsig + 1/2 * sum(z^2)

	return(out)
}

jWarp = 1:5
init <- c(rep(0, 3), rep(0, 2 * length(jWarp)))
fit1 <- optim(par = init, fn = NEGLOGLIK_DEFORM, control = list(trace = 5, maxit = 3000)) #
#324
p <- fit1$par

theta <- exp(p[1:3])

beta1 <- p[3 + 1:length(jWarp)]
beta2 <- p[3 + length(jWarp) + 1:length(jWarp)]

parWarpsSum <- cbind(rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T)),
	       rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T)))

Y <- locs
Y[, 1] <- Y[, 1] + c(parWarpsSum[, 1] %*% sigma)
Y[, 2] <- Y[, 2] + c(parWarpsSum[, 2] %*% sigma)

pdf(file = paste(root, 'Figures/2-scratch-cov2-heatmap.pdf', sep = ''), width = 10, height = 10)

#plot(Y, pch = 3, cex = 2)
plot(cbind(c(parWarpsSum[, 1] %*% sigma), c(parWarpsSum[, 2] %*% sigma)), pch = 3, cex = 2)

dev.off()

####################################################################################################



NEGLOGLIK_DEFORM <- function(p){

	theta <- exp(p[1:3])

	beta1 <- p[3 + 1:length(jWarp)]
	beta2 <- p[3 + length(jWarp) + 1:length(jWarp)]
  
	parWarpsSum <- cbind(rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T)),
		       rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T)))

	Y1 <- Xtarg
	Y1[, 1] <- Y1[, 1] + c(parWarpsSum[, 1] %*% sigma)
	Y1[, 2] <- Y1[, 2] + c(parWarpsSum[, 2] %*% sigma)

	beta1 <- p[3 + 2 * length(jWarp) + 1:length(jWarp)]
	beta2 <- p[3 + 3 * length(jWarp) + 1:length(jWarp)]
  
	parWarpsSum <- cbind(rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T)),
		       rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T)))

	Y2 <- Xtarg
	Y2[, 1] <- Y2[, 1] + c(parWarpsSum[, 1] %*% sigma)
	Y2[, 2] <- Y2[, 2] + c(parWarpsSum[, 2] %*% sigma)

	Y <- rbind(Y1[1:550, ], Y2[550 + 1:550, ])

	dist0 <- as.matrix(dist(Y, diag = TRUE, upper = TRUE))	

	Sigma <- theta[1] * Matern(dist0, range = theta[2], nu = theta[3])
	cholmat <- t(cholesky(Sigma, parallel = TRUE))
	z <- forwardsolve(cholmat, t(Z_rand_sample))
	logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
	out  <- 1/2 * logsig + 1/2 * sum(z^2)

	return(out)
}

jWarp = 1:5
init <- c(rep(0, 3), rep(0, 4 * length(jWarp)))
fit1 <- optim(par = init, fn = NEGLOGLIK_DEFORM, control = list(trace = 5, maxit = 3000)) #

####################################################################################################


NEGLOGLIK_SPATIALLY_VARYING <- function(p){

	theta <- exp(p[1:3])

	beta1 <- p[3 + 1:length(jWarp)]
	beta2 <- p[3 + length(jWarp) + 1:length(jWarp)]
	beta3 <- p[3 + 2 * length(jWarp) + 1:length(jWarp)]
  
	parWarpsSum <- cbind(rowSums( g[,3+jWarp] * matrix(beta3, ncol=length(beta3), nrow=nrow(X), byrow=T)),
		       rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T)),
			rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T)))

	PARAMETER_NONSTAT <- t(sigma) %*% parWarpsSum

	Sigma <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(theta, 0.1001, 0.1001, 0.001, 0, 0, 0.001), LOCATION = locs, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT, FITTING = T)
	cholmat <- t(cholesky(Sigma[['covariance']], parallel = TRUE))
	z <- forwardsolve(cholmat, t(Z_rand_sample))
	logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
	out  <- 1/2 * logsig + 1/2 * sum(z^2)

	return(out)
}

jWarp = 1:10
init <- c(rep(0, 3), rep(0, 3 * length(jWarp)))
fit1 <- optim(par = init, fn = NEGLOGLIK_SPATIALLY_VARYING, control = list(trace = 5, maxit = 3000)) #

p <- fit1$par
len_old <- (length(p) - 3) / 3
jWarp = 1:15
new_init <- c(p[1:3], p[3 + 1:len_old], rep(0, max(jWarp) - len_old), p[3 + len_old + 1:len_old], rep(0, max(jWarp) - len_old), p[3 + 2 * len_old + 1:len_old], rep(0, max(jWarp) - len_old))
fit1 <- optim(par = new_init, fn = NEGLOGLIK_SPATIALLY_VARYING, control = list(trace = 5, maxit = 3000)) #

p <- fit1$par
theta <- exp(p[1:3])

beta1 <- p[3 + 1:length(jWarp)]
beta2 <- p[3 + length(jWarp) + 1:length(jWarp)]
beta3 <- p[3 + 2 * length(jWarp) + 1:length(jWarp)]

parWarpsSum <- cbind(rowSums( g[,3+jWarp] * matrix(beta3, ncol=length(beta3), nrow=nrow(X), byrow=T)),
	       rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T)),
		rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T)))

PARAMETER_NONSTAT <- t(sigma) %*% parWarpsSum

pdf(file = paste(root, 'Figures/2-scratch-cov1-heatmap.pdf', sep = ''), width = 30, height = 11)

par(mfrow = c(1, 3))
for(i in 1:3){
	quilt.plot(locs[, 1], locs[, 2], PARAMETER_NONSTAT[, i], nx = 25, ny = 25)
}
dev.off()



####################################################################################################

####################################################################################################

Z_rand_sample <- matrix(c(DATA[[1]][1, ], DATA[[1]][2, ], DATA[[1]][8760 + 1, ], DATA[[1]][8760 + 2, ]), nrow = 1)

NEGLOGLIK_MULTI_SPATIALLY_VARYING <- function(p){

	theta <- exp(p[1:5])

	beta1 <- p[5 + 1:length(jWarp)]
	beta2 <- p[5 + length(jWarp) + 1:length(jWarp)]
	beta3 <- p[5 + 2 * length(jWarp) + 1:length(jWarp)]
  
	parWarpsSum <- cbind(rowSums( g[,3+jWarp] * matrix(beta3, ncol=length(beta3), nrow=nrow(X), byrow=T)),
		       rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T)),
			rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T)))

	PARAMETER_NONSTAT <- t(sigma) %*% parWarpsSum

	beta4 <- p[5 + 3 * length(jWarp) + 1:length(jWarp)]
	beta5 <- p[5 + 3 * length(jWarp) + length(jWarp) + 1:length(jWarp)]
	beta6 <- p[5 + 3 * length(jWarp) + 2 * length(jWarp) + 1:length(jWarp)]
  
	parWarpsSum <- cbind(rowSums( g[,3+jWarp] * matrix(beta6, ncol=length(beta6), nrow=nrow(X), byrow=T)),
		       rowSums( g[,3+jWarp] * matrix(beta5, ncol=length(beta5), nrow=nrow(X), byrow=T)),
			rowSums( g[,3+jWarp] * matrix(beta4, ncol=length(beta4), nrow=nrow(X), byrow=T)))

	PARAMETER_NONSTAT2 <- t(sigma) %*% parWarpsSum

	Sigma <- MULTIVARIATE_MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(theta, 0.5, WIND, 0.001, 0, 0, 0.001), LOCATION = locs, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT, PARAMETER_NONSTAT2 = PARAMETER_NONSTAT2, FITTING = T)
	cholmat <- t(cholesky(Sigma[['covariance']], parallel = TRUE))
	z <- forwardsolve(cholmat, t(Z_rand_sample))
	logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
	out  <- 1/2 * logsig + 1/2 * sum(z^2)

	return(out)
}

jWarp = 1:5
init <- c(rep(0, 5), rep(0, 2 * 3 * length(jWarp)))
fit1 <- optim(par = init, fn = NEGLOGLIK_MULTI_SPATIALLY_VARYING, control = list(trace = 5, maxit = 3000)) #

