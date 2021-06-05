## INPUT: N x T matrix of log PM2.5 concentrations, where N is the number of spatial locations and T is the number of temporal locations 
## INPUT: N x 2 matrix of locations containing longitude and latitude

## OUTPUT: textfile of large training and testing datasets (measurements, spatial locations, temporal locations)

directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))

sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))


AREA <- 'SAUDI'

DATA <- data_format(aggregate = 1, area = AREA)

locs <- as.matrix(DATA[[2]])
locs <- cbind((locs[, 1] - mean(locs[, 1])) / sd(locs[, 1]), (locs[, 2] - mean(locs[, 2])) / sd(locs[, 2]))


####################################################################################################

plot_realdata_for_checking_stationarity(data_list = DATA, file_name = '2-application-saudi-data-scratch.jpg', start_hr = 8760 + 1, saudi = T)
plot_realdata_for_manuscript(data_list = DATA, file_name = '2-application-saudi-data-scratch.jpg', start_hr = 8760 + 1, saudi = T)

####################################################################################################



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

htarg <- as.matrix(dist(rbind(X, locs),diag=TRUE,upper=TRUE))
htarg <- htarg[1:nn, (nn + 1):(nn + nrow(locs))]
sigma <- htarg^2 * log(htarg)
sigma[htarg == 0] <- 0

Z_rand_sample <- matrix(DATA[[1]][8760 * 3 + 1:10, ], nrow = 1)



####################################################################################################



NEGLOGLIK_DEFORM <- function(p){

	theta <- exp(p[1:3])

	beta1 <- p[3 + 1:length(jWarp)]
	beta2 <- p[3 + length(jWarp) + 1:length(jWarp)]
  
	parWarpsSum <- cbind(rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T)),
		       rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T)))

	Y <- locs
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

jWarp = 1:15
init <- c(rep(0, 3), rep(0, 2 * length(jWarp)))
fit1 <- optim(par = init, fn = NEGLOGLIK_DEFORM, control = list(trace = 5, maxit = 3000)) #

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


NEGLOGLIK_SPATIALLY_VARYING <- function(p){

	theta <- exp(p[1:3])

	beta1 <- p[3 + 1:length(jWarp)]
	beta2 <- p[3 + length(jWarp) + 1:length(jWarp)]
	beta3 <- p[3 + 2 * length(jWarp) + 1:length(jWarp)]
  
	parWarpsSum <- cbind(rowSums( g[,3+jWarp] * matrix(beta3, ncol=length(beta3), nrow=nrow(X), byrow=T)),
		       rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T)),
			rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T)))

	PARAMETER_NONSTAT <- t(sigma) %*% parWarpsSum

	Sigma <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(theta, 0.1001, 0.1001, 0.001, 0, 0, 0.001), LOCATION = locs, TIME = 1, PARAMETER_NONSTAT = PARAMETER_NONSTAT, FITTING = T)
	cholmat <- t(cholesky(Sigma[['covariance']], parallel = TRUE))
	z <- forwardsolve(cholmat, t(Z_rand_sample))
	logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
	out  <- 1/2 * logsig + 1/2 * sum(z^2)

	return(out)
}

jWarp = 1:15
init <- c(rep(0, 3), rep(0, 3 * length(jWarp)))
fit1 <- optim(par = init, fn = NEGLOGLIK_SPATIALLY_VARYING, control = list(trace = 5, maxit = 10000)) #

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

