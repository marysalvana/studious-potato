


source("./pkg-config.R")



AREA <- 'SAUDI'

DATA <- data_format(aggregate = 1, area = AREA)

locs <- as.matrix(DATA[[2]])

#normalize the location for unit square simulation later

locs <- cbind((locs[, 1] - min(locs[, 1])) / (max(locs[, 1]) - min(locs[, 1])), (locs[, 2] - min(locs[, 2])) / (max(locs[, 2]) - min(locs[, 2])))



NN <- 50
grid_x <- seq(from = -1, to = 2, length.out = NN)
grid_y <- seq(from = -1, to = 2, length.out = NN)
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

Xtarg <- locs



htarg <- as.matrix(dist(rbind(X, Xtarg),diag=TRUE,upper=TRUE))
htarg <- htarg[1:nn, (nn + 1):(nn + nrow(Xtarg))]
sigma <- htarg^2 * log(htarg)
sigma[htarg == 0] <- 0



Z_rand_sample <- matrix(DATA[[1]][1, ], nrow = 1)



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

p <- fit1$par

theta <- exp(p[1:3])

beta1 <- p[3 + 1:length(jWarp)]
beta2 <- p[3 + length(jWarp) + 1:length(jWarp)]

parWarpsSum <- cbind(rowSums( g[,3+jWarp] * matrix(beta2, ncol=length(beta2), nrow=nrow(X), byrow=T)),
	       rowSums( g[,3+jWarp] * matrix(beta1, ncol=length(beta1), nrow=nrow(X), byrow=T)))

Y <- locs
Y[, 1] <- Y[, 1] + c(parWarpsSum[, 1] %*% sigma)
Y[, 2] <- Y[, 2] + c(parWarpsSum[, 2] %*% sigma)

pdf(file = paste(root, 'Figures/extra-real-data-estimated-deformation.pdf', sep = ''), width = 10, height = 5)

par(mfrow = c(1, 2))
plot(locs, pch = 3, cex = 2)
plot(Y, pch = 3, cex = 2)

dev.off()

