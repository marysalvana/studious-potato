
#The data should reflect the true unknown covariance model which generated it.
#Meaning: if the smoothness parameter is high, the data should show a smooth random field.
#Meaning: if the spatial correlation is weak, it should show as well.
#So, when we compute the sample covariance from the data, it should be a rough estimate of the true unknown covariance.

#Step 1: compute the sample covariance

#Step 2: Plot dist vs. covariance for easier analysis

#Step 3: Try to fit a curve on the computed sample covariance based on a theoretical covariance model.
#The parameters you used should be a good guess to the MLE parameters.

library(MASS)
library(dplyr)

N <- 40
n <- N^2
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

cat('Computing covariances...', '\n')

MATERN_UNI_STATIONARY <- function(PARAMETER, LOCATION){
  
  nu <- PARAMETER[3]
  beta <- PARAMETER[2]
  sigma <- PARAMETER[1]
  
  dist0 <- dist(x = LOCATION, diag = TRUE, upper = TRUE) %>% as.matrix()
  
  S <- ifelse(dist0 != 0, (dist0 / beta)^nu * besselK(dist0 / beta, nu) / (2^(nu - 1) * gamma(nu)) * sigma, sigma)
  
  return(S)
}

cov1 <- MATERN_UNI_STATIONARY(PARAMETER = c(1, 0.23, 1), LOCATION = locations)

cat('Generating realizations...', '\n')

set.seed(1234)
DATA <- mvrnorm(1, rep(0, ncol(cov1)), cov1)

cat('Remove the mean from the observations...', '\n')

Z <- matrix(DATA - mean(DATA), nrow = 1)

cat('Compute the empirical covariance between all possible location pairs...', '\n')

emp_covariance_temp <- t(Z) %*% Z

cat('Arrange the empirical covariance values in one dataframe/matrix...', '\n')

xlag <- ylag <- emp_vals <- NULL
for(loc in 1:nrow(locations)){
	xlag <- c(xlag, locations[, 1] - locations[loc, 1])
	ylag <- c(ylag, locations[, 2] - locations[loc, 2])
	emp_vals <- c(emp_vals, emp_covariance_temp[, loc])
}

empirical <- data.frame(xlag, ylag, sqrt(xlag^2 + ylag^2), emp_vals)
colnames(empirical) <- c('xlag', 'ylag', 'dist', 'covariance_value')

cat('Take the average of the empirical covariance values for each spatial lag...', '\n')

binned_temp <- empirical %>% group_by(dist) %>% summarize(avg1 = mean(covariance_value))
binned <- cbind(binned_temp$dist, binned_temp$avg1)

cat('Bin the spatial lags into a few representative spatial lags...', '\n')

nbins <- 20
dist.bin <- seq(0, max(sqrt(xlag^2 + ylag^2)), length = nbins)

EMP_VEC <- EMP_COUNT <- rep(0, nbins)

for(gg in 1:nrow(binned)){

	IND <- findInterval(binned[gg, 1], dist.bin)
	EMP_VEC[IND] <- EMP_VEC[IND] + binned[gg, 2] 
	EMP_COUNT[IND] <- EMP_COUNT[IND] + 1
}

EMP_COV <- EMP_VEC / EMP_COUNT

cat('Plot the empirical covariance...', '\n')

cat('Guess estimate of the theoretical covariance...', '\n')

sigma <- 1
beta <- 0.15
nu <- 1
dist0 <- seq(0, max(sqrt(xlag^2 + ylag^2)), length = 500)
theoretical <- ifelse(dist0 != 0, (dist0 / beta)^nu * besselK(dist0 / beta, nu) / (2^(nu - 1) * gamma(nu)) * sigma, sigma)

plot(dist.bin, EMP_COV, pch = 20, xlab = 'Distance', ylab = 'Empirical Covariance')
lines(dist0, theoretical, col = 'red', lwd = 2)

