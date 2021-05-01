
directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))

N <- 11
n <- N^2
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

cat('Computing covariances...', '\n')

cov1 <- MATERN_UNI_STATIONARY(PARAMETER = c(1, 0.23, 1), LOCATION = locations)

cat('Generating realizations...', '\n')

set.seed(1)
r1 <- rmvn(1, rep(0, ncol(cov1)), cov1, ncores = 25)

cat('Remove the mean from the observations...', '\n')

Z <- r1 - mean(r1)

cat('Compute the empirical covariance between all possible location pairs...', '\n')

emp_covariance_temp <- t(Z) %*% Z

xlag <- ylag <- emp_vals <- NULL  

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

nbins <- 10
dist.bin <- seq(0, max(sqrt(xlag^2 + ylag^2)), length = nbins)

EMP_VEC <- EMP_COUNT <- rep(0, nbins)

for(gg in 1:nrow(binned)){

	IND <- findInterval(binned[gg, 3], dist.bin)
	EMP_VEC[IND] <- EMP_VEC[IND] + binned[gg, 4] 
	EMP_COUNT[IND] <- EMP_COUNT[IND] + 1
}

EMP_COV <- EMP_VEC / EMP_COUNT

cat('Plot the empirical covariance...', '\n')

pdf(file = paste(root, 'Figures/empirical-covariance.pdf', sep = ''), width = 11, height = 6.5)

plot(dist.bin, EMP_VEC, pch = 20)

dev.off()
