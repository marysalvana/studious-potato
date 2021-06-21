
workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/' 		else 		directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))

sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))

N <- 70
n <- N^2
TT <- 1
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

ind <- 14 
set.seed(ind)
PARAMETER_NONSTAT <- runif(15, -3, 3)

set.seed(3)
PARAMETER_DEFORMATION <- c(runif(2, -8, 8), runif(1, 0, 8), runif(2, -1, 1))

cat('Computing covariances...', '\n')

cov1 <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(1, 0.23, 1, 0.1001, 0.1001, 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT)

set.seed(1)
r1 <- rmvn(500, rep(0, ncol(cov1)), cov1, ncores = 25)

emp_cov <- cov(r1)

reference_locations <- c(1, ceiling(n / 2) + ceiling(N / 2))

Z <- r1[1, 1:n] - mean(r1[1, 1:n])

Z_dot_product <- Z %*% t(Z)

kernel <- MATERN_UNI_STATIONARY(PARAMETER = c(1, 0.001, 0.5), LOCATION = sim_grid_locations)

vals <- c()
for(k in 1:n){

	cat(k, '\n')

	loc1 <- reference_locations[2]
	loc2 <- k
	
	KERNEL <- kernel[loc1, ] %*% t(kernel[loc2, ])
	emp_cov_est <- sum(KERNEL * Z_dot_product) / sum(KERNEL)

	vals[k] <- emp_cov_est
}

pdf(file = paste(root, 'Figures/4-empirical-covariance-estimation.pdf', sep = ''), width = 15, height = 15)

par(mfrow = c(2, 2))

#image.plot(matrix(emp_cov[reference_locations[1], 1:n], N, N))
#image.plot(matrix(cov1[reference_locations[1], 1:n], N, N))

#image.plot(matrix(emp_cov[reference_locations[2], 1:n], N, N))
image.plot(matrix(cov1[reference_locations[2], 1:n], N, N))

image.plot(matrix(vals, N, N))
#image.plot(matrix(r1[3, 1:n], N, N))

dev.off()

