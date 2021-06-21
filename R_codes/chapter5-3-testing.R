
workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/' 		else 		directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
library(future.apply)
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))

sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))
sourceCpp(file = paste(root, "R_codes/Functions/distR.cpp",sep=''))

N <- 10
n <- N^2
TT <- 10
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

ind <- 14 
set.seed(ind)
PARAMETER_NONSTAT <- runif(15, -3, 3)

set.seed(3)
PARAMETER_DEFORMATION <- c(runif(2, -8, 8), runif(1, 0, 8), runif(2, -1, 1))

cat('Computing covariances...', '\n')

#cov1 <- MATERN_UNI_DEFORMATION(PARAMETER = c(1, 0.23, 1, 0.3001, 0.3001, 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_DEFORMATION = PARAMETER_DEFORMATION)

loc2 <- sim_grid_locations
for(tt in 1:(TT - 1)){
	loc2 <- rbind(loc2, cbind(sim_grid_locations[, 1] - 0.2001 * tt, sim_grid_locations[, 2] - 0.2001 * tt))
}

cov1 <- MATERN_UNI_STATIONARY(PARAMETER = c(1, 0.23, 1), LOCATION = loc2)

set.seed(1)
r1 <- rmvn(500, rep(0, ncol(cov1)), cov1, ncores = 25)

Y <- matrix(r1[1, ], nrow = n, ncol = TT, byrow = F)

pdf(file = paste(root, 'Figures/4-simulated-asymmetric.pdf', sep = ''), width = 15, height = 15)

par(mfrow = c(2, 2))

for(tt in 1:4){
	par(mai=c(1.5, 1.5, 1.5, 1.5))
	par(pty="s") 
	image.plot(matrix(Y[, tt], N, N))
}

dev.off()

Y <- matrix(simulateSymmData(size.dim = N, p = TT, blocksize = 4, lambda = 0.5), nrow = n, ncol = TT, byrow = F)

