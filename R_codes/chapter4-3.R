
workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/'          else            directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))

sourceCpp(file=paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))
sourceCpp(file=paste(root, "R_codes/Functions/distR.cpp",sep=''))

N <- 20
n <- N^2
TT <- 5
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

set.seed(14)
PARAMETER_NONSTAT <- runif(15, -3, 3)

set.seed(3)
PARAMETER_DEFORMATION <- c(runif(2, -8, 8), runif(1, 0, 8), runif(2, -1, 1))

cat('Computing covariances...', '\n')

wind_mu <- c(0.05263158, 0.05263158)
wind_sigma <- c(0.01, 0, 0, 0.01)
n_sim = 500


cov1 <- MATERN_UNI_DEFORMATION(PARAMETER = c(1, 0.23, 1, wind_mu, wind_sigma), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_DEFORMATION = PARAMETER_DEFORMATION, N_SIM = n_sim)

cov2 <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(1, 0.23, 1, wind_mu, wind_sigma), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT, N_SIM = n_sim)

set.seed(1234)
r1 <- mvrnorm(500, rep(0, ncol(cov1)), cov1)

set.seed(1234)
r2 <- mvrnorm(500, rep(0, ncol(cov2)), cov2)

EMPCOV <- array(, dim = c(dim(cov1), 2))

EMPCOV[, , 1] <- cov(r1)
EMPCOV[, , 2] <- cov(r2)

adj_mu <- c(0, 0, 0, 0.05)
adj_sig <- c(1, 0.1, 5, 1)

DIFF_ARRAY <- array(, dim = c(dim(diff_cov), 4, 2))

for(model in 1:2){

	empcov <- EMPCOV[, , model]
	for(m in 1:length(adj_mu)){

		set.seed(1234)
		WIND_SIMULATED <- matrix(mvrnorm(n_sim, mu = wind_mu + adj_mu[m], Sigma = matrix(wind_sigma * adj_sig[m], ncol = 2, nrow = 2)), ncol = 2, byrow = T)

		locs_sub_index <- which(sim_grid_locations[, 1] >= 0.25 & sim_grid_locations[, 1] <= 0.75 & sim_grid_locations[, 2] >= 0.25 & sim_grid_locations[, 2] <= 0.75)

		locs_sub_length <- length(locs_sub_index)

		count <- 1
		diff_cov <- matrix(, ncol = 4, nrow = locs_sub_length)
		for(l2 in locs_sub_index){
			diff_cov_temp <- NULL
			for(t2 in 1:(TT - 1)){
				cov_purely_time <- empcov[l2, l2 + n * t2]
				cov_purely_space_temp <- 0
				for(k in 1:n_sim){
					wind <- WIND_SIMULATED[k, ]
					new_loc <- matrix(c(sim_grid_locations[l2, 1] - wind[1] * t2, sim_grid_locations[l2, 2] - wind[2] * t2), ncol = 2)
					find_new_loc_index <- which.min(distR_C(sim_grid_locations, new_loc))[1]

					cov_purely_space_temp <- cov_purely_space_temp + empcov[l2, find_new_loc_index]
				}
				cov_purely_space <- cov_purely_space_temp / n_sim
				diff_cov_temp <- c(diff_cov_temp, cov_purely_time - cov_purely_space)
			}
			diff_cov[count, ] <- diff_cov_temp
			count <- count + 1
		}

		DIFF_ARRAY[, , m, model] <- diff_cov
	}
}

pdf(file = paste(root, 'Figures/4-taylors-hypothesis-example.pdf', sep = ''), width = 15, height = 8)

par(mfrow = c(2, 4))

for(m in 1:4){

	if(m == 1){
		fbplot(t(DIFF_ARRAY[, , m, 1]), method='MBD', ylab = '', xlab = '', xaxt = 'n', ylim = c(-0.5, 0.5), cex.axis = 2)
		mtext(expression(hat(f)), side = 2, cex = 1.5, line = 2.5)
	}else{
		fbplot(t(DIFF_ARRAY[, , m, 1]), method='MBD', ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', ylim = c(-0.5, 0.5))
	}
	abline(h = 0, col = "#39FF14", lty = 2, lwd = 1.5)
	axis(1, at = seq(1, 4, by = 1), cex.axis = 2)
	mtext(expression(t[2]), side = 1, cex = 1.5, line = 2.5)
}

dev.off()

pdf(file = paste(root, 'Figures/4-taylors-hypothesis-example.pdf', sep = ''), width = 25, height = 10)

split.screen( rbind(c(0.08,0.98,0.1,0.95), c(0.98,0.99,0.1,0.95)))
split.screen( figs = c( 2, 4 ), screen = 1 )

hr_label <- c('i', 'ii', 'iii', 'iv')
mod_label <- c('A', 'B')

for(model in 1:2){

	for(m in 1:(TT - 1)){
	
		screen((model - 1) * 4 + 2 + m)
		par(mai=c(0.2,0.2,0.2,0.2))
		
		fbplot(t(DIFF_ARRAY[, , m, model]), method='MBD', ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', ylim = c(-0.5, 0.5), col = "#FFAE42", barcol = "#0074B7")
		abline(h = 0, col = "red", lty = 2, lwd = 3)

		if(m == 1){
			mtext(expression(hat(f)), side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
			text(-0.275, 0, mod_label[model], col = 'blue', xpd = NA, cex = 4, font = 2)
			#mtext(mod_label[model], side = 2, line = 8, adj = 0.5, cex = 4, font = 2, col = 'blue')
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


