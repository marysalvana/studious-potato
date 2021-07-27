


source("./pkg-config.R")



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

wind_mu <- c(0.05, 0.05)
wind_sigma <- c(0.01, 0, 0, 0.01)
n_sim = 500

cov1 <- MATERN_UNI_DEFORMATION(PARAMETER = c(1, 0.23, 1, wind_mu, wind_sigma), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_DEFORMATION = PARAMETER_DEFORMATION, N_SIM = n_sim, DIST = "NORMAL")

cov2 <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(1, 0.23, 1, wind_mu, wind_sigma), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT, N_SIM = n_sim, DIST = "NORMAL")

set.seed(1234)
r1 <- mvrnorm(1000, rep(0, ncol(cov1[["covariance"]])), cov1[["covariance"]])

set.seed(1234)
r2 <- mvrnorm(1000, rep(0, ncol(cov2[["covariance"]])), cov2[["covariance"]])

EMPCOV <- THEOCOV <- array(, dim = c(dim(cov1[["covariance"]]), 2))

EMPCOV[, , 1] <- cov(r1)
EMPCOV[, , 2] <- cov(r2)

THEOCOV[, , 1] <- cov1[["covariance"]]
THEOCOV[, , 2] <- cov2[["covariance"]]

#adj_mu <- c(0, 0, 0, 0.05)
#adj_sig <- c(1, 0.1, 5, 1)

adj_mu <- c(0, 0, 0, 0)
adj_sig <- c(1, 1, 1, 1)

locs_sub_index <- which(sim_grid_locations[, 1] >= 0.25 & sim_grid_locations[, 1] <= 0.75 & sim_grid_locations[, 2] >= 0.25 & sim_grid_locations[, 2] <= 0.75)
locs_sub_length <- length(locs_sub_index)

DIFF_ARRAY_EMP <- DIFF_ARRAY_THEO <- array(, dim = c(locs_sub_length, TT - 1, length(adj_mu), 2))

for(model in 1:2){

	empcov <- EMPCOV[, , model]
	theocov <- THEOCOV[, , model]

	for(m in 1:length(adj_mu)){

		set.seed(1234)
		WIND_SIMULATED <- mvrnorm(n_sim, mu = wind_mu + adj_mu[m], Sigma = matrix(wind_sigma * adj_sig[m], ncol = 2, nrow = 2))

		count <- 1
		diff_cov_emp <- diff_cov_theo <- matrix(, ncol = 4, nrow = locs_sub_length)
		for(l2 in locs_sub_index){
			diff_cov_emp_temp <- diff_cov_theo_temp <- NULL
			for(t2 in 1:(TT - 1)){
				cov_purely_time_emp <- empcov[l2, l2 + n * t2]
				cov_purely_time_theo <- theocov[l2, l2 + n * t2]
				cov_purely_space_emp_temp <- cov_purely_space_theo_temp <- 0
				for(k in 1:n_sim){
					wind <- WIND_SIMULATED[k, ]
					new_loc <- matrix(c(sim_grid_locations[l2, 1] - wind[1] * t2, sim_grid_locations[l2, 2] - wind[2] * t2), ncol = 2)
					find_new_loc_index <- which.min(distR_C(sim_grid_locations, new_loc))[1]

					cov_purely_space_emp_temp <- cov_purely_space_emp_temp + empcov[l2, find_new_loc_index]
					cov_purely_space_theo_temp <- cov_purely_space_theo_temp + theocov[l2, find_new_loc_index]
				}
				cov_purely_space_emp <- cov_purely_space_emp_temp / n_sim
				cov_purely_space_theo <- cov_purely_space_theo_temp / n_sim
				diff_cov_emp_temp <- c(diff_cov_emp_temp, cov_purely_time_emp - cov_purely_space_emp)
				diff_cov_theo_temp <- c(diff_cov_theo_temp, cov_purely_time_theo - cov_purely_space_theo)
			}
			diff_cov_emp[count, ] <- diff_cov_emp_temp
			diff_cov_theo[count, ] <- diff_cov_theo_temp
			count <- count + 1
		}

		DIFF_ARRAY_EMP[, , m, model] <- diff_cov_emp
		DIFF_ARRAY_THEO[, , m, model] <- diff_cov_theo
	}
}

DIFF_ARRAY_THEO[, , 1, 1] <- DIFF_ARRAY_THEO[, , 1, 2] <- 0

pdf(file = paste(root, 'Figures/5-taylors-hypothesis-theoretical-test-functions-TEST.pdf', sep = ''), width = 25, height = 10)

split.screen( rbind(c(0.08,0.98,0.1,0.95), c(0.98,0.99,0.1,0.95)))
split.screen( figs = c( 2, 4 ), screen = 1 )

hr_label <- c('i', 'ii', 'iii', 'iv')
mod_label <- c('A', 'B')

for(model in 1:2){

	for(m in 1:length(adj_mu)){
	
		screen((model - 1) * 4 + 2 + m)
		par(mai=c(0.2,0.2,0.2,0.2))
		
		fbplot(t(DIFF_ARRAY_THEO[, , m, model]), method='MBD', ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', ylim = c(-0.25, 0.25))
		abline(h = 0, col = 3, lty = 2, lwd = 5)

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

pdf(file = paste(root, 'Figures/5-taylors-hypothesis-empirical-test-functions-TEST.pdf', sep = ''), width = 25, height = 10)

split.screen( rbind(c(0.08,0.98,0.1,0.95), c(0.98,0.99,0.1,0.95)))
split.screen( figs = c( 2, 4 ), screen = 1 )

hr_label <- c('i', 'ii', 'iii', 'iv')
mod_label <- c('A', 'B')

for(model in 1:2){

	for(m in 1:length(adj_mu)){
	
		screen((model - 1) * 4 + 2 + m)
		par(mai=c(0.2,0.2,0.2,0.2))
		
		fbplot(t(DIFF_ARRAY_EMP[, , m, model]), method='MBD', ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', ylim = c(-0.35, 0.35))
		abline(h = 0, col = 3, lty = 2, lwd = 5)

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


