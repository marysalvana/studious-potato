directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))

sourceCpp(file = paste(root, "R_codes/Functions/nonfrozen_integration.cpp",sep=''))




sec2_examples <- function(load_default_data = T, plot_covariances = T, plot_realizations = T){
 
	N <- 30
	n <- N^2
	TT <- 3
	grid_x <- seq(from = 0, to = 1, length.out = N)
	grid_y <- seq(from = 0, to = 1, length.out = N)
	locs <- sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

	reference_locations <- c(1, ceiling(n / 2) + ceiling(N / 2))

	if(load_default_data){
		cat('Loading covariance model 1 values...', '\n')
		cov1 <- read.table(paste(root, "Data/sec2_example1_frozen_covariance", sep = ''), header = FALSE, sep = " ") %>% as.matrix() 
		cat('Loading covariance model 2 values...', '\n')
		cov2 <- read.table(paste(root, "Data/sec2_example2_frozen_covariance", sep = ''), header = FALSE, sep = " ") %>% as.matrix() 

		cat('Loading realizations from model 1...', '\n')
		r1 <- read.table(paste(root, "Data/sec2_example1_frozen_realizations", sep = ''), header = FALSE, sep = " ") %>% as.matrix() 
		cat('Loading realizations from model 2...', '\n')
		r2 <- read.table(paste(root, "Data/sec2_example2_frozen_realizations", sep = ''), header = FALSE, sep = " ") %>% as.matrix() 

	}else{

		cov1 <- nonfrozen_matern_uni_cov_distributed_nonstationary(theta = c(1, 1, 1, 0.1001, 0.1001, 0.001, 0.001, 0), locs = sim_grid_locations, TT = 3)
		cov2 <- nonfrozen_matern_uni_cov_distributed_deform(theta = c(1, 1, 1, 0.1001, 0.1001, c(0.001, 0.001, 0)), locs = sim_grid_locations, TT = 3)
		#save(cov1, cov2, file=paste(root, "Results/sec2_covariances.RData", sep = ''))

		set.seed(1234)
		r1 <- mvrnorm(100, rep(0, ncol(cov1)), cov1 )
		set.seed(1234)
		r2 <- rmvn(100, rep(0, ncol(cov2)), cov2 + 0.009 * diag(n * TT), ncores = 25)
		#save(r1, r2, file=paste(root, "Results/sec2_realizations.RData", sep = ''))

	}

	COV_LIST <- list()
	COV_LIST[[1]] <- cov1[reference_locations, ]
	COV_LIST[[2]] <- cov2[reference_locations,]

	## Plot the heatmaps of the covariances at different times and different reference locations

	new_grid_x <- matrix(locs[1:n, 1], N, N)
	new_grid_y <- matrix(locs[1:n, 2], N, N, byrow = F)

	mod_labels <- c('(a)', '(b)', '(c)', '(d)')

	if(plot_covariances){

		cat('Plotting covariances...', '\n')

		point_label <- c('  Ref Loc 1', '  Ref Loc 2')


		pdf(file = paste(root, 'Figures/sec2_example_frozen_covariance.pdf', sep = ''), width = 30, height = 11)

		split.screen( rbind(c(0.05,0.48,0.06,0.88), c(0.52,0.95,0.06,0.88), c(0.96,0.99,0.06,0.88)))
		split.screen( figs = c( 2, 3 ), screen = 1 )
		split.screen( figs = c( 2, 3 ), screen = 2 ) 

		for(cov_mod in 1:length(COV_LIST)){

			for(ll in 1:2){

				for(tt in 1:3){

					matrix_vals <- matrix(COV_LIST[[cov_mod]][ll, (tt - 1) * n + 1:n], N, N)

					screen((cov_mod - 1) * 6 + (ll - 1) * 3 + tt + 3)	

					par(pty="s")
					if(cov_mod == 1) par(mai=c(0.1,0.03,0.03,0.03)) else par(mai=c(0.03,0.03,0.03,0.03))

					
					poly.image(new_grid_x, new_grid_y, matrix_vals, xlab = " ", ylab = " ", cex.axis = 2, yaxt = 'n', xaxt = 'n', zlim = range(COV_LIST[[cov_mod]][ll, 1:(n * 3)]))	
					points(matrix(sim_grid_locations[reference_locations[ll], ], ncol = 2), col = 'black', pch = 4, cex = 4, lwd = 4)
  					text(matrix(sim_grid_locations[reference_locations[ll], ], ncol = 2), labels = point_label[ll], col = 'black', pch = 17, cex = 2, adj = c(1, 1), pos = 4)
					if(tt == 1 & (cov_mod == 1 | cov_mod == 3)){
						mtext(expression(s[y]), side = 2, line = 3, cex = 2.5)
						axis(2, at = seq(0, 1, length.out = 6), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis = 2)
						mtext(paste("Reference Location ", ll, sep = ""), side = 2, line = 6, adj = 0.5, cex = 2, font = 2, col = 4)
					}
					if(ll == 2){
						mtext(expression(s[x]), side = 1, line = 3, cex = 2.5)
						axis(1, at = seq(0, 1, length.out = 6), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis = 2)
					}
					if(ll == 1 & cov_mod <= 2){
						label1 <- bquote(paste("C(", bold(s)[1], "=", bold(s)[ref], ",", bold(s)[2], "=\U2022;", t[1], "=", 0, ",", t[2], "=", .(tt - 1), ")", sep=""))
						mtext(label1, side = 3, line = 1, cex = 2)
					}
					if(tt == 2 & ll == 1){
						if(cov_mod <= 2){
							mtext(mod_labels[cov_mod], side = 3, line = 4, cex = 3, font = 2, col=4)
						}else{
							mtext(mod_labels[cov_mod], side = 3, line = 2, cex = 3, font = 2, col=4)
						}
					}
				}
			}
		}

		screen(3)
		x1 <- c(0.025,0.12,0.12,0.025)
		y1 <- c(0.3,0.3,0.7,0.7)
		legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(0, 1, length.out = 3), 1), CEX = 2)

		close.screen( all=TRUE)
		dev.off()


	}

	## Plot the realizations

	REALIZATIONS_LIST <- array(, dim = c(nrow(r1), ncol(r1), 2))
	REALIZATIONS_LIST[,, 1] <- r1
	REALIZATIONS_LIST[,, 2] <- r2

	zlim_range1 <- range(REALIZATIONS_LIST[8,,])

	if(plot_realizations){

		cat('Plotting realizations...', '\n')

		pdf(file = paste(root, 'Figures/sec2_example_frozen_realizations.pdf', sep = ''), width = 30, height = 6.5)

		split.screen( rbind(c(0.05,0.48,0.06,0.88), c(0.52,0.95,0.06,0.88), c(0.95,0.99,0.06,0.88)))
		split.screen( figs = c( 1, 3 ), screen = 1 )
		split.screen( figs = c( 1, 3 ), screen = 2 ) 

		for(cov_mod in 1:2){

			for(tt in 1:3){

				matrix_vals <- matrix(REALIZATIONS_LIST[8, (tt - 1) * n + 1:n, cov_mod], N, N)

				screen((cov_mod - 1) * 3 + tt + 3)	

				par(pty="s")
				par(mai=c(0.03,0.03,0.03,0.03))

				
				poly.image(new_grid_x, new_grid_y, matrix_vals, zlim = zlim_range1, xlab = " ", ylab = " ", cex.axis = 2, yaxt = 'n', xaxt = 'n')	
				points(sim_grid_locations[reference_locations, ], col = 'black', pch = 4, cex = 4, lwd = 4)
  				text(sim_grid_locations[reference_locations, ], labels = c('  Ref Loc 1', '  Ref Loc 2'), col = 'black', pch = 17, cex = 2, adj = c(1, 1), pos = 4)
				if(tt == 1 & (cov_mod == 1 | cov_mod == 3)){
					mtext(expression(s[y]), side = 2, line = 3, cex = 2.5)
					axis(2, at = seq(0, 1, length.out = 6), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis = 2)
				}
					mtext(expression(s[x]), side = 1, line = 3, cex = 2.5)
					axis(1, at = seq(0, 1, length.out = 6), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis = 2)
				if(cov_mod <= 2){
					label1 <- bquote(paste(t, "=", .(tt - 1), sep=""))
					mtext(label1, side = 3, line = 1, cex = 2.5, font = 2)
				}
				if(tt == 2){
					if(cov_mod <= 2){
						mtext(mod_labels[cov_mod], side = 3, line = 3.5, cex = 3, font = 2, col=4)
					}else{
						mtext(mod_labels[cov_mod], side = 3, line = 2, cex = 3, font = 2, col=4)
					}
				}
			}
		}

		screen(3)
		x1 <- c(0.025,0.12,0.12,0.025) + 0.1
		y1 <- c(0.3,0.3,0.7,0.7)
		legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(-3, 3, length.out = 5), 1), CEX = 2)

		close.screen( all=TRUE)
		dev.off()
	}
}






sec2_examples(load_default_data = F, plot_covariances = T, plot_realizations = T)


