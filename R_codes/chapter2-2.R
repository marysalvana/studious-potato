directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))

sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))


N <- 50
n <- N^2
TT <- 3
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

reference_locations <- c(1, ceiling(n / 2) + ceiling(N / 2))

ind <- 6 #3, 4, 6, 8, 9
set.seed(ind)
PARAMETER_NONSTAT <- runif(15, -3, 3)

set.seed(3)
PARAMETER_DEFORMATION <- c(runif(2, -8, 8), runif(1, 0, 8), runif(2, -1, 1))

cat('Computing covariances...', '\n')

cov1 <- MATERN_UNI_SPATIALLY_VARYING_PARAMETERS(PARAMETER = c(1, 0.23, 1, 0.1001, 0.1001, 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_NONSTAT = PARAMETER_NONSTAT)

#chol(cov1)
#ind <- ind + 1

cov2 <- MATERN_UNI_DEFORMATION(PARAMETER = c(1, 0.23, 1, 0.1001, 0.1001, 0.001, 0, 0, 0.001), LOCATION = sim_grid_locations, TIME = TT, PARAMETER_DEFORMATION = PARAMETER_DEFORMATION)

cat('Generating realizations...', '\n')

set.seed(1)
r1 <- rmvn(100, rep(0, ncol(cov1)), cov1, ncores = 25)

set.seed(1)
r2 <- rmvn(100, rep(0, ncol(cov2)), cov2, ncores = 25)

plot_covariances <- T
plot_realizations <- T

COV_LIST <- list()
COV_LIST[[1]] <- cov1[reference_locations, ]
COV_LIST[[2]] <- cov2[reference_locations,]

## Plot the heatmaps of the covariances at different times and different reference locations

new_grid_x <- matrix(sim_grid_locations[1:n, 1], N, N)
new_grid_y <- matrix(sim_grid_locations[1:n, 2], N, N, byrow = F)

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

warnings()

## Plot the realizations

if(plot_realizations){

	cat('Plotting realizations...', '\n')

	REALIZATIONS_LIST <- array(, dim = c(nrow(r1), ncol(r1), 2))
	REALIZATIONS_LIST[,, 1] <- r1
	REALIZATIONS_LIST[,, 2] <- r2

	zlim_range1 <- range(REALIZATIONS_LIST[8,,])

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

