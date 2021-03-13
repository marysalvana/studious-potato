directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))
source(file = paste(root, "R_codes/Functions/load_packages.R",sep=''))

yr <- 2017
DAT <- DAT_temp <- read.table(paste('/home/salvanmo/Desktop/ipdps/Data/pm_', yr, sep = ''), header = FALSE, sep = " ")
locs <- read.table(paste('/home/salvanmo/Desktop/ipdps/Data/locations_550', sep = ''), header = FALSE, sep = ",")

Yhat1 <- res_mat1 <- obs_mat_standardized1 <- t(DAT)

obs_mat_SVD1 <- svd(obs_mat_standardized1)

variance.explained1 = prop.table(obs_mat_SVD1$d^2)

percent_sum_squared_variation1 <- cumsum(variance.explained1)[cumsum(variance.explained1) >= 0.9]
min_percent_sum_squared_variation1 <- min(percent_sum_squared_variation1)
num_singular_vec1 <- which(cumsum(variance.explained1) == min_percent_sum_squared_variation1) 

X1 <- cbind(rep(1, nrow(obs_mat_SVD1$u)), obs_mat_SVD1$u)

for(nn in 1:ncol(obs_mat_standardized1)){

	beta_hat <- solve(t(X1[, 1:num_singular_vec1]) %*% X1[, 1:num_singular_vec1]) %*% t(X1[, 1:num_singular_vec1]) %*% obs_mat_standardized1[, nn]

	Yhat1[, nn] <- Yhat_temp <- X1[, 1:num_singular_vec1] %*% beta_hat

	err <- Yhat_temp - obs_mat_standardized1[, nn]

	res_mat1[, nn] <- err/sd(err)

}
Z <- res_mat1

hr_ind <- seq(1, 4 * 5, by = 4)			

reference_locations <- c(188, 448)

zlim_range1 <- range(res_mat1[hr_ind,])
zlim_range1 <- c(sign(min(zlim_range1)) * round(abs(min(zlim_range1)), 1) - 0.1, sign(max(zlim_range1)) * round(abs(max(zlim_range1)), 1) + 0.1)
jpeg(file = paste(root, 'Figures/application_data_residuals_sec5.jpg', sep = ''), width = 1800, height = 500)

split.screen( rbind(c(0.08,0.95,0.1,0.95), c(0.95,0.99,0.1,0.95)))
split.screen( figs = c( 1, 5 ), screen = 1 )

hr_count <- 0
for(hr in hr_ind){
	
	hr_count <- hr_count + 1
	
	screen(2 + hr_count)

	par(pty = 's')
	par(mai=c(0.2,0.2,0.2,0.2))
	
	quilt.plot(locs[, 1], locs[, 2], res_mat1[hr, ], zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, axes = F)
	points(locs[reference_locations, ], col = 'black', pch = 4, cex = 4, lwd = 4)
	text(locs[reference_locations, ], labels = c('  Ref Loc 1', '  Ref Loc 2'), col = 'black', cex = 2, adj = c(1, 1), pos = 1, offset = 1)

	if(hr_count == 1){
		axis(1, cex.axis = 2)
		axis(2, cex.axis = 2)
		mtext('log PM 2.5', side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 'blue')
	}else{
		axis(1, cex.axis = 2)
	}
	map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
	
	if(hr_count == 1){
		mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
	}
	mtext(paste(hr - 1, ':00', sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
	mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
}

screen(2)

x1 <- c(0.025,0.12,0.12,0.025) + 0.1
y1 <- c(0.3,0.3,0.7,0.7)
legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(zlim_range1[1], zlim_range1[2], length.out = 5), 1), CEX = 2)

close.screen( all=TRUE)
dev.off()


