


source("./pkg-config.R")



radian <- seq(-2 * pi, 2 * pi, length = 2000)
radius <- 1
hx = radius *  cos(radian)
hy = radius *  sin(radian)

C_G_list <- list()
C_LGR_list <- list()

for(u in 1:3){
	sigma_sq <- 1
	C_G <- 1 / (1 + sigma_sq * u^2) * exp( - ((hx)^2 + (hy)^2) / (1 + sigma_sq * u^2))

	sigma_11_sq <- 1
	sigma_22_sq <- 1
	sigma_12_sq_vals <- c(0.1, 0.5, 0.9)

	C_LGR <- matrix(, ncol = length(sigma_12_sq_vals), nrow = length(C_G))
	for(cc in 1:length(sigma_12_sq_vals)){
	
		sigma_12_sq <- sigma_12_sq_vals[cc]
		det_term <- (1 + sigma_11_sq * u^2) * (1 + sigma_22_sq * u^2) - sigma_12_sq^2 * u^4
		Sigma_inv_numerator <- hx^2 * (1 + sigma_22_sq * u^2) + hy^2 * (1 + sigma_11_sq * u^2) - hx * hy * sigma_12_sq * u^2  - hx * hy * sigma_12_sq * u^2 
		Sigma_inv <- Sigma_inv_numerator / det_term

		C_LGR[, cc] <- 1 / sqrt(det_term) * exp(- Sigma_inv)

	}
	C_G_list[[u]] <- C_G
	C_LGR_list[[u]] <- C_LGR
}

label1 <- c(expression(paste(rho, " = 0", sep = '')), expression(paste(rho, " = 0.1", sep = '')), expression(paste(rho, " = 0.5", sep = '')), expression(paste(rho, " = 0.9", sep = '')))

pdf(file = paste(root, 'Figures/simulation_stationary.pdf', sep = ''), width = 15, height = 6)

split.screen( rbind(c(0.05,0.99,0.1,0.95), c(0.99,0.99,0.1,0.95)))
split.screen( figs = c( 1, 3 ), screen = 1 ) 

lims <- max(abs(c(C_LGR_list[[1]], C_G_list[[1]])))

for(u in 1:3){

	screen(u + 2)

	par(pty="s")
	par(mai=c(0.3, 0.3, 0.3, 0.3))

	if(u == 1)	plot(0,0, 'n', ylab = '', xlab = '', ylim = c(-lims, lims), xlim = c(-lims, lims), cex.axis = 2)
	else plot(0,0, 'n', ylab = '', xlab = '', ylim = c(-lims, lims), xlim = c(-lims, lims), cex.axis = 2, yaxt = 'n')
	mtext(paste('u = ', u, sep = ''), line = 1, cex = 2, side = 3)
	if(u == 1)	mtext(expression(h[y]), line = 3, cex = 2, side = 2)
	mtext(expression(h[x]), line = 4, cex = 2, side = 1)

	if(u == 3){
		legend(0.15, 0.4, legend = label1, lty = rep(1, 4), col = 1:4, cex=1, box.lty=0)
	}

	abline(h = 0, v = 0, lty = 3)
	lines(x = C_G_list[[u]] * cos(radian), y = C_G_list[[u]] * sin(radian), col = 1)

	for(cc in 1:length(sigma_12_sq_vals)){
		lines(x = C_LGR_list[[u]][, cc] * cos(radian), y = C_LGR_list[[u]][, cc] * sin(radian), col = cc + 1)
	}
}
close.screen( all=TRUE)
dev.off()

#############################################  PARAMETER ESTIMATES ANALYSIS  ###########################################

load(paste(root, "Results/est_params_C_LGR_D_LGR.RData", sep = ''))
load(paste(root, "Results/est_params_C_LGR_D_LGR_B.RData", sep = ''))
load(paste(root, "Results/est_params_C_LGR_D_LGR_C.RData", sep = ''))

est_sigma_vals <- matrix(est_params_C_LGR_D_LGR_final[, 3, 1], ncol = 1)
est_sigma_vals2 <- matrix(est_params_C_LGR_D_LGR_B_final[, 3, 1], ncol = 1)
est_sigma_vals3 <- matrix(est_params_C_LGR_D_LGR_C_final[, 3, 1], ncol = 1)

for(aa in 2:9){
	est_sigma_vals <- cbind(est_sigma_vals, matrix(est_params_C_LGR_D_LGR_final[, 3, aa], ncol = 1))
	est_sigma_vals2 <- cbind(est_sigma_vals2, matrix(est_params_C_LGR_D_LGR_B_final[, 3, aa], ncol = 1))
	est_sigma_vals3 <- cbind(est_sigma_vals3, matrix(est_params_C_LGR_D_LGR_C_final[, 3, aa], ncol = 1))
}

pdf(file = paste(root, 'Figures/simulation_stationary_sigma_estimates.pdf', sep = ''), width = 15, height = 6)

split.screen( rbind(c(0.05,0.99,0.1,0.95), c(0.99,0.99,0.1,0.95)))
split.screen( figs = c( 1, 3 ), screen = 1 ) 


screen(3)
par(mai=c(0.4, 0.4, 0.4, 0.4))
boxplot(est_sigma_vals, ylab = '', xlab = '', xaxt = 'n')
axis(1, at = 1:9, labels = (1:9)/10 )
mtext(expression(rho), side = 1, line = 4, cex = 2)
mtext(expression(hat(sigma)[bold(V)]^2), side = 2, line = 3, cex = 2)
mtext('(a)', line = 1, cex = 2, side = 3)

screen(4)
par(mai=c(0.4, 0.4, 0.4, 0.4))
boxplot(est_sigma_vals2, ylab = '', xlab = '', xaxt = 'n')
axis(1, at = 1:9, labels = (1:9)/10 )
mtext(expression(sigma[y]^2), side = 1, line = 4, cex = 2)
mtext('(b)', line = 1, cex = 2, side = 3)

screen(5)
par(mai=c(0.4, 0.4, 0.4, 0.4))
boxplot(est_sigma_vals3, ylab = '', xlab = '', xaxt = 'n')
axis(1, at = 1:9, labels = (1:9)/10 )
mtext(expression(mu), side = 1, line = 4, cex = 2)
mtext('(c)', line = 1, cex = 2, side = 3)

close.screen( all=TRUE)
dev.off()


#############################################  MSPE ANALYSIS  ###########################################

N <- 10
n <- N^2
TT <- 10
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
locs <- sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

insample_index <- 1:900
outsample_index <- 901:1000

rho_vals <- -0.1 * seq(1, 9, by = 1) / 10 

C_G_stationary_est <- gneiting_st_uni_cov(theta = c(1, 0.23, 0.1), max_time_lag = TT - 1, LOCS = sim_grid_locations)
C_G_insample <- C_G_stationary_est[insample_index, insample_index] + 0.00001 * diag(length(insample_index))
C_G_outsample <- C_G_stationary_est[outsample_index, insample_index]		

MSE_C_G <- MSE_C_LGR <- matrix(, ncol = length(rho_vals), nrow = 100)

for(cc in 1:length(rho_vals)){

	C_LGR_stationary <- nonfrozen_schlather_uni_cov_distributed(theta = c(1, 0.23, 0, 0, 0.1, 0.1, rho_vals[cc]), locs = sim_grid_locations, TT = TT)
	C_LGR_insample <- C_LGR_stationary[insample_index, insample_index] + 0.00001 * diag(length(insample_index))
	C_LGR_outsample <- C_LGR_stationary[outsample_index, insample_index]		

	MSE_C_G_temp <- 1 - diag(2 * C_LGR_outsample %*% solve(C_G_insample) %*% t(C_G_outsample) - C_G_outsample %*% solve(C_G_insample) %*% C_LGR_insample %*% solve(C_G_insample) %*% t(C_G_outsample))

	MSE_C_G[, cc] <- MSE_C_G_temp

	MSE_C_LGR_temp <- 1 - diag(C_LGR_outsample %*% solve(C_LGR_insample) %*% t(C_LGR_outsample))
	MSE_C_LGR[, cc] <- MSE_C_LGR_temp

}

MSE_C_G_OLD <- 1 - diag(C_G_outsample %*% solve(C_G_insample) %*% t(C_G_outsample))

MLOE <- MSE_C_G / MSE_C_LGR - 1
MMOM <- MSE_C_G_OLD / MSE_C_G - 1

pdf(file = paste(root, 'Figures/simulation_mse_stationary.pdf', sep = ''), width = 35, height = 10)

split.screen( rbind(c(0.05,0.95,0.526,0.95), c(0.05,0.95,0.1,0.524), c(0.95,0.99,0.526,0.95), c(0.95,0.99,0.1,0.524)))
split.screen( figs = c( 1, 9 ), screen = 1 ) 
split.screen( figs = c( 1, 9 ), screen = 2 ) 


for(zz in 1:ncol(MLOE)){

	RHO <- zz / 10
	screen(4 + zz)
	par(pty = 's')
	par(mai=c(0.05, 0.05, 0.05, 0.05))
	quilt.plot(locs[, 1], locs[, 2], MLOE[, zz], zlim =  range(MLOE), nx = 10, ny = 10, ylab = '', xlab = '', cex.axis = 2, add.legend = F, axes = F)
	mtext(bquote(rho == .(RHO)), line = 1, cex = 3, side = 3)

	if(zz == 1){
		 axis(2); mtext(expression(s[y]), line = 3, cex = 2, side = 2); mtext('LOE', side = 2, line = 6, adj = 0.5, cex = 3, font = 2, col = 'blue')
	}
}

for(zz in 1:ncol(MLOE)){
	screen(13 + zz)
	par(pty = 's')
	par(mai=c(0.05, 0.05, 0.05, 0.05))
	quilt.plot(locs[, 1], locs[, 2], MMOM[, zz], zlim =  range(MMOM), nx = 10, ny = 10, ylab = '', xlab = '', cex.axis = 2, add.legend = F, axes = F)
	mtext(expression(s[x]), line = 3, cex = 2, side = 1)
	axis(1)
	if(zz == 1){
		 axis(2); mtext(expression(s[y]), line = 3, cex = 2, side = 2); mtext('MOM', side = 2, line = 6, adj = 0.5, cex = 3, font = 2, col = 'blue')
	}

}

screen(3)
x1 <- c(0.025,0.12,0.12,0.025) + 0.05
y1 <- c(0.25,0.25,0.7,0.7)
legend.gradient2(cbind(x1,y1), title = "", limits = seq(min(MLOE), max(MLOE), length.out = 3), CEX = 2)

screen(4)
x1 <- c(0.025,0.12,0.12,0.025) + 0.05
y1 <- c(0.25,0.25,0.7,0.7)
legend.gradient2(cbind(x1,y1), title = "", limits = seq(min(MMOM), max(MMOM), length.out = 3), CEX = 2)

close.screen( all=TRUE)
dev.off()

################## COMPARE C_G_NS WITH C_LGR_matern_nonstationary VIA LOE AND MOM ###################

LOE_mean <- MOM_mean <- LOE_mean_2 <- MOM_mean_2 <- c()

for(cc in 1:length(rho_vals)){
	cat(cc, '\n')

	C_LGR_matern_nonstationary <- nonfrozen_matern_uni_cov_distributed_nonstationary(theta = c(1, 0.23, 1, 0, 0, c(0.1, 0.1, rho_vals[cc])), locs = sim_grid_locations, TT = TT)
	C_G_NS <- matern_st_uni_cov_NS(theta = c(1, 1, 0.1), max_time_lag = 9, LOCS = sim_grid_locations)


	C_G_insample <- C_G_NS[insample_index, insample_index] + 0.00001 * diag(length(insample_index))
	C_G_outsample <- C_G_NS[outsample_index, insample_index]		

	C_LGR_insample <- C_LGR_matern_nonstationary[insample_index, insample_index] + 0.00001 * diag(length(insample_index))
	C_LGR_outsample <- C_LGR_matern_nonstationary[outsample_index, insample_index]		

	MSE_C_G_temp <- 1 - diag(2 * C_LGR_outsample %*% solve(C_G_insample) %*% t(C_G_outsample) - C_G_outsample %*% solve(C_G_insample) %*% C_LGR_insample %*% solve(C_G_insample) %*% t(C_G_outsample))

	MSE_C_G <- MSE_C_G_temp

	MSE_C_LGR_temp <- 1 - diag(C_LGR_outsample %*% solve(C_LGR_insample) %*% t(C_LGR_outsample))
	MSE_C_LGR <- MSE_C_LGR_temp

	MSE_C_G_OLD <- 1 - diag(C_G_outsample %*% solve(C_G_insample) %*% t(C_G_outsample))

	MLOE <- MSE_C_G / MSE_C_LGR - 1
	MMOM <- MSE_C_G_OLD / MSE_C_G - 1

	MSE_C_G_mis <- 1 - diag(2 * C_G_outsample %*% solve(C_LGR_insample) %*% t(C_LGR_outsample) - C_LGR_outsample %*% solve(C_LGR_insample) %*% C_G_insample %*% solve(C_LGR_insample) %*% t(C_LGR_outsample))
	MLOE_CLGR_INSTEAD_OF_CG <- MSE_C_G_mis / MSE_C_G_OLD - 1
	MMOM_CLGR_INSTEAD_OF_CG <- MSE_C_LGR / MSE_C_G_mis - 1

	LOE_mean <- c(LOE_mean, mean(MLOE))
	MOM_mean <- c(MOM_mean, mean(MMOM))
	LOE_mean_2 <- c(LOE_mean_2, mean(MLOE_CLGR_INSTEAD_OF_CG))
	MOM_mean_2 <- c(MOM_mean_2, mean(MMOM_CLGR_INSTEAD_OF_CG))
}


pdf(file = paste(root, 'Figures/simulation_mse_NS.pdf', sep = ''), width = 15, height = 5)

par(mfrow = c(1,2))
plot(LOE_mean, type = 'l', col = 1, lwd = 3, xlab = expression(rho), ylim = c(0, max(LOE_mean)), ylab = 'MLOE', xaxt = 'n')
lines(LOE_mean_2, col = 2, lwd = 3)
axis(1, at = 1:6, labels = (1:6)/10)

plot(MOM_mean, type = 'l', col = 1, lwd = 3, xlab = expression(rho), ylim = c(min(MOM_mean), max(MOM_mean_2)), ylab = 'MMOM', xaxt = 'n')
lines(MOM_mean_2, col = 2, lwd = 3)
axis(1, at = 1:6, labels = (1:6)/10)

dev.off()
