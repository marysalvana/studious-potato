


source("./pkg-config.R")


PLOT_REALIZATIONS = F
PLOT_PARAMETER_ESTIMATION = T

if(PLOT_REALIZATIONS){



	################################################                                      			      	################################################
	################################################   PLOT REALIZATION FROM EXAGEOSTAT LAGRANGIAN  MODEL   	################################################
	################################################                                      				################################################

	obs1 <- read.table(paste(root, 'Data/synthetic/Z_mu_0.1_var_0.001', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	locs1 <- read.table(paste(root, 'Data/synthetic/LOC_mu_0.1_var_0.001', sep = ''), header = FALSE, sep = ",") %>% as.matrix()

	obs2 <- read.table(paste(root, 'Data/synthetic/Z_mu_0.1_var_0.1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	locs2 <- read.table(paste(root, 'Data/synthetic/LOC_mu_0.1_var_0.1', sep = ''), header = FALSE, sep = ",") %>% as.matrix()

	zlim_range <- range(c(obs1, obs2))

	jpeg(file = paste(root, 'Figures/8-synthetic_data_lagrangian1.jpg', sep = ''), width = 1500, height = 700)

	split.screen( rbind(c(0.08,0.94,0.12,0.92), c(0.94,0.98,0.12,0.92)))
	split.screen( figs = c( 2, 5 ), screen = 1 )

	hr_count <- 0

	for(hr in 1:5){
		
		hr_count <- hr_count + 1
		
		screen(2 + hr_count)

		ind <- which(locs1[, 3] == hr)

		par(pty = 's')
		par(mai=c(0.2,0.2,0.2,0.2))
		quilt.plot(locs1[ind, 1], locs1[ind, 2], obs1[ind, 1], ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, nx = 50, ny = 50, zlim = zlim_range, axes = F)

		if(hr == 1){
			mtext(expression(s[y]), side = 2, line = 4, adj = 0.5, cex = 2.5)
			mtext('Strong', side = 2, line = 7, adj = 0.5, cex = 2.5, col = 4, font = 2)
			axis(2, cex.axis = 1)
		}
		mtext(paste('t = ', hr, sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
	}

	for(hr in 1:5){
		
		hr_count <- hr_count + 1
		
		screen(2 + hr_count)
		
		ind <- which(locs2[, 3] == hr)

		par(pty = 's')
		par(mai=c(0.2,0.2,0.2,0.2))
		quilt.plot(locs2[ind, 1], locs2[ind, 2], obs2[ind, 1], ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, nx = 50, ny = 50, zlim = zlim_range, axes = F)

		if(hr == 1){
			mtext(expression(s[y]), side = 2, line = 4, adj = 0.5, cex = 2.5)
			mtext('Weak', side = 2, line = 7, adj = 0.5, cex = 2.5, col = 4, font = 2)
			axis(2, cex.axis = 1)
		}
		mtext(expression(s[x]), side = 1, line = 4, adj = 0.5,  cex = 2.5)
		axis(1, at = seq(0, 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"), cex.axis = 1)
	}

	screen(2)

	x1 <- c(0.025,0.1,0.1,0.025) + 0.1
	y1 <- c(0.3,0.3,0.7,0.7)
	legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(-3, 3, length.out = 3), 1), CEX = 2)

	close.screen( all=TRUE)
	dev.off()


	obs1 <- read.table(paste(root, 'Data/synthetic/Z_mu_0_var_0.001', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	locs1 <- read.table(paste(root, 'Data/synthetic/LOC_mu_0_var_0.001', sep = ''), header = FALSE, sep = ",") %>% as.matrix()

	obs2 <- read.table(paste(root, 'Data/synthetic/Z_mu_0_var_0.1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	locs2 <- read.table(paste(root, 'Data/synthetic/LOC_mu_0_var_0.1', sep = ''), header = FALSE, sep = ",") %>% as.matrix()

	zlim_range <- range(c(obs1, obs2))

	jpeg(file = paste(root, 'Figures/8-synthetic_data_lagrangian2.jpg', sep = ''), width = 1500, height = 700)

	split.screen( rbind(c(0.08,0.94,0.12,0.92), c(0.94,0.98,0.12,0.92)))
	split.screen( figs = c( 2, 5 ), screen = 1 )

	hr_count <- 0

	for(hr in 1:5){
		
		hr_count <- hr_count + 1
		
		screen(2 + hr_count)

		ind <- which(locs1[, 3] == hr)

		par(pty = 's')
		par(mai=c(0.2,0.2,0.2,0.2))
		quilt.plot(locs1[ind, 1], locs1[ind, 2], obs1[ind, 1], ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, nx = 50, ny = 50, zlim = zlim_range, axes = F)

		if(hr == 1){
			mtext(expression(s[y]), side = 2, line = 4, adj = 0.5, cex = 2.5)
			mtext('Strong', side = 2, line = 7, adj = 0.5, cex = 2.5, col = 4, font = 2)
			axis(2, cex.axis = 1)
		}
		mtext(paste('t = ', hr, sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
	}

	for(hr in 1:5){
		
		hr_count <- hr_count + 1
		
		screen(2 + hr_count)
		
		ind <- which(locs2[, 3] == hr)

		par(pty = 's')
		par(mai=c(0.2,0.2,0.2,0.2))
		quilt.plot(locs2[ind, 1], locs2[ind, 2], obs2[ind, 1], ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, nx = 50, ny = 50, zlim = zlim_range, axes = F)

		if(hr == 1){
			mtext(expression(s[y]), side = 2, line = 4, adj = 0.5, cex = 2.5)
			mtext('Weak', side = 2, line = 7, adj = 0.5, cex = 2.5, col = 4, font = 2)
			axis(2, cex.axis = 1)
		}
		mtext(expression(s[x]), side = 1, line = 4, adj = 0.5,  cex = 2.5)
		axis(1, at = seq(0, 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"), cex.axis = 1)
	}

	screen(2)

	x1 <- c(0.025,0.1,0.1,0.025) + 0.1
	y1 <- c(0.3,0.3,0.7,0.7)
	legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(-3, 3, length.out = 3), 1), CEX = 2)

	close.screen( all=TRUE)
	dev.off()

}



if(PLOT_PARAMETER_ESTIMATION){
	################################################                                 ################################################
	################################################      PARAMETER ESTIMATION       ################################################
	################################################                                 ################################################

	DAT <- read.table(paste(root, 'Results/parameter-estimates-lagrangian.log', sep = ''), header = FALSE, sep = " ")

	n <- 8
	# width of each boxplot is 0.8
	x0s <- 1:n - 0.4
	x1s <- 1:n + 0.4
	# these are the y-coordinates for the horizontal lines
	# that you need to set to the desired values.
	y0s <- matrix(c(1, 0.03, 1, 0.1, 0.1, 0.1, 0, 0.1,
			1, 0.09, 1, 0.1, 0.1, 0.1, 0, 0.1,
			1, 0.23, 1, 0.1, 0.1, 0.1, 0, 0.1,
			1, 0.03, 1, 0.1, 0.1, 0.1, 0, 0.1,
			1, 0.09, 1, 0.1, 0.1, 0.1, 0, 0.1,
			1, 0.23, 1, 0.1, 0.1, 0.1, 0, 0.1,
			1, 0.03, 1, 0.1, 0.1, 0.1, 0, 0.1,
			1, 0.09, 1, 0.1, 0.1, 0.1, 0, 0.1,
			1, 0.23, 1, 0.1, 0.1, 0.1, 0, 0.1), nrow = 9, ncol = 8, byrow = T)

	ER_lables <- c('Weak', 'Moderate', 'Strong')

	pdf(file = paste(root, 'Figures/8-parameter-estimates-lagrangian.pdf', sep = ''), width = 8, height = 9)

	split.screen( rbind(c(0.12,0.98,0.1,0.95), c(0.93,0.99,0.1,0.95)))
	split.screen( figs = c( 3, 3 ), screen = 1 )

	bb <- 1
	ff <- 1

	for(dd in 1:3){
		screen(dd + 2)
		par(pty = 's')
		par(mai = c(0.1, 0.1, 0.1, 0.1))
	  
		if(dd %% 3 == 1){
			boxplot(DAT[which(DAT[, 16] == y0s[dd, 2] & DAT[, 20] == y0s[dd, 6] & DAT[, 21] == y0s[dd, 7] & DAT[, 22] == y0s[dd, 8]), 31:38], ylim = c(0, 1.5), xaxt = 'n')
			mtext(ER_lables[bb], side = 2, line = 2.5, cex = 0.8, col = '#0993f3', font = 2)
			if(bb == 2){
				mtext('T    I    M    E', side = 2, line = 3.5, cex = 1.5, col = 'blue', font = 2)
			}
			bb <- bb + 1
		}else{
			boxplot(DAT[which(DAT[, 16] == y0s[dd, 2] & DAT[, 20] == y0s[dd, 6] & DAT[, 21] == y0s[dd, 7] & DAT[, 22] == y0s[dd, 8]), 31:38], ylim = c(0, 1.5), xaxt = 'n', yaxt = 'n')
		}
		segments(x0 = x0s, x1 = x1s, y0 = y0s[dd, ], col = "red", lwd = 2)
	  
		if(dd > 6){
			axis(1, at = 1:6, labels = c(expression(sigma^2), expression(a[s]), expression(nu), expression(a[t]), expression(alpha), expression(beta)))
		}
		if(dd < 4){
			mtext(ER_lables[ff], side = 3, line = 0.5, cex = 0.8, col = '#0993f3', font = 2)
			if(ff == 2){
				mtext('S    P    A    C    E', side = 3, line = 1.5, cex = 1.5, col = 'blue', font = 2)
			}
			ff <- ff + 1
		}  
	}

	close.screen( all=TRUE)
	dev.off()
}
