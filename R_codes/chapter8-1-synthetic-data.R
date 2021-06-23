directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))

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
		axis(2, cex.axis = 2)
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
		axis(2, cex.axis = 2)
	}
	mtext(expression(s[x]), side = 1, line = 4, adj = 0.5,  cex = 2.5)
	axis(1, at = seq(0, 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"), cex.axis = 2, mgp=c(3, 3, 0))
}

screen(2)

x1 <- c(0.025,0.1,0.1,0.025) + 0.1
y1 <- c(0.3,0.3,0.7,0.7)
legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(-3, 3, length.out = 3), 1), CEX = 2)

close.screen( all=TRUE)
dev.off()


