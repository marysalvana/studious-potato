directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))


################################################                                      ################################################
################################################   NONSEPARABLE VS SEPARABLE MODELS   ################################################
################################################                                      ################################################

# simulate 100 datasets from exageostat with space-time interaction parameter beta = 0.1, 0.5, 0.9 at N = 400 spatial locations and T = 100 temporal locations

# save the measurements and space time locations in folder Data/synthetic/beta-0.1/

# split simulated data into 60% training and 40% testing

N = 40000

# if the subset for training is chosen randomly

for(set in 1:100){
	obs <- read.table(paste(root, 'Data/synthetic/beta-0.1/Z_' , N, '_', set, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	locs <- read.table(paste(root, 'Data/synthetic/beta-0.1/LOC_', N, '_', set, sep = ''), header = FALSE, sep = ",") %>% as.matrix()

	n <- nrow(locs)
	set.seed(1234)
	subset_index <- sample(1:n, floor(0.4 * n))

	# SAVING TRAINING AND TESTING DATASETS training and testing datasets

	write.table(locs[-subset_index, 1:2], file = paste(root, "Data/synthetic/beta-0.1/LOC_SPACE_", N, "_", set, "_training", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
	write.table(locs[-subset_index, 3], file = paste(root, "Data/synthetic/beta-0.1/LOC_TIME_", N, "_", set, "_training", sep = ''), sep = ",",row.names = FALSE, col.names = FALSE)
	write.table(obs[-subset_index] - mean(obs[-subset_index]), file = paste(root, "Data/synthetic/beta-0.1/Z_", N, "_", set, "_training", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)

	write.table(locs[subset_index, 1:2], file = paste(root, "Data/synthetic/beta-0.1/LOC_SPACE_", N, "_", set, "_testing", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
	write.table(locs[subset_index, 3], file = paste(root, "Data/synthetic/beta-0.1/LOC_TIME_", N, "_", set, "_testing", sep = ''), sep = ",",row.names = FALSE, col.names = FALSE)
	write.table(obs[subset_index] - mean(obs[subset_index]), file = paste(root, "Data/synthetic/beta-0.1/Z_", N, "_", set, "_testing", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)

}

# if the subset for training is chosen using the 10 days and testing is the future 10 days
# this is for N = 400 and T = 20

in_sample <- 1:(400 * 10)
out_sample <- max(in_sample) + 1:(400 * 1)

for(set in c(2:5, 7:67, 69:100)){
	obs <- read.table(paste(root, 'Data/synthetic/beta-1/Z_', N, '_', set, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	locs <- read.table(paste(root, 'Data/synthetic/beta-1/LOC_', N, '_', set, sep = ''), header = FALSE, sep = ",") %>% as.matrix()

	n <- nrow(locs)
	write.table(locs[in_sample, 1:2], file = paste(root, "Data/synthetic/beta-1/LOC_SPACE_", N, "_", set, "_training", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
	write.table(locs[in_sample, 3], file = paste(root, "Data/synthetic/beta-1/LOC_TIME_", N, "_", set, "_training", sep = ''), sep = ",",row.names = FALSE, col.names = FALSE)
	write.table(obs[in_sample] - mean(obs[in_sample]), file = paste(root, "Data/synthetic/beta-1/Z_", N, "_", set, "_training", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)

	write.table(locs[out_sample, 1:2], file = paste(root, "Data/synthetic/beta-1/LOC_SPACE_", N, "_", set, "_testing", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
	write.table(locs[out_sample, 3], file = paste(root, "Data/synthetic/beta-1/LOC_TIME_", N, "_", set, "_testing", sep = ''), sep = ",",row.names = FALSE, col.names = FALSE)
	write.table(obs[out_sample] - mean(obs[out_sample]), file = paste(root, "Data/synthetic/beta-1/Z_", N, "_", set, "_testing", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
}


# RETRIEVE THE RESULTS FROM EXAGEOSTAT AND FORMAT THE MSPE AND TOTAL KRIGING VALUES INTO A MATRIX OF 6 COLUMNS AND SAVE THEM IN THE TEXT FILE nonsep_vs_sep_mspe AND nonsep_vs_sep_kriging_variance

# PLOT THE VALUES

DAT_MSPE <- read.table(paste(root, 'Results/nonsep_vs_sep_mspe', sep = ''), header = FALSE, sep = ",") %>% as.matrix()
DAT_KRIGING_VARIANCE <- read.table(paste(root, 'Results/nonsep_vs_sep_kriging_variance', sep = ''), header = FALSE, sep = ",") %>% as.matrix()

pdf(file = paste(root, 'Figures/6-total-kriging-variance.pdf', sep = ''), width = 10, height = 10)

boxplot(DAT_KRIGING_VARIANCE, xlab = "beta", ylab = "Total Kriging Variance", xaxt = 'n', col = rep(c("#28908C", "#F28F20"), 3))

points(1:6, c(rep(312.03, 2), rep(356.29, 2), rep(387.24, 2)), pch = 18, cex = 2, lwd = 2)#, col = '#BEB200')

axis(1, at = c(1.5, 3.5, 5.5), labels = c('0.1', '0.5', '1'))

legend(0.5, 395, legend = c("Separable", "Nonseparable", "True"), pt.cex = c(2, 2, 2), box.lty = 0, inset = .02, col = c("#28908C", "#F28F20", 'black'), pch = c(15, 15, 18))

dev.off()



pdf(file = paste(root, 'Figures/6-MSPE.pdf', sep = ''), width = 10, height = 10)

par(mai=c(1,1.2,1,1))

boxplot(DAT_MSPE, xlab = "", ylab = "", xaxt = 'n', col = rep(c("#28908C", "#F28F20"), 3), yaxt = 'n')
axis(2, cex.axis = 2)
mtext('MSPE', side = 2, line = 3, adj = 0.5, cex = 2.5, col = 4, font = 2)

legend(0.5, 0.0258, legend = c("Separable", "Nonseparable"), pt.cex = c(3, 3), box.lty = 0, inset = .02, col = c("#28908C", "#F28F20"), pch = 15, cex = 2)
axis(1, at = c(1.5, 3.5, 5.5), labels = c(expression(beta==0.1), expression(beta==0.5), expression(beta==1)), cex.axis = 2, tick = F)
dev.off()

################################################                                      			      	################################################
################################################   PLOT REALIZATION FROM A NONSEPARABLE AND SEPARABLE MODEL   	################################################
################################################                                      				################################################

obs1 <- read.table(paste(root, 'Data/synthetic/Z_nonseparable', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
locs1 <- read.table(paste(root, 'Data/synthetic/LOC_nonseparable', sep = ''), header = FALSE, sep = ",") %>% as.matrix()

obs2 <- read.table(paste(root, 'Data/synthetic/Z_separable', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
locs2 <- read.table(paste(root, 'Data/synthetic/LOC_separable', sep = ''), header = FALSE, sep = ",") %>% as.matrix()

zlim_range <- range(c(obs1, obs2))

jpeg(file = paste(root, 'Figures/6-synthetic_data_matern_nonsep_vs_sep.jpg', sep = ''), width = 1500, height = 700)

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
		mtext('Nonseparable', side = 2, line = 7, adj = 0.5, cex = 2.5, col = 4, font = 2)
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
		mtext('Separable', side = 2, line = 7, adj = 0.5, cex = 2.5, col = 4, font = 2)
		axis(2, cex.axis = 2)
	}
	mtext(expression(s[x]), side = 1, line = 4, adj = 0.5,  cex = 2.5)
	axis(1, at = seq(0, 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"), cex.axis = 3, mgp=c(3, 3, 0))
}

screen(2)

x1 <- c(0.025,0.1,0.1,0.025) + 0.1
y1 <- c(0.3,0.3,0.7,0.7)
legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(-3, 3, length.out = 3), 1), CEX = 2)

close.screen( all=TRUE)
dev.off()



################################################                                 ################################################
################################################   PARAMETER ESTIMATION PSWARM   ################################################
################################################                                 ################################################

DAT <- read.table(paste(root, 'Results/parameter-estimates-pswarm.log', sep = ''), header = FALSE, sep = ",") %>% as.matrix()

n <- 6
# width of each boxplot is 0.8
x0s <- 1:n - 0.4
x1s <- 1:n + 0.4
# these are the y-coordinates for the horizontal lines
# that you need to set to the desired values.
y0s <- matrix(c(1, 0.03, 1, 1, 0.5, 0.5,
		1, 0.09, 1, 1, 0.5, 0.5,
		1, 0.23, 1, 1, 0.5, 0.5,
		1, 0.03, 1, 0.48, 0.5, 0.5,
		1, 0.09, 1, 0.48, 0.5, 0.5,
		1, 0.23, 1, 0.48, 0.5, 0.5,
		1, 0.03, 1, 0.24, 0.5, 0.5,
		1, 0.09, 1, 0.24, 0.5, 0.5,
		1, 0.23, 1, 0.24, 0.5, 0.5), nrow = 9, ncol = 6, byrow = T)

ER_lables <- c('Weak', 'Moderate', 'Strong')

pdf(file = paste(root, 'Figures/6-parameter-estimates-pswarm.pdf', sep = ''), width = 8, height = 9)

split.screen( rbind(c(0.12,0.98,0.1,0.95), c(0.93,0.99,0.1,0.95)))
split.screen( figs = c( 3, 3 ), screen = 1 )

bb <- 1
ff <- 1

for(dd in 1:9){
	screen(dd + 2)
	par(pty = 's')
	par(mai = c(0.1, 0.1, 0.1, 0.1))
  
  	if(mod(dd, 3) == 1){
    		boxplot(DAT[(dd - 1) * 20 + 1:20, 1:6], ylim = c(0, 1.5), xaxt = 'n', col = 'white')
    		mtext(ER_lables[bb], side = 2, line = 2.2, cex = 1.2, col = '#0993f3', font = 2)
    		if(bb == 2){
      			mtext('T    I    M    E', side = 2, line = 3.5, cex = 1.5, col = 'blue', font = 2)
    		}
    		bb <- bb + 1
  	}else{
    		boxplot(DAT[(dd - 1) * 20 + 1:20, 1:6], ylim = c(0, 1.5), xaxt = 'n', yaxt = 'n', col = 'white')
  	}
  	segments(x0 = x0s, x1 = x1s, y0 = y0s[dd, ], col = "red", lwd = 2)
  
  	if(dd > 6){
    		axis(1, at = 1:6, labels = c(expression(sigma^2), expression(a[s]), expression(nu), expression(a[t]), expression(alpha), expression(beta)), cex.axis = 1.2)
  	}
  	if(dd < 4){
    		mtext(ER_lables[ff], side = 3, line = 0.2, cex = 1.2, col = '#0993f3', font = 2)
    		if(ff == 2){
      			mtext('S    P    A    C    E', side = 3, line = 1.5, cex = 1.5, col = 'blue', font = 2)
    		}
    		ff <- ff + 1
  	}  
}

#screen(1)

#mtext(expression(bold('N = 400, T = 100, nreps = 30, opt_tol = 5, ' * delta *' = 0, exact matern st kernel')), side = 1, line = 8, adj = 1, cex = 0.6, col = "#FF7F50", font = 4)

close.screen( all=TRUE)
dev.off()








################################################                                        ################################################
################################################   HOW NUMBER OF SWARMS REACH THE MLE   ################################################
################################################                                        ################################################


DAT <- NULL

for(swarm in 1:6){

	# LOAD DATA FOR EACH YEAR, EACH YEAR HAS A 550 x 8760 MATRIX, WHERE COLUMNS ARE MEASUREMENTS IN TIME AND ROWS ARE MEASUREMENTS IN SPACE
	cat("LOADING PSWARM ", swarm * 5, " LIKELIHOOD DATA . . .", '\n')
	DAT_temp <- read.table(paste(root, 'Results/pswarm_', swarm * 5, '_lower_bound', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	DAT <- cbind(DAT, DAT_temp[1:672, ])
}

for(swarm in 1:6){

	# LOAD DATA FOR EACH YEAR, EACH YEAR HAS A 550 x 8760 MATRIX, WHERE COLUMNS ARE MEASUREMENTS IN TIME AND ROWS ARE MEASUREMENTS IN SPACE
	cat("LOADING PSWARM ", swarm * 5, " LIKELIHOOD DATA . . .", '\n')
	DAT_temp <- read.table(paste(root, 'Results/pswarm_', swarm * 5, '_initial_theta', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	DAT <- cbind(DAT, DAT_temp[1:672, ])
}

pdf(file = paste(root, 'Figures/6-pswarm-likelihood-values.pdf', sep = ''), width = 25, height = 10)

par(mfrow = c(1, 2))

plot(DAT[, 1], type = 'l', ylim = c(-100000, -40000))
for(aa in 2:5){
	lines(DAT[, aa], col = aa)
}

plot(DAT[, 6], type = 'l', ylim = c(-100000, -40000))
for(aa in 2:5){
	lines(DAT[, 5 + aa], col = aa)
}

dev.off()


