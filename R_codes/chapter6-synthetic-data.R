directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))


################################################                                      ################################################
################################################   NONSEPARABLE VS SEPARABLE MODELS   ################################################
################################################                                      ################################################

DAT_MSPE <- read.table(paste(root, 'Results/nonsep_vs_sep_mspe', sep = ''), header = FALSE, sep = ",") %>% as.matrix()
DAT_KRIGING_VARIANCE <- read.table(paste(root, 'Results/nonsep_vs_sep_kriging_variance', sep = ''), header = FALSE, sep = ",") %>% as.matrix()

pdf(file = paste(root, 'Figures/6-total-kriging-variance.pdf', sep = ''), width = 10, height = 10)

boxplot(DAT_KRIGING_VARIANCE, xlab = "beta", ylab = "Total Kriging Variance", xaxt = 'n', col = rep(c("#28908C", "#F28F20"), 3))

points(1:6, c(rep(true_mspe_var_0.1[2], 2), rep(true_mspe_var_0.5[2], 2), rep(true_mspe_var_1[2], 2)), pch = 18, cex = 2, lwd = 2)#, col = '#BEB200')

axis(1, at = c(1.5, 3.5, 5.5), labels = c('0.1', '0.5', '1'))

legend(0.5, 395, legend = c("Separable", "Nonseparable", "True"), pt.cex = c(2, 2, 2), box.lty = 0, inset = .02, col = c("#28908C", "#F28F20", 'black'), pch = c(15, 15, 18))

dev.off()



pdf(file = paste(root, 'Figures/6-MSPE.pdf', sep = ''), width = 10, height = 10)

boxplot(DAT_MSPE, xlab = expression(beta), ylab = "MSPE", xaxt = 'n', col = rep(c("#28908C", "#F28F20"), 3))

legend(0.5, 0.0258, legend = c("Separable", "Nonseparable"), pt.cex = c(2, 2), box.lty = 0, inset = .02, col = c("#28908C", "#F28F20"), pch = 15)
axis(1, at = c(1.5, 3.5, 5.5), labels = c('0.1', '0.5', '1'))

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
    		boxplot(DAT[(dd - 1) * 20 + 1:20, 1:6], ylim = c(0, 1.5), xaxt = 'n')
    		mtext(ER_lables[bb], side = 2, line = 2.5, cex = 0.8, col = '#0993f3', font = 2)
    		if(bb == 2){
      			mtext('T    I    M    E', side = 2, line = 3.5, cex = 1.5, col = 'blue', font = 2)
    		}
    		bb <- bb + 1
  	}else{
    		boxplot(DAT[(dd - 1) * 20 + 1:20, 1:6], ylim = c(0, 1.5), xaxt = 'n', yaxt = 'n')
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


