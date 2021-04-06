directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))

# LOAD HOURLY DATA FROM 2016 TO 2019

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


