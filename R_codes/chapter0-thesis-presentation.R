
directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R",sep=''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R",sep=''))


yr <- 2018

for(variable in 1:2){

	if(variable == 1){
		dat <- read.table(paste(root, 'Data/motivation/DUSMASS25_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	}else{
		dat <- read.table(paste(root, 'Data/motivation/BCSMASS_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	}

	dat3 <- read.table(paste(root, 'Data/motivation/LOCS', sep = ''), header = FALSE, sep = " ") %>% as.matrix()


	start_hr <- 1
	subset_ind <- start_hr:(start_hr + 23)

	zlim_range1 <- range(dat[subset_ind,])

	for(hr in subset_ind){

		jpeg(file = paste(root, 'Figures/0-spacetime-maps_variable_', variable, '_t', hr, '.jpg', sep = ''), width = 1000, height = 1000)
		
		split.screen( rbind(c(0.05,0.95,0.1,0.95), c(0.90,0.99,0.1,0.95)))

		screen(1)

		par(pty = 's')
		par(mai=c(0.5, 0.5, 0.5, 0.5))
		
		quilt.plot(dat3[, 1], dat3[, 2], dat[hr, ], zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2)
		map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
		
		mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
		mtext(paste(hr - 1, ':00', sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
		mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)

		screen(2)

		x1 <- c(0.01,0.1,0.1,0.01)
		y1 <- c(0.2,0.2,0.8,0.8)
		legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(round(zlim_range1[1], 0), round(zlim_range1[2], 0), length.out = 5), 1), CEX = 2)

		close.screen( all=TRUE)
		dev.off()

	}

}

#################################################################################################################

cov_example_1 <- read.table(paste(root, 'Data/univariate-nonstationary/cov-example-1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
realizations_example_1 <- read.table(paste(root, 'Data/univariate-nonstationary/realizations-example-1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()


zlim_range1 <- c(0, 1)
zlim_range2 <- range(realizations_example_1[1, 1:(n * 5)])

jpeg(file = paste(root, 'Figures/0-scratch-cov1-heatmap.jpg', sep = ''), width = 1200, height = 600)

split.screen( rbind(c(0.06,0.94,0.08,0.93), c(0.94,0.98,0.08,0.93)))
split.screen( figs = c( 3, 6 ), screen = 1 )


hr_count <- 0
for(tt in 1:5){
	
	hr_count <- hr_count + 1
	
	screen(3 + hr_count)

	par(pty = 's')
	par(mai=c(0.2,0.2,0.2,0.2))
	
	if(tt == 1){
	quilt.plot(sim_grid_locations[, 1], sim_grid_locations[, 2], realizations_example_1[1, (tt - 1) * n + 1:n], zlim = zlim_range2, nx = N, ny = N, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 1, xaxt = 'n')
	}else{
	quilt.plot(sim_grid_locations[, 1], sim_grid_locations[, 2], realizations_example_1[1, (tt - 1) * n + 1:n], zlim = zlim_range2, nx = N, ny = N, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 1, xaxt = 'n', yaxt = 'n')
	}
	mtext(paste('t = ', tt, sep = ''), side = 3, line = 1, adj = 0.5, cex = 2, font = 2)
}	

hr_count <- hr_count + 1
screen(3 + hr_count)

par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

quilt.plot(sim_grid_locations[, 1], sim_grid_locations[, 2], realizations_example_1[1, 1:n], zlim = zlim_range2, nx = N, ny = N, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 1)
points(matrix(sim_grid_locations[reference_locations[1], ], ncol = 2), col = 'black', pch = 4, cex = 3, lwd = 4)
mtext(paste('Ref Loc 1', sep = ''), side = 2, line = 4, adj = 0.5, cex = 2, font = 2, col = 'blue')

for(tt in 1:5){
	
	hr_count <- hr_count + 1
	
	screen(3 + hr_count)

	par(pty = 's')
	par(mai=c(0.2,0.2,0.2,0.2))
	
	quilt.plot(sim_grid_locations[, 1], sim_grid_locations[, 2], cov_example_1[1, (tt - 1) * n + 1:n], zlim = zlim_range1, nx = N, ny = N, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 1, yaxt = 'n', xaxt = 'n')

}	

hr_count <- hr_count + 1
screen(3 + hr_count)

par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

quilt.plot(sim_grid_locations[, 1], sim_grid_locations[, 2], realizations_example_1[1, 1:n], zlim = zlim_range2, nx = N, ny = N, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 1)
points(matrix(sim_grid_locations[reference_locations[2], ], ncol = 2), col = 'black', pch = 4, cex = 3, lwd = 4)
mtext(paste('Ref Loc 2', sep = ''), side = 2, line = 4, adj = 0.5, cex = 2, font = 2, col = 'blue')

mtext(paste('t = ', 1, sep = ''), side = 2, line = 2, adj = 1.5, cex = 2, font = 2)

for(tt in 1:5){
	
	hr_count <- hr_count + 1
	
	screen(3 + hr_count)

	par(pty = 's')
	par(mai=c(0.2,0.2,0.2,0.2))
	
	quilt.plot(sim_grid_locations[, 1], sim_grid_locations[, 2], cov_example_1[2, (tt - 1) * n + 1:n], zlim = zlim_range1, nx = N, ny = N, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 1, yaxt = 'n')

}	

screen(2)

x1 <- c(0.025,0.1,0.1,0.025) + 0.1
y1 <- c(0.17,0.17,0.45,0.45)
legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(0, 1, length.out = 3), 1), CEX = 1.5)

close.screen( all=TRUE)
dev.off()


