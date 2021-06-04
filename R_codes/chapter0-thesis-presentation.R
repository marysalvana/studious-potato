
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
