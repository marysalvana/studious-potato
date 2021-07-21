


source("./pkg-config.R")



AREA <- 'SAUDI'

DATA <- data_format(aggregate = 1, area = AREA, start_yr = 2016, end_yr = 2019)

locs <- as.matrix(DATA[[2]])
locs <- cbind((locs[, 1] - mean(locs[, 1])) / sd(locs[, 1]), (locs[, 2] - mean(locs[, 2])) / sd(locs[, 2]))




hr_ind <- 1:10
label <- c('00:00', '00:05', '00:10', '00:15', '00:20', '00:25', '00:30', '00:35', '00:40', '00:45')


zlim_range1 <- range(DATA[[1]][hr_ind, ], na.rm = T)

jpeg(file = paste(root, 'Figures/5-application.jpg', sep = ''), width = 1600, height = 800)

split.screen( rbind(c(0.05,0.95,0.12,0.92), c(0.95,0.99,0.12,0.92)))
split.screen( figs = c( 2, 5 ), screen = 1 )

hr_count <- 0
for(hr in hr_ind){
	
	hr_count <- hr_count + 1
	
	screen(2 + hr_count)

	par(pty = 's')
	#par(pin=c(6, 1.5))
	par(mai=c(0.3,0.3,0.3,0.6))
	
	if(hr_count %in% c(1, 6)){
		quilt.plot(DATA[["locations"]][, 1], DATA[["locations"]][, 2], DATA[["data_matrix"]][hr, ], zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, xaxt = 'n')
	}else{
		quilt.plot(DATA[["locations"]][, 1], DATA[["locations"]][, 2], DATA[["data_matrix"]][hr, ], zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, xaxt = 'n', yaxt = 'n')
	}

	map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
	
	if(hr_count %in% c(1, 6)){
		mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
		axis(2, cex.axis = 2)
	}

	mtext(paste(hr - 1, ":00", sep = ""), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)

	if(hr > 5){
		mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
		axis(1, cex.axis = 2)
	}
}

screen(2)

x1 <- c(0.025,0.1,0.1,0.025) + 0.1
y1 <- c(0.3,0.3,0.7,0.7)
legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(-3, 3, length.out = 3), 1), CEX = 2)

close.screen( all=TRUE)
dev.off()


