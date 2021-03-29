
directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))

ncname <- paste("/home/salvanmo/Downloads/unidata_TP.nc4", sep='')   

ncin <- nc_open(ncname)

dname1 <- "Base_reflectivity_surface_layer"

u_array <- ncvar_get(ncin, dname1)

# get longitude and latitude
lon <- ncvar_get(ncin, "lon")
lat <- ncvar_get(ncin, "lat")
time <- ncvar_get(ncin, "reftime")

nc_close(ncin)

TT <- dim(u_array)[3]

zlim_range <- range(u_array, na.rm = T)

cat("PLOTTING...", '\n')

for(tt in 1:3){
	jpeg(file = paste(root, 'Figures/4_application_', tt,'.jpg', sep = ''), width = 1500, height = 900)

	split.screen( rbind(c(0.05,0.96,0.1,0.93), c(0.93,0.99,0.1,0.95)))

	screen(1)
	par(mai=c(0.5,0.5,0.5,0.5))
	quilt.plot(c(lon), c(lat), u_array[, , tt], zlim = zlim_range, nx = 1000, ny = 1000, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 3, tick = F, line = 0.5)
	map("state", xlim =  c(-120, -70), ylim = c(30, 50), lwd = 0.75, add = T)

	mtext('Longitude', side = 1, line = 4.5, adj = 0.5,  cex = 3.5, font = 2)
	mtext('Latitude', side = 2, line = 4.5, adj = 0.5, cex = 3.5, font = 2)

	#screen(2)
	#x1 <- c(0.025,0.12,0.12,0.025)
	#y1 <- c(0.25,0.25,0.7,0.7)
	#legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(-zlim_range, zlim_range,length.out = 5), 1), cex = 3)

	close.screen( all=TRUE)
	dev.off()
}

saudi <- map("state", "Virginia", fill = TRUE)

IDs <- sapply(strsplit(saudi$names, ":"), function(x) x[1])
saudi <- map2SpatialPolygons(saudi, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

test1 <- data.frame(c(lon), c(lat), c(u_array[,, 1]), c(u_array[,, 2]), c(u_array[,, 3]), c(u_array[,, 4]), c(u_array[,, 5]))
colnames(test1) <- c('lon', 'lat', 'Z1', 'Z2', 'Z3', 'Z4', 'Z5')
spdf <- SpatialPointsDataFrame(coords = test1[, c("lon", "lat")], data = test1,
			       proj4string = CRS("+proj=longlat +datum=WGS84"))

DAT_SUB <- data.frame(spdf[!is.na(over(spdf, as(saudi, "SpatialPolygons"))), ])

zlim_range <- range(DAT_SUB[, 3:7], na.rm = T)

hr_ind <- 1:5
label <- c('00:00', '00:05', '00:10', '00:15', '00:20')

asratio = 1

jpeg(file = paste(root, 'Figures/4_application.jpg', sep = ''), width = 2000, height = 500)

split.screen( rbind(c(0.08,0.95,0.1,0.95), c(0.95,0.99,0.1,0.95)))
split.screen( figs = c( 1, 5 ), screen = 1 )

hr_count <- 0
for(hr in hr_ind){
	
	hr_count <- hr_count + 1
	
	screen(2 + hr_count)

	par(pty = 's')
	par(mai=c(0.2,0.2,0.2,0.2))
	
	quilt.plot(DAT_SUB[, 1], DAT_SUB[, 2], DAT_SUB[, 2 + hr], zlim = zlim_range, nx = 600, ny = 100, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, axes = F, asp = asratio)

	if(hr_count == 1){
		axis(1, cex.axis = 2)
		axis(2, cex.axis = 2)
		mtext('Base Reflectivity (dBZ)', side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 'blue')
	}else{
		axis(1, cex.axis = 2)
	}
	map("state", xlim =  c(-120, -70), ylim = c(30, 50), lwd = 0.75, add = T)
	
	if(hr_count == 1){
		mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
	}
	mtext(label[hr], side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
	mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
}

screen(2)

x1 <- c(0.025,0.12,0.12,0.025) + 0.1
y1 <- c(0.3,0.3,0.7,0.7)
legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(zlim_range[1], zlim_range[2], length.out = 5), 1), CEX = 2)

close.screen( all=TRUE)
dev.off()



for(tt in 1:3){

	jpeg(file = paste(root, 'Figures/4_application_FULL_', tt,'.jpg', sep = ''), width = 1500, height = 900)

	split.screen( rbind(c(0.05,0.96,0.1,0.93), c(0.93,0.99,0.1,0.95)))

	screen(1)
	par(mai=c(0.5,0.5,0.5,0.5))
	quilt.plot(DAT_SUB[, 1], DAT_SUB[, 2], DAT_SUB[, 3], zlim = zlim_range, nx = 200, ny = 200, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 3, tick = F, line = 0.5)
	map("state", xlim =  c(-120, -70), ylim = c(30, 50), lwd = 0.75, add = T)

	mtext('Longitude', side = 1, line = 4.5, adj = 0.5,  cex = 3.5, font = 2)
	mtext('Latitude', side = 2, line = 4.5, adj = 0.5, cex = 3.5, font = 2)

	#screen(2)
	#x1 <- c(0.025,0.12,0.12,0.025)
	#y1 <- c(0.25,0.25,0.7,0.7)
	#legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(-zlim_range, zlim_range,length.out = 5), 1), cex = 3)

	close.screen( all=TRUE)
	dev.off()

}

FULL <- T

if(!FULL){
	saudi <- map("state", "Virginia", fill = TRUE)

	IDs <- sapply(strsplit(saudi$names, ":"), function(x) x[1])
	saudi <- map2SpatialPolygons(saudi, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))


	for(tt in 1:3){

		test1 <- data.frame(c(lon), c(lat), c(u_array[,, tt]))
		colnames(test1) <- c('lon', 'lat', 'Z1')
		spdf <- SpatialPointsDataFrame(coords = test1[, c("lon", "lat")], data = test1,
					       proj4string = CRS("+proj=longlat +datum=WGS84"))

		DAT_SUB <- data.frame(spdf[!is.na(over(spdf, as(saudi, "SpatialPolygons"))), ])

		zlim_range <- range(DAT_SUB[, 3], na.rm = T)

		jpeg(file = paste(root, 'Figures/4_application_FULL_', tt,'.jpg', sep = ''), width = 1500, height = 900)

		split.screen( rbind(c(0.05,0.96,0.1,0.93), c(0.93,0.99,0.1,0.95)))

		screen(1)
		par(mai=c(0.5,0.5,0.5,0.5))
		quilt.plot(DAT_SUB[, 1], DAT_SUB[, 2], DAT_SUB[, 3], zlim = zlim_range, nx = 200, ny = 200, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 3, tick = F, line = 0.5)
		map("state", xlim =  c(-120, -70), ylim = c(30, 50), lwd = 0.75, add = T)

		mtext('Longitude', side = 1, line = 4.5, adj = 0.5,  cex = 3.5, font = 2)
		mtext('Latitude', side = 2, line = 4.5, adj = 0.5, cex = 3.5, font = 2)

		#screen(2)
		#x1 <- c(0.025,0.12,0.12,0.025)
		#y1 <- c(0.25,0.25,0.7,0.7)
		#legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(-zlim_range, zlim_range,length.out = 5), 1), cex = 3)

		close.screen( all=TRUE)
		dev.off()

	}
}

#Virginia, 2021-02-27T00:00:00Z 
#2021-03-15T00:00:00Z
#2021-03-18T00:00:00Z
