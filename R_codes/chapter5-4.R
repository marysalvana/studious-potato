
source("./pkg-config.R")

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

cat("PLOTTING...", '\n')

saudi <- map("state", "Virginia", fill = TRUE)

IDs <- sapply(strsplit(saudi$names, ":"), function(x) x[1])
saudi <- map2SpatialPolygons(saudi, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

test1 <- data.frame(c(lon), c(lat), c(u_array[,, 1]), c(u_array[,, 2]), c(u_array[,, 3]), c(u_array[,, 4]), c(u_array[,, 5]), c(u_array[,, 6]), c(u_array[,, 7]), c(u_array[,, 8]), c(u_array[,, 9]), c(u_array[,, 10]))
colnames(test1) <- c('lon', 'lat', 'Z1', 'Z2', 'Z3', 'Z4', 'Z5', 'Z6', 'Z7', 'Z8', 'Z9', 'Z10')
spdf <- SpatialPointsDataFrame(coords = test1[, c("lon", "lat")], data = test1,
			       proj4string = CRS("+proj=longlat +datum=WGS84"))

DAT_SUB <- data.frame(spdf[!is.na(over(spdf, as(saudi, "SpatialPolygons"))), ])

zlim_range <- range(DAT_SUB[, 3:12], na.rm = T)

hr_ind <- 1:10
label <- c('00:00', '00:05', '00:10', '00:15', '00:20', '00:25', '00:30', '00:35', '00:40', '00:45')

jpeg(file = paste(root, 'Figures/4_application.jpg', sep = ''), width = 1800, height = 500)
#pdf(file = paste(root, 'Figures/4_application.pdf', sep = ''), width = 27, height = 7)

split.screen( rbind(c(0.05,0.95,0.12,0.92), c(0.95,0.99,0.12,0.92)))
split.screen( figs = c( 2, 5 ), screen = 1 )

hr_count <- 0
for(hr in hr_ind){
	
	hr_count <- hr_count + 1
	
	screen(2 + hr_count)

	#par(pty = 's')
	par(pin=c(6, 1.5))
	par(mai=c(0.5,0.5,0.5,0.6))
	
	quilt.plot(DAT_SUB[, 1], DAT_SUB[, 2], DAT_SUB[, 2 + hr], zlim = zlim_range, nx = 900, ny = 150, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, axes = F, xlim = c(-82, -76))

	#mtext('Base Reflectivity (dBZ)', side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 'blue')
	map("state", xlim =  c(-120, -70), ylim = c(30, 50), lwd = 0.75, add = T)
	
	if(hr_count %in% c(1, 6)){
		mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
		axis(2, at = c(37, 38, 39), cex.axis = 2)
	}
	mtext(label[hr], side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
	if(hr > 5){
		mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
		axis(1, cex.axis = 2)
	}
}

screen(2)

x1 <- c(0.025,0.1,0.1,0.025) + 0.1
y1 <- c(0.3,0.3,0.7,0.7)
legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(zlim_range[1], zlim_range[2], length.out = 3), 1), CEX = 2)

close.screen( all=TRUE)
dev.off()


#Virginia, 2021-02-27T00:00:00Z 
#2021-03-15T00:00:00Z
#2021-03-18T00:00:00Z
