
directory <- '/home/salvanmo/Desktop/'

data_directory <- '/media/salvanmo/yourlainess/phd/data/motivation/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))

dname1 <- "DUSMASS25"
dname2 <- "BCSMASS"

saudi<- map("world", "Saudi", fill = TRUE)
IDs <- sapply(strsplit(saudi$names, ":"), function(x) x[1])
saudi <- map2SpatialPolygons(saudi, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

for(yr in 1980:1980){

	data_array1 <- data_array2 <- data_array3 <- NULL

	if(yr < 1992){
		merra_ind <- 100
	}else if(yr >= 1992 & yr < 2001){
		merra_ind <- 200
	}else if (yr >= 2001 & yr < 2011){
		merra_ind <- 300
	}else{
		merra_ind <- 400
	}

	for(mnth in 1:12){
		if(mnth == 2){
			mnth_end <- 28
		}else if(mnth %in% c(1, 3, 5, 7, 8, 10, 12)){
			mnth_end <- 31
		}else{
			mnth_end <- 30
		}

		if(mnth < 10){
			mo <- paste("0", mnth, sep='')
		}else{
			mo <- mnth
		}
		ncname <- paste(data_directory, "MERRA2_", merra_ind, ".tavg1_2d_aer_Nx.", yr, mo, "01.SUB.nc", sep='')

		ncin <- nc_open(ncname)

		u_array <- ncvar_get(ncin,dname1)
		v_array <- ncvar_get(ncin,dname2)

		# get longitude and latitude
		lon <- ncvar_get(ncin,"lon")
		lat <- ncvar_get(ncin,"lat")
		lon.lat <- expand.grid(lon, lat)

		nc_close(ncin)

		U <- u_array
		V <- v_array

		test1 <- data.frame(rep(lon.lat[,1], 24), rep(lon.lat[,2], 24), c(U), c(V))

		for(day in 2:mnth_end){
			cat('READING NETCDF DATA ===> year: ', yr, 'month: ', mnth, 'day: ', day, '\n')
			if(day > 9){
				ncname <- paste(data_directory, "MERRA2_", merra_ind, ".tavg1_2d_aer_Nx.", yr, mo, day,".SUB.nc", sep='')
			}else{
				ncname <- paste(data_directory, "MERRA2_", merra_ind, ".tavg1_2d_aer_Nx.", yr, mo, "0",day,".SUB.nc", sep='')
			}
			ncin <- nc_open(ncname)

			u_array <- ncvar_get(ncin,dname1)
			v_array <- ncvar_get(ncin,dname2)

			nc_close(ncin)

			U <- u_array
			V <- v_array

			test1 <- rbind(test1, data.frame(rep(lon.lat[,1], 24), rep(lon.lat[,2], 24),  c(U),  c(V)))
		}
		colnames(test1) <- c('lon', 'lat', 'Y1', 'Y2')
		spdf <- SpatialPointsDataFrame(coords = test1[, c("lon", "lat")], data = test1, proj4string = CRS("+proj=longlat +datum=WGS84"))
		saudi_data_orig <- data.frame(spdf[!is.na(over(spdf, as(saudi, "SpatialPolygons"))), ])

		data_temp <- matrix(log(saudi_data_orig[, 3]), ncol = 550, byrow = T)
		data_array1 <- rbind(data_array1, data_temp)

		data_temp <- matrix(log(saudi_data_orig[, 4]), ncol = 550, byrow = T)
		data_array2 <- rbind(data_array2, data_temp)


	}

	data_array3 <- saudi_data_orig[1:550, 1:2]
		
	write.table(data_array1, file = paste(root, "Data/ncdf/DUSMASS25_", yr, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
	write.table(data_array2, file = paste(root, "Data/ncdf/BCSMASS_", yr, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
	write.table(data_array3, file = paste(root, "Data/ncdf/LOCS", sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)

}
