## INPUT: netcdf dataset downloaded from https://disc.gsfc.nasa.gov/datasets/M2T1NXAER_5.12.4/summary

## OUTPUT: N x T matrix of log PM2.5 concentrations, where N is the number of spatial locations and T is the number of temporal locations 
## OUTPUT: N x 2 matrix of locations containing longitude and latitude

directory <- '/home/salvanmo/Desktop/'

area <- 'SAUDI'

if(area == 'US'){
	data_directory <- '/media/salvanmo/yourlainess/phd/data/sc21/US/'
	VARIABLE_NAME = "DUSMASS25"
	SUBSET = NULL
	REGION = NULL
	GET_SUBSET = F
	NCDF_EXTENSION = ".SUB.nc"
}else if(area == 'SAUDI'){
	data_directory <- '/media/salvanmo/yourlainess/phd/data/sc21/SAUDI/'
	VARIABLE_NAME = "DUCMASS25"
	SUBSET = "Saudi"
	REGION = "world"
	GET_SUBSET = T
	NCDF_EXTENSION = ".nc4.nc4"
}



source("./pkg-config.R")



extract_data <- function(yr, variable_name, get_subset = F, subset = NULL, region = NULL){

	#subset = "Saudi", "Colorado", any other country name when region == world, any US state name when region == state
	#region = "world", "state"

	if(yr < 1992){
		merra_ind <- 100
	}else if(yr >= 1992 & yr < 2001){
		merra_ind <- 200
	}else if (yr >= 2001 & yr < 2011){
		merra_ind <- 300
	}else{
		merra_ind <- 400
	}

	CONSO_DATA <- NULL

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

		for(day in 1:mnth_end){

			cat('READING NETCDF DATA ===>', '\t', 'year: ', yr, '\t', ' month: ', mnth, '\t', 'day: ', day, '\n')

			if(day > 9){
				ncname <- paste(data_directory, "MERRA2_", merra_ind, ".tavg1_2d_aer_Nx.", yr, mo, day, NCDF_EXTENSION, sep='')
			}else{
				ncname <- paste(data_directory, "MERRA2_", merra_ind, ".tavg1_2d_aer_Nx.", yr, mo, "0",day, NCDF_EXTENSION, sep='')
			}
			
			ncin <- nc_open(ncname)

			u_array <- ncvar_get(ncin, variable_name)

			# get longitude and latitude
			lon <- ncvar_get(ncin, "lon")
			lat <- ncvar_get(ncin, "lat")
			lon.lat <- expand.grid(lon, lat)

			nc_close(ncin)

			MEAN <- mean(log(u_array))

			for(tt in 1:dim(u_array)[3]){

				###########################   GET ONLY DATA OVER SUBSET   ###########################

				if(area == 'US'){
					Y <- cbind(lon.lat, c(log(u_array[, , tt])))
				}else if(area == 'SAUDI'){
					Y <- cbind(lon.lat, c(log(u_array[, , tt])) - MEAN)
				}
				colnames(Y) <- c('lon', 'lat', 'Y1')
				if(get_subset){
					subregion <- map(region, subset, fill = TRUE)
					IDs <- sapply(strsplit(subregion$names, ":"), function(x) x[1])
					subregion <- map2SpatialPolygons(subregion, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

					spdf <- SpatialPointsDataFrame(coords = Y[, c("lon", "lat")], data = Y, proj4string = CRS("+proj=longlat +datum=WGS84"))

					subregion_final <- data.frame(spdf[!is.na(over(spdf, as(subregion, "SpatialPolygons"))), ])
					CONSO_DATA <- cbind(CONSO_DATA, subregion_final$Y1)
					locations <- subregion_final[, 1:2]
				}else{
					CONSO_DATA <- cbind(CONSO_DATA, Y$Y1)
					locations <- cbind(Y$lon, Y$lat)
				}
				N <- nrow(locations)
			}
		}	
	}	

	FINAL_DATA <- list("log_measurements" = CONSO_DATA, "locations" = locations)
	
	return(FINAL_DATA)
}

## Indicate (1) the variable name and (2) the year to which you want to get the data.

for(YEAR in 1980:2009){

	data_matrix <- extract_data(yr = YEAR, variable_name = VARIABLE_NAME, get_subset = GET_SUBSET, subset = SUBSET, region = REGION)

	if(area == 'US'){
		write.table(data_matrix[["locations"]], file = paste(root, "Data/sc21/locations_US_", YEAR, sep = ""), sep = ",", row.names = FALSE, col.names = FALSE)
		write.table(data_matrix[["log_measurements"]], file = paste(root, "Data/sc21/pm_US_", YEAR, sep = ""), sep = ",", row.names = FALSE, col.names = FALSE)
	}else if(area == 'SAUDI'){
		write.table(data_matrix[["locations"]], file = paste(root, "Data/sc21/locations_", YEAR, sep = ""), sep = ",", row.names = FALSE, col.names = FALSE)
		write.table(data_matrix[["log_measurements"]], file = paste(root, "Data/sc21/pm_", YEAR, sep = ""), sep = ",", row.names = FALSE, col.names = FALSE)
	}
}

cat(paste("Data is saved in: ", root, "Data/sc21/", sep = ""), '\n')



