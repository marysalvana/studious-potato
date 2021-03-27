
directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))

ncname <- paste("/home/salvanmo/Downloads/unidata_TP (7).nc4", sep='')   

ncin <- nc_open(ncname)

dname1 <- "Base_reflectivity_surface_layer"

u_array <- ncvar_get(ncin, dname1)

# get longitude and latitude
lon <- ncvar_get(ncin, "lon")
lat <- ncvar_get(ncin, "lat")
time <- ncvar_get(ncin, "reftime")

nc_close(ncin)

zlim_range <- range(u_array, na.rm = T)

for(tt in 1:11){
	jpeg(file = paste(root, 'Figures/application_', tt,'.jpg', sep = ''), width = 1950, height = 1000)

	quilt.plot(c(lon), c(lat), u_array[,, tt], zlim = zlim_range, nx = 300, ny = 300, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 3, yaxt = 'n', tick = F, line = 0.5)

	dev.off()
}

