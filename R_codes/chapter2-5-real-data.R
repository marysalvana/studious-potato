## INPUT: N x T matrix of log PM2.5 concentrations, where N is the number of spatial locations and T is the number of temporal locations 
## INPUT: N x 2 matrix of locations containing longitude and latitude

## OUTPUT: textfile of large training and testing datasets (measurements, spatial locations, temporal locations)

directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))


AREA <- 'SAUDI'

DATA <- data_format(aggregate = 1, area = AREA)

####################################################################################################

plot_realdata_for_checking_stationarity(data_list = DATA, file_name = '2-application-saudi-data-scratch.jpg', start_hr = 8760 + 1, saudi = T)
plot_realdata_for_manuscript(data_list = DATA, file_name = '2-application-saudi-data-scratch.jpg', start_hr = 8760 + 1, saudi = T)

####################################################################################################


