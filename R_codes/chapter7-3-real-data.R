## INPUT: N x T matrix of log PM2.5 concentrations, where N is the number of spatial locations and T is the number of temporal locations 
## INPUT: N x 2 matrix of locations containing longitude and latitude

## OUTPUT: textfile of large training and testing datasets (measurements, spatial locations, temporal locations)

directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))


AREA <- 'SAUDI'

DATA <- data_format(aggregate = 48, area = AREA)


# Save the new testing file for Saudi from January 2020.

area <- 'SAUDI'

data_format_saudi_testing_only(aggregate = 48, forward_time_predict = 1)

# Save the new testing file for US from 2016.

AREA <- 'US'

DATA <- data_format(aggregate = 1, area = AREA)
data_format_US_testing_only(data_list = DATA, temporal_length = 500, forward_time_predict = 2, file_name = '_US')


####################################################################################################

plot_realdata_for_checking_stationarity(data_list = DATA, file_name = '6-application-US-data-scratch.jpg', start_hr = 1, saudi = F)

data_split(data_list = DATA, training_percent = 0.9, temporal_length = 500, file_name = '_US')

plot_realdata_for_manuscript(data_list = DATA, file_name = '6-application-US-data.jpg', start_hr = 1, saudi = F, aggregate = 'hourly')

####################################################################################################

results_list <- list()

#Saudi insample
sep_pred <- read.table(paste(root, "Data/sc21/predicted_values_May_29_2021_04:14:28", sep = ''), header = FALSE, sep = " ") %>% as.matrix()
nonsep_pred <- read.table(paste(root, "Data/sc21/predicted_values_May_29_2021_04:20:21", sep = ''), header = FALSE, sep = " ") %>% as.matrix()
nonsep_crippled_pred <- read.table(paste(root, "Data/sc21/predicted_values_May_29_2021_04:20:02", sep = ''), header = FALSE, sep = " ") %>% as.matrix()
obs_test <- read.table(paste(root, "Data/sc21/data_st_testing_FULL", sep = ''), header = FALSE, sep = " ") %>% as.matrix()

results_list[[1]] <- cbind(obs_test, sep_pred, nonsep_pred, nonsep_crippled_pred)

#Saudi forward
sep_pred <- read.table(paste(root, "Data/sc21/predicted_values_May_29_2021_04:25:20", sep = ''), header = FALSE, sep = " ") %>% as.matrix()
nonsep_pred <- read.table(paste(root, "Data/sc21/predicted_values_May_29_2021_04:26:11", sep = ''), header = FALSE, sep = " ") %>% as.matrix()
nonsep_crippled_pred <- read.table(paste(root, "Data/sc21/predicted_values_May_29_2021_04:26:52", sep = ''), header = FALSE, sep = " ") %>% as.matrix()
obs_test <- read.table(paste(root, "Data/sc21/data_st_testing_NEW", sep = ''), header = FALSE, sep = " ") %>% as.matrix()

results_list[[2]] <- cbind(obs_test, sep_pred, nonsep_pred, nonsep_crippled_pred)

#US insample
sep_pred <- read.table(paste(root, "Data/sc21/predicted_values_May_29_2021_03:53:26", sep = ''), header = FALSE, sep = " ") %>% as.matrix()
nonsep_pred <- read.table(paste(root, "Data/sc21/predicted_values_May_29_2021_03:54:20", sep = ''), header = FALSE, sep = " ") %>% as.matrix()
nonsep_crippled_pred <- read.table(paste(root, "Data/sc21/predicted_values_May_29_2021_04:01:12", sep = ''), header = FALSE, sep = " ") %>% as.matrix()
obs_test <- read.table(paste(root, "Data/sc21/data_st_testing_FULL_US", sep = ''), header = FALSE, sep = " ") %>% as.matrix()

results_list[[3]] <- cbind(obs_test, sep_pred, nonsep_pred, nonsep_crippled_pred)

#US forward
sep_pred <- read.table(paste(root, "Data/sc21/predicted_values_May_29_2021_04:03:17", sep = ''), header = FALSE, sep = " ") %>% as.matrix()
nonsep_pred <- read.table(paste(root, "Data/sc21/predicted_values_May_29_2021_04:07:35", sep = ''), header = FALSE, sep = " ") %>% as.matrix()
nonsep_crippled_pred <- read.table(paste(root, "Data/sc21/predicted_values_May_29_2021_04:13:33", sep = ''), header = FALSE, sep = " ") %>% as.matrix()
obs_test <- read.table(paste(root, "Data/sc21/data_st_testing_NEW_US", sep = ''), header = FALSE, sep = " ") %>% as.matrix()

results_list[[4]] <- cbind(obs_test, sep_pred, nonsep_pred, nonsep_crippled_pred)


MSPE_BOXPLOTS <- list()

for(jj in 1:length(results_list)){
	mspe_vals_sep <- mspe_vals_nonsep <- c()
	for(ii in 1:100){
		set.seed(ii)
		subset <- sample(1:nrow(results_list[[jj]]), nrow(results_list[[jj]]) * 0.25)
		mspe_vals_sep[ii] <- mean((results_list[[jj]][subset, 1] - results_list[[jj]][subset, 2])^2)
		mspe_vals_nonsep[ii] <- mean((results_list[[jj]][subset, 1] - results_list[[jj]][subset, 3])^2)
	}
	MSPE_BOXPLOTS[[jj]] <- cbind(mspe_vals_sep, mspe_vals_nonsep)
}

dataset_label <- rep(c('Testing Dataset 1', 'Testing Dataset 2'), 2)
country_label <- c('Saudi Arabia', 'Saudi Arabia', 'US', 'US')
label_adj <- c(2.5, 2.5, 1.25, 1.25)

pdf(file = paste(root, 'Figures/6-MSPE-real-data.pdf', sep = ''), width = 20, height = 8)

split.screen( rbind(c(0.04,0.99,0.06,0.85), c(0.99,0.99,0.1,0.95)))
split.screen( figs = c( 1, 4 ), screen = 1 )

for(jj in 1:length(MSPE_BOXPLOTS)){
	screen(jj + 2)
	par(mai=c(0.5, 0.5, 0.5, 0.5))

	boxplot(MSPE_BOXPLOTS[[jj]], xlab = "", ylab = "", xaxt = 'n', col = rep(c("#28908C", "#F28F20"), 3), yaxt = 'n')
	axis(2, cex.axis = 2)
	mtext(dataset_label[jj], side = 3, line = 1, adj = 0.5, cex = 2)
	if(jj %in% c(1, 3)){
		mtext(country_label[jj], side = 3, line = 4, adj = label_adj[jj], cex = 2.5, col = 4, font = 2)
	}
	if(jj == 1){
		mtext('MSPE', side = 2, line = 3.5, adj = 0.5, cex = 2.5, col = 4, font = 2)
	}

	axis(1, at = 1:2, labels = c("Separable", "Nonseparable"), cex.axis = 1.5, tick = F)
}

close.screen( all=TRUE)
dev.off()

################################################################################################

obs_training <- read.table(paste(root, "Data/sc21/data_st_training_FULL", sep = ''), header = FALSE, sep = " ") %>% as.matrix()
obs_training <- read.table(paste(root, "Data/sc21/data_st_training_FULL_US", sep = ''), header = FALSE, sep = " ") %>% as.matrix()

var(obs_training)



obs_test <- read.table(paste(root, "Data/sc21/data_st_testing_FULL", file_name, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
locs_s_test <- read.table(paste(root, "Data/sc21/locations_space_testing_FULL", file_name, sep = ''), header = FALSE, sep = ",") %>% as.matrix()
locs_t_test <- read.table(paste(root, "Data/sc21/locations_time_testing_FULL", file_name, sep = ''), header = FALSE, sep = " ") %>% as.matrix()

obs_train <- read.table(paste(root, "Data/sc21/data_st_training_FULL", file_name, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
locs_s_train <- read.table(paste(root, "Data/sc21/locations_space_training_FULL", file_name, sep = ''), header = FALSE, sep = ",") %>% as.matrix()
locs_t_train <- read.table(paste(root, "Data/sc21/locations_time_training_FULL", file_name, sep = ''), header = FALSE, sep = " ") %>% as.matrix()



# PLOTTING REAL DATA WITH PREDICTIONS

obs_sub_pred <- read.table(paste(root, 'Results/6-predicted_values_FULL', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

start_hr <- 1
zlim_range1 <- range(c(res_mat1[start_hr:(start_hr + 4),], obs_sub_pred))

jpeg(file = paste(root, 'Figures/6-application_data.jpg', sep = ''), width = 2000, height = 1000)

split.screen( rbind(c(0.08,0.94,0.08,0.88), c(0.94,0.98,0.08,0.88)))
split.screen( figs = c( 2, 5 ), screen = 1 )

hr_count <- 0
for(hr in start_hr:(start_hr + 4)){
	
	ind_train <- which(locs_t_sub_train == hr)
	ind_test <- which(locs_t_sub_test == hr)

	hr_count <- hr_count + 1
	
	screen(2 + hr_count)

	par(pty = 's')
# PLOTTING REAL DATA WITH PREDICTIONS

obs_sub_pred <- read.table(paste(root, 'Results/6-predicted_values_FULL', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

start_hr <- 1
zlim_range1 <- range(c(res_mat1[start_hr:(start_hr + 4),], obs_sub_pred))

jpeg(file = paste(root, 'Figures/6-application_data.jpg', sep = ''), width = 2000, height = 1000)

split.screen( rbind(c(0.08,0.94,0.08,0.88), c(0.94,0.98,0.08,0.88)))
split.screen( figs = c( 2, 5 ), screen = 1 )

hr_count <- 0
for(hr in start_hr:(start_hr + 4)){
	
	ind_train <- which(locs_t_sub_train == hr)
	ind_test <- which(locs_t_sub_test == hr)

	hr_count <- hr_count + 1
	
	screen(2 + hr_count)

	par(pty = 's')
	par(mai=c(0.2,0.2,0.2,0.2))
	
	if(hr == 1){
	quilt.plot(c(locs_s_sub_train[ind_train, 1], locs_s_sub_test[ind_test, 1]), c(locs_s_sub_train[ind_train, 2], locs_s_sub_test[ind_test, 2]), c(obs_sub_train[ind_train], obs_sub_test[ind_test]), zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, xaxt = 'n')
	mtext('Real Data', side = 2, line = 7, adj = 0.5, cex = 2.5, font = 2, col = 'blue')
	}else{
	quilt.plot(c(locs_s_sub_train[ind_train, 1], locs_s_sub_test[ind_test, 1]), c(locs_s_sub_train[ind_train, 2], locs_s_sub_test[ind_test, 2]), c(obs_sub_train[ind_train], obs_sub_test[ind_test]), zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2, xaxt = 'n')
	}
	points(locs_s_sub_test[ind_test, ], col = 'black', pch = 4, cex = 1, lwd = 2)
	map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
	
	if(hr == 1){
		mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
	}

	if(aggregate == 'hourly'){
		mtext(paste(hr - 1 - floor(hr/24) * 24, ':00', sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
	}else if(aggregate == 'daily'){
		mtext(paste('January ', hr, ', ', yr, sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
	}else{
		mtext(paste('January ', 2 * hr_count - 1, '-', 2 * hr_count, sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
	}
	if(hr == 3) mtext('Mean log PM 2.5 Concentration for the Period', side = 3, line = 7, cex = 3, font = 2, col = 4)
}

for(hr in start_hr:(start_hr + 4)){
	
	hr_count <- hr_count + 1

	ind_train <- which(locs_t_sub_train == hr)
	ind_test <- which(locs_t_sub_test == hr)

	screen(2 + hr_count)

	par(pty = 's')
	par(mai=c(0.2,0.2,0.2,0.2))
	
	if(hr == 1){
	quilt.plot(c(locs_s_sub_train[ind_train, 1], locs_s_sub_test[ind_test, 1]), c(locs_s_sub_train[ind_train, 2], locs_s_sub_test[ind_test, 2]), c(obs_sub_train[ind_train], obs_sub_pred[ind_test]), zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2)
	mtext('Training & Predictions', side = 2, line = 7, adj = 0.5, cex = 2.5, font = 2, col = 'blue')
	}else{
	quilt.plot(c(locs_s_sub_train[ind_train, 1], locs_s_sub_test[ind_test, 1]), c(locs_s_sub_train[ind_train, 2], locs_s_sub_test[ind_test, 2]), c(obs_sub_train[ind_train], obs_sub_pred[ind_test]), zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
	}
	points(locs_s_sub_test[ind_test, ], col = 'black', pch = 4, cex = 1, lwd = 2)
	map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
	
	if(hr == 1){
		mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
	}

	mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
}
screen(2)

x1 <- c(0.025,0.1,0.1,0.025) + 0.1
y1 <- c(0.3,0.3,0.7,0.7)
legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(-5, 5, length.out = 3), 1), CEX = 3)

close.screen( all=TRUE)
dev.off()

################################################                                      ################################################
################################################        SUBSET OF FULL DATASET        ################################################
################################################                                      ################################################


# CHOOSE ONLY A SMALL SAMPLE WITH 100,000 SPACE-TIME MEASUREMENTS

small_scale_n <- 100000

set.seed(1234)
subset_index <- sample(1:nrow(Z), small_scale_n)

locs_sub <- Z[subset_index, 2:3]
time_sub <- Z[subset_index, 4]
Z_sub <- Z[subset_index, 1]

# DIVIDE THE 100,000 DATA INTO 10 SETS OF 10, FOR 10-FOLD CROSS VALIDATION

cv_n_test = 10000	# set the number or testing locations for the k-fold cross validation

for(set in 1:10){

	set_ind <- (set - 1) * cv_n_test + 1:cv_n_test
	
	# SAVING TESTING DATA

	cat("SAVING TESTING DATA FOR SET ", set, "\n")

	write.table(locs_sub[set_ind, 1:2], file = paste(root, "Data/locations_space_testing_set_", set, "_", cv_n_test, sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
	write.table(time_sub[set_ind], file = paste(root, "Data/locations_time_testing_set_", set, "_", cv_n_test, sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
	write.table(Z_sub[set_ind] - mean(Z_sub[-set_ind]), file = paste(root, "Data/data_st_testing_set_", set, "_", cv_n_test, sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)

	# SAVING TRAINING DATA

	cat("SAVING TRAINING DATA FOR SET ", set, "\n")

	write.table(locs_sub[-set_ind, 1:2], file = paste(root, "Data/locations_space_training_set_", set, "_", small_scale_n - cv_n_test, sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
	write.table(time_sub[-set_ind], file = paste(root, "Data/locations_time_training_set_", set, "_", small_scale_n - cv_n_test, sep = ''), sep = ",",row.names = FALSE, col.names = FALSE)
	write.table(Z_sub[-set_ind] - mean(Z_sub[-set_ind]), file = paste(root, "Data/data_st_training_set_", set, "_", small_scale_n - cv_n_test, sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)

}

#################################################### 		FIT THE MODEL IN EXAGEOSTAT 			####################################################
#################################################### 		RETRIEVE THE PREDICTIONS IN EXAGEOSTAT 		####################################################

# LOAD PREDICTIONS TXT FILE IMPORTED FROM EXAGEOSTAT

preds <- read.table(paste('/home/salvanmo/Desktop/ipdps/Predictions/predicted_values_September_08_2020_15:27:29', sep = ''), header = FALSE, sep = " ")

# LOAD THE CORRESPONDING TRAINING AND TESTING DATASETS

obs_TEST <- read.table(paste('/home/salvanmo/Desktop/ipdps/Data/data_st_testing_1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
locs_s_TEST <- read.table(paste('/home/salvanmo/Desktop/ipdps/Data/locations_space_testing_1', sep = ''), header = FALSE, sep = ",") %>% as.matrix()
locs_t_TEST <- read.table(paste('/home/salvanmo/Desktop/ipdps/Data/locations_time_testing_1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

obs_TRAIN <- read.table(paste('/home/salvanmo/Desktop/ipdps/Data/data_st_training_1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()
locs_s_TRAIN <- read.table(paste('/home/salvanmo/Desktop/ipdps/Data/locations_space_training_1', sep = ''), header = FALSE, sep = ",") %>% as.matrix()
locs_t_TRAIN <- read.table(paste('/home/salvanmo/Desktop/ipdps/Data/locations_time_training_1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

# CHOOSE 1 TEMPORAL LOCATION TO PLOT THE RESULTS. CHOOSE THE TEMPORAL LOCATION WITH THE MOST NUMBER OF SPATIAL MEASUREMENTS IN THE TRAINING SET

#time_index_for_plotting = which.max(count(locs_t_TRAIN)[, 2])

time_index_for_plotting = 413

ind_sub_test <- which(locs_t_TEST == time_index_for_plotting)		#indices of measurements in testing set at time point time_index_for_plotting
ind_sub_train <- which(locs_t_TRAIN == time_index_for_plotting)		#indices of measurements in training set at time point time_index_for_plotting

locs_s_sub_train <- locs_s_TRAIN[ind_sub_train, ]			#extract all spatial locations in training set at time point time_index_for_plotting
locs_t_sub_train <- locs_t_TRAIN[ind_sub_train, ]			#extract all temporal locations in training set at time point time_index_for_plotting
obs_sub_train <- obs_TRAIN[ind_sub_train]				#extract all spatial measurements in training set at time point time_index_for_plotting

locs_FULL <- rbind(locs_s_sub_train, locs_s_TEST[ind_sub_test, ])	#concatenate the spatial locations from the chosen training and testing sets 
obs_FULL <- c(obs_sub_train, obs_TEST[ind_sub_test])			#concatenate the spatial measurements from the chosen training and testing sets


locs_s_sub_test <- locs_s_TEST[ind_sub_test, ]			#extract all spatial locations in testing set at time point time_index_for_plotting
locs_t_sub_test <- locs_t_TEST[ind_sub_test, ]			#extract all temporal locations in testing set at time point time_index_for_plotting
obs_sub_test <- obs_TEST[ind_sub_test]				#extract all spatial measurements in testing set at time point time_index_for_plotting

obs_sub_pred <- preds[ind_sub_test, ]				#extract all spatial measurements in predictions file at time point time_index_for_plotting

zlim_range <- range(c(obs_sub_pred, obs_sub_pred, obs_sub_train))

jpeg(file = paste(root, 'Figures/6-application_data_predictions.jpg', sep = ''), width = 1800, height = 700)

split.screen( rbind(c(0.08,0.95,0.1,0.95), c(0.95,0.99,0.1,0.95)))
split.screen( figs = c( 1, 3 ), screen = 1 )

screen(3)
par(pty = 's')

quilt.plot(locs_s_sub_train[, 1], locs_s_sub_train[, 2], obs_sub_train, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, nx = 25, ny = 25, zlim = zlim_range)
map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)

mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)

mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)

mtext('Training', side = 3, line = 1, adj = 0.5, cex = 2.5, font = 2)

screen(4)
par(pty = 's')

quilt.plot(c(locs_s_sub_train[, 1], locs_s_sub_test[, 1]), c(locs_s_sub_train[, 2], locs_s_sub_test[, 2]), c(obs_sub_train, obs_sub_pred), ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, nx = 25, ny = 25, zlim = zlim_range)
map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)

mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)

mtext('Training & Predictions', side = 3, line = 1, adj = 0.5, cex = 2.5, font = 2)

screen(5)
par(pty = 's')

err <- (obs_sub_pred - obs_sub_test)^2

quilt.plot(locs_s_sub_test[, 1], locs_s_sub_test[, 2], err, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, nx = 25, ny = 25)
map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)

mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)

mtext('(Testing - Prediction)^2', side = 3, line = 1, adj = 0.5, cex = 2.5, font = 2)

screen(2)

x1 <- c(0.025,0.12,0.12,0.025) + 0.1
y1 <- c(0.3,0.3,0.7,0.7)
legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(min(err), max(err), length.out = 5), 1), CEX = 2)

close.screen( all=TRUE)
dev.off()


# FOR FULL MAPPING OF TRAINING AND PREDICTION, CHOOSE ALL THE REMAINING POINTS AT TIME time_index_for_plotting IN TESTING DATASET 

set <- 1

ind_sub_train <- which(locs_t_TRAIN == time_index_for_plotting) 	#find the index of spatial measurements in the training dataset at time time_index_for_plotting

locs_s_sub_train <- locs_s_TRAIN[ind_sub_train, ]			#extract the spatial locations at the training dataset at time time_index_for_plotting

# RETRIEVE THE INDEX OF THE SPATIAL LOCATIONS AT THE TRAINING DATASET AT TIME time_index_for_plotting IN THE ORIGINAL 550 x 2 LOCATIONS MATRIX

new_test_ind <- NULL

for(qq in 1:nrow(locs_s_sub_train)){
        new_test_ind <- c(new_test_ind, which(locs[, 1] == locs_s_sub_train[qq, 1] & locs[, 2] == locs_s_sub_train[qq, 2]))
}

# EXTRACT ONLY THE SPATIAL MEASUREMENTS AT TIME time_index_for_plotting AT LOCATIONS NOT INCLUDED IN THE SUBSET OF LOCATIONS ALREADY IN THE TRAINING SET

obs_sub_test <- res_mat1[time_index_for_plotting, -new_test_ind]
locs_s_sub_test <- locs[-new_test_ind, ]
locs_t_sub_test <- rep(time_index_for_plotting, nrow(locs_s_sub))

# SAVE THE NEW TESTING DATASET ALL FROM THE SAME TIME

write.table(locs_s_sub_test, file = paste(root, "Data/locations_space_testing_set_", set, "_", cv_n_test, "-FULL", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
write.table(locs_t_sub_test, file = paste(root, "Data/locations_time_testing_set_", set, "_", cv_n_test, "-FULL", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
write.table(obs_sub_test - mean(obs_sub_test), file = paste(root, "Data/data_st_testing_set_", set, "_", cv_n_test, "-FULL", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)


# LOAD PREDICTIONS TXT FILE IMPORTED FROM EXAGEOSTAT

obs_sub_pred <- read.table(paste('/home/salvanmo/Desktop/ipdps/Predictions/predicted_values_September_09_2020_16:01:44', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

zlim_range <- range(c(obs_sub_pred, obs_sub_pred, obs_sub_train))

jpeg(file = paste(root, 'Figures/6-application_data_predictions.jpg', sep = ''), width = 1800, height = 700)

split.screen( rbind(c(0.08,0.95,0.1,0.95), c(0.95,0.99,0.1,0.95)))
split.screen( figs = c( 1, 3 ), screen = 1 )

screen(3)
par(pty = 's')

quilt.plot(locs_s_sub_train[, 1], locs_s_sub_train[, 2], obs_sub_train, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, nx = 25, ny = 25, zlim = zlim_range)
map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)

mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)

mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)

mtext('Training', side = 3, line = 1, adj = 0.5, cex = 2.5, font = 2)

screen(4)
par(pty = 's')

quilt.plot(c(locs_s_sub_train[, 1], locs_s_sub_test[, 1]), c(locs_s_sub_train[, 2], locs_s_sub_test[, 2]), c(obs_sub_train, obs_sub_pred), ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, nx = 25, ny = 25, zlim = zlim_range)
map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)

mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
