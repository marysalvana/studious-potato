directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'studious-potato/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))

# LOAD HOURLY DATA FROM 2016 TO 2019

start_yr <- 2016

DAT <- NULL

for(yr in start_yr:(start_yr + 3)){

	# LOAD DATA FOR EACH YEAR, EACH YEAR HAS A 550 x 8760 MATRIX, WHERE COLUMNS ARE MEASUREMENTS IN TIME AND ROWS ARE MEASUREMENTS IN SPACE
	cat("LOADING HOURLY DATA IN ", yr, ". . .", '\n')
	DAT_temp <- read.table(paste('/home/salvanmo/Desktop/ipdps/Data/pm_', yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	DAT <- cbind(DAT, DAT_temp)
}

# LOAD LOCATION COORDINATES

locs <- read.table(paste('/home/salvanmo/Desktop/ipdps/Data/locations_550', sep = ''), header = FALSE, sep = ",")

N <- nrow(locs)

# TAKE THE AVERAGE OF THE MEASUREMENTS FOR TWO DAYS BECAUSE THEY YIELD APPROXIMATELY SECOND-ORDER STATIONARY MEASUREMENTS

aggregate <- 48
TT <- floor(ncol(DAT) / aggregate)		#total number of temporal locations to be analyzed

DAT2 <- matrix(, ncol = TT, nrow = N)

for(aa in 1:TT){
	DAT2[, aa] <- apply(DAT[, (aa - 1) * aggregate + 1:aggregate], 1, mean)
}

# REMOVE THE SPATIO-TEMPORALLY VARYING TREND USING EMPIRICAL ORTHOGONAL FUNCTIONS

Yhat1 <- res_mat1 <- obs_mat_standardized1 <- t(DAT2)

obs_mat_SVD1 <- svd(obs_mat_standardized1)

variance.explained1 = prop.table(obs_mat_SVD1$d^2)

percent_sum_squared_variation1 <- cumsum(variance.explained1)[cumsum(variance.explained1) >= 0.9]
min_percent_sum_squared_variation1 <- min(percent_sum_squared_variation1)
num_singular_vec1 <- which(cumsum(variance.explained1) == min_percent_sum_squared_variation1) 

X1 <- cbind(rep(1, nrow(obs_mat_SVD1$u)), obs_mat_SVD1$u)

for(nn in 1:ncol(obs_mat_standardized1)){

	beta_hat <- solve(t(X1[, 1:num_singular_vec1]) %*% X1[, 1:num_singular_vec1]) %*% t(X1[, 1:num_singular_vec1]) %*% obs_mat_standardized1[, nn]

	Yhat1[, nn] <- Yhat_temp <- X1[, 1:num_singular_vec1] %*% beta_hat

	err <- Yhat_temp - obs_mat_standardized1[, nn]

	res_mat1[, nn] <- err/sd(err)

}

# CREATE A MATRIX OF MEASUREMENTS AND THEIR CORRESPONDING SPACE-TIME LOCATIONS: (RESIDUALS, LONGITUDE, LATITUDE, TIME) OF DIMENSION 401500 x 4 

Z <- NULL

for(tt in 1:TT){
	Z <- rbind(Z, cbind(res_mat1[tt, ], cbind(locs, rep(tt, nrow(locs)))))	
}

# SPLIT THE DATA INTO 90% TRAINING AND 10% TESTING

set.seed(1234)
subset_index <- sample(1:nrow(Z), 0.9 * N * TT)

# SAVE THE TRAINING AND TESTING DATASETS IN A TXT FILE

write.table(Z[-subset_index, 2:3], file = paste(root, "Data/locations_space_testing_FULL", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
write.table(Z[-subset_index, 4], file = paste(root, "Data/locations_time_testing_FULL", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
write.table(Z[-subset_index, 1] - mean(Z[-subset_index, 1]), file = paste(root, "Data/data_st_testing_FULL", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)

write.table(Z[subset_index, 2:3], file = paste(root, "Data/locations_space_training_FULL", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
write.table(Z[subset_index, 4], file = paste(root, "Data/locations_time_training_FULL", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
write.table(Z[subset_index, 1] - mean(Z[subset_index, 1]), file = paste(root, "Data/data_st_training_FULL", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)

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
