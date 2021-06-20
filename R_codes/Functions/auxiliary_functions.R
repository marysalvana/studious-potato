simulate_locations <- function(N, grid = T){

	grid_x <- seq(from = 0, to = 1, length.out = N)
	grid_y <- seq(from = 0, to = 1, length.out = N)
	sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

	return(sim_grid_locations)
}

write_to_txt <- function(data, file_name){

	#comma-separated if you save the bivariate realizations data.
	#space-separated if you save the cross-covariance matrix.	
	if(!is.matrix(data)){
		write.table(data, file = file_name, sep = ",", row.names = FALSE, col.names = FALSE)
	}else{	
		if(ncol(data) > 2 ){
			write.table(data, file = file_name, sep = " ", row.names = FALSE, col.names = FALSE)
		}else{	
			write.table(data, file = file_name, sep = ",", row.names = FALSE, col.names = FALSE)
		}
	}
}

colors=c("blue","yellow","red")
colsteps=100

legend.gradient2 = function(pnts,cols=tim.colors(64),limits=c(0,1), title='Legend', CEX = 1, ...){
  	pnts = try(as.matrix(pnts),silent=T)
  	if(!is.matrix(pnts)) stop("you must have a 4x2 matrix")
  	if(dim(pnts)[1]!=4 || dim (pnts)[2]!=2) stop ("Matrix must have dimensions of 4 rows and 2 columms")
  	if(length(cols)<2) stop("You must have 2 or more colors")
  	#break up the min and max into a number of values == length(cols)
  	yvals = seq(min(pnts[, 2]), max(pnts[, 2]), length=length(cols)+1)
  	#cycle through each of the yvals and create polygons
  	for (i in 1:length(cols)){  #create the polygon for that color
    		polygon(x=pnts[,1],y=c(yvals[i],yvals[i],yvals[i+1],yvals[i+1]),col=cols[i],border=F)
  	}
  	#add the text
	if(length(limits) == 5){
  		locationn <- seq(min(pnts[,2]),max(pnts[,2]),length.out = 5)
  		text(max(pnts[,1]),locationn[1],labels=round(limits[1],2),pos=4, cex = CEX)
  		text(max(pnts[,1]),locationn[2],labels=round(limits[2],2),pos=4, cex = CEX)
  		text(max(pnts[,1]),locationn[3],labels=round(limits[3],2),pos=4, cex = CEX)
  		text(max(pnts[,1]),locationn[4],labels=round(limits[4],2),pos=4, cex = CEX)
  		text(max(pnts[,1]),locationn[5],labels=round(limits[5],2),pos=4, cex = CEX)
	}else if(length(limits) == 3){
  		locationn <- seq(min(pnts[,2]),max(pnts[,2]),length.out = 3)
  		text(max(pnts[,1]),locationn[1],labels=round(limits[1],2),pos=4, cex = CEX)
  		text(max(pnts[,1]),locationn[2],labels=round(limits[2],2),pos=4, cex = CEX)
  		text(max(pnts[,1]),locationn[3],labels=round(limits[3],2),pos=4, cex = CEX)
	}
}

legend.gradient6 = function(pnts,cols=tim.colors(64),limits=c(0,1), title='Legend', cex = 1, ...){
  pnts = try(as.matrix(pnts),silent=T)
  if(!is.matrix(pnts)) stop("you must have a 4x2 matrix")
  if(dim(pnts)[1]!=4 || dim (pnts)[2]!=2) stop ("Matrix must have dimensions of 4 rows and 2 columms")
  if(length(cols)<2) stop("You must have 2 or more colors")
  #break up the min and max into a number of values == length(cols)
  yvals = seq(min(pnts[, 2]), max(pnts[, 2]), length=length(cols)+1)
  #cycle through each of the yvals and create polygons
  for (i in 1:length(cols)){  #create the polygon for that color
    polygon(x=pnts[,1],y=c(yvals[i],yvals[i],yvals[i+1],yvals[i+1]),col=cols[i],border=F)
  }
  #add the text
  locationn <- seq(min(pnts[,2]),max(pnts[,2]),length.out = 3)
  text(max(pnts[,1]) - 0.1,locationn[1],labels=round(limits[1],0),pos=4, cex = cex)
  text(max(pnts[,1]) - 0.1,locationn[2],labels=0.5 * (round(limits[1],0) + round(limits[5],0)),pos=4, cex = cex)
  text(max(pnts[,1]) - 0.1,locationn[3],labels=round(limits[5],0),pos=4, cex = cex)

}

toeplitz_mat <- function(S_list){
	k <- min(unlist(lapply(S_list, dim)))
	n <- length(S_list)
	#
	# Create the "strip".
	#
	strip <- array(NA, dim=c(k,k,2*n-1))
	for (i in 1:n) strip[,,i] <- S_list[[n+1-i]]
	if (n > 1) for (i in 2:n) strip[,,n+i-1] <- t(S_list[[i]])
	#
	# Assemble into "block-Toeplitz" form.
	#
	X <- array(NA, dim=c(k,k,n,n))
	# Blast the strip across X.
	#
	for (i in 1:n) X[,,,i] <- strip[,,(n+1-i):(2*n-i)]
	X <- matrix(aperm(X, c(1,3,2,4)), n*k)
}


generate_h_matrix <- function(LOCS){
	
	n <- nrow(LOCS)	

	H_MAT <- matrix(, ncol = ncol(LOCS), nrow = n^2)

        for(rr in 1:n){
                for(ss in 1:n){
			for(cc in 1:ncol(LOCS)){
                        		H_MAT[ (rr - 1) * n + ss, cc] <- LOCS[rr, cc] - LOCS[ss, cc]
			}
                }
        }
	
	return(H_MAT)

}

mahalanobis_dist <- function(x, SIGS) {
	SIGS_INV <- solve(SIGS)
	u <- apply(x, 1, function(y) y %*% SIGS_INV %*% y)
        d <- outer(u, u, `+`) - 2 * x %*% SIGS_INV %*% t(x)
        return(d)
}


data_format <- function(aggregate = 1, area){

	## LOAD HOURLY DATA FROM 2016 TO 2019

	DAT <- NULL

	for(YEAR in 2016:2019){

		cat("LOADING MEASUREMENTS FROM ", YEAR, ". . .", '\n')
		if(area == 'US'){
			DAT_temp <- read.table(paste(root, "Data/sc21/pm_US_", YEAR, sep = ""), header = FALSE, sep = ",") %>% as.matrix()
		}else if(area == 'SAUDI'){
			DAT_temp <- read.table(paste(root, "Data/sc21/pm_", YEAR, sep = ""), header = FALSE, sep = ",") %>% as.matrix()
		}
		DAT <- cbind(DAT, DAT_temp)

	}

	cat("LOADING LOCATIONS . . .", '\n')

	if(area == 'US'){
		LOCS <- read.table(paste(root, "Data/sc21/locations_US_", YEAR, sep = ""), header = FALSE, sep = ",")
	}else if(area == 'SAUDI'){
		LOCS <- read.table(paste(root, "Data/sc21/locations_", YEAR, sep = ""), header = FALSE, sep = ",")
	}

	N <- nrow(LOCS)

	## Indicate the degree of aggregation you want.

	#### FOR PM2.5 DATA OVER SAUDI, TAKE THE AVERAGE OF THE MEASUREMENTS FOR TWO DAYS (aggregate = 48) BECAUSE THEY YIELD APPROXIMATELY SECOND-ORDER STATIONARY MEASUREMENTS

	TT <- floor(ncol(DAT) / aggregate)		#total number of temporal locations to be analyzed

	if(aggregate > 1){
		DAT2 <- matrix(, ncol = TT, nrow = N)

		cat("AGGREGATING DATA BY TAKING THE MEAN OF ", aggregate, "CONSECUTIVE MEASUREMENTS", '\n')

		for(aa in 1:TT){
			DAT2[, aa] <- apply(DAT[, (aa - 1) * aggregate + 1:aggregate], 1, mean)
		}
	}else{
		DAT2 <- DAT
	}

	#Remove the nonstationarity in the mean
	#### REMOVE THE SPATIO-TEMPORALLY VARYING TREND USING EMPIRICAL ORTHOGONAL FUNCTIONS

	cat('Remove the nonstationarity in the mean', '\n')
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
	FINAL_DATA <- list("data_matrix" = res_mat1, "locations" = LOCS)
	
	return(FINAL_DATA)
}


# For saudi dataset only

data_format_saudi_testing_only <- function(aggregate = 48, forward_time_predict = 1){

	## LOAD HOURLY DATA FROM 2016 TO 2019

	DAT <- NULL

	for(YEAR in 2020:2020){

		cat("LOADING MEASUREMENTS FROM ", YEAR, ". . .", '\n')
		DAT_temp <- read.table(paste(root, "Data/sc21/pm_", YEAR, sep = ""), header = FALSE, sep = ",") %>% as.matrix()
		DAT <- cbind(DAT, DAT_temp)

	}

	cat("LOADING LOCATIONS . . .", '\n')

	LOCS <- read.table(paste(root, "Data/sc21/locations_", YEAR, sep = ""), header = FALSE, sep = ",")

	N <- nrow(LOCS)

	## Indicate the degree of aggregation you want.

	#### FOR PM2.5 DATA OVER SAUDI, TAKE THE AVERAGE OF THE MEASUREMENTS FOR TWO DAYS (aggregate = 48) BECAUSE THEY YIELD APPROXIMATELY SECOND-ORDER STATIONARY MEASUREMENTS

	TT <- floor(ncol(DAT) / aggregate)		#total number of temporal locations to be analyzed

	if(aggregate > 1){
		DAT2 <- matrix(, ncol = TT, nrow = N)

		cat("AGGREGATING DATA BY TAKING THE MEAN OF ", aggregate, "CONSECUTIVE MEASUREMENTS", '\n')

		for(aa in 1:TT){
			DAT2[, aa] <- apply(DAT[, (aa - 1) * aggregate + 1:aggregate], 1, mean)
		}
	}else{
		DAT2 <- DAT
	}

	#Remove the nonstationarity in the mean
	#### REMOVE THE SPATIO-TEMPORALLY VARYING TREND USING EMPIRICAL ORTHOGONAL FUNCTIONS

	cat('Remove the nonstationarity in the mean', '\n')

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
	FINAL_DATA <- list("data_matrix" = res_mat1, "locations" = LOCS)

	N <- dim(FINAL_DATA[["data_matrix"]])[2]
	TT <- forward_time_predict

	Z <- NULL

	for(tt in 1:TT){
		Z <- rbind(Z, cbind(FINAL_DATA[["data_matrix"]][tt, ], FINAL_DATA[["locations"]], rep(tt + 730, N)))	
	}

	#### SPLIT THE DATA INTO 90% TRAINING AND 10% TESTING

	locs_s_test <- Z[, 2:3]                
	locs_t_test <- Z[, 4]              
	obs_test <- Z[, 1]

	#### SAVE THE TESTING DATASETS IN A TXT FILE

	write.table(locs_s_test, file = paste(root, "Data/sc21/locations_space_testing_NEW", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
	write.table(locs_t_test, file = paste(root, "Data/sc21/locations_time_testing_NEW", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
	write.table(obs_test - mean(obs_test), file = paste(root, "Data/sc21/data_st_testing_NEW", sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)

	cat("Textfiles are saved in ", paste(root, 'Data/sc21/', sep = ''), '\n')
	
}

#Function to check visually if the plot is second-order stationary in its covariance structure

plot_realdata_for_checking_stationarity <- function(data_list, file_name, start_hr = 1, saudi = F, aggregate = 'hourly'){

	cat('PLOTTING THE FIRST FIVE SPACETIME IMAGES TO CHECK STATIONARITY', '\n')

	zlim_range1 <- range(data_list[["data_matrix"]][start_hr:(start_hr + 4),])
	zlim_range1 <- c(sign(min(zlim_range1)) * round(abs(min(zlim_range1)), 1) - 0.1, sign(max(zlim_range1)) * round(abs(max(zlim_range1)), 1) + 0.1)

	jpeg(file = paste(root, 'Figures/', file_name, sep = ''), width = 1800, height = 600)

	split.screen( rbind(c(0.05,0.95,0.1,0.85), c(0.95,0.99,0.1,0.95)))
	split.screen( figs = c( 1, 5 ), screen = 1 )

	hr_count <- 0
	for(hr in start_hr:(start_hr + 4)){
		
		hr_count <- hr_count + 1
		
		screen(2 + hr_count)

		par(pty = 's')
		par(mai=c(0.2,0.2,0.2,0.2))
		
		if(hr_count == 1){
		quilt.plot(data_list[["locations"]][, 1], data_list[["locations"]][, 2], data_list[["data_matrix"]][hr, ], zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2)
		#mtext('log PM 2.5', side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 'blue')
		}else{
		quilt.plot(data_list[["locations"]][, 1], data_list[["locations"]][, 2], data_list[["data_matrix"]][hr, ], zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
		}

		if(saudi){
			map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
		}else{
			map("state", xlim =  c(-120, -70), ylim = c(30, 50), lwd = 0.75, add = T)
		}
		
		if(hr_count == 1){
			mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
		}

		if(aggregate == 'hourly'){
			mtext(paste(hr - 1 - floor(hr/24) * 24, ':00', sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
		}else if(aggregate == 'daily'){
			mtext(paste('January ', hr, ', ', yr, sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
		}else{
			mtext(paste('January ', 2 * hr_count - 1, '-', 2 * hr_count, sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
			#mtext(paste((hr - 1) * aggregate, ':00', sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
		}
		mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
		if(hr == 3) mtext('Mean log PM2.5 Concentration for the Period', side = 3, line = 7, cex = 3, font = 2, col = 4)
	}

	screen(2)

	x1 <- c(0.025,0.12,0.12,0.025) + 0.1
	y1 <- c(0.3,0.3,0.7,0.7)
	legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(-5, 5, length.out = 5), 1), CEX = 2)

	close.screen( all=TRUE)
	dev.off()

	cat("Check image in ", paste(root, 'Figures/', file_name, sep = ''), '\n')
}

#Function for manuscript ready real data plots

plot_realdata_for_manuscript <- function(data_list, file_name, start_hr = 1, saudi = F, aggregate = 'hourly'){

	cat('PLOTTING THE FIRST SIX SPACETIME IMAGES FOR MANUSCRIPT', '\n')

	step_size_for_display <- 4
	subset_for_display <- seq(start_hr, start_hr + 5 * step_size_for_display, by = step_size_for_display)

	zlim_range1 <- range(data_list[["data_matrix"]][subset_for_display, ])
	zlim_range1 <- c(sign(min(zlim_range1)) * round(abs(min(zlim_range1)), 1) - 0.1, sign(max(zlim_range1)) * round(abs(max(zlim_range1)), 1) + 0.1)

	jpeg(file = paste(root, 'Figures/', file_name, sep = ''), width = 1200, height = 1000)

	split.screen( rbind(c(0.06,0.94,0.08,0.93), c(0.94,0.98,0.08,0.93)))
	split.screen( figs = c( 2, 3 ), screen = 1 )


	hr_count <- 0
	for(hr in subset_for_display){
		
		hr_count <- hr_count + 1
		
		screen(2 + hr_count)

		par(pty = 's')
		par(mai=c(0.2,0.2,1,0.2))
		
		if(hr_count %in% c(1, 4)){
		quilt.plot(data_list[["locations"]][, 1], data_list[["locations"]][, 2], data_list[["data_matrix"]][hr, ], zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, xaxt = 'n')
		}else{
		quilt.plot(data_list[["locations"]][, 1], data_list[["locations"]][, 2], data_list[["data_matrix"]][hr, ], zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2, xaxt = 'n')
		}
		if(saudi){
			map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
		}else{
			map("state", xlim =  c(-120, -70), ylim = c(30, 50), lwd = 0.75, add = T)
		}
		
		if(hr_count %in% c(1, 4)){
			mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
		}
		if(hr_count >= 4){	
			axis(1, cex.axis = 2)
			mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
		}

		if(aggregate == 'hourly'){
			mtext(paste(hr - 1 - floor(hr/24) * 24, ':00', sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
		}else if(aggregate == 'daily'){
			mtext(paste('January ', hr, ', ', yr, sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
		}else{
			mtext(paste('January ', 2 * hr_count - 1, '-', 2 * hr_count, sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
		}
		if(hr_count == 2) mtext('Mean log PM 2.5 Concentration for the Period', side = 3, line = 6, cex = 3, font = 2, col = 4)
	}

	screen(2)

	x1 <- c(0.025,0.1,0.1,0.025) + 0.1
	y1 <- c(0.27,0.27,0.65,0.65)
	legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(-5, 5, length.out = 3), 1), CEX = 3)

	close.screen( all=TRUE)
	dev.off()

	cat("Check image in ", paste(root, 'Figures/', file_name, sep = ''), '\n')
}


#Function to create training and testing datasets

data_split <- function(data_list, training_percent = 0.9, temporal_length = NULL, file_name){

	cat('SPLITTING THE DATA INTO TRAINING AND TESTING DATASETS', '\n')

	#### CREATE A MATRIX OF MEASUREMENTS AND THEIR CORRESPONDING SPACE-TIME LOCATIONS: (RESIDUALS, LONGITUDE, LATITUDE, TIME) OF DIMENSION (NT) x 4 
	#### IF TT is very big, you can just get a subset 1:temporal_length

	N <- dim(data_list[["data_matrix"]])[2]
	if(is.null(temporal_length)){
		TT <- dim(data_list[["data_matrix"]])[1]
	}else{
		TT <- temporal_length
	}

	Z <- NULL

	for(tt in 1:TT){
		Z <- rbind(Z, cbind(data_list[["data_matrix"]][tt, ], data_list[["locations"]], rep(tt, N)))	
	}

	#### SPLIT THE DATA INTO 90% TRAINING AND 10% TESTING

	set.seed(1234)
	subset_index <- sample(1:nrow(Z), training_percent * N * TT)

	locs_s_train <- Z[subset_index, 2:3]                
	locs_t_train <- Z[subset_index, 4]              
	obs_train <- Z[subset_index, 1]

	locs_s_test <- Z[-subset_index, 2:3]                
	locs_t_test <- Z[-subset_index, 4]              
	obs_test <- Z[-subset_index, 1]

	#### SAVE THE TRAINING AND TESTING DATASETS IN A TXT FILE

	write.table(locs_s_test, file = paste(root, "Data/sc21/locations_space_testing_FULL", file_name, sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
	write.table(locs_t_test, file = paste(root, "Data/sc21/locations_time_testing_FULL", file_name, sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
	write.table(obs_test - mean(obs_test), file = paste(root, "Data/sc21/data_st_testing_FULL", file_name, sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)

	write.table(locs_s_train, file = paste(root, "Data/sc21/locations_space_training_FULL", file_name, sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
	write.table(locs_t_train, file = paste(root, "Data/sc21/locations_time_training_FULL", file_name, sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
	write.table(obs_train - mean(obs_train), file = paste(root, "Data/sc21/data_st_training_FULL", file_name, sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)

	cat("Textfiles are saved in ", paste(root, 'Data/sc21/', sep = ''), '\n')
}

# For US dataset only

data_format_US_testing_only <- function(data_list, temporal_length = NULL, forward_time_predict = 1, file_name){

	cat('SPLITTING THE DATA INTO TRAINING AND TESTING DATASETS', '\n')

	#### CREATE A MATRIX OF MEASUREMENTS AND THEIR CORRESPONDING SPACE-TIME LOCATIONS: (RESIDUALS, LONGITUDE, LATITUDE, TIME) OF DIMENSION (NT) x 4 
	#### IF TT is very big, you can just get a subset 1:temporal_length

	N <- dim(data_list[["data_matrix"]])[2]
	if(is.null(temporal_length)){
		TT <- dim(data_list[["data_matrix"]])[1]
	}else{
		TT <- temporal_length
	}

	Z <- NULL

	for(tt in 1:(TT + forward_time_predict)){
		Z <- rbind(Z, cbind(data_list[["data_matrix"]][tt, ], data_list[["locations"]], rep(tt, N)))	
	}

	locs_s_forward <- Z[-(1:(TT * N)), 2:3]                
	locs_t_forward <- Z[-(1:(TT * N)), 4]              
	obs_test_forward <- Z[-(1:(TT * N)), 1]

	#### SAVE THE TRAINING AND TESTING DATASETS IN A TXT FILE

	write.table(locs_s_forward, file = paste(root, "Data/sc21/locations_space_testing_NEW", file_name, sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
	write.table(locs_t_forward, file = paste(root, "Data/sc21/locations_time_testing_NEW", file_name, sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)
	write.table(obs_test_forward - mean(obs_test_forward), file = paste(root, "Data/sc21/data_st_testing_NEW", file_name, sep = ''), sep = ",", row.names = FALSE, col.names = FALSE)

	cat("Textfiles are saved in ", paste(root, 'Data/sc21/', sep = ''), '\n')
}

plot_simulated_data_for_beamer <- function(covariance, realizations, locations, file_name, reference_locations){

	n <- nrow(locations)
	N <- sqrt(n)

	zlim_range1 <- c(0, 1)
	zlim_range2 <- range(realizations[1, 1:(n * 5)])

	jpeg(file = paste(root, 'Figures/', file_name, sep = ''), width = 1200, height = 600)

	split.screen( rbind(c(0.06,0.94,0.08,0.93), c(0.94,0.98,0.08,0.93)))
	split.screen( figs = c( 3, 6 ), screen = 1 )


	hr_count <- 0
	for(tt in 1:5){
		
		hr_count <- hr_count + 1
		
		screen(3 + hr_count)

		par(pty = 's')
		par(mai=c(0.2,0.2,0.2,0.2))
		
		quilt.plot(locations[, 1], locations[, 2], realizations[1, (tt - 1) * n + 1:n], zlim = zlim_range2, nx = N, ny = N, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 1, xaxt = 'n', yaxt = 'n')
		if(tt == 1){
			axis(2, at = seq(min(locations[, 2]), max(locations[, 2]), length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 1)
		}
		mtext(paste('t = ', tt, sep = ''), side = 3, line = 1, adj = 0.5, cex = 2, font = 2)
	}	

	hr_count <- hr_count + 1
	screen(3 + hr_count)

	par(pty = 's')
	par(mai=c(0.2,0.2,0.2,0.2))

	quilt.plot(locations[, 1], locations[, 2], realizations[1, 1:n], zlim = zlim_range2, nx = N, ny = N, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 1, xaxt = 'n', yaxt = 'n')
	axis(2, at = seq(min(locations[, 2]), max(locations[, 2]), length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 1)
	points(matrix(locations[reference_locations[1], ], ncol = 2), col = 'black', pch = 4, cex = 3, lwd = 4)
	mtext(paste('Ref Loc 1', sep = ''), side = 2, line = 4, adj = 0.5, cex = 2, font = 2, col = 'blue')

	for(tt in 1:5){
		
		hr_count <- hr_count + 1
		
		screen(3 + hr_count)

		par(pty = 's')
		par(mai=c(0.2,0.2,0.2,0.2))
		
		quilt.plot(locations[, 1], locations[, 2], covariance[1, (tt - 1) * n + 1:n], zlim = zlim_range1, nx = N, ny = N, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 1, yaxt = 'n', xaxt = 'n')

	}	

	hr_count <- hr_count + 1
	screen(3 + hr_count)

	par(pty = 's')
	par(mai=c(0.2,0.2,0.2,0.2))

	quilt.plot(locations[, 1], locations[, 2], realizations[1, 1:n], zlim = zlim_range2, nx = N, ny = N, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 1, xaxt = 'n', yaxt = 'n')
	axis(1, at = seq(min(locations[, 1]), max(locations[, 1]), length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 1)
	axis(2, at = seq(min(locations[, 2]), max(locations[, 2]), length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 1)
	points(matrix(locations[reference_locations[2], ], ncol = 2), col = 'black', pch = 4, cex = 3, lwd = 4)
	mtext(paste('Ref Loc 2', sep = ''), side = 2, line = 4, adj = 0.5, cex = 2, font = 2, col = 'blue')

	mtext(paste('t = ', 1, sep = ''), side = 2, line = 2, adj = 1.5, cex = 2, font = 2)

	for(tt in 1:5){
		
		hr_count <- hr_count + 1
		
		screen(3 + hr_count)

		par(pty = 's')
		par(mai=c(0.2,0.2,0.2,0.2))
		
		quilt.plot(locations[, 1], locations[, 2], covariance[2, (tt - 1) * n + 1:n], zlim = zlim_range1, nx = N, ny = N, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 1, yaxt = 'n', xaxt = 'n')
		axis(1, at = seq(min(locations[, 1]), max(locations[, 1]), length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 1)

	}	

	screen(2)

	x1 <- c(0.025,0.1,0.1,0.025) + 0.1
	y1 <- c(0.17,0.17,0.45,0.45)
	legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(0, 1, length.out = 3), 1), CEX = 1.5)

	close.screen( all=TRUE)
	dev.off()

	cat("Check image in ", paste(root, 'Figures/', file_name, sep = ''), '\n')

}


movie_simulated_data_for_beamer <- function(realizations, locations, file_name, col_labels = c("I", "II", "III"), row_labels = c("A", "B")){

	n <- nrow(locations)
	N <- sqrt(n)

	zlim_range2 <- range(realizations[, 1:(n * 5)])

	mod_labels <- c(col_labels, row_labels)

	for(tt in 1:5){

		jpeg(file = paste(root, 'Figures/', file_name, '_t', tt, '.jpg', sep = ''), width = 1300, height = 900)

		split.screen( rbind(c(0.1,0.92,0.06,0.93), c(0.92,0.98,0.06,0.93)))
		split.screen( figs = c( 2, 3 ), screen = 1 )

		hr_count <- 0

		for(mod in 1:6){
			
			hr_count <- hr_count + 1
			
			screen(2 + hr_count)

			par(pty = 's')
			par(mai=c(0.3,0.3,0.3,0.3))
			
			quilt.plot(locations[, 1], locations[, 2], realizations[mod, (tt - 1) * n + 1:n], zlim = zlim_range2, nx = N, ny = N, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 1, xaxt = 'n', yaxt = 'n')
			if(mod <= 3){
				mtext(paste('t = ', tt, sep = ''), side = 3, line = 1, adj = 0.5, cex = 2, font = 2)
				mtext(mod_labels[mod], side = 3, line = 4, adj = 0.5, cex = 3, font = 2, col = 4)
			}
			if(mod %in% c(1, 4)){
				mtext(expression(s[y]), side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
				axis(2, at = seq(min(locations[, 2]), max(locations[, 2]), length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 2)
			}
			if(mod >= 4){
				mtext(expression(s[x]), side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
				axis(1, at = seq(min(locations[, 1]), max(locations[, 1]), length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 2)
			}
			if(mod == 1){
				text(min(locations[, 1]) - 2.3, 0, mod_labels[4], pos = 4, cex = 3, font = 2, col = 4, xpd = NA)
				#mtext(mod_labels[4], side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 4)
			}
			if(mod == 4){
				text(min(locations[, 1]) - 2.3, 0, mod_labels[5], pos = 4, cex = 3, font = 2, col = 4, xpd = NA)
			}
		}	

		screen(2)

		x1 <- c(0.01,0.1,0.1,0.01)
		y1 <- c(0.3,0.3,0.7,0.7)
		legend.gradient2(cbind(x1,y1), title = "", limits = seq(-5, 5, length.out = 5), CEX = 2)

		close.screen( all=TRUE)
		dev.off()

	}
}


plot_univariate_nonstationary_covariance_heatmap <- function(covariance, locations, file_name, reference_locations){

	n <- nrow(locations)
	N <- sqrt(n)

	zlim_range1 <- c(0, 1)

	pdf(file = paste(root, 'Figures/', file_name, sep = ''), width = 30, height = 11)

	split.screen( rbind(c(0.05,0.48,0.06,0.88), c(0.52,0.95,0.06,0.88), c(0.96,0.99,0.06,0.88)))
	split.screen( figs = c( 2, 3 ), screen = 1 )
	split.screen( figs = c( 2, 3 ), screen = 2 ) 


	hr_count <- 0

	for(tt in 1:3){
		
		hr_count <- hr_count + 1
		
		screen(3 + hr_count)

		par(pty = 's')
		par(mai=c(0.2,0.2,0.2,0.2))
		
		quilt.plot(locations[, 1], locations[, 2], covariance[2, (tt - 1) * n + 1:n], zlim = zlim_range1, nx = N, ny = N, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 1, yaxt = 'n', xaxt = 'n')
		points(matrix(locations[reference_locations[2], ], ncol = 2), col = 'black', pch = 4, cex = 3, lwd = 4)
		mtext(paste('t = ', tt, sep = ''), side = 3, line = 1, adj = 0.5, cex = 2, font = 2)
		
		if(tt == 1){
			axis(2, at = seq(min(locations[, 2]), max(locations[, 2]), length.out = 5), labels = seq(0, 1, length.out = 5), cex.axis = 1.5)
		}

	}	

	screen(3)

	x1 <- c(0.025,0.1,0.1,0.025) + 0.1
	y1 <- c(0.17,0.17,0.45,0.45)
	legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(0, 1, length.out = 3), 1), CEX = 1.5)

	close.screen( all=TRUE)
	dev.off()

	cat("Check image in ", paste(root, 'Figures/', file_name, sep = ''), '\n')

}


