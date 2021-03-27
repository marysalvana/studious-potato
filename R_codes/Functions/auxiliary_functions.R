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
