generate_params <- function (beta0, beta1, beta2, beta3, beta4, w, locs, max_T = 1, generate_new_locs = F) {
  
  	theta <- beta0 + beta1 * (locs[, 1] - .5) + beta2 * (locs[, 2] - .5) + 
    		beta3 * (locs[, 1] - .5)^2 + beta4 * (locs[, 2] - .5)^2

	if(generate_new_locs)     locs_new <- locs
 
	if(max_T > 1){

  		for(tt in 1:(max_T - 1)){
    			theta.temp <- beta0 + beta1 * (locs[, 1] - tt * w[1] - .5) + beta2 * (locs[, 2]  - tt * w[2] - .5) + 
				beta3 * (locs[, 1] - tt * w[1] - .5)^2 + beta4 * (locs[, 2] - tt * w[2] - .5)^2
    

			if(generate_new_locs){
				locs_new <- rbind(locs_new, cbind(locs[, 1] - tt * w[1], locs[, 2] - tt * w[2]))
				theta <- c(theta, theta.temp - tt * w[3])
			}else{
				theta <- c(theta, theta.temp)
			}			
  		}
	}

	if(generate_new_locs)     return(cbind(locs_new, theta))      else return(theta)
}
