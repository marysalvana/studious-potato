// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <Rcpp.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <RcppNumerical.h>

using namespace Rcpp;
using namespace boost::math;

NumericMatrix rbind_cpp(NumericMatrix a,NumericMatrix b) {
//the colnumber of first matrix
int acoln=a.nrow();
//the colnumber of second matrix
int bcoln=b.nrow();
//build a new matrix, the dim is a.nrow() and acoln+bcoln
NumericMatrix out(acoln+bcoln, a.ncol()) ;
for (int j = 0; j < acoln + bcoln; j++) {
if (j < acoln) {
out(j,_) = a(j,_);
} else {
//put the context in the second matrix to the new matrix
out(j,_) = b(j-acoln,_);
}
}
return out ;
}

NumericMatrix cbind_cpp(NumericMatrix a, NumericMatrix b) {
//the colnumber of first matrix
int acoln=a.ncol();
//the colnumber of second matrix
int bcoln=b.ncol();
//build a new matrix, the dim is a.nrow() and acoln+bcoln
NumericMatrix out(a.nrow(),acoln+bcoln) ;
for (int j = 0; j < acoln + bcoln; j++) {
if (j < acoln) {
out(_,j) = a(_,j);
} else {
//put the context in the second matrix to the new matrix
out(_,j) = b(_,j-acoln);
}
}
return out ;
}

// [[Rcpp::export]]

List SPATIALLY_VARYING_PARAMETERS(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind, NumericVector & param_nonstat, int & time) {

  	const int m = Loc.nrow(), n_wind = wind.nrow();
  	float sigma2 = param(0), beta = param(1), nu = param(2);
	NumericMatrix cor(m, m), param_matrix(m, 6);
	List conso(2);

	for (int i = 0; i < m; ++i) {
    		for (int j = 0; j <= i; ++j) {

        		double temp_val = 0.0;

        		for (int k = 0; k < n_wind; ++k) {

			      	double new_loc1_x = Loc(i, 0) - wind(k, 0) * Loc(i, 2);
			      	double new_loc1_y = Loc(i, 1) - wind(k, 1) * Loc(i, 2);
			      	double new_loc2_x = Loc(j, 0) - wind(k, 0) * Loc(j, 2);
			      	double new_loc2_y = Loc(j, 1) - wind(k, 1) * Loc(j, 2);

			      	double x1 = new_loc1_x - new_loc2_x;
			      	double x2 = new_loc1_y - new_loc2_y;
			      	NumericVector z2 = {x1, x2};

			      	double omega1 = param_nonstat(0) + param_nonstat(1) * (new_loc1_x - .5) + param_nonstat(2) * (new_loc1_y - .5) + param_nonstat(3) * pow(new_loc1_x - .5, 2) + param_nonstat(4) * pow(new_loc1_y - .5, 2);
			      	double log_lam1_1 = param_nonstat(5) + param_nonstat(6) * (new_loc1_x - .5) + param_nonstat(7) * (new_loc1_y - .5) + param_nonstat(8) * pow(new_loc1_x - .5, 2) + param_nonstat(9) * pow(new_loc1_y - .5, 2);
			      	double log_lam1_2 = param_nonstat(10) + param_nonstat(11) * (new_loc1_x - .5) + param_nonstat(12) * (new_loc1_y - .5) + param_nonstat(13) * pow(new_loc1_x - .5, 2) + param_nonstat(14) * pow(new_loc1_y - .5, 2);

			      	double omega2 = param_nonstat(0) + param_nonstat(1) * (new_loc2_x - .5) + param_nonstat(2) * (new_loc2_y - .5) + param_nonstat(3) * pow(new_loc2_x - .5, 2) + param_nonstat(4) * pow(new_loc2_y - .5, 2);
			      	double log_lam2_1 = param_nonstat(5) + param_nonstat(6) * (new_loc2_x - .5) + param_nonstat(7) * (new_loc2_y - .5) + param_nonstat(8) * pow(new_loc2_x - .5, 2) + param_nonstat(9) * pow(new_loc2_y - .5, 2);
			      	double log_lam2_2 = param_nonstat(10) + param_nonstat(11) * (new_loc2_x - .5) + param_nonstat(12) * (new_loc2_y - .5) + param_nonstat(13) * pow(new_loc2_x - .5, 2) + param_nonstat(14) * pow(new_loc2_y - .5, 2);

				double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
				double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
				double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

				double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
				double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
				double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

				if(i == j){
					param_matrix(i, 0) = Loc(i, 0);
					param_matrix(i, 1) = Loc(i, 1);
					param_matrix(i, 2) = Loc(i, 2);
					param_matrix(i, 3) = Sigma1_11;
					param_matrix(i, 4) = Sigma1_12;
					param_matrix(i, 5) = Sigma1_22;
				}

				double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
				double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

				double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
				double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
				double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

				double Inv_ij_11 = Kernel_ij_22;
				double Inv_ij_22 = Kernel_ij_11;
				double Inv_ij_12 = - Kernel_ij_12;
				double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

				double sigma = sqrt(sqrt(det_i * det_j)/det_ij);

				double dist = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22);

				if (dist == 0) {
					temp_val = temp_val + sigma2 * pow(sigma, 2);
				}else {
					temp_val = temp_val + sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
				}
        		}
        		cor(i, j) = temp_val / n_wind;
        		cor(j, i) = cor(i, j);
    		}
  	}
	conso(0) = cor;
	conso(1) = param_matrix;
  	return conso;
}

// [[Rcpp::export]]

List SPATIALLY_VARYING_PARAMETERS_FOR_FITTING(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind, NumericMatrix & param_nonstat, int & time) {

  	const int m = Loc.nrow(), n_wind = wind.nrow();
  	float sigma2 = param(0), beta = param(1), nu = param(2);
	NumericMatrix cor(m, m);
	List conso(1);

	for (int i = 0; i < m; ++i) {
    		for (int j = 0; j <= i; ++j) {

        		double temp_val = 0.0;

        		for (int k = 0; k < n_wind; ++k) {

			      	double new_loc1_x = Loc(i, 0) - wind(k, 0) * Loc(i, 2);
			      	double new_loc1_y = Loc(i, 1) - wind(k, 1) * Loc(i, 2);
			      	double new_loc2_x = Loc(j, 0) - wind(k, 0) * Loc(j, 2);
			      	double new_loc2_y = Loc(j, 1) - wind(k, 1) * Loc(j, 2);

			      	double x1 = new_loc1_x - new_loc2_x;
			      	double x2 = new_loc1_y - new_loc2_y;
			      	NumericVector z2 = {x1, x2};

			      	double omega1 = param_nonstat(i, 0);
			      	double log_lam1_1 = param_nonstat(i, 1);
			      	double log_lam1_2 = param_nonstat(i, 2);

			      	double omega2 = param_nonstat(j, 0);
			      	double log_lam2_1 = param_nonstat(j, 1);
			      	double log_lam2_2 = param_nonstat(j, 2);

				double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
				double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
				double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

				double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
				double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
				double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

				double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
				double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

				double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
				double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
				double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

				double Inv_ij_11 = Kernel_ij_22;
				double Inv_ij_22 = Kernel_ij_11;
				double Inv_ij_12 = - Kernel_ij_12;
				double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

				double sigma = sqrt(sqrt(det_i * det_j)/det_ij);

				double dist = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22);

				if (dist == 0) {
					temp_val = temp_val + sigma2 * pow(sigma, 2);
				}else {
					temp_val = temp_val + sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
				}
        		}
        		cor(i, j) = temp_val / n_wind;
        		cor(j, i) = cor(i, j);
    		}
  	}
	conso(0) = cor;
  	return conso;
}

// [[Rcpp::export]]

NumericMatrix SPATIALLY_VARYING_PARAMETERS_PARALLEL(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind, NumericVector & param_nonstat, int & time) {

  	const int m = Loc.nrow(), n_wind = wind.nrow();
  	float sigma2 = param(0), beta = param(1), nu = param(2);
	NumericMatrix cor(m, m);

	for (int i = 0; i < m; ++i) {
    		for (int j = 0; j <= i; ++j) {

        		double temp_val = 0.0;

        		for (int k = 0; k < n_wind; ++k) {

			      	double new_loc1_x = Loc(i, 0) - wind(k, 0) * Loc(i, 2);
			      	double new_loc1_y = Loc(i, 1) - wind(k, 1) * Loc(i, 2);
			      	double new_loc2_x = Loc(j, 0) - wind(k, 0) * Loc(j, 2);
			      	double new_loc2_y = Loc(j, 1) - wind(k, 1) * Loc(j, 2);

			      	double x1 = new_loc1_x - new_loc2_x;
			      	double x2 = new_loc1_y - new_loc2_y;
			      	NumericVector z2 = {x1, x2};

			      	double omega1 = param_nonstat(0) + param_nonstat(1) * (new_loc1_x - .5) + param_nonstat(2) * (new_loc1_y - .5) + param_nonstat(3) * pow(new_loc1_x - .5, 2) + param_nonstat(4) * pow(new_loc1_y - .5, 2);
			      	double log_lam1_1 = param_nonstat(5) + param_nonstat(6) * (new_loc1_x - .5) + param_nonstat(7) * (new_loc1_y - .5) + param_nonstat(8) * pow(new_loc1_x - .5, 2) + param_nonstat(9) * pow(new_loc1_y - .5, 2);
			      	double log_lam1_2 = param_nonstat(10) + param_nonstat(11) * (new_loc1_x - .5) + param_nonstat(12) * (new_loc1_y - .5) + param_nonstat(13) * pow(new_loc1_x - .5, 2) + param_nonstat(14) * pow(new_loc1_y - .5, 2);

			      	double omega2 = param_nonstat(0) + param_nonstat(1) * (new_loc2_x - .5) + param_nonstat(2) * (new_loc2_y - .5) + param_nonstat(3) * pow(new_loc2_x - .5, 2) + param_nonstat(4) * pow(new_loc2_y - .5, 2);
			      	double log_lam2_1 = param_nonstat(5) + param_nonstat(6) * (new_loc2_x - .5) + param_nonstat(7) * (new_loc2_y - .5) + param_nonstat(8) * pow(new_loc2_x - .5, 2) + param_nonstat(9) * pow(new_loc2_y - .5, 2);
			      	double log_lam2_2 = param_nonstat(10) + param_nonstat(11) * (new_loc2_x - .5) + param_nonstat(12) * (new_loc2_y - .5) + param_nonstat(13) * pow(new_loc2_x - .5, 2) + param_nonstat(14) * pow(new_loc2_y - .5, 2);

				double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
				double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
				double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

				double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
				double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
				double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

				double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
				double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

				double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
				double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
				double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

				double Inv_ij_11 = Kernel_ij_22;
				double Inv_ij_22 = Kernel_ij_11;
				double Inv_ij_12 = - Kernel_ij_12;
				double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

				double sigma = sqrt(sqrt(det_i * det_j)/det_ij);

				double dist = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22);

				if (dist == 0) {
					temp_val = temp_val + sigma2 * pow(sigma, 2);
				}else {
					temp_val = temp_val + sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
				}
        		}
        		cor(i, j) = temp_val / n_wind;
        		cor(j, i) = cor(i, j);
    		}
  	}
  	return cor;
}

// [[Rcpp::export]]

NumericMatrix SPATIALLY_VARYING_PARAMETERS_FOR_FITTING_PARALLEL(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind, NumericMatrix & param_nonstat, int & time) {

  	const int m = Loc.nrow(), n_wind = wind.nrow();
  	float sigma2 = param(0), beta = param(1), nu = param(2);
	NumericMatrix cor(m, m);

	for (int i = 0; i < m; ++i) {
    		for (int j = 0; j <= i; ++j) {

        		double temp_val = 0.0;

        		for (int k = 0; k < n_wind; ++k) {

			      	double new_loc1_x = Loc(i, 0) - wind(k, 0) * Loc(i, 2);
			      	double new_loc1_y = Loc(i, 1) - wind(k, 1) * Loc(i, 2);
			      	double new_loc2_x = Loc(j, 0) - wind(k, 0) * Loc(j, 2);
			      	double new_loc2_y = Loc(j, 1) - wind(k, 1) * Loc(j, 2);

			      	double x1 = new_loc1_x - new_loc2_x;
			      	double x2 = new_loc1_y - new_loc2_y;
			      	NumericVector z2 = {x1, x2};

			      	double omega1 = param_nonstat(i, 0);
			      	double log_lam1_1 = param_nonstat(i, 1);
			      	double log_lam1_2 = param_nonstat(i, 2);

			      	double omega2 = param_nonstat(j, 0);
			      	double log_lam2_1 = param_nonstat(j, 1);
			      	double log_lam2_2 = param_nonstat(j, 2);

				double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
				double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
				double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

				double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
				double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
				double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

				double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
				double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

				double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
				double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
				double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

				double Inv_ij_11 = Kernel_ij_22;
				double Inv_ij_22 = Kernel_ij_11;
				double Inv_ij_12 = - Kernel_ij_12;
				double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

				double sigma = sqrt(sqrt(det_i * det_j)/det_ij);

				double dist = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22);

				if (dist == 0) {
					temp_val = temp_val + sigma2 * pow(sigma, 2);
				}else {
					temp_val = temp_val + sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
				}
        		}
        		cor(i, j) = temp_val / n_wind;
        		cor(j, i) = cor(i, j);
    		}
  	}
  	return cor;
}

// [[Rcpp::export]]

NumericMatrix ORIG_SPATIALLY_VARYING_PARAMETERS(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind, NumericVector & param_nonstat) {
  
  const int m = Loc.nrow(), n_wind = wind.nrow();
  
  float sigma2 = param(0), beta = param(1), nu = param(2);
  
  NumericMatrix cor(m, m);
  
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
 
 	double temp_val = 0.0;

    	for (int k = 0; k < n_wind; ++k) {
      
	      double new_loc1_x = Loc(i, 0) - wind(k, 0) * Loc(i, 2);
	      double new_loc1_y = Loc(i, 1) - wind(k, 1) * Loc(i, 2);
	      double new_loc2_x = Loc(j, 0) - wind(k, 0) * Loc(j, 2);
	      double new_loc2_y = Loc(j, 1) - wind(k, 1) * Loc(j, 2);

	      double x1 = new_loc1_x - new_loc2_x;
	      double x2 = new_loc1_y - new_loc2_y;
	      NumericVector z2 = {x1, x2};
		
	      double omega1 = param_nonstat(0) + param_nonstat(1) * (new_loc1_x - .5) + param_nonstat(2) * (new_loc1_y - .5) + param_nonstat(3) * pow(new_loc1_x - .5, 2) + param_nonstat(4) * pow(new_loc1_y - .5, 2);
	      double log_lam1_1 = param_nonstat(5) + param_nonstat(6) * (new_loc1_x - .5) + param_nonstat(7) * (new_loc1_y - .5) + param_nonstat(8) * pow(new_loc1_x - .5, 2) + param_nonstat(9) * pow(new_loc1_y - .5, 2);
	      double log_lam1_2 = param_nonstat(10) + param_nonstat(11) * (new_loc1_x - .5) + param_nonstat(12) * (new_loc1_y - .5) + param_nonstat(13) * pow(new_loc1_x - .5, 2) + param_nonstat(14) * pow(new_loc1_y - .5, 2);
		
	      double omega2 = param_nonstat(0) + param_nonstat(1) * (new_loc2_x - .5) + param_nonstat(2) * (new_loc2_y - .5) + param_nonstat(3) * pow(new_loc2_x - .5, 2) + param_nonstat(4) * pow(new_loc2_y - .5, 2);
	      double log_lam2_1 = param_nonstat(5) + param_nonstat(6) * (new_loc2_x - .5) + param_nonstat(7) * (new_loc2_y - .5) + param_nonstat(8) * pow(new_loc2_x - .5, 2) + param_nonstat(9) * pow(new_loc2_y - .5, 2);
	      double log_lam2_2 = param_nonstat(10) + param_nonstat(11) * (new_loc2_x - .5) + param_nonstat(12) * (new_loc2_y - .5) + param_nonstat(13) * pow(new_loc2_x - .5, 2) + param_nonstat(14) * pow(new_loc2_y - .5, 2);

	      double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
	      double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
	      double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);
	      
	      double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
	      double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
	      double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);
	      
	      double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
	      double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;
	      
	      double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
	      double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
	      double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);
	      
	      double Inv_ij_11 = Kernel_ij_22; 
	      double Inv_ij_22 = Kernel_ij_11;
	      double Inv_ij_12 = - Kernel_ij_12; 
	      double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;
	      
	      double sigma = sqrt(sqrt(det_i * det_j)/det_ij);
	      
	      double dist = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22);


	      if (dist == 0) {
		temp_val = temp_val + sigma2 * pow(sigma, 2);
	      } else {
		temp_val = temp_val + sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
	      }
//	      Rprintf("dist : %f, sigma : %f, sigma2 : %f, nu : %f, cov value : %f, tgamma(nu) : %f \n", dist, sigma, sigma2, nu, temp_val, tgamma(nu));
	}
	cor(i, j) = temp_val / n_wind;
      	cor(j, i) = cor(i, j);
    }
  }
  return cor;
}

// [[Rcpp::export]]

NumericMatrix DEFORMATION(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind, NumericVector & param_nonstat) {
  
  const int m = Loc.nrow(), n_wind = wind.nrow();
  
  float sigma2 = param(0), beta = param(1), nu = param(2);
  
  NumericMatrix cor(m, m);
  
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
 
 	double temp_val = 0.0;

    	for (int k = 0; k < n_wind; ++k) {
      
	      double new_loc1_x = Loc(i, 0) - wind(k, 0) * Loc(i, 2);
	      double new_loc1_y = Loc(i, 1) - wind(k, 1) * Loc(i, 2);
	      double new_loc2_x = Loc(j, 0) - wind(k, 0) * Loc(j, 2);
	      double new_loc2_y = Loc(j, 1) - wind(k, 1) * Loc(j, 2);

		double dist_source1 = sqrt(pow(new_loc1_x - param_nonstat(3), 2) + pow(new_loc1_y - param_nonstat(4), 2)); 
		double dist_source2 = sqrt(pow(new_loc2_x - param_nonstat(3), 2) + pow(new_loc2_y - param_nonstat(4), 2)); 
		double deform_loc1_x = param_nonstat(3) + (new_loc1_x - param_nonstat(3)) * ( param_nonstat(0) + param_nonstat(1) * exp(-param_nonstat(2) * pow(dist_source1, 2 )));
		double deform_loc1_y = param_nonstat(4) + (new_loc1_y - param_nonstat(4)) * ( param_nonstat(0) + param_nonstat(1) * exp(-param_nonstat(2) * pow(dist_source1, 2 )));
		double deform_loc2_x = param_nonstat(3) + (new_loc2_x - param_nonstat(3)) * ( param_nonstat(0) + param_nonstat(1) * exp(-param_nonstat(2) * pow(dist_source2, 2 )));
		double deform_loc2_y = param_nonstat(4) + (new_loc2_y - param_nonstat(4)) * ( param_nonstat(0) + param_nonstat(1) * exp(-param_nonstat(2) * pow(dist_source2, 2 )));

		double dist = sqrt(pow(deform_loc2_x - deform_loc1_x, 2) + pow(deform_loc2_y - deform_loc1_y, 2));

	      if (dist == 0) {
		temp_val = temp_val + sigma2;
	      } else {
		temp_val = temp_val + sigma2 * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
	      }
	}
	cor(i, j) = temp_val / n_wind;
      	cor(j, i) = cor(i, j);
    }
  }
  return cor;
}

// [[Rcpp::export]]

NumericMatrix spatially_varying_parameters(NumericMatrix & Loc, NumericVector & param, NumericVector & v, int & max_T) {

	const int m = Loc.nrow();

	double sigma2 = param(0), nu = param(1);

	NumericMatrix cor(m, m);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j <= i; ++j) {

		double s1x = Loc(i, 0);
		double s2x = Loc(j, 0);
		double s1y = Loc(i, 1);
		double s2y = Loc(j, 1);
		int t1 = Loc(i, 2);
		int t2 = Loc(j, 2);

		double new_s1x = s1x - v[0] * t1;
		double new_s1y = s1y - v[1] * t1;
		double new_s2x = s2x - v[0] * t2;
		double new_s2y = s2y - v[1] * t2;
		double hx = new_s1x - new_s2x;
		double hy = new_s1y - new_s2y;

		double omega1 = (1 / pow(max_T, 2)) * (new_s1x - .5 * max_T) - (2 / pow(max_T, 2)) * (new_s1y - .5 * max_T) + (0 / pow(max_T, 2)) * pow(new_s1x - .5 * max_T, 2) + (1 / pow(max_T, 2)) * pow(new_s1y - .5 * max_T, 2);
		double log_lam1_1 = (-3 / pow(max_T, 2)) + (0 / pow(max_T, 2)) * (new_s1x - .5 * max_T) + (0 / pow(max_T, 2)) * (new_s1y - .5 * max_T) - (6 / pow(max_T, 2)) * pow(new_s1x - .5 * max_T, 2) - (7 / pow(max_T, 2)) * pow(new_s1y - .5 * max_T, 2);
		double log_lam1_2 = (-5 / pow(max_T, 2)) + (0 / pow(max_T, 2)) * (new_s1x - .5 * max_T) + (0 / pow(max_T, 2)) * (new_s1y - .5 * max_T) + (1 / pow(max_T, 2)) * pow(new_s1x - .5 * max_T, 2) - (4 / pow(max_T, 2)) * pow(new_s1y - .5 * max_T, 2);

		double omega2 = (1 / pow(max_T, 2)) * (new_s2x - .5 * max_T) - (2 / pow(max_T, 2)) * (new_s2y - .5 * max_T) + (0 / pow(max_T, 2)) * pow(new_s2x - .5 * max_T, 2) + (1 / pow(max_T, 2)) * pow(new_s2y - .5 * max_T, 2);
		double log_lam2_1 = (-3 / pow(max_T, 2)) + (0 / pow(max_T, 2)) * (new_s2x - .5 * max_T) + (0 / pow(max_T, 2)) * (new_s2y - .5 * max_T) - (6 / pow(max_T, 2)) * pow(new_s2x - .5 * max_T, 2) - (7 / pow(max_T, 2)) * pow(new_s2y - .5 * max_T, 2);
		double log_lam2_2 = (-5 / pow(max_T, 2)) + (0 / pow(max_T, 2)) * (new_s2x - .5 * max_T) + (0 / pow(max_T, 2)) * (new_s2y - .5 * max_T) + (1 / pow(max_T, 2)) * pow(new_s2x - .5 * max_T, 2) - (4 / pow(max_T, 2)) * pow(new_s2y - .5 * max_T, 2);

		//double omega1 = 1 * (new_s1x - .5) + 1 * (new_s1y - .5) + 1 * pow(new_s1x - .5, 2) + 1 * pow(new_s1y - .5, 2);
		//double log_lam1_1 = 3 * (new_s1x - .5) + - 3 * (new_s1y - .5) + 0 * pow(new_s1x - .5, 2) + -3 * pow(new_s1y - .5, 2);
		//double log_lam1_2 = - 3 * (new_s1x - .5) + 0 * (new_s1y - .5) + - 1 * pow(new_s1x - .5, 2) + 3 * pow(new_s1y - .5, 2);

		//double omega2 = 1 * (new_s2x - .5) + 1 * (new_s2y - .5) + 1 * pow(new_s2x - .5, 2) + 1 * pow(new_s2y - .5, 2);
		//double log_lam2_1 = 3 * (new_s2x - .5) + - 3 * (new_s2y - .5) + 0 * pow(new_s2x - .5, 2) + -3 * pow(new_s2y - .5, 2);
		//double log_lam2_2 = - 3 * (new_s2x - .5) + 0 * (new_s2y - .5) + - 1 * pow(new_s2x - .5, 2) + 3 * pow(new_s2y - .5, 2);

		double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
		double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
		double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

		double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
		double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
		double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

		double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
		double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

		double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
		double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
		double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

		double Inv_ij_11 = Kernel_ij_22; 
		double Inv_ij_22 = Kernel_ij_11;
		double Inv_ij_12 = - Kernel_ij_12; 
		double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

		double sigma = sqrt(sqrt(det_i * det_j)/det_ij);

		double dist = sqrt((hx * hx * Inv_ij_11 + 2 * hx * hy * Inv_ij_12 + hy * hy * Inv_ij_22) / det_ij);

		if (dist == 0) {
			cor(i, j) = sigma2 * pow(sigma, 2);
		} else {
			cor(i, j) = sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
			//cor(i, j) = sigma2 * pow(sigma, 2) * exp(-pow(dist, 2));
		}
		cor(j, i) = cor(i, j);
		}
	}
	return cor;
}

// [[Rcpp::export]]

NumericMatrix spatially_varying_parameters_for_fitting(NumericMatrix & Loc, NumericVector & param, NumericVector & v, NumericMatrix & Nonstat_params) {

	const int m = Loc.nrow();

	double sigma2 = param(0), nu = param(1);

	NumericMatrix cor(m, m);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j <= i; ++j) {

		double s1x = Loc(i, 0);
		double s2x = Loc(j, 0);
		double s1y = Loc(i, 1);
		double s2y = Loc(j, 1);
		int t1 = Loc(i, 2);
		int t2 = Loc(j, 2);

		double new_s1x = s1x - v[0] * t1;
		double new_s1y = s1y - v[1] * t1;
		double new_s2x = s2x - v[0] * t2;
		double new_s2y = s2y - v[1] * t2;
		double hx = new_s1x - new_s2x;
		double hy = new_s1y - new_s2y;

		double omega1 = Nonstat_params(i, 0);
		double log_lam1_1 = Nonstat_params(i, 1);
		double log_lam1_2 = Nonstat_params(i, 2);

		double omega2 = Nonstat_params(j, 0);
		double log_lam2_1 = Nonstat_params(j, 1);
		double log_lam2_2 = Nonstat_params(j, 2);

		double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
		double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
		double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

		double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
		double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
		double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

		double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
		double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

		double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
		double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
		double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

		double Inv_ij_11 = Kernel_ij_22; 
		double Inv_ij_22 = Kernel_ij_11;
		double Inv_ij_12 = - Kernel_ij_12; 
		double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

		double sigma = sqrt(sqrt(det_i * det_j)/det_ij);

		double dist = sqrt((hx * hx * Inv_ij_11 + 2 * hx * hy * Inv_ij_12 + hy * hy * Inv_ij_22) / det_ij);

		if (dist == 0) {
			cor(i, j) = sigma2 * pow(sigma, 2);
		} else {
			cor(i, j) = sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
		}
		cor(j, i) = cor(i, j);
		}
	}
	return cor;
}

// [[Rcpp::export]]

NumericMatrix spatially_varying_parameters_for_fitting2(NumericMatrix & Loc, NumericVector & param, NumericVector & v, NumericMatrix & Nonstat_params) {

	const int m = Loc.nrow();

	double sigma2 = param(0), range = param(1), nu = param(2);

	NumericMatrix cor(m, m);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j <= i; ++j) {

		double s1x = Loc(i, 0);
		double s2x = Loc(j, 0);
		double s1y = Loc(i, 1);
		double s2y = Loc(j, 1);
		int t1 = Loc(i, 2);
		int t2 = Loc(j, 2);

		double new_s1x = s1x - v[0] * t1;
		double new_s1y = s1y - v[1] * t1;
		double new_s2x = s2x - v[0] * t2;
		double new_s2y = s2y - v[1] * t2;
		double hx = new_s1x - new_s2x;
		double hy = new_s1y - new_s2y;

		double omega1 = Nonstat_params(i, 0);
		double log_lam1_1 = Nonstat_params(i, 1);
		double log_lam1_2 = Nonstat_params(i, 2);

		double omega2 = Nonstat_params(j, 0);
		double log_lam2_1 = Nonstat_params(j, 1);
		double log_lam2_2 = Nonstat_params(j, 2);

		double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
		double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
		double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

		double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
		double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
		double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

		double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
		double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

		double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
		double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
		double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

		double Inv_ij_11 = Kernel_ij_22; 
		double Inv_ij_22 = Kernel_ij_11;
		double Inv_ij_12 = - Kernel_ij_12; 
		double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

		double sigma = sqrt(sqrt(det_i * det_j)/det_ij);

		double dist = sqrt((hx * hx * Inv_ij_11 + 2 * hx * hy * Inv_ij_12 + hy * hy * Inv_ij_22) / det_ij) / range;

		if (dist == 0) {
			cor(i, j) = sigma2 * pow(sigma, 2);
		} else {
			cor(i, j) = sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
		}
		cor(j, i) = cor(i, j);
		}
	}
	return cor;
}

// [[Rcpp::export]]

NumericMatrix spatially_varying_parameters_for_fitting3(NumericMatrix & Loc1, NumericMatrix & Loc2, NumericVector & param, NumericVector & v, NumericMatrix & Nonstat_params1, NumericMatrix & Nonstat_params2) {

	const int m = Loc1.nrow();

	double sigma2 = param(0), range = param(1), nu = param(2);

	NumericMatrix cor(m, m);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {

		double s1x = Loc1(i, 0);
		double s2x = Loc2(j, 0);
		double s1y = Loc1(i, 1);
		double s2y = Loc2(j, 1);
		int t1 = Loc1(i, 2);
		int t2 = Loc2(j, 2);

		double new_s1x = s1x - v[0] * t1;
		double new_s1y = s1y - v[1] * t1;
		double new_s2x = s2x - v[0] * t2;
		double new_s2y = s2y - v[1] * t2;
		double hx = new_s1x - new_s2x;
		double hy = new_s1y - new_s2y;

		double omega1 = Nonstat_params1(i, 0);
		double log_lam1_1 = Nonstat_params1(i, 1);
		double log_lam1_2 = Nonstat_params1(i, 2);

		double omega2 = Nonstat_params2(j, 0);
		double log_lam2_1 = Nonstat_params2(j, 1);
		double log_lam2_2 = Nonstat_params2(j, 2);

		double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
		double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
		double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

		double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
		double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
		double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

		double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
		double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

		double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
		double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
		double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

		double Inv_ij_11 = Kernel_ij_22; 
		double Inv_ij_22 = Kernel_ij_11;
		double Inv_ij_12 = - Kernel_ij_12; 
		double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

		double sigma = sqrt(sqrt(det_i * det_j)/det_ij);

		double dist = sqrt((hx * hx * Inv_ij_11 + 2 * hx * hy * Inv_ij_12 + hy * hy * Inv_ij_22) / det_ij) / range;

		if (dist == 0) {
			cor(i, j) = sigma2 * pow(sigma, 2);
		} else {
			cor(i, j) = sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
		}
		//cor(j, i) = cor(i, j);
		}
	}
	return cor;
}

// [[Rcpp::export]]

NumericMatrix nonfrozen_rcpp(NumericMatrix & Loc1, NumericMatrix & Loc2, NumericVector & param, NumericVector & v_mean, NumericMatrix & v_var) {

	const int m = Loc1.nrow();

	double sigma2 = param(0), range = param(1), nu = param(2);

	NumericMatrix cor(m, m);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {

		double s1x = Loc1(i, 0);
		double s2x = Loc2(j, 0);
		double s1y = Loc1(i, 1);
		double s2y = Loc2(j, 1);
		int t1 = Loc1(i, 2);
		int t2 = Loc2(j, 2);
		int u_2 = pow(t1 - t2, 2);

		double new_s1x = s1x - v_mean[0] * t1;
		double new_s1y = s1y - v_mean[1] * t1;
		double new_s2x = s2x - v_mean[0] * t2;
		double new_s2y = s2y - v_mean[1] * t2;
		double hx = new_s1x - new_s2x;
		double hy = new_s1y - new_s2y;

		double Kernel_ij_11 = 1 + u_2 * v_var(0, 0);
		double Kernel_ij_12 = u_2 * v_var(0, 1);
		double Kernel_ij_22 = 1 + u_2 * v_var(1, 1);

		double Inv_ij_11 = Kernel_ij_22; 
		double Inv_ij_22 = Kernel_ij_11;
		double Inv_ij_12 = - Kernel_ij_12; 
		double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

		double dist = sqrt((hx * hx * Inv_ij_11 + 2 * hx * hy * Inv_ij_12 + hy * hy * Inv_ij_22) / det_ij) / range;
	
		//Rprintf("u : %d, dist : %f \n", u_2, dist);

		if (dist == 0) {
			cor(i, j) = sigma2 / pow(det_ij, 0.5);
		} else {
			cor(i, j) = sigma2 * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu) * pow(det_ij, 0.5));
		}
		}
	}
	return cor;
}


// [[Rcpp::export]]

NumericMatrix nonfrozen_rcpp_LS(NumericMatrix & Loc1, NumericMatrix & Loc2, NumericVector & param, NumericVector & v_mean, NumericMatrix & v_var) {

	const int m = Loc1.nrow();
	const int n = Loc2.nrow();

	double sigma2 = param(0), range = param(1), nu = param(2);

	NumericMatrix cor(m, n);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {

		double s1x = Loc1(i, 0);
		double s2x = Loc2(j, 0);
		double s1y = Loc1(i, 1);
		double s2y = Loc2(j, 1);
		int t1 = Loc1(i, 2);
		int t2 = Loc2(j, 2);
		int u_2 = pow(t1 - t2, 2);

		double new_s1x = s1x - v_mean[0] * t1;
		double new_s1y = s1y - v_mean[1] * t1;
		double new_s2x = s2x - v_mean[0] * t2;
		double new_s2y = s2y - v_mean[1] * t2;
		double hx = new_s1x - new_s2x;
		double hy = new_s1y - new_s2y;

		double Kernel_ij_11 = 1 + u_2 * v_var(0, 0);
		double Kernel_ij_12 = u_2 * v_var(0, 1);
		double Kernel_ij_22 = 1 + u_2 * v_var(1, 1);

		double Inv_ij_11 = Kernel_ij_22; 
		double Inv_ij_22 = Kernel_ij_11;
		double Inv_ij_12 = - Kernel_ij_12; 
		double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

		double dist = sqrt((hx * hx * Inv_ij_11 + 2 * hx * hy * Inv_ij_12 + hy * hy * Inv_ij_22) / det_ij) / range;
	
		//Rprintf("u : %d, dist : %f \n", u_2, dist);

		if (dist == 0) {
			cor(i, j) = sigma2 / pow(det_ij, 0.5);
		} else {
			cor(i, j) = sigma2 * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu) * pow(det_ij, 0.5));
		}
		}
	}
	return cor;
}

// [[Rcpp::export]]

NumericMatrix nonfrozen_rcpp_added_dimensions(NumericMatrix & Loc1, NumericMatrix & Loc2, NumericVector & param, NumericVector & v_mean, NumericMatrix & v_var) {

	const int m = Loc1.nrow();

	double sigma2 = param(0), range = param(1), nu = param(2);

	Eigen::MatrixXd Sigma_wind = as<Eigen::MatrixXd>(v_var);
	NumericMatrix cor(m, m);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {

		double s1x = Loc1(i, 0);
		double s2x = Loc2(j, 0);
		double s1y = Loc1(i, 1);
		double s2y = Loc2(j, 1);
		double s1z = Loc1(i, 2);
		double s2z = Loc2(j, 2);

		int t1 = Loc1(i, 3);
		int t2 = Loc2(j, 3);
		int u_2 = pow(t1 - t2, 2);

		double new_s1x = s1x - v_mean[0] * t1;
		double new_s1y = s1y - v_mean[1] * t1;
		double new_s1z = s1z - v_mean[2] * t1;
		double new_s2x = s2x - v_mean[0] * t2;
		double new_s2y = s2y - v_mean[1] * t2;
		double new_s2z = s2z - v_mean[2] * t2;

		Eigen::VectorXd h(3);		

		h[0] = new_s1x - new_s2x;
		h[1] = new_s1y - new_s2y;
		h[2] = new_s1z - new_s2z;

		Eigen::MatrixXd Sigma_new = Sigma_wind * u_2 + Eigen::MatrixXd::Identity(3, 3);
		Eigen::MatrixXd Sigma_tilde = Sigma_new.inverse();

		double det_ij = (Sigma_new).determinant();
		double dist = sqrt(h.transpose() * Sigma_tilde * h) / range;

		//Rprintf("u : %d, dist : %f \n", u_2, dist);

		if (dist == 0) {
			cor(i, j) = sigma2 / pow(det_ij, 0.5);
		} else {
			cor(i, j) = sigma2 * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu) * pow(det_ij, 0.5));
		}
		}
	}
	return cor;
}

// [[Rcpp::export]]

NumericMatrix nonfrozen_rcpp_multi_cross_added_dimensions(NumericMatrix & Loc1, NumericMatrix & Loc2, NumericVector & param, NumericVector & v_mean, NumericMatrix & v_var) {

	const int m = Loc1.nrow();

	double sigma2 = param(0), range = param(1), nu = param(2);

	Eigen::MatrixXd Sigma(3, 3);
	Eigen::MatrixXd Sigma_wind = as<Eigen::MatrixXd>(v_var);

	NumericMatrix cor(m, m);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {

		double s1x = Loc1(i, 0);
		double s2x = Loc2(j, 0);
		double s1y = Loc1(i, 1);
		double s2y = Loc2(j, 1);
		double s1z = Loc1(i, 2);
		double s2z = Loc2(j, 2);

		int t1 = Loc1(i, 3);
		int t2 = Loc2(j, 3);

		Eigen::MatrixXd t_mat(3, 6);
		t_mat.leftCols(3) = -t1 * Eigen::MatrixXd::Identity(3, 3);
		t_mat.rightCols(3) = t2 * Eigen::MatrixXd::Identity(3, 3);

		Eigen::MatrixXd Sigma_new = t_mat.transpose() * t_mat + Sigma_wind.inverse(); 
		Eigen::MatrixXd Sigma_new_inv = Sigma_new.inverse(); 
		Eigen::MatrixXd Sigma_tilde = Eigen::MatrixXd::Identity(3, 3) - t_mat * Sigma_new_inv * t_mat.transpose();

		Eigen::MatrixXd Sigma_for_determinant = Sigma_wind * t_mat.transpose() * t_mat + Eigen::MatrixXd::Identity(6, 6);
		//double det_ij = (2 * Sigma_tilde).inverse().determinant();
		double det_ij = (Sigma_for_determinant).determinant();

		//Rprintf("det : %f \n", det_ij);

		double new_s1x = s1x - v_mean[0] * t1;
		double new_s1y = s1y - v_mean[1] * t1;
		double new_s1z = s1z - v_mean[2] * t1;
		double new_s2x = s2x - v_mean[3] * t2;
		double new_s2y = s2y - v_mean[4] * t2;
		double new_s2z = s2z - v_mean[5] * t2;

		Eigen::VectorXd h(3);		

		h[0] = new_s1x - new_s2x;
		h[1] = new_s1y - new_s2y;
		h[2] = new_s1z - new_s2z;
			
		double dist = sqrt(h.transpose() * Sigma_tilde * h) / range;
		//Rprintf("det : %f \n", det_ij);
		//std::cout << " s1z: " << s1z << " v_mean[2] * t1: " << v_mean[2] * t1 << "\n" << std::endl;
		//std::cout << " s2z: " << s2z << " v_mean[5] * t2: " << v_mean[5] * t2 << "\n" << std::endl;
		//std::cout << " h: " << h  << "\n" <<  std::endl;
	
		if (dist == 0) {
			cor(i, j) = sigma2 / pow(det_ij, 0.5);
		} else {
			cor(i, j) = sigma2 * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu) * pow(det_ij, 0.5));
		}

		}
	}
	return cor;
}

// [[Rcpp::export]]

NumericMatrix nonfrozen_rcpp_multi_cross_LS(NumericMatrix & Loc1, NumericMatrix & Loc2, NumericVector & param, NumericVector & v_mean, NumericMatrix & v_var) {

	const int m = Loc1.nrow();
	const int n = Loc2.nrow();

	double sigma2 = param(0), range = param(1), nu = param(2);

	Eigen::MatrixXd Sigma(2, 2);
	Eigen::MatrixXd Sigma_wind = as<Eigen::MatrixXd>(v_var);

	NumericMatrix cor(m, n);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {

		double s1x = Loc1(i, 0);
		double s2x = Loc2(j, 0);
		double s1y = Loc1(i, 1);
		double s2y = Loc2(j, 1);
		int t1 = Loc1(i, 2);
		int t2 = Loc2(j, 2);

		Eigen::MatrixXd t_mat(2, 4);
		t_mat.leftCols(2) = -t1 * Eigen::MatrixXd::Identity(2, 2);
		t_mat.rightCols(2) = t2 * Eigen::MatrixXd::Identity(2, 2);

		Eigen::MatrixXd Sigma_new = t_mat.transpose() * t_mat + Sigma_wind.inverse(); 
		Eigen::MatrixXd Sigma_new_inv = Sigma_new.inverse(); 
		Eigen::MatrixXd Sigma_tilde = Eigen::MatrixXd::Identity(2, 2) - t_mat * Sigma_new_inv * t_mat.transpose();

		Eigen::MatrixXd Sigma_for_determinant = Sigma_wind * t_mat.transpose() * t_mat + Eigen::MatrixXd::Identity(4, 4);
		//double det_ij = (2 * Sigma_tilde).inverse().determinant();
		double det_ij = (Sigma_for_determinant).determinant();

		//Rprintf("det : %f \n", det_ij);

		double new_s1x = s1x - v_mean[0] * t1;
		double new_s1y = s1y - v_mean[1] * t1;
		double new_s2x = s2x - v_mean[2] * t2;
		double new_s2y = s2y - v_mean[3] * t2;

		Eigen::VectorXd h(2);		

		h[0] = new_s1x - new_s2x;
		h[1] = new_s1y - new_s2y;
			
		double dist = sqrt(h.transpose() * Sigma_tilde * h) / range;
	
		if (dist == 0) {
			cor(i, j) = sigma2 / pow(det_ij, 0.5);
		} else {
			cor(i, j) = sigma2 * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu) * pow(det_ij, 0.5));
		}
		//std::cout << "Sigma_for_determinant: \n" << Sigma_for_determinant << std::endl;
		//std::cout << "det_ij: \n" << det_ij << std::endl;
		}
	}
	return cor;
}

// [[Rcpp::export]]

NumericMatrix nonfrozen_rcpp_multi_cross(NumericMatrix & Loc1, NumericMatrix & Loc2, NumericVector & param, NumericVector & v_mean, NumericMatrix & v_var) {

	const int m = Loc1.nrow();

	double sigma2 = param(0), range = param(1), nu = param(2);

	Eigen::MatrixXd Sigma(2, 2);
	Eigen::MatrixXd Sigma_wind = as<Eigen::MatrixXd>(v_var);

	NumericMatrix cor(m, m);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {

		double s1x = Loc1(i, 0);
		double s2x = Loc2(j, 0);
		double s1y = Loc1(i, 1);
		double s2y = Loc2(j, 1);
		int t1 = Loc1(i, 2);
		int t2 = Loc2(j, 2);

		Eigen::MatrixXd t_mat(2, 4);
		t_mat.leftCols(2) = -t1 * Eigen::MatrixXd::Identity(2, 2);
		t_mat.rightCols(2) = t2 * Eigen::MatrixXd::Identity(2, 2);

		Eigen::MatrixXd Sigma_new = t_mat.transpose() * t_mat + Sigma_wind.inverse(); 
		Eigen::MatrixXd Sigma_new_inv = Sigma_new.inverse(); 
		Eigen::MatrixXd Sigma_tilde = Eigen::MatrixXd::Identity(2, 2) - t_mat * Sigma_new_inv * t_mat.transpose();

		Eigen::MatrixXd Sigma_for_determinant = Sigma_wind * t_mat.transpose() * t_mat + Eigen::MatrixXd::Identity(4, 4);
		//double det_ij = (2 * Sigma_tilde).inverse().determinant();
		double det_ij = (Sigma_for_determinant).determinant();

		//Rprintf("det : %f \n", det_ij);

		double new_s1x = s1x - v_mean[0] * t1;
		double new_s1y = s1y - v_mean[1] * t1;
		double new_s2x = s2x - v_mean[2] * t2;
		double new_s2y = s2y - v_mean[3] * t2;

		Eigen::VectorXd h(2);		

		h[0] = new_s1x - new_s2x;
		h[1] = new_s1y - new_s2y;
			
		double dist = sqrt(h.transpose() * Sigma_tilde * h) / range;
	
		if (dist == 0) {
			cor(i, j) = sigma2 / pow(det_ij, 0.5);
		} else {
			cor(i, j) = sigma2 * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu) * pow(det_ij, 0.5));
		}
		//std::cout << "Sigma_for_determinant: \n" << Sigma_for_determinant << std::endl;
		//std::cout << "det_ij: \n" << det_ij << std::endl;
		}
	}
	return cor;
}

// [[Rcpp::export]]

NumericMatrix nonfrozen_rcpp_NEW(NumericMatrix & Loc1, NumericMatrix & Loc2, NumericVector & param, NumericVector & v_mean1, NumericVector & v_mean2, NumericMatrix & v_var1, NumericMatrix & v_var2, NumericMatrix & v_var12) {

	const int m = Loc1.nrow();

	double sigma2 = param(0), range = param(1), nu = param(2);

	Eigen::VectorXd mean_wind1 = as<Eigen::VectorXd>(v_mean1);
	Eigen::VectorXd mean_wind2 = as<Eigen::VectorXd>(v_mean2);

	Eigen::MatrixXd Sigma_wind1 = as<Eigen::MatrixXd>(v_var1);
	Eigen::MatrixXd Sigma_wind2 = as<Eigen::MatrixXd>(v_var2);
	Eigen::MatrixXd Sigma_wind12 = as<Eigen::MatrixXd>(v_var12);

	NumericMatrix cor(m, m);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {

		double s1x = Loc1(i, 0);
		double s2x = Loc2(j, 0);
		double s1y = Loc1(i, 1);
		double s2y = Loc2(j, 1);
		int t1 = Loc1(i, 2);
		int t2 = Loc2(j, 2);


		double new_s1x = s1x - mean_wind1[0] * t1;
		double new_s1y = s1y - mean_wind1[1] * t1;
		double new_s2x = s2x - mean_wind1[0] * t2;
		double new_s2y = s2y - mean_wind1[1] * t2;

		Eigen::VectorXd h(2), h_tilde(2), h_new(2);		
		Eigen::MatrixXd u_tilde(2, 2), Sigma_check(2, 2), Sigma_check1(2, 2), Sigma_tilde(2, 2), Sigma_temp(2, 2), Sigma_new(2, 2), det_temp(2, 2);	

		h[0] = new_s1x - new_s2x;
		h[1] = new_s1y - new_s2y;

		h_tilde = h + Sigma_wind12 * Sigma_wind2.inverse() * mean_wind2 * (t1 - t2);
		u_tilde = Sigma_wind12 * Sigma_wind2.inverse() * (t1 - t2);
		h_new = h_tilde - u_tilde * mean_wind2;

		Sigma_tilde = Sigma_wind1 - Sigma_wind12 * Sigma_wind2.inverse() * Sigma_wind12.transpose();
		Sigma_check = Sigma_tilde * pow(t1 - t2, 2) + Eigen::MatrixXd::Identity(2, 2);
		Sigma_check1 = Sigma_check.inverse();
		Sigma_temp = u_tilde.transpose() * Sigma_check1 * u_tilde + Sigma_wind2.inverse();
		Sigma_new = Sigma_check1 - Sigma_check1 * u_tilde * Sigma_temp.inverse() * u_tilde.transpose() * Sigma_check1;

		double det1 = Sigma_check.determinant();
		det_temp = Sigma_wind2 * u_tilde.transpose() * Sigma_check1 * u_tilde + Eigen::MatrixXd::Identity(2, 2);
		double det2 = det_temp.determinant();
			
		double dist = sqrt(h_new.transpose() * Sigma_new * h_new) / range;
	
		if (dist == 0) {
			cor(i, j) = sigma2 / (pow(det1, 0.5) * pow(det2, 0.5));
		} else {
			cor(i, j) = sigma2 * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu) * pow(det1, 0.5) * pow(det2, 0.5));
		}

		}
	}
	return cor;
}

// [[Rcpp::export]]

NumericMatrix spatially_varying_parameters_for_fitting2_VELOCITY_FIELD(NumericMatrix & Loc, NumericVector & param, NumericMatrix & v, NumericMatrix & Nonstat_params) {

	const int m = Loc.nrow();

	double sigma2 = param(0), range = param(1), nu = param(2);

	NumericMatrix cor(m, m);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j <= i; ++j) {

		double s1x = Loc(i, 0);
		double s2x = Loc(j, 0);
		double s1y = Loc(i, 1);
		double s2y = Loc(j, 1);
		int t1 = Loc(i, 2);
		int t2 = Loc(j, 2);

		double new_s1x = s1x - v(i, 0) * t1;
		double new_s1y = s1y - v(i, 1) * t1;
		double new_s2x = s2x - v(j, 0) * t2;
		double new_s2y = s2y - v(j, 1) * t2;
		double hx = new_s1x - new_s2x;
		double hy = new_s1y - new_s2y;

		double omega1 = Nonstat_params(i, 0);
		double log_lam1_1 = Nonstat_params(i, 1);
		double log_lam1_2 = Nonstat_params(i, 2);

		double omega2 = Nonstat_params(j, 0);
		double log_lam2_1 = Nonstat_params(j, 1);
		double log_lam2_2 = Nonstat_params(j, 2);

		double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
		double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
		double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

		double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
		double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
		double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

		double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
		double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

		double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
		double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
		double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

		double Inv_ij_11 = Kernel_ij_22; 
		double Inv_ij_22 = Kernel_ij_11;
		double Inv_ij_12 = - Kernel_ij_12; 
		double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

		double sigma = sqrt(sqrt(det_i * det_j)/det_ij);

		double dist = sqrt((hx * hx * Inv_ij_11 + 2 * hx * hy * Inv_ij_12 + hy * hy * Inv_ij_22) / det_ij) / range;

		if (dist == 0) {
			cor(i, j) = sigma2 * pow(sigma, 2);
		} else {
			cor(i, j) = sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
		}
		cor(j, i) = cor(i, j);
		}
	}
	return cor;
}


// [[Rcpp::export]]

NumericMatrix Gneiting_spatially_varying_parameters_for_fitting(NumericMatrix & Loc1, NumericMatrix & Loc2, NumericVector & param, NumericMatrix & Nonstat_params1, NumericMatrix & Nonstat_params2) {

	const int m = Loc1.nrow();

	double sigma2 = param(0), range = param(1), nu = param(2), range_time = param(3), nu_time = param(4), beta = param(5);

	NumericMatrix cor(m, m);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {

		double s1x = Loc1(i, 0);
		double s2x = Loc2(j, 0);
		double s1y = Loc1(i, 1);
		double s2y = Loc2(j, 1);
		int t1 = Loc1(i, 2);
		int t2 = Loc2(j, 2);

		double hx = s1x - s2x;
		double hy = s1y - s2y;
		int u = pow(abs(t1 - t2), 2 * nu_time) / range_time + 1;

		double omega1 = Nonstat_params1(i, 0);
		double log_lam1_1 = Nonstat_params1(i, 1);
		double log_lam1_2 = Nonstat_params1(i, 2);

		double omega2 = Nonstat_params2(j, 0);
		double log_lam2_1 = Nonstat_params2(j, 1);
		double log_lam2_2 = Nonstat_params2(j, 2);

		double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
		double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
		double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

		double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
		double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
		double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

		double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
		double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

		double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
		double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
		double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

		double Inv_ij_11 = Kernel_ij_22; 
		double Inv_ij_22 = Kernel_ij_11;
		double Inv_ij_12 = - Kernel_ij_12; 
		double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

		double sigma = sqrt(sqrt(det_i * det_j)/det_ij);

		double dist_s = sqrt((hx * hx * Inv_ij_11 + 2 * hx * hy * Inv_ij_12 + hy * hy * Inv_ij_22) / det_ij) / range;
		double dist = dist_s / pow(u, beta / 2);

		if (dist == 0) {
			cor(i, j) = sigma2 * pow(sigma, 2) / pow(u, beta);
		} else {
			cor(i, j) = sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu) * pow(u, beta));
		}
		}
	}
	return cor;
}

// [[Rcpp::export]]

NumericMatrix deformation(NumericMatrix & Loc, NumericVector & param, NumericVector & v) {

	const int m = Loc.nrow();

	double sigma2 = param(0), range = param(1), nu = param(2);

	NumericMatrix cor(m, m);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j <= i; ++j) {

		double s1x = Loc(i, 0);
		double s2x = Loc(j, 0);
		double s1y = Loc(i, 1);
		double s2y = Loc(j, 1);
		int t1 = Loc(i, 2);
		int t2 = Loc(j, 2);

		double new_s1x = s1x - v[0] * t1;
		double new_s1y = s1y - v[1] * t1;
		double new_s2x = s2x - v[0] * t2;
		double new_s2y = s2y - v[1] * t2;
		double hx = new_s1x - new_s2x;
		double hy = new_s1y - new_s2y;
	
		double temp_loc1 = sqrt(pow(new_s1x - 0.5, 2) + pow(new_s1y - 0.5, 2));
		double temp_loc2 = sqrt(pow(new_s2x - 0.5, 2) + pow(new_s2y - 0.5, 2));

		double new_new_s1x = 0.5 + (new_s1x - 0.5) * temp_loc1;
		double new_new_s1y = 0.5 + (new_s1y - 0.5) * temp_loc1;
		double new_new_s2x = 0.5 + (new_s2x - 0.5) * temp_loc2;
		double new_new_s2y = 0.5 + (new_s2y - 0.5) * temp_loc2;

		double dist = sqrt(pow(new_new_s1x - new_new_s2x, 2) + pow(new_new_s1y - new_new_s2y, 2)) / range;

		if (dist == 0) {
			cor(i, j) = sigma2;
		} else {
			cor(i, j) = sigma2 * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
		}
		cor(j, i) = cor(i, j);
		}
	}
	return cor;
}

// [[Rcpp::export]]

NumericMatrix MATERN(NumericMatrix & Loc, NumericVector & param) {

	const int m = Loc.nrow();

	double sigma2 = param(0), range = param(1), nu = param(2);

	NumericMatrix cor(m, m);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j <= i; ++j) {

		double s1x = Loc(i, 0);
		double s2x = Loc(j, 0);
		double s1y = Loc(i, 1);
		double s2y = Loc(j, 1);

		double dist = sqrt(pow(s1x - s2x, 2) + pow(s1y - s2y, 2)) / range;

		if (dist == 0) {
			cor(i, j) = sigma2;
		} else {
			cor(i, j) = sigma2 * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
		}
		cor(j, i) = cor(i, j);
		}
	}
	return cor;
}


// [[Rcpp::export]]

NumericMatrix spatially_varying_parameters_multi(NumericMatrix & Loc, NumericVector & param, NumericVector & v) {

	//supports same spatial nonstationarity for both variables.
	//TO DO different nonstationarity
  
	const int m = Loc.nrow();
  
	double sigma11 = param(0), sigma22 = param(1), nu1 = param(2), nu2 = param(3), rho = param(4);
  
	NumericMatrix cor11(m, m), cor22(m, m), cor12(m, m), cor21(m, m);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {

		double s1x = Loc(i, 0);
		double s2x = Loc(j, 0);
		double s1y = Loc(i, 1);
		double s2y = Loc(j, 1);
		int t1 = Loc(i, 2);
		int t2 = Loc(j, 2);

		double new_s1x = s1x - v[0] * t1;
		double new_s1y = s1y - v[1] * t1;
		double new_s2x = s2x - v[0] * t2;
		double new_s2y = s2y - v[1] * t2;
		double hx = new_s1x - new_s2x;
		double hy = new_s1y - new_s2y;

		double omega1 = 1 * (new_s1x - .5) + 1 * (new_s1y - .5) + 1 * pow(new_s1x - .5, 2) + 1 * pow(new_s1y - .5, 2);
		double log_lam1_1 = 3 * (new_s1x - .5) + - 3 * (new_s1y - .5) + 0 * pow(new_s1x - .5, 2) + -3 * pow(new_s1y - .5, 2);
		double log_lam1_2 = - 3 * (new_s1x - .5) + 0 * (new_s1y - .5) + - 1 * pow(new_s1x - .5, 2) + 3 * pow(new_s1y - .5, 2);

		double omega2 = 1 * (new_s2x - .5) + 1 * (new_s2y - .5) + 1 * pow(new_s2x - .5, 2) + 1 * pow(new_s2y - .5, 2);
		double log_lam2_1 = 3 * (new_s2x - .5) + - 3 * (new_s2y - .5) + 0 * pow(new_s2x - .5, 2) + -3 * pow(new_s2y - .5, 2);
		double log_lam2_2 = - 3 * (new_s2x - .5) + 0 * (new_s2y - .5) + - 1 * pow(new_s2x - .5, 2) + 3 * pow(new_s2y - .5, 2);

		double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
		double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
		double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

		double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
		double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
		double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

		double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
		double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

		double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
		double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
		double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

		double Inv_ij_11 = Kernel_ij_22; 
		double Inv_ij_22 = Kernel_ij_11;
		double Inv_ij_12 = - Kernel_ij_12; 
		double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

		double sigma = sqrt(sqrt(det_i * det_j)/det_ij);

		double dist = sqrt((hx * hx * Inv_ij_11 + 2 * hx * hy * Inv_ij_12 + hy * hy * Inv_ij_22) / det_ij);

		if (dist == 0) {
			cor11(i, j) = sigma11 * pow(sigma, 2);
			cor22(i, j) = sigma22 * pow(sigma, 2);
			cor12(i, j) = rho * sigma11 * sigma22 * pow(sigma, 2);
		} else {
			cor11(i, j) = sigma11 * pow(sigma, 2) * pow(dist, nu1) * cyl_bessel_k(nu1, dist) / (pow(2, nu1 - 1) * tgamma(nu1));
			cor22(i, j) = sigma22 * pow(sigma, 2) * pow(dist, nu2) * cyl_bessel_k(nu2, dist) / (pow(2, nu2 - 1) * tgamma(nu2));
			cor12(i, j) = rho * sigma11 * sigma22 * pow(sigma, 2) * pow(dist, 0.5 * (nu1 + nu2)) * cyl_bessel_k(0.5 * (nu1 + nu2), dist)/(pow(2, 0.5 * (nu1 + nu2) - 1) * tgamma(0.5 * (nu1 + nu2)));
		}
      		cor21(j, i) = cor12(i, j);
		}
	}
  	return rbind_cpp(cbind_cpp(cor11, cor12), cbind_cpp(cor21, cor22));
}


// [[Rcpp::export]]

NumericMatrix spatially_varying_parameters3_for_fitting(NumericMatrix & Loc, NumericVector & param, NumericMatrix & Nonstat_params) {

const int m = Loc.nrow();

float sigma2 = param(0), beta = param(1), time_param = param(2), nu = 1;
//float sigma2 = param(0), nu = param(1), time_param = param(2);

NumericMatrix cor(m, m);

for (int i = 0; i < m; ++i) {
for (int j = 0; j <= i; ++j) {

double x1 = Loc(i, 0) - Loc(j, 0);
double x2 = Loc(i, 1) - Loc(j, 1);
double tloc = Loc(i, 2) - Loc(j, 2);
NumericVector z2 = {x1, x2};

double omega1 = Nonstat_params(i, 0);
double log_lam1_1 = Nonstat_params(i, 1);
double log_lam1_2 = Nonstat_params(i, 2);

double omega2 = Nonstat_params(j, 0);
double log_lam2_1 = Nonstat_params(j, 1);
double log_lam2_2 = Nonstat_params(j, 2);

double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

double Inv_ij_11 = Kernel_ij_22; 
double Inv_ij_22 = Kernel_ij_11;
double Inv_ij_12 = - Kernel_ij_12; 
double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

double sigma = sqrt(sqrt(det_i * det_j)/det_ij);

double dist = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22) / det_ij / time_param;

if (dist == 0) {
cor(i, j) = sigma2 * pow(sigma, 2);
} else {
cor(i, j) = sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * gamma(nu));
}
cor(j, i) = cor(i, j);
}
}
return cor;
}

// [[Rcpp::export]]

NumericMatrix spatially_varying_parameters3(NumericMatrix & Loc, NumericVector & param) {

const int m = Loc.nrow();

float sigma2 = param(0), beta = param(1), time_param = param(2), nu = 1;

NumericMatrix cor(m, m);

for (int i = 0; i < m; ++i) {
for (int j = 0; j <= i; ++j) {

double new_loc1_x = Loc(i, 0);
double new_loc1_y = Loc(i, 1);
double new_loc2_x = Loc(j, 0);
double new_loc2_y = Loc(j, 1);

double x1 = Loc(i, 0) - Loc(j, 0);
double x2 = Loc(i, 1) - Loc(j, 1);
NumericVector z2 = {x1, x2};

double omega1 = 0 + 1 * (new_loc1_x - .5) + 1 * (new_loc1_y - .5) + 3 * pow(new_loc1_x - .5, 2) + -1 * pow(new_loc1_y - .5, 2);
double log_lam1_1 = exp(- pow(new_loc1_x, 2)) - exp(- pow(new_loc1_y, 2));
double log_lam1_2 = sin(new_loc1_x) * exp(- pow(new_loc1_y, 2)) - sin(new_loc1_y) * exp(- pow(new_loc1_x, 2));

double omega2 = 0 + 1 * (new_loc2_x - .5) + 1 * (new_loc2_y - .5) + 3 * pow(new_loc2_x - .5, 2) + -1 * pow(new_loc2_y - .5, 2);
double log_lam2_1 = exp(- pow(new_loc2_x, 2)) - exp(- pow(new_loc2_y, 2));
double log_lam2_2 = sin(new_loc2_x) * exp(- pow(new_loc2_y, 2)) - sin(new_loc2_y) * exp(- pow(new_loc2_x, 2));

double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

double Inv_ij_11 = Kernel_ij_22; 
double Inv_ij_22 = Kernel_ij_11;
double Inv_ij_12 = - Kernel_ij_12; 
double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

double sigma = sqrt(sqrt(det_i * det_j)/det_ij);

double dist = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22) / det_ij / time_param;

if (dist == 0) {
cor(i, j) = sigma2 * pow(sigma, 2);
} else {
cor(i, j) = sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * gamma(nu));
}
cor(j, i) = cor(i, j);
}
}
return cor;
}

// [[Rcpp::export]]

NumericMatrix spatially_varying_parameters3_montecarlo_for_fitting(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind1, NumericMatrix & wind2, NumericMatrix & Nonstat_params) {
  
  const int m = Loc.nrow(), n = wind1.ncol();
  
  float sigma2 = param(0), beta = param(1), nu = 1;
  
  NumericMatrix cor(m, m);
  
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
  	double temp_val = 0.0;
    	for (int k = 0; k < n; ++k) {
      
	      double new_loc1_x = Loc(i, 0) - wind1(i, k) * Loc(i, 2);
	      double new_loc1_y = Loc(i, 1) - wind2(i, k) * Loc(i, 2);
	      double new_loc2_x = Loc(j, 0) - wind1(j, k) * Loc(j, 2);
	      double new_loc2_y = Loc(j, 1) - wind2(j, k) * Loc(j, 2);

	      double x1 = new_loc1_x - new_loc2_x;
	      double x2 = new_loc1_y - new_loc2_y;
	      NumericVector z2 = {x1, x2};

		double omega1 = Nonstat_params(i, 0);
		double log_lam1_1 = Nonstat_params(i, 1);
		double log_lam1_2 = Nonstat_params(i, 2);

		double omega2 = Nonstat_params(j, 0);
		double log_lam2_1 = Nonstat_params(j, 1);
		double log_lam2_2 = Nonstat_params(j, 2);
			
	      double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
	      double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
	      double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);
	      
	      double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
	      double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
	      double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);
	      
	      double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
	      double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;
	      
	      double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
	      double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
	      double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);
	      
	      double Inv_ij_11 = Kernel_ij_22; 
	      double Inv_ij_22 = Kernel_ij_11;
	      double Inv_ij_12 = - Kernel_ij_12; 
	      double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;
	      
	      double sigma = sqrt(sqrt(det_i * det_j)/det_ij);
	      
	      double dist = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22);

	      if (dist == 0) {
		temp_val = temp_val + sigma2 * pow(sigma, 2);
	      } else {
		temp_val = temp_val + sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * gamma(nu));
	      }
	}
	cor(i, j) = temp_val;
      cor(j, i) = cor(i, j);
    }
  }
  return cor;
}


// [[Rcpp::export]]

NumericMatrix spatially_varying_parameters3_montecarlo(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind1, NumericMatrix & wind2) {
  
  const int m = Loc.nrow(), n = wind1.ncol();
  
  float sigma2 = param(0), beta = param(1), nu = 1;
  
  NumericMatrix cor(m, m);
  
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
  	double temp_val = 0.0;
    	for (int k = 0; k < n; ++k) {
      
	      double new_loc1_x = Loc(i, 0) - wind1(i, k) * Loc(i, 2);
	      double new_loc1_y = Loc(i, 1) - wind2(i, k) * Loc(i, 2);
	      double new_loc2_x = Loc(j, 0) - wind1(j, k) * Loc(j, 2);
	      double new_loc2_y = Loc(j, 1) - wind2(j, k) * Loc(j, 2);

	      double x1 = new_loc1_x - new_loc2_x;
	      double x2 = new_loc1_y - new_loc2_y;
	      NumericVector z2 = {x1, x2};
		
	      double omega1 = 0 + 1 * (new_loc1_x - .5) + 1 * (new_loc1_y - .5) + 3 * pow(new_loc1_x - .5, 2) + -1 * pow(new_loc1_y - .5, 2);
	      double log_lam1_1 = exp(- pow(new_loc1_x, 2)) - exp(- pow(new_loc1_y, 2));
	      double log_lam1_2 = sin(new_loc1_x) * exp(- pow(new_loc1_y, 2)) - sin(new_loc1_y) * exp(- pow(new_loc1_x, 2));

	      double omega2 = 0 + 1 * (new_loc2_x - .5) + 1 * (new_loc2_y - .5) + 3 * pow(new_loc2_x - .5, 2) + -1 * pow(new_loc2_y - .5, 2);
	      double log_lam2_1 = exp(- pow(new_loc2_x, 2)) - exp(- pow(new_loc2_y, 2));
	      double log_lam2_2 = sin(new_loc2_x) * exp(- pow(new_loc2_y, 2)) - sin(new_loc2_y) * exp(- pow(new_loc2_x, 2));

	      double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
	      double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
	      double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);
	      
	      double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
	      double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
	      double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);
	      
	      double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
	      double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;
	      
	      double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
	      double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
	      double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);
	      
	      double Inv_ij_11 = Kernel_ij_22; 
	      double Inv_ij_22 = Kernel_ij_11;
	      double Inv_ij_12 = - Kernel_ij_12; 
	      double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;
	      
	      double sigma = sqrt(sqrt(det_i * det_j)/det_ij);
	      
	      double dist = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22);

	      if (dist == 0) {
		temp_val = temp_val + sigma2 * pow(sigma, 2);
	      } else {
		temp_val = temp_val + sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * gamma(nu));
	      }
	}
	cor(i, j) = temp_val;
      cor(j, i) = cor(i, j);
    }
  }
  return cor;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]

NumericMatrix stationary_montecarlo(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind1, NumericMatrix & wind2) {
  
  const int m = Loc.nrow(), n = wind1.ncol();
  
  float sigma2 = param(0), beta = param(1), nu = 1;
  
  NumericMatrix cor(m, m);
  
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m ; ++j) {
  	double temp_val = 0.0;
    	for (int k = 0; k < n; ++k) {
      
	      double new_loc1_x = Loc(i, 0) - wind1(i, k) * Loc(i, 2);
	      double new_loc1_y = Loc(i, 1) - wind2(i, k) * Loc(i, 2);
	      double new_loc2_x = Loc(j, 0) - wind1(j, k) * Loc(j, 2);
	      double new_loc2_y = Loc(j, 1) - wind2(j, k) * Loc(j, 2);

	      double x1 = new_loc1_x - new_loc2_x;
	      double x2 = new_loc1_y - new_loc2_y;
	      double dist = sqrt(pow(x1, 2) + pow(x2, 2)) / beta;

	      if (dist == 0) {
		temp_val = temp_val + sigma2;
	      } else {
		temp_val = temp_val + sigma2 * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * gamma(nu));
	      }
	}
	cor(i, j) = temp_val;
      cor(j, i) = cor(i, j);
    }
  }
  return cor;
}


// [[Rcpp::export]]

NumericVector spatially_varying_parameters_TODELETE(NumericMatrix & Loc, NumericMatrix & Nonstat_params) {
  
  int nrow_x = Loc.nrow();
  NumericVector out(19322436);
  int count = 0;
  
  for (int i = 0; i < nrow_x; ++i) {
    for(int j = 0; j <= i; ++j){
      
      double x1 = Loc(i, 0) - Loc(j, 0);
      double x2 = Loc(i, 1) - Loc(j, 1);
      //double tloc = Loc(i, 2) - Loc(j, 2);
      //NumericVector z2 = {x1 - param(2) * tloc, x2 - param(3) * tloc};
      NumericVector z2 = {x1 , x2};
      
      double omega1 = Nonstat_params(i, 0);
      double log_lam1_1 = Nonstat_params(i, 1);
      double log_lam1_2 = Nonstat_params(i, 2);
      
      double omega2 = Nonstat_params(j, 0);
      double log_lam2_1 = Nonstat_params(j, 1);
      double log_lam2_2 = Nonstat_params(j, 2);
      
      double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
      double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
      double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);
      
      double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
      double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
      double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);
      
      double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
      double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;
      
      double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
      double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
      double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);
      
      double Inv_ij_11 = Kernel_ij_22; 
      double Inv_ij_22 = Kernel_ij_11;
      double Inv_ij_12 = - Kernel_ij_12; 
      double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;
      
      double sigma = sqrt(sqrt(det_i * det_j)/det_ij);
      
      double dist = sqrt((z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22) / det_ij);
      
      out(count) = exp( -dist);
      count = count + 1;
    }
  }
  
  return out;
}

// [[Rcpp::export]]

NumericMatrix spatially_varying_parametersBEST(NumericMatrix & Loc, NumericVector & param1, NumericVector & param2) {
  
  const int m = Loc.nrow();
  
  NumericMatrix cor(m, m);
  
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
      
      double x1 = Loc(i, 0) - Loc(j, 0);
      double x2 = Loc(i, 1) - Loc(j, 1);
      double tloc = Loc(i, 2) - Loc(j, 2);
      NumericVector z2 = {x1, x2};
      
      double omega1 = param1(i);
      double log_lam1_1 = param2(i);
      double log_lam1_2 = log_lam1_1;
      
      double omega2 = param1(j);
      double log_lam2_1 = param2(j);
      double log_lam2_2 = log_lam2_1;
      
      double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
      double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
      double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);
      
      double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
      double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
      double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);
      
      double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
      double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;
      
      double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
      double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
      double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);
      
      double Inv_ij_11 = Kernel_ij_22; 
      double Inv_ij_22 = Kernel_ij_11;
      double Inv_ij_12 = - Kernel_ij_12; 
      double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;
      
      double sigma = sqrt(sqrt(det_i * det_j)/det_ij);
      
      double dist = sqrt((z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22) / det_ij);
      //double dist = sqrt((z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22));
      //cor(i, j) = sigma + dist;
      //cor(i, j) = sigma * pow(dist / param(1), param(0)) * cyl_bessel_k(1, dist / param(1))/(pow(2, param(0) - 1));
      if (dist == 0) {
        cor(i, j) = sigma;
      } else {
        cor(i, j) = sigma * pow(dist / 0.23, 1) * cyl_bessel_k(1, dist / 0.23)/(pow(2, 1 - 1) * gamma(1));
      }
      cor(j, i) = cor(i, j);
    }
  }
  return cor;
}

// [[Rcpp::export]]

NumericMatrix spatially_varying_parameters_garg(NumericMatrix & Loc, NumericVector & param) {
  
  const int m = Loc.nrow();
  
  NumericMatrix cor(m, m);
  
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
      
      double x1 = Loc(i, 0) - Loc(j, 0);
      double x2 = Loc(i, 1) - Loc(j, 1);
      double tloc = Loc(i, 2) - Loc(j, 2);
      NumericVector z2 = {x1 - param(2) * tloc, x2 - param(3) * tloc};
      
      double omega1 = (Loc(i, 0) - .5) - 2*(Loc(i, 1) - .5) + pow(Loc(i, 1) - .5, 2);
      double log_lam1_1 = -1 - 6*pow(Loc(i, 0) - .5, 2) - 7*pow(Loc(i, 1) - .5, 2);
      double log_lam1_2 = log_lam1_1;
      
      double omega2 = (Loc(j, 0) - .5) - 2*(Loc(j, 1) - .5) + pow(Loc(j, 1) - .5, 2);
      double log_lam2_1 = -1 - 6*pow(Loc(j, 0) - .5, 2) - 7*pow(Loc(j, 1) - .5, 2);
      double log_lam2_2 = log_lam2_1;
      
      double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
      double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
      double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);
      
      double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
      double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
      double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);
      
      double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
      double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;
      
      double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
      double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
      double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);
      
      double Inv_ij_11 = Kernel_ij_22; 
      double Inv_ij_22 = Kernel_ij_11;
      double Inv_ij_12 = - Kernel_ij_12; 
      double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;
      
      double sigma = sqrt(sqrt(det_i * det_j)/det_ij);
      
      double dist = sqrt((z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22) / det_ij);
      if (dist == 0) {
        cor(i, j) = sigma;
      } else {
        cor(i, j) = sigma * pow(dist / param(1), param(0)) * cyl_bessel_k(1, dist / param(1))/(pow(2, param(0) - 1) * gamma(param(0)));
      }
      cor(j, i) = cor(i, j);
    }
  }
  return cor;
}

// [[Rcpp::export]]

NumericVector spatially_varying_parameters_garg_WLS(NumericVector & Loc0, NumericMatrix & Loc, NumericVector & param) {
  
  const int m = Loc.nrow();
  
  NumericVector cor(m);
  
  for (int j = 0; j < m/2; ++j) {
    
    double x1 = Loc0(0) - Loc(j, 0);
    double x2 = Loc0(1) - Loc(j, 1);
    double tloc = Loc0(2) - Loc(j, 2);
    NumericVector z2 = {x1 + param(2), x2 + param(3)};
    
    double omega1 = (Loc0(0) - .5) - 2*(Loc0(1) - .5) + pow(Loc0(1) - .5, 2);
    double log_lam1_1 = -1 - 6*pow(Loc0(0) - .5, 2) - 7*pow(Loc0(1) - .5, 2);
    double log_lam1_2 = log_lam1_1;
    
    double omega2 = (Loc(j, 0) - .5) - 2*(Loc(j, 1) - .5) + pow(Loc(j, 1) - .5, 2);
    double log_lam2_1 = -1 - 6*pow(Loc(j, 0) - .5, 2) - 7*pow(Loc(j, 1) - .5, 2);
    double log_lam2_2 = log_lam2_1;
    
    double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
    double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
    double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);
    
    double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
    double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
    double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);
    
    double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
    double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;
    
    double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
    double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
    double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);
    
    double Inv_ij_11 = Kernel_ij_22; 
    double Inv_ij_22 = Kernel_ij_11;
    double Inv_ij_12 = - Kernel_ij_12; 
    double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;
    
    double sigma = sqrt(sqrt(det_i * det_j)/det_ij);
    
    double dist = sqrt((z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22) / det_ij);
    if (dist == 0) {
      cor(j) = sigma;
    } else {
      cor(j) = sigma * pow(dist / param(1), param(0)) * cyl_bessel_k(1, dist / param(1))/(pow(2, param(0) - 1) * gamma(param(0)));
    }
  }
  
  for (int j = m/2; j < m; ++j) {
    
    double x1 = Loc0(0) - Loc(j, 0);
    double x2 = Loc0(1) - Loc(j, 1);
    double tloc = Loc0(2) - Loc(j, 2);
    NumericVector z2 = {x1 + (param(2) + param(4)), x2 + (param(3) + param(5))};
    
    double omega1 = (Loc0(0) - .5) - 2*(Loc0(1) - .5) + pow(Loc0(1) - .5, 2);
    double log_lam1_1 = -1 - 6*pow(Loc0(0) - .5, 2) - 7*pow(Loc0(1) - .5, 2);
    double log_lam1_2 = log_lam1_1;
    
    double omega2 = (Loc(j, 0) - .5) - 2*(Loc(j, 1) - .5) + pow(Loc(j, 1) - .5, 2);
    double log_lam2_1 = -1 - 6*pow(Loc(j, 0) - .5, 2) - 7*pow(Loc(j, 1) - .5, 2);
    double log_lam2_2 = log_lam2_1;
    
    double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
    double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
    double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);
    
    double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
    double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
    double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);
    
    double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
    double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;
    
    double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
    double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
    double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);
    
    double Inv_ij_11 = Kernel_ij_22; 
    double Inv_ij_22 = Kernel_ij_11;
    double Inv_ij_12 = - Kernel_ij_12; 
    double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;
    
    double sigma = sqrt(sqrt(det_i * det_j)/det_ij);
    
    double dist = sqrt((z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22) / det_ij);
    if (dist == 0) {
      cor(j) = sigma;
    } else {
      cor(j) = sigma * pow(dist / param(1), param(0)) * cyl_bessel_k(1, dist / param(1))/(pow(2, param(0) - 1) * gamma(param(0)));
    }
  }
  return cor;
}

// [[Rcpp::export]]

NumericVector spatially_varying_parameters2_WLS(NumericVector & Loc0, NumericMatrix & Loc, NumericVector & param) {
  
  const int m = Loc.nrow();
  
  NumericVector cor(m);
  
  for (int j = 0; j < m/2; ++j) {
    
    double x1 = Loc0(0) - Loc(j, 0);
    double x2 = Loc0(1) - Loc(j, 1);
    double tloc = Loc0(2) - Loc(j, 2);
    NumericVector z2 = {x1 - param(2) * tloc, x2 - param(3) * tloc};
    
    double omega1 = (Loc0(0) - param(2) * Loc0(2) - .5) - 2*(Loc0(1) - param(3) * Loc0(2) - .5) + pow(Loc0(1) - param(3) * Loc0(2) - .5, 2);
    double log_lam1_1 = -1 - 6*pow(Loc0(0) - param(2) * Loc0(2) - .5, 2) - 7*pow(Loc0(1) - param(3) * Loc0(2) - .5, 2);
    double log_lam1_2 = log_lam1_1;
    
    double omega2 = (Loc(j, 0) - param(2) * Loc(j, 2) - .5) - 2*(Loc(j, 1) - param(3) * Loc(j, 2) - .5) + pow(Loc(j, 1) - param(3) * Loc(j, 2) - .5, 2);
    double log_lam2_1 = -1 - 6*pow(Loc(j, 0) - param(2) * Loc(j, 2) - .5, 2) - 7*pow(Loc(j, 1) - param(3) * Loc(j, 2) - .5, 2);
    //double log_lam2_2 = -1 + 6*pow(Loc2(0) - x(0) * Loc2(2) - .5, 2) - 4*pow(Loc2(1) - x(1) * Loc2(2) - .5, 2);
    //double log_lam2_2 = log(0.37 - exp(log_lam2_1));
    double log_lam2_2 = log_lam2_1;
    
    double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
    double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
    double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);
    
    double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
    double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
    double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);
    
    double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
    double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;
    
    double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
    double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
    double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);
    
    double Inv_ij_11 = Kernel_ij_22; 
    double Inv_ij_22 = Kernel_ij_11;
    double Inv_ij_12 = - Kernel_ij_12; 
    double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;
    
    double sigma = sqrt(sqrt(det_i * det_j)/det_ij);
    
    double dist = sqrt((z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22) / det_ij);
    //double dist = sqrt((z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22));
    //cor(i, j) = sigma + dist;
    //cor(i, j) = sigma * pow(dist / param(1), param(0)) * cyl_bessel_k(1, dist / param(1))/(pow(2, param(0) - 1));
    if (dist == 0) {
      cor(j) = sigma;
    } else {
      cor(j) = sigma * pow(dist / param(1), param(0)) * cyl_bessel_k(1, dist / param(1))/(pow(2, param(0) - 1) * gamma(param(0)));
    }
  }
  for (int j = m/2; j < m; ++j) {
    
    double x1 = Loc0(0) - Loc(j, 0);
    double x2 = Loc0(1) - Loc(j, 1);
    double tloc = Loc0(2) - Loc(j, 2);
    NumericVector z2 = {x1 - param(2) * tloc, x2 - param(3) * tloc};
    
    double omega1 = (Loc0(0) - param(2) * Loc0(2) - .5) - 2*(Loc0(1) - param(3) * Loc0(2) - .5) + pow(Loc0(1) - param(3) * Loc0(2) - .5, 2);
    double log_lam1_1 = -1 - 6*pow(Loc0(0) - param(2) * Loc0(2) - .5, 2) - 7*pow(Loc0(1) - param(3) * Loc0(2) - .5, 2);
    double log_lam1_2 = log_lam1_1;
    
    double omega2 = (Loc(j, 0) - param(2) * Loc(j, 2) - .5) - 2*(Loc(j, 1) - param(3) * Loc(j, 2) - .5) + pow(Loc(j, 1) - param(3) * Loc(j, 2) - .5, 2);
    double log_lam2_1 = -1 - 6*pow(Loc(j, 0) - param(2) * Loc(j, 2) - .5, 2) - 7*pow(Loc(j, 1) - param(3) * Loc(j, 2) - .5, 2);
    //double log_lam2_2 = -1 + 6*pow(Loc2(0) - x(0) * Loc2(2) - .5, 2) - 4*pow(Loc2(1) - x(1) * Loc2(2) - .5, 2);
    //double log_lam2_2 = log(0.37 - exp(log_lam2_1));
    double log_lam2_2 = log_lam2_1;
    
    double Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
    double Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
    double Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);
    
    double Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
    double Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
    double Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);
    
    double det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
    double det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;
    
    double Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
    double Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
    double Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);
    
    double Inv_ij_11 = Kernel_ij_22; 
    double Inv_ij_22 = Kernel_ij_11;
    double Inv_ij_12 = - Kernel_ij_12; 
    double det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;
    
    double sigma = sqrt(sqrt(det_i * det_j)/det_ij);
    
    double dist = sqrt((z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22) / det_ij);
    //double dist = sqrt((z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22));
    //cor(i, j) = sigma + dist;
    //cor(i, j) = sigma * pow(dist / param(1), param(0)) * cyl_bessel_k(1, dist / param(1))/(pow(2, param(0) - 1));
    if (dist == 0) {
      cor(j) = sigma;
    } else {
      cor(j) = sigma * pow(dist / param(1), param(0)) * cyl_bessel_k(1, dist / param(1))/(pow(2, param(0) - 1) * gamma(param(0)));
    }
  }
  return cor;
}
