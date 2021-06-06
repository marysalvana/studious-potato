// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <RcppEigen.h>
//#include <gsl/gsl_sf_bessel.h>

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

List DEFORMATION(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind, NumericVector & param_nonstat) {
  
  	const int m = Loc.nrow(), n_wind = wind.nrow();
  	float sigma2 = param(0), beta = param(1), nu = param(2);
  	NumericMatrix cor(m, m), param_matrix(m, 5);
	List conso(2);
  
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

				if(i == j){
					param_matrix(i, 0) = Loc(i, 0);
					param_matrix(i, 1) = Loc(i, 1);
					param_matrix(i, 2) = Loc(i, 2);
					param_matrix(i, 3) = deform_loc1_x;
					param_matrix(i, 4) = deform_loc1_y;
				}

				double dist = sqrt(pow(deform_loc2_x - deform_loc1_x, 2) + pow(deform_loc2_y - deform_loc1_y, 2));

	      			if (dist == 0) {
					temp_val = temp_val + sigma2;
	      			}else {
					temp_val = temp_val + sigma2 * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
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

List DEFORMATION_FOR_FITTING(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind, NumericVector & param_nonstat) {
  
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

				double deform_loc1_x = new_loc1_x + param_nonstat(i, 0);
				double deform_loc1_y = new_loc1_y + param_nonstat(i, 1);
				double deform_loc2_x = new_loc2_x + param_nonstat(j, 0);
				double deform_loc2_y = new_loc2_y + param_nonstat(j, 1);

				double dist = sqrt(pow(deform_loc2_x - deform_loc1_x, 2) + pow(deform_loc2_y - deform_loc1_y, 2));

	      			if (dist == 0) {
					temp_val = temp_val + sigma2;
	      			}else {
					temp_val = temp_val + sigma2 * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
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

NumericMatrix DEFORMATION_PARALLEL(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind, NumericVector & param_nonstat) {
  
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
	      			}else {
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

NumericMatrix DEFORMATION_FOR_FITTING_PARALLEL(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind, NumericVector & param_nonstat) {
  
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

				double deform_loc1_x = new_loc1_x + param_nonstat(i, 0);
				double deform_loc1_y = new_loc1_y + param_nonstat(i, 1);
				double deform_loc2_x = new_loc2_x + param_nonstat(j, 0);
				double deform_loc2_y = new_loc2_y + param_nonstat(j, 1);

				double dist = sqrt(pow(deform_loc2_x - deform_loc1_x, 2) + pow(deform_loc2_y - deform_loc1_y, 2));

	      			if (dist == 0) {
					temp_val = temp_val + sigma2;
	      			}else {
					temp_val = temp_val + sigma2 * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
	      			}
			}
		cor(i, j) = temp_val / n_wind;
      		cor(j, i) = cor(i, j);
    		}
  	}
  	return cor;
}
///////////////////////////////////////////////   MULTIVARIATE   ///////////////////////////////////////////////

// [[Rcpp::export]]

List MULTIVARIATE_SPATIALLY_VARYING_PARAMETERS_FOR_FITTING(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind, NumericMatrix & param_nonstat, NumericMatrix & param_nonstat2, int & time) {

  	const int m = Loc.nrow(), n_wind = wind.nrow();
  	float sigma2_1 = param(0), sigma2_2 = param(1), beta = param(2), nu1 = param(3), nu2 = param(4), rho = param(5);
	float nu12 = 0.5 * (nu1 + nu2);
	NumericMatrix cor11(m, m), cor22(m, m), cor12(m, m), cor21(m, m);
	List conso(1);

	for (int i = 0; i < m; ++i) {
    		for (int j = 0; j <m; ++j) {
    		//for (int j = 0; j <= i; ++j) {

        		double temp_val11 = 0.0, temp_val22 = 0.0, temp_val12 = 0.0;

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

				double sigma11 = sqrt(sqrt(det_i * det_j)/det_ij);

				double dist11 = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22);

			      	omega1 = param_nonstat2(i, 0);
			      	log_lam1_1 = param_nonstat2(i, 1);
			      	log_lam1_2 = param_nonstat2(i, 2);

			      	omega2 = param_nonstat2(j, 0);
			      	log_lam2_1 = param_nonstat2(j, 1);
			      	log_lam2_2 = param_nonstat2(j, 2);

				Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
				Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
				Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

				Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
				Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
				Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

				det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
				det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

				Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
				Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
				Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

				Inv_ij_11 = Kernel_ij_22;
				Inv_ij_22 = Kernel_ij_11;
				Inv_ij_12 = - Kernel_ij_12;
				det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

				double sigma22 = sqrt(sqrt(det_i * det_j)/det_ij);

				double dist22 = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22);

			      	omega1 = param_nonstat(i, 0);
			      	log_lam1_1 = param_nonstat(i, 1);
			      	log_lam1_2 = param_nonstat(i, 2);

			      	omega2 = param_nonstat2(j, 0);
			      	log_lam2_1 = param_nonstat2(j, 1);
			      	log_lam2_2 = param_nonstat2(j, 2);

				Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
				Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
				Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

				Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
				Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
				Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

				det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
				det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

				Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
				Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
				Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

				Inv_ij_11 = Kernel_ij_22;
				Inv_ij_22 = Kernel_ij_11;
				Inv_ij_12 = - Kernel_ij_12;
				det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

				double sigma12 = sqrt(sqrt(det_i * det_j)/det_ij);

				double dist12 = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22);

				if (dist11 == 0) {
					temp_val11 = temp_val11 + sigma2_1 * pow(sigma11, 2);
					temp_val22 = temp_val22 + sigma2_2 * pow(sigma22, 2);
					temp_val12 = temp_val12 + rho * sqrt(sigma2_1 * sigma2_2) * pow(sigma12, 2);
				}else {
					temp_val11 = temp_val11 + sigma2_1 * pow(sigma11, 2) * pow(dist11, nu1) * cyl_bessel_k(nu1, dist11) / (pow(2, nu1 - 1) * tgamma(nu1));
					temp_val22 = temp_val22 + sigma2_2 * pow(sigma22, 2) * pow(dist22, nu2) * cyl_bessel_k(nu2, dist22) / (pow(2, nu2 - 1) * tgamma(nu2));
					temp_val12 = temp_val12 + rho * sqrt(sigma2_1 * sigma2_2) * pow(sigma12, 2) * pow(dist12, nu12) * cyl_bessel_k(nu12, dist12) / (pow(2, nu12 - 1) * tgamma(nu12));
				}

        		}
        		cor11(i, j) = temp_val11 / n_wind;
        		//cor11(j, i) = cor11(i, j);

        		cor22(i, j) = temp_val22 / n_wind;
        		//cor22(j, i) = cor22(i, j);

        		cor12(i, j) = temp_val12 / n_wind;
      			cor21(j, i) = cor12(i, j);

			//cor12(j, i) = cor12(i, j);
      			//cor21(i, j) = cor12(j, i);
		}
	}
	//return rbind_cpp(cbind_cpp(cor11, cor12), cbind_cpp(cor21, cor22));
	conso(0) = rbind_cpp(cbind_cpp(cor11, cor12), cbind_cpp(cor21, cor22));
  	return conso;
}

// [[Rcpp::export]]

NumericMatrix MULTIVARIATE_SPATIALLY_VARYING_PARAMETERS_FOR_FITTING_PARALLEL(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind, NumericMatrix & param_nonstat, NumericMatrix & param_nonstat2, int & time) {

  	const int m = Loc.nrow(), n_wind = wind.nrow();
  	float sigma2_1 = param(0), sigma2_2 = param(1), beta = param(2), nu1 = param(3), nu2 = param(4), rho = param(5);
	float nu12 = 0.5 * (nu1 + nu2);
	NumericMatrix cor11(m, m), cor22(m, m), cor12(m, m), cor21(m, m);

	for (int i = 0; i < m; ++i) {
    		for (int j = 0; j <m; ++j) {
    		//for (int j = 0; j <= i; ++j) {

        		double temp_val11 = 0.0, temp_val22 = 0.0, temp_val12 = 0.0;

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

				double sigma11 = sqrt(sqrt(det_i * det_j)/det_ij);

				double dist11 = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22);

			      	omega1 = param_nonstat2(i, 0);
			      	log_lam1_1 = param_nonstat2(i, 1);
			      	log_lam1_2 = param_nonstat2(i, 2);

			      	omega2 = param_nonstat2(j, 0);
			      	log_lam2_1 = param_nonstat2(j, 1);
			      	log_lam2_2 = param_nonstat2(j, 2);

				Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
				Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
				Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

				Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
				Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
				Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

				det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
				det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

				Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
				Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
				Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

				Inv_ij_11 = Kernel_ij_22;
				Inv_ij_22 = Kernel_ij_11;
				Inv_ij_12 = - Kernel_ij_12;
				det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

				double sigma22 = sqrt(sqrt(det_i * det_j)/det_ij);

				double dist22 = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22);

			      	omega1 = param_nonstat(i, 0);
			      	log_lam1_1 = param_nonstat(i, 1);
			      	log_lam1_2 = param_nonstat(i, 2);

			      	omega2 = param_nonstat2(j, 0);
			      	log_lam2_1 = param_nonstat2(j, 1);
			      	log_lam2_2 = param_nonstat2(j, 2);

				Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
				Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
				Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

				Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
				Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
				Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

				det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
				det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

				Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
				Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
				Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

				Inv_ij_11 = Kernel_ij_22;
				Inv_ij_22 = Kernel_ij_11;
				Inv_ij_12 = - Kernel_ij_12;
				det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

				double sigma12 = sqrt(sqrt(det_i * det_j)/det_ij);

				double dist12 = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22);

				if (dist11 == 0) {
					temp_val11 = temp_val11 + sigma2_1 * pow(sigma11, 2);
					temp_val22 = temp_val22 + sigma2_2 * pow(sigma22, 2);
					temp_val12 = temp_val12 + rho * sqrt(sigma2_1 * sigma2_2) * pow(sigma12, 2);
				}else {
					temp_val11 = temp_val11 + sigma2_1 * pow(sigma11, 2) * pow(dist11, nu1) * cyl_bessel_k(nu1, dist11) / (pow(2, nu1 - 1) * tgamma(nu1));
					temp_val22 = temp_val22 + sigma2_2 * pow(sigma22, 2) * pow(dist22, nu2) * cyl_bessel_k(nu2, dist22) / (pow(2, nu2 - 1) * tgamma(nu2));
					temp_val12 = temp_val12 + rho * sqrt(sigma2_1 * sigma2_2) * pow(sigma12, 2) * pow(dist12, nu12) * cyl_bessel_k(nu12, dist12) / (pow(2, nu12 - 1) * tgamma(nu12));
				}

        		}
        		cor11(i, j) = temp_val11 / n_wind;
        		//cor11(j, i) = cor11(i, j);

        		cor22(i, j) = temp_val22 / n_wind;
        		//cor22(j, i) = cor22(i, j);

        		cor12(i, j) = temp_val12 / n_wind;
      			cor21(j, i) = cor12(i, j);

			//cor12(j, i) = cor12(i, j);
      			//cor21(i, j) = cor12(j, i);
		}
	}
	return rbind_cpp(cbind_cpp(cor11, cor12), cbind_cpp(cor21, cor22));
}

// [[Rcpp::export]]

NumericMatrix MULTIPLE_ADVEC_MULTIVARIATE_SPATIALLY_VARYING_PARAMETERS_FOR_FITTING_PARALLEL(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind, NumericMatrix & wind2, NumericMatrix & param_nonstat, NumericMatrix & param_nonstat2, int & time) {

  	const int m = Loc.nrow(), n_wind = wind.nrow();
  	float sigma2_1 = param(0), sigma2_2 = param(1), beta = param(2), nu1 = param(3), nu2 = param(4), rho = param(5);
	float nu12 = 0.5 * (nu1 + nu2);
	NumericMatrix cor11(m, m), cor22(m, m), cor12(m, m), cor21(m, m);

	for (int i = 0; i < m; ++i) {
    		for (int j = 0; j <m; ++j) {
    		//for (int j = 0; j <= i; ++j) {

        		double temp_val11 = 0.0, temp_val22 = 0.0, temp_val12 = 0.0;

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

				double sigma11 = sqrt(sqrt(det_i * det_j)/det_ij);

				double dist11 = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22);
			      	new_loc1_x = Loc(i, 0) - wind2(k, 0) * Loc(i, 2);
			      	new_loc1_y = Loc(i, 1) - wind2(k, 1) * Loc(i, 2);
			      	new_loc2_x = Loc(j, 0) - wind2(k, 0) * Loc(j, 2);
			      	new_loc2_y = Loc(j, 1) - wind2(k, 1) * Loc(j, 2);

			      	x1 = new_loc1_x - new_loc2_x;
			      	x2 = new_loc1_y - new_loc2_y;
			      	z2 = {x1, x2};

			      	omega1 = param_nonstat2(i, 0);
			      	log_lam1_1 = param_nonstat2(i, 1);
			      	log_lam1_2 = param_nonstat2(i, 2);

			      	omega2 = param_nonstat2(j, 0);
			      	log_lam2_1 = param_nonstat2(j, 1);
			      	log_lam2_2 = param_nonstat2(j, 2);

				Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
				Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
				Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

				Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
				Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
				Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

				det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
				det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

				Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
				Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
				Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

				Inv_ij_11 = Kernel_ij_22;
				Inv_ij_22 = Kernel_ij_11;
				Inv_ij_12 = - Kernel_ij_12;
				det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

				double sigma22 = sqrt(sqrt(det_i * det_j)/det_ij);

				double dist22 = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22);

			      	new_loc1_x = Loc(i, 0) - wind(k, 0) * Loc(i, 2);
			      	new_loc1_y = Loc(i, 1) - wind(k, 1) * Loc(i, 2);
			      	new_loc2_x = Loc(j, 0) - wind2(k, 0) * Loc(j, 2);
			      	new_loc2_y = Loc(j, 1) - wind2(k, 1) * Loc(j, 2);

			      	x1 = new_loc1_x - new_loc2_x;
			      	x2 = new_loc1_y - new_loc2_y;
			      	z2 = {x1, x2};

			      	omega1 = param_nonstat(i, 0);
			      	log_lam1_1 = param_nonstat(i, 1);
			      	log_lam1_2 = param_nonstat(i, 2);

			      	omega2 = param_nonstat2(j, 0);
			      	log_lam2_1 = param_nonstat2(j, 1);
			      	log_lam2_2 = param_nonstat2(j, 2);

				Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
				Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
				Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);

				Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
				Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
				Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

				det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
				det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;

				Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
				Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
				Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);

				Inv_ij_11 = Kernel_ij_22;
				Inv_ij_22 = Kernel_ij_11;
				Inv_ij_12 = - Kernel_ij_12;
				det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;

				double sigma12 = sqrt(sqrt(det_i * det_j)/det_ij);

				double dist12 = sqrt(z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22);

				if (dist11 == 0) {
					temp_val11 = temp_val11 + sigma2_1 * pow(sigma11, 2);
					temp_val22 = temp_val22 + sigma2_2 * pow(sigma22, 2);
					temp_val12 = temp_val12 + rho * sqrt(sigma2_1 * sigma2_2) * pow(sigma12, 2);
				}else {
					temp_val11 = temp_val11 + sigma2_1 * pow(sigma11, 2) * pow(dist11, nu1) * cyl_bessel_k(nu1, dist11) / (pow(2, nu1 - 1) * tgamma(nu1));
					temp_val22 = temp_val22 + sigma2_2 * pow(sigma22, 2) * pow(dist22, nu2) * cyl_bessel_k(nu2, dist22) / (pow(2, nu2 - 1) * tgamma(nu2));
					temp_val12 = temp_val12 + rho * sqrt(sigma2_1 * sigma2_2) * pow(sigma12, 2) * pow(dist12, nu12) * cyl_bessel_k(nu12, dist12) / (pow(2, nu12 - 1) * tgamma(nu12));
				}

        		}
        		cor11(i, j) = temp_val11 / n_wind;
        		//cor11(j, i) = cor11(i, j);

        		cor22(i, j) = temp_val22 / n_wind;
        		//cor22(j, i) = cor22(i, j);

        		cor12(i, j) = temp_val12 / n_wind;
      			cor21(j, i) = cor12(i, j);

			//cor12(j, i) = cor12(i, j);
      			//cor21(i, j) = cor12(j, i);
		}
	}
	return rbind_cpp(cbind_cpp(cor11, cor12), cbind_cpp(cor21, cor22));
}

// [[Rcpp::export]]

NumericMatrix MULTIVARIATE_DEFORMATION_FOR_FITTING_PARALLEL(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind, NumericVector & param_nonstat, NumericVector & param_nonstat2) {
  
  	const int m = Loc.nrow(), n_wind = wind.nrow();
  	float sigma2_1 = param(0), sigma2_2 = param(1), beta = param(2), nu1 = param(3), nu2 = param(4), rho = param(5);
	float nu12 = 0.5 * (nu1 + nu2);
	NumericMatrix cor11(m, m), cor22(m, m), cor12(m, m), cor21(m, m);
  
  	for (int i = 0; i < m; ++i) {
    		for (int j = 0; j < m; ++j) {
 
        		double temp_val11 = 0.0, temp_val22 = 0.0, temp_val12 = 0.0;

    			for (int k = 0; k < n_wind; ++k) {
      
	      			double new_loc1_x = Loc(i, 0) - wind(k, 0) * Loc(i, 2);
	      			double new_loc1_y = Loc(i, 1) - wind(k, 1) * Loc(i, 2);
	      			double new_loc2_x = Loc(j, 0) - wind(k, 0) * Loc(j, 2);
	      			double new_loc2_y = Loc(j, 1) - wind(k, 1) * Loc(j, 2);

				double deform_loc1_x = new_loc1_x + param_nonstat(i, 0);
				double deform_loc1_y = new_loc1_y + param_nonstat(i, 1);
				double deform_loc2_x = new_loc2_x + param_nonstat(j, 0);
				double deform_loc2_y = new_loc2_y + param_nonstat(j, 1);

				double dist11 = sqrt(pow(deform_loc2_x - deform_loc1_x, 2) + pow(deform_loc2_y - deform_loc1_y, 2));

				deform_loc1_x = new_loc1_x + param_nonstat2(i, 0);
				deform_loc1_y = new_loc1_y + param_nonstat2(i, 1);
				deform_loc2_x = new_loc2_x + param_nonstat2(j, 0);
				deform_loc2_y = new_loc2_y + param_nonstat2(j, 1);

				double dist22 = sqrt(pow(deform_loc2_x - deform_loc1_x, 2) + pow(deform_loc2_y - deform_loc1_y, 2));

				deform_loc1_x = new_loc1_x + param_nonstat(i, 0);
				deform_loc1_y = new_loc1_y + param_nonstat(i, 1);
				deform_loc2_x = new_loc2_x + param_nonstat2(j, 0);
				deform_loc2_y = new_loc2_y + param_nonstat2(j, 1);

				double dist12 = sqrt(pow(deform_loc2_x - deform_loc1_x, 2) + pow(deform_loc2_y - deform_loc1_y, 2));

				if (dist11 == 0) {
					temp_val11 = temp_val11 + sigma2_1;
					temp_val22 = temp_val22 + sigma2_2;
					temp_val12 = temp_val12 + rho * sqrt(sigma2_1 * sigma2_2);
				}else {
					temp_val11 = temp_val11 + sigma2_1 * pow(dist11, nu1) * cyl_bessel_k(nu1, dist11) / (pow(2, nu1 - 1) * tgamma(nu1));
					temp_val22 = temp_val22 + sigma2_2 * pow(dist22, nu2) * cyl_bessel_k(nu2, dist22) / (pow(2, nu2 - 1) * tgamma(nu2));
					temp_val12 = temp_val12 + rho * sqrt(sigma2_1 * sigma2_2) * pow(dist12, nu12) * cyl_bessel_k(nu12, dist12) / (pow(2, nu12 - 1) * tgamma(nu12));
				}

			}
        		cor11(i, j) = temp_val11 / n_wind;
        		//cor11(j, i) = cor11(i, j);

        		cor22(i, j) = temp_val22 / n_wind;
        		//cor22(j, i) = cor22(i, j);

        		cor12(i, j) = temp_val12 / n_wind;
      			cor21(j, i) = cor12(i, j);

			//cor12(j, i) = cor12(i, j);
      			//cor21(i, j) = cor12(j, i);
    		}
  	}
	return rbind_cpp(cbind_cpp(cor11, cor12), cbind_cpp(cor21, cor22));
}
		
NumericMatrix ORIG_MULTIVARIATE_SPATIALLY_VARYING_PARAMETERS_FOR_FITTING_PARALLEL(NumericMatrix & Loc, NumericVector & param, NumericMatrix & Nonstat_params, NumericMatrix & Nonstat_params2) {
  
  const int m = Loc.nrow();
  
  float nu1 = param(0), nu2 = param(1);
  
  NumericMatrix cor11(m, m), cor22(m, m), cor12(m, m), cor21(m, m);
  
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      
      double x1 = Loc(i, 0) - Loc(j, 0);
      double x2 = Loc(i, 1) - Loc(j, 1);
      double tloc = Loc(i, 2) - Loc(j, 2);
      NumericVector z2 = {x1 - param(4) * tloc, x2 - param(5) * tloc};
      
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
      
      double sigma11 = sqrt(sqrt(det_i * det_j)/det_ij);
      
      double dist11 = sqrt((z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22) / det_ij);
      
      omega1 = Nonstat_params2(i, 0);
      log_lam1_1 = Nonstat_params2(i, 1);
      log_lam1_2 = Nonstat_params2(i, 2);
      
      omega2 = Nonstat_params2(j, 0);
      log_lam2_1 = Nonstat_params2(j, 1);
      log_lam2_2 = Nonstat_params2(j, 2);
      
      Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
      Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
      Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);
      
      Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
      Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
      Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);
      
      det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
      det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;
      
      Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
      Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
      Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);
      
      Inv_ij_11 = Kernel_ij_22; 
      Inv_ij_22 = Kernel_ij_11;
      Inv_ij_12 = - Kernel_ij_12; 
      det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;
      
      double sigma22 = sqrt(sqrt(det_i * det_j)/det_ij);
      
      double dist22 = sqrt((z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22) / det_ij);
      
      omega1 = Nonstat_params(i, 0);
      log_lam1_1 = Nonstat_params(i, 1);
      log_lam1_2 = Nonstat_params(i, 2);
      
      omega2 = Nonstat_params2(j, 0);
      log_lam2_1 = Nonstat_params2(j, 1);
      log_lam2_2 = Nonstat_params2(j, 2);
      
      Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
      Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
      Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);
      
      Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
      Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
      Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);
      
      det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
      det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;
      
      Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
      Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
      Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);
      
      Inv_ij_11 = Kernel_ij_22; 
      Inv_ij_22 = Kernel_ij_11;
      Inv_ij_12 = - Kernel_ij_12; 
      det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;
      
      double sigma12 = sqrt(sqrt(det_i * det_j)/det_ij);
      
      double dist12 = sqrt((z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22) / det_ij);
      
      if (dist11 == 0) {
        cor11(i, j) = pow(sigma11, 2);
        cor22(i, j) = pow(sigma22, 2);
        cor12(i, j) = param(3) * sigma11 * sigma22;
      } else {
        cor11(i, j) = pow(sigma11, 2) * pow(dist11 / param(2), nu1) * cyl_bessel_k(nu1, dist11 / param(2))/(pow(2, nu1 - 1) * tgamma(nu1));
        cor22(i, j) = pow(sigma22, 2) * pow(dist22 / param(2), nu2) * cyl_bessel_k(nu2, dist22 / param(2))/(pow(2, nu2 - 1) * tgamma(nu2));
        cor12(i, j) = param(3) * sigma11 * sigma22 * pow(dist12 / param(2), 0.5 * (nu1 + nu2)) * cyl_bessel_k(0.5 * (nu1 + nu2), dist12 / param(2))/(pow(2, 0.5 * (nu1 + nu2) - 1) * tgamma(0.5 * (nu1 + nu2)));
      }
      //cor11(j, i) = cor11(i, j);
      //cor22(j, i) = cor22(i, j);
      //cor12(j, i) = cor12(i, j);
      cor21(j, i) = cor12(i, j);
      //cor21(i, j) = cor12(j, i);
    }
  }
  return rbind_cpp(cbind_cpp(cor11, cor12), cbind_cpp(cor21, cor22));
}

///////////////////////////////////////////////   SCRATCH   ///////////////////////////////////////////////  

// [[Rcpp::export]]

NumericMatrix POINT_SOURCE_DEFORMATION(NumericMatrix & Loc1, NumericMatrix & Loc2, NumericVector & param, NumericVector & param_nonstat) {
  
  const int m = Loc1.nrow(), n = Loc2.nrow();
  
  float sigma2 = param(0), beta = param(1), nu = param(2);
  
  NumericMatrix cor(m, n);
  
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
 
 	double temp_val = 0.0;

	      double new_loc1_x = Loc1(i, 0);
	      double new_loc1_y = Loc1(i, 1);
	      double new_loc2_x = Loc2(j, 0);
	      double new_loc2_y = Loc2(j, 1);

		double dist_source1 = sqrt(pow(new_loc1_x - param_nonstat(0), 2) + pow(new_loc1_y - param_nonstat(1), 2)); 
		double dist_source2 = sqrt(pow(new_loc2_x - param_nonstat(0), 2) + pow(new_loc2_y - param_nonstat(1), 2)); 
		double deform_loc1_x = param_nonstat(0) + (new_loc1_x - param_nonstat(0)) * dist_source1;
		double deform_loc1_y = param_nonstat(1) + (new_loc1_y - param_nonstat(1)) * dist_source1;
		double deform_loc2_x = param_nonstat(0) + (new_loc2_x - param_nonstat(0)) * dist_source2;
		double deform_loc2_y = param_nonstat(1) + (new_loc2_y - param_nonstat(1)) * dist_source2;

		double dist = sqrt(pow(deform_loc2_x - deform_loc1_x, 2) + pow(deform_loc2_y - deform_loc1_y, 2));

	      if (dist == 0) {
		temp_val = temp_val + sigma2;
	      } else {
		temp_val = temp_val + sigma2 * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
	      }
	cor(i, j) = temp_val;
    }
  }
  return cor;
}

// [[Rcpp::export]]

NumericMatrix spatially_varying_parameters2(NumericMatrix & Loc, NumericVector & param, NumericMatrix & Nonstat_params) {

const int m = Loc.nrow();

double sigma2 = param(0), beta = param(1), nu = param(2);

NumericMatrix cor(m, m);

for (int i = 0; i < m; ++i) {
for (int j = 0; j <= i; ++j) {

double x1 = Loc(i, 0) - Loc(j, 0);
double x2 = Loc(i, 1) - Loc(j, 1);
double tloc = Loc(i, 2) - Loc(j, 2);
//NumericVector z2 = {x1 - param(2) * tloc, x2 - param(3) * tloc};
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

double dist = sqrt((z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22) / det_ij);

//Rprintf("dist : %f, sigma : %f, sigma2 : %f, nu : %f \n", dist, sigma, sigma2, nu);
//Rprintf("the value of v[%i] : %f \n", i, v[i]);

if (dist == 0) {
cor(i, j) = sigma2 * pow(sigma, 2);
} else {
//cor(i, j) = sigma2 * pow(sigma, 2) * pow(dist, nu) * gsl_sf_bessel_Knu(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
//cor(i, j) = sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
cor(i, j) = sigma2 * pow(sigma, 2) * exp(-pow(dist, 2) / beta);
}
cor(j, i) = cor(i, j);
}
}
return cor;
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
cor(i, j) = sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
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
cor(i, j) = sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
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
		temp_val = temp_val + sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
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
		temp_val = temp_val + sigma2 * pow(sigma, 2) * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
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
		temp_val = temp_val + sigma2 * pow(dist, nu) * cyl_bessel_k(nu, dist) / (pow(2, nu - 1) * tgamma(nu));
	      }
	}
	cor(i, j) = temp_val;
      cor(j, i) = cor(i, j);
    }
  }
  return cor;
}


// [[Rcpp::export]]

NumericVector spatially_varying_parameters(NumericMatrix & Loc, NumericMatrix & Nonstat_params) {
  
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

NumericMatrix multivariate_spatially_varying_orig(NumericMatrix & Loc, NumericVector & param, NumericMatrix & Nonstat_params, NumericMatrix & Nonstat_params2) {
  
  const int m = Loc.nrow();
  
  float nu1 = param(0), nu2 = param(1);
  
  NumericMatrix cor11(m, m), cor22(m, m), cor12(m, m), cor21(m, m);
  
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      
      double x1 = Loc(i, 0) - Loc(j, 0);
      double x2 = Loc(i, 1) - Loc(j, 1);
      double tloc = Loc(i, 2) - Loc(j, 2);
      NumericVector z2 = {x1 - param(4) * tloc, x2 - param(5) * tloc};
      
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
      
      double sigma11 = sqrt(sqrt(det_i * det_j)/det_ij);
      
      double dist11 = sqrt((z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22) / det_ij);
      
      omega1 = Nonstat_params2(i, 0);
      log_lam1_1 = Nonstat_params2(i, 1);
      log_lam1_2 = Nonstat_params2(i, 2);
      
      omega2 = Nonstat_params2(j, 0);
      log_lam2_1 = Nonstat_params2(j, 1);
      log_lam2_2 = Nonstat_params2(j, 2);
      
      Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
      Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
      Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);
      
      Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
      Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
      Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);
      
      det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
      det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;
      
      Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
      Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
      Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);
      
      Inv_ij_11 = Kernel_ij_22; 
      Inv_ij_22 = Kernel_ij_11;
      Inv_ij_12 = - Kernel_ij_12; 
      det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;
      
      double sigma22 = sqrt(sqrt(det_i * det_j)/det_ij);
      
      double dist22 = sqrt((z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22) / det_ij);
      
      omega1 = Nonstat_params(i, 0);
      log_lam1_1 = Nonstat_params(i, 1);
      log_lam1_2 = Nonstat_params(i, 2);
      
      omega2 = Nonstat_params2(j, 0);
      log_lam2_1 = Nonstat_params2(j, 1);
      log_lam2_2 = Nonstat_params2(j, 2);
      
      Sigma1_11 = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
      Sigma1_12 = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
      Sigma1_22 = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);
      
      Sigma2_11 = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
      Sigma2_12 = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
      Sigma2_22 = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);
      
      det_i = Sigma1_11 * Sigma1_22 - Sigma1_12 * Sigma1_12;
      det_j = Sigma2_11 * Sigma2_22 - Sigma2_12 * Sigma2_12;
      
      Kernel_ij_11 = 0.5 * (Sigma1_11 + Sigma2_11);
      Kernel_ij_12 = 0.5 * (Sigma1_12 + Sigma2_12);
      Kernel_ij_22 = 0.5 * (Sigma1_22 + Sigma2_22);
      
      Inv_ij_11 = Kernel_ij_22; 
      Inv_ij_22 = Kernel_ij_11;
      Inv_ij_12 = - Kernel_ij_12; 
      det_ij = Kernel_ij_11 * Kernel_ij_22 - Kernel_ij_12 * Kernel_ij_12;
      
      double sigma12 = sqrt(sqrt(det_i * det_j)/det_ij);
      
      double dist12 = sqrt((z2(0) * z2(0) * Inv_ij_11 + z2(0) * z2(1) * Inv_ij_12 + z2(0) * z2(1) * Inv_ij_12 + z2(1) * z2(1) * Inv_ij_22) / det_ij);
      
      if (dist11 == 0) {
        cor11(i, j) = pow(sigma11, 2);
        cor22(i, j) = pow(sigma22, 2);
        cor12(i, j) = param(3) * sigma11 * sigma22;
      } else {
        cor11(i, j) = pow(sigma11, 2) * pow(dist11 / param(2), nu1) * cyl_bessel_k(nu1, dist11 / param(2))/(pow(2, nu1 - 1) * tgamma(nu1));
        cor22(i, j) = pow(sigma22, 2) * pow(dist22 / param(2), nu2) * cyl_bessel_k(nu2, dist22 / param(2))/(pow(2, nu2 - 1) * tgamma(nu2));
        cor12(i, j) = param(3) * sigma11 * sigma22 * pow(dist12 / param(2), 0.5 * (nu1 + nu2)) * cyl_bessel_k(0.5 * (nu1 + nu2), dist12 / param(2))/(pow(2, 0.5 * (nu1 + nu2) - 1) * tgamma(0.5 * (nu1 + nu2)));
      }
      //cor11(j, i) = cor11(i, j);
      //cor22(j, i) = cor22(i, j);
      //cor12(j, i) = cor12(i, j);
      cor21(j, i) = cor12(i, j);
      //cor21(i, j) = cor12(j, i);
    }
  }
  return rbind_cpp(cbind_cpp(cor11, cor12), cbind_cpp(cor21, cor22));
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
        cor(i, j) = sigma * pow(dist / 0.23, 1) * cyl_bessel_k(1, dist / 0.23)/(pow(2, 1 - 1) * tgamma(1));
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
        cor(i, j) = sigma * pow(dist / param(1), param(0)) * cyl_bessel_k(1, dist / param(1))/(pow(2, param(0) - 1) * tgamma(param(0)));
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
      cor(j) = sigma * pow(dist / param(1), param(0)) * cyl_bessel_k(1, dist / param(1))/(pow(2, param(0) - 1) * tgamma(param(0)));
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
      cor(j) = sigma * pow(dist / param(1), param(0)) * cyl_bessel_k(1, dist / param(1))/(pow(2, param(0) - 1) * tgamma(param(0)));
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
      cor(j) = sigma * pow(dist / param(1), param(0)) * cyl_bessel_k(1, dist / param(1))/(pow(2, param(0) - 1) * tgamma(param(0)));
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
      cor(j) = sigma * pow(dist / param(1), param(0)) * cyl_bessel_k(1, dist / param(1))/(pow(2, param(0) - 1) * tgamma(param(0)));
    }
  }
  return cor;
}



// [[Rcpp::export]]

List NEW_SPATIALLY_VARYING_PARAMETERS(NumericMatrix & Loc, NumericVector & param, NumericMatrix & wind, NumericVector & param_nonstat, int & time) {

	Eigen::MatrixXd SIG1(2, 2), SIG2(2, 2), SIG12(2, 2);
	Eigen::VectorXd z(2);
 	const int m = Loc.nrow(), n_wind = wind.nrow();
  	float sigma2 = param(0), beta = param(1), nu = param(2);
  	NumericMatrix cor(m, m), param_matrix(m / time, 5);
	List conso(2);
  
  	for (int i = 0; i < m; ++i) {
    		for (int j = 0; j <= i; ++j) {
 
 			double temp_val = 0.0;

			for (int k = 0; k < n_wind; ++k) {
      
				double new_loc1_x = Loc(i, 0) - wind(k, 0) * Loc(i, 2);
				double new_loc1_y = Loc(i, 1) - wind(k, 1) * Loc(i, 2);
				double new_loc2_x = Loc(j, 0) - wind(k, 0) * Loc(j, 2);
				double new_loc2_y = Loc(j, 1) - wind(k, 1) * Loc(j, 2);
				z[0] = new_loc1_x - new_loc2_x;
				z[1] = new_loc1_y - new_loc2_y;

				double omega1 = param_nonstat(0) + param_nonstat(1) * (new_loc1_x - .5) + param_nonstat(2) * (new_loc1_y - .5) + param_nonstat(3) * pow(new_loc1_x - .5, 2) + param_nonstat(4) * pow(new_loc1_y - .5, 2);
				double log_lam1_1 = param_nonstat(5) + param_nonstat(6) * (new_loc1_x - .5) + param_nonstat(7) * (new_loc1_y - .5) + param_nonstat(8) * pow(new_loc1_x - .5, 2) + param_nonstat(9) * pow(new_loc1_y - .5, 2);
				double log_lam1_2 = param_nonstat(10) + param_nonstat(11) * (new_loc1_x - .5) + param_nonstat(12) * (new_loc1_y - .5) + param_nonstat(13) * pow(new_loc1_x - .5, 2) + param_nonstat(14) * pow(new_loc1_y - .5, 2);
				
				double omega2 = param_nonstat(0) + param_nonstat(1) * (new_loc2_x - .5) + param_nonstat(2) * (new_loc2_y - .5) + param_nonstat(3) * pow(new_loc2_x - .5, 2) + param_nonstat(4) * pow(new_loc2_y - .5, 2);
				double log_lam2_1 = param_nonstat(5) + param_nonstat(6) * (new_loc2_x - .5) + param_nonstat(7) * (new_loc2_y - .5) + param_nonstat(8) * pow(new_loc2_x - .5, 2) + param_nonstat(9) * pow(new_loc2_y - .5, 2);
				double log_lam2_2 = param_nonstat(10) + param_nonstat(11) * (new_loc2_x - .5) + param_nonstat(12) * (new_loc2_y - .5) + param_nonstat(13) * pow(new_loc2_x - .5, 2) + param_nonstat(14) * pow(new_loc2_y - .5, 2);
			
				SIG1(0, 0) = exp(log_lam1_1) * cos(omega1) * cos(omega1) + exp(log_lam1_2) * sin(omega1) * sin(omega1);
				SIG1(0, 1) = SIG1(1, 0) = exp(log_lam1_1) * cos(omega1) * sin(omega1) - exp(log_lam1_2) * sin(omega1) * cos(omega1);
				SIG1(1, 1) = exp(log_lam1_1) * sin(omega1) * sin(omega1) + exp(log_lam1_2) * cos(omega1) * cos(omega1);
			      
				SIG2(0, 0) = exp(log_lam2_1) * cos(omega2) * cos(omega2) + exp(log_lam2_2) * sin(omega2) * sin(omega2);
				SIG2(0, 1) = SIG2(1, 0) = exp(log_lam2_1) * cos(omega2) * sin(omega2) - exp(log_lam2_2) * sin(omega2) * cos(omega2);
				SIG2(1, 1) = exp(log_lam2_1) * sin(omega2) * sin(omega2) + exp(log_lam2_2) * cos(omega2) * cos(omega2);

				SIG12 = 0.5 * (SIG1 + SIG2);

				if(Loc(i, 2) == 0){
					param_matrix(i, 0) = Loc(i, 0);
					param_matrix(i, 1) = Loc(i, 1);
					param_matrix(i, 2) = SIG1(0, 0);
					param_matrix(i, 3) = SIG1(0, 1);
					param_matrix(i, 4) = SIG1(1, 1);
				}
			      
				double det_i = SIG1.determinant();
				double det_j = SIG2.determinant();
				double det_ij = SIG12.determinant();

				//Rprintf("dist : %f, sigma : %f, sigma2 : %f, nu : %f \n", dist, sigma, sigma2, nu);
			      
				double sigma = sqrt(sqrt(det_i * det_j)/det_ij);
				double dist  = z.transpose() * SIG12.inverse() * z;

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

