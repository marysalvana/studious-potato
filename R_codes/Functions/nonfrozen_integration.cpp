// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <Rcpp.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <RcppNumerical.h>
#include <iostream>

using namespace Rcpp;
using namespace boost::math;
using namespace std;
using namespace Numer;

class Mvn
{
public:
  Mvn(const Eigen::VectorXd& mu,
      const Eigen::MatrixXd& s);
  ~Mvn();
  double pdf(const Eigen::VectorXd& x) const;
  Eigen::VectorXd sample(unsigned int nr_iterations = 20) const;
  Eigen::VectorXd mean;
  Eigen::MatrixXd sigma;
};

double Mvn::pdf(const Eigen::VectorXd& x) const
{
  double n = x.rows();
  double sqrt2pi = std::sqrt(2 * M_PI);
  double quadform  = (x - mean).transpose() * sigma.inverse() * (x - mean);
  double norm = std::pow(sqrt2pi, - n) *
                std::pow(sigma.determinant(), - 0.5);

  return norm * exp(-0.5 * quadform);
}

class PDF: public MFunc
{
private:
  double mu1, mu2;
  double sigma11, sigma22, sigma12;
  double h1, h2;
  int u;
  double var, nu, beta;

  double denom = 0.0, logretval = 0.0, determinant = 0.0, dist0 = 0.0, val = 0.0;

public:
  PDF(double h1_, double h2_, int u_, double mu1_, double mu2_, double sigma11_, double sigma22_, double sigma12_, double var_, double nu_, double beta_): h1(h1_), h2(h2_), u(u_), mu1(mu1_), mu2(mu2_), sigma11(sigma11_), sigma22(sigma22_), sigma12(sigma12_), var(var_), nu(nu_), beta(beta_) {};

  double operator()(Constvec& x)
  { 
	Eigen::MatrixXd sigma(2, 2);
	sigma(0, 0) = sigma11;
	sigma(0, 1) = sigma12;
	sigma(1, 0) = sigma12;
	sigma(1, 1) = sigma22;
	
	Eigen::VectorXd mean(2);
	mean[0] = mu1;
	mean[1] = mu2;

	Eigen::VectorXd x_new(2);
	x_new[0] = x[0];
	x_new[1] = x[1];
	
	double n = x_new.rows();
  	double sqrt2pi = std::sqrt(2 * M_PI);
  	double quadform  = (x_new - mean).transpose() * sigma.inverse() * (x_new - mean);
  	double norm = std::pow(sqrt2pi, - n) * std::pow(sigma.determinant(), - 0.5);
  
  	double temp = norm * exp(-0.5 * quadform);
	//determinant = sigma11 * sigma22 - pow(sigma12, 2);
	//denom = 2 * M_PI * sqrt(determinant);
	//logretval = pow(x[0] - mu1, 2) * sigma22 - 2 * (x[0] - mu1) * (x[1] - mu2) * sigma12 + pow(x[1] - mu2, 2) * sigma11;
	//dist0 = (h1 - x[0] * u) * (h1 - x[0] * u) + (h2 - x[1] * u) * (h2 - x[1] * u);
	dist0 = sqrt((h1 - x[0] * u) * (h1 - x[0] * u) + (h2 - x[1] * u) * (h2 - x[1] * u));
	
	//cout << x[0] << ' ' << x[1] << endl;	

	if (dist0 == 0){
        	//val = var * exp(-0.5 * logretval / determinant) / denom;
        	val = var * temp;
	} else{
        	//val = var * exp(- dist0 / beta) * temp;
        	//val = var * exp(- dist0 / beta) * exp(-0.5 * logretval / determinant) / denom;
        	val = var * pow(dist0, nu) * cyl_bessel_k(nu, dist0) / (pow(2, nu - 1) * tgamma(nu)) * temp;
        	//val = var * pow(dist0 / beta, nu) * cyl_bessel_k(nu, dist0 / beta) / (pow(2, nu - 1) * tgamma(nu)) * temp;
        	//val = var * pow(dist0 / beta, nu) * cyl_bessel_k(nu, dist0 / beta) / (pow(2, nu - 1) * tgamma(nu)) * exp(-0.5 * logretval / determinant) / denom;
	}
	return val;
  }
};

class PDF_nonstat: public MFunc
{
private:
  double mu1, mu2;
  double sigma11, sigma22, sigma12;
  double x1, x2, y1, y2, new_x1, new_x2, new_y1, new_y2;
  int t1, t2;
  double var, nu, beta;

  double denom = 0.0, logretval = 0.0, determinant = 0.0, dist0 = 0.0, val = 0.0;

public:
  PDF_nonstat(double x1_, double x2_, double y1_, double y2_, int t1_, int t2_, double mu1_, double mu2_, double sigma11_, double sigma22_, double sigma12_, double var_, double nu_, double beta_): x1(x1_), x2(x2_), y1(y1_), y2(y2_), t1(t1_), t2(t2_), mu1(mu1_), mu2(mu2_), sigma11(sigma11_), sigma22(sigma22_), sigma12(sigma12_), var(var_), nu(nu_), beta(beta_) {};

  double operator()(Constvec& x)
  { 
	Eigen::MatrixXd sigma(2, 2);
	sigma(0, 0) = sigma11;
	sigma(0, 1) = sigma12;
	sigma(1, 0) = sigma12;
	sigma(1, 1) = sigma22;
	
	Eigen::VectorXd mean(2);
	mean[0] = mu1;
	mean[1] = mu2;

	Eigen::VectorXd x_new(2);
	x_new[0] = x[0];
	x_new[1] = x[1];
	
	double n = x_new.rows();
  	double sqrt2pi = std::sqrt(2 * M_PI);
  	double quadform  = (x_new - mean).transpose() * sigma.inverse() * (x_new - mean);
  	double norm = std::pow(sqrt2pi, - n) * std::pow(sigma.determinant(), - 0.5);
  
  	double temp = norm * exp(-0.5 * quadform);

      new_x1 = x1 - x[0] * (t1 - 1);
      new_x2 = x2 - x[1] * (t1 - 1);
      new_y1 = y1 - x[0] * (t2 - 1);
      new_y2 = y2 - x[1] * (t2 - 1);
      double z1 = new_x1 - new_y1;
      double z2 = new_x2 - new_y2;

      double omega1 = 1 * (new_x1 - .5) + 1 * (new_x2 - .5) + 3 * pow(new_x1 - .5, 2) + -1 * pow(new_x2 - .5, 2);
      double log_lam1_1 = exp(- pow(new_x1, 2)) - exp(- pow(new_x2, 2));
      double log_lam1_2 = sin(new_x1) * exp(- pow(new_x2, 2)) - sin(new_x2) * exp(- pow(new_x1, 2));

      double omega2 = 1 * (new_y1 - .5) + 1 * (new_y2 - .5) + 3 * pow(new_y1 - .5, 2) + -1 * pow(new_y2 - .5, 2);
      double log_lam2_1 = exp(- pow(new_y1, 2)) - exp(- pow(new_y2, 2));
      double log_lam2_2 = sin(new_y1) * exp(- pow(new_y2, 2)) - sin(new_y2) * exp(- pow(new_y1, 2));

      //double omega1 = sin(4 * new_x1) * sin(4 * new_x2);
      //double log_lam1_1 = sin(6 * new_x1) + sin(6 * new_x2);
      //double log_lam1_2 = cos(4 * new_x1) * cos(4 * new_x2);
      
      //double omega2 = sin(4 * new_y1) * sin(4 * new_y2);
      //double log_lam2_1 = sin(6 * new_y1) + sin(6 * new_y2);
      //double log_lam2_2 = cos(4 * new_y1) * cos(4 * new_y2);

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
      
      double VAR = sqrt(sqrt(det_i * det_j)/det_ij);
      
      dist0 = sqrt((z1 * z1 * Inv_ij_11 + z1 * z2 * Inv_ij_12 + z1 * z2 * Inv_ij_12 + z2 * z2 * Inv_ij_22) / det_ij);
	//cout << det_ij << ' ' << VAR << endl;	
	//cout << i << ' ' << j << endl;	
	
	if (dist0 == 0){
        	val = var * pow(VAR, 2) * temp;
	} else{
        	//val = VAR * exp(- dist0 / beta) * temp;
        	val = var * pow(VAR, 2) * pow(dist0, nu) * cyl_bessel_k(nu, dist0) / (pow(2, nu - 1) * tgamma(nu)) * temp;
	}
	return val;
  }
};

class PDF_deform: public MFunc
{
private:
  double mu1, mu2;
  double sigma11, sigma22, sigma12;
  double x1, x2, y1, y2, new_x1, new_x2, new_y1, new_y2;
  int t1, t2;
  double var, nu, beta;

  double denom = 0.0, logretval = 0.0, determinant = 0.0, dist0 = 0.0, val = 0.0;

public:
  PDF_deform(double x1_, double x2_, double y1_, double y2_, int t1_, int t2_, double mu1_, double mu2_, double sigma11_, double sigma22_, double sigma12_, double var_, double nu_, double beta_): x1(x1_), x2(x2_), y1(y1_), y2(y2_), t1(t1_), t2(t2_), mu1(mu1_), mu2(mu2_), sigma11(sigma11_), sigma22(sigma22_), sigma12(sigma12_), var(var_), nu(nu_), beta(beta_) {};

  double operator()(Constvec& x)
  { 
	Eigen::MatrixXd sigma(2, 2);
	sigma(0, 0) = sigma11;
	sigma(0, 1) = sigma12;
	sigma(1, 0) = sigma12;
	sigma(1, 1) = sigma22;
	
	Eigen::VectorXd mean(2);
	mean[0] = mu1;
	mean[1] = mu2;

	Eigen::VectorXd x_new(2);
	x_new[0] = x[0];
	x_new[1] = x[1];
	
	double n = x_new.rows();
  	double sqrt2pi = std::sqrt(2 * M_PI);
  	double quadform  = (x_new - mean).transpose() * sigma.inverse() * (x_new - mean);
  	double norm = std::pow(sqrt2pi, - n) * std::pow(sigma.determinant(), - 0.5);
  
  	double temp = norm * exp(-0.5 * quadform);

      new_x1 = x1 - x[0] * (t1 - 1);
      new_x2 = x2 - x[1] * (t1 - 1);
      new_y1 = y1 - x[0] * (t2 - 1);
      new_y2 = y2 - x[1] * (t2 - 1);

	double dist_source1 = sqrt(pow(new_x1 - 0.15, 2) + pow(new_x2 - 0.15, 2)); 
	double dist_source2 = sqrt(pow(new_y1 - 0.15, 2) + pow(new_y2 - 0.15, 2)); 
	double deform_x1 = 0.15 + (new_x1 - 0.15) * ( 1 + 2 * exp(-0.5 * pow(dist_source1, 2 )));
	double deform_x2 = 0.15 + (new_x2 - 0.15) * ( 1 + 2 * exp(-0.5 * pow(dist_source1, 2 )));
	double deform_y1 = 0.15 + (new_y1 - 0.15) * ( 1 + 2 * exp(-0.5 * pow(dist_source2, 2 )));
	double deform_y2 = 0.15 + (new_y2 - 0.15) * ( 1 + 2 * exp(-0.5 * pow(dist_source2, 2 )));

      dist0 = sqrt(pow(deform_x1 - deform_y1, 2) + pow(deform_x2 - deform_y2, 2));
	//cout << det_ij << ' ' << VAR << endl;	
	//cout << i << ' ' << j << endl;	
	
	if (dist0 == 0){
        	val = var * temp;
	} else{
        	val = var * pow(dist0, nu) * cyl_bessel_k(nu, dist0) / (pow(2, nu - 1) * tgamma(nu)) * temp;
	}
	return val;
  }
};

class PDF3D: public MFunc
{
private:
  double mu1, mu2, mu3;
  double sigma11, sigma22, sigma33, sigma12, sigma13, sigma23;
  double x1, x2, y1, y2, z1, z2, h1 = 0.0, h2 = 0.0, h3 = 0.0;
  int t1, t2;
  double var, beta;

  double dist0 = 0.0, val = 0.0;

public:
  PDF3D(double x1_, double x2_, double y1_, double y2_, int t1_, int t2_, double mu1_, double mu2_, double mu3_, double sigma11_, double sigma22_, double sigma33_, double sigma12_, double sigma13_, double sigma23_, double var_, double beta_): x1(x1_), x2(x2_), y1(y1_), y2(y2_), t1(t1_), t2(t2_), mu1(mu1_), mu2(mu2_), mu3(mu3_), sigma11(sigma11_), sigma22(sigma22_), sigma33(sigma33_), sigma12(sigma12_), sigma13(sigma13_), sigma23(sigma23_), var(var_), beta(beta_) {};

  double operator()(Constvec& x)
  { 
	Eigen::MatrixXd sigma(3, 3);
	sigma(0, 0) = sigma11;
	sigma(0, 1) = sigma12;
	sigma(0, 2) = sigma13;
	sigma(1, 0) = sigma12;
	sigma(1, 1) = sigma22;
	sigma(1, 2) = sigma23;
	sigma(2, 0) = sigma13;
	sigma(2, 1) = sigma23;
	sigma(2, 2) = sigma33;
	
	Eigen::VectorXd mean(3);
	mean[0] = mu1;
	mean[1] = mu2;
	mean[2] = mu3;

	Eigen::VectorXd x_new(3);
	x_new[0] = x[0];
	x_new[1] = x[1];
	x_new[2] = x[2];
	
	double n = x_new.rows();
  	double sqrt2pi = std::sqrt(2 * M_PI);
  	double quadform  = (x_new - mean).transpose() * sigma.inverse() * (x_new - mean);
  	double norm = std::pow(sqrt2pi, - n) * std::pow(sigma.determinant(), - 0.5);
  
  	double temp = norm * exp(-0.5 * quadform);

	z1 = (x1 - x[0] * t1 - 0.5) * sqrt(pow(x1 - x[0] * t1 - 0.5, 2) + pow(y1 - x[1] * t1 - 0.5, 2));
	z2 = (x2 - x[0] * t2 - 0.5) * sqrt(pow(x2 - x[0] * t1 - 0.5, 2) + pow(y2 - x[1] * t2 - 0.5, 2));

	//z1 = 1 + (x1 - x[0] * t1 - 0.5) + (y1 - x[1] * t1 - 0.5) + 3 * pow(x1 - x[0] * t1 - 0.5, 2) - pow(y1 - x[1] * t1 - 0.5, 2);
	//z2 = 1 + (x2 - x[0] * t2 - 0.5) + (y2 - x[1] * t2 - 0.5) + 3 * pow(x2 - x[0] * t2 - 0.5, 2) - pow(y2 - x[1] * t2 - 0.5, 2);

	h1 = x2 - x1 - (x[0] * t2 - x[0] * t1);
	h2 = y2 - y1 - (x[1] * t2 - x[1] * t1);
	h3 = z2 - z1 - (x[2] * t2 - x[2] * t1);

	dist0 = pow(h1, 2) + pow(h2, 2) + pow(h3, 2);
	
	if (dist0 == 0){
        	val = var * temp;
	} else{
        	val = var * exp(- dist0 / beta) * temp;
	}
	return val;
  }
};

// [[Rcpp::export]]

double integrate_velocity(double h1, double h2, int u, double mu1, double mu2, double sigma11, double sigma22, double sigma12, double var, double nu, double beta, Eigen::VectorXd  lower, Eigen::VectorXd upper)
{
  PDF f(h1, h2, u, mu1, mu2, sigma11, sigma22, sigma12, var, nu, beta);
  double err_est;
  int err_code;
  const int maxeval = 5000;
  double res = integrate ( f, lower, upper,err_est,err_code, maxeval);
  return res;
    
}
// [[Rcpp::export]]

double integrate_velocity_nonstat(double x1, double x2, double y1, double y2, int t1, int t2, double mu1, double mu2, double sigma11, double sigma22, double sigma12, double var, double nu, double beta, Eigen::VectorXd  lower, Eigen::VectorXd upper)
{
  PDF_nonstat f(x1, x2, y1, y2, t1, t2, mu1, mu2, sigma11, sigma22, sigma12, var, nu, beta);
  double err_est;
  int err_code;
  const int maxeval = 5000;
  double res = integrate ( f, lower, upper,err_est,err_code, maxeval);
  return res;
    
}

// [[Rcpp::export]]

double integrate_velocity_deform(double x1, double x2, double y1, double y2, int t1, int t2, double mu1, double mu2, double sigma11, double sigma22, double sigma12, double var, double nu, double beta, Eigen::VectorXd  lower, Eigen::VectorXd upper)
{
  PDF_deform f(x1, x2, y1, y2, t1, t2, mu1, mu2, sigma11, sigma22, sigma12, var, nu, beta);
  double err_est;
  int err_code;
  const int maxeval = 5000;
  double res = integrate ( f, lower, upper,err_est,err_code, maxeval);
  return res;
    
}



// [[Rcpp::export]]

double integrate_velocity_3D(double x1, double x2, double y1, double y2, int t1, int t2, double mu1, double mu2, double mu3, double sigma11, double sigma22, double sigma33, double sigma12, double sigma13, double sigma23, double var, double beta, Eigen::VectorXd lower, Eigen::VectorXd upper)
{
  PDF3D f(x1, x2, y1, y2, t1, t2, mu1, mu2, mu3, sigma11, sigma22, sigma33, sigma12, sigma13, sigma23, var, beta);
  double err_est;
  int err_code;
  const int maxeval = 100000;
  double res = integrate(f, lower, upper, err_est, err_code, maxeval);
  return res;
}
// [[Rcpp::export]]

NumericMatrix nonfrozen_matern_uni_cpp(NumericMatrix & Loc, NumericVector & params) {
  const int m = Loc.nrow();
  
  NumericMatrix cor(m, m);
  Eigen::VectorXd upper(2);
  Eigen::VectorXd lower(2);

  double var = params(0), nu = params(2), beta = params(1);
  double mu1 = params(3), mu2 = params(4), sigma11 = params(5), sigma22 = params(6), sigma12 = params(7);

  upper[0] = 5, upper[1] = 5, lower[0] = -5, lower[1] = -5;
  
//  for (int i = 0; i < m; ++i) {
//    for (int j = 0; j < m; ++j) {
for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
      
	cout << i << ' ' << j << endl;	
      double x1 = Loc(i, 0) - Loc(j, 0);
      double x2 = Loc(i, 1) - Loc(j, 1);
      int tloc = Loc(i, 2) - Loc(j, 2);

      double val = integrate_velocity(x1, x2, tloc, mu1, mu2, sigma11, sigma22, sigma12, var, nu, beta,  lower,  upper);
      cor(i, j) = val;
      cor(j, i) = cor(i, j);
    }
  }
  return cor;
}

// [[Rcpp::export]]

NumericMatrix nonfrozen_matern_uni_cpp_distributed(NumericMatrix & Loc1, NumericMatrix & Loc2, NumericVector & params) {
  const int m1 = Loc1.nrow();
  const int m2 = Loc2.nrow();
  
  NumericMatrix cor(m1, m2);
  Eigen::VectorXd upper(2);
  Eigen::VectorXd lower(2);

  upper[0] = 1, upper[1] = 1, lower[0] = -1, lower[1] = -1;

  double var = params(0), nu = params(2), beta = params(1);
  double mu1 = params(3), mu2 = params(4), sigma11 = params(5), sigma22 = params(6), sigma12 = params(7);
  
  for (int i = 0; i < m1; ++i) {
    for (int j = 0; j < m2; ++j) {
      
      double x1 = Loc1(i, 0) - Loc2(j, 0);
      double x2 = Loc1(i, 1) - Loc2(j, 1);
      int tloc = Loc1(i, 2) - Loc2(j, 2);

      double val = integrate_velocity(x1, x2, tloc, mu1, mu2, sigma11, sigma22, sigma12, var, nu, beta,  lower,  upper);
      cor(i, j) = val;
    }
  }
  return cor;
}

// [[Rcpp::export]]

NumericMatrix nonfrozen_matern_uni_cpp_distributed_nonstat(NumericMatrix & Loc1, NumericMatrix & Loc2, NumericVector & params) {
  const int m1 = Loc1.nrow();
  const int m2 = Loc2.nrow();
  
  NumericMatrix cor(m1, m2);
  Eigen::VectorXd upper(2);
  Eigen::VectorXd lower(2);

  upper[0] = 5, upper[1] = 5, lower[0] = -5, lower[1] = -5;

  double var = params(0), nu = params(2), beta = params(1);
  double mu1 = params(3), mu2 = params(4), sigma11 = params(5), sigma22 = params(6), sigma12 = params(7);
  
  for (int i = 0; i < m1; ++i) {
    for (int j = 0; j < m2; ++j) {
      
      double x1 = Loc1(i, 0);
      double y1 = Loc2(j, 0);
      double x2 = Loc1(i, 1);
      double y2 = Loc2(j, 1);
      int t1 = Loc1(i, 2);
      int t2 = Loc2(j, 2);

      double val = integrate_velocity_nonstat(x1, x2, y1, y2, t1, t2, mu1, mu2, sigma11, sigma22, sigma12, var, nu, beta,  lower,  upper);
      cor(i, j) = val;
    }
  }
  return cor;
}

// [[Rcpp::export]]

NumericMatrix nonfrozen_matern_uni_cpp_distributed_deform(NumericMatrix & Loc1, NumericMatrix & Loc2, NumericVector & params) {
  const int m1 = Loc1.nrow();
  const int m2 = Loc2.nrow();
  
  NumericMatrix cor(m1, m2);
  Eigen::VectorXd upper(2);
  Eigen::VectorXd lower(2);

  upper[0] = 5, upper[1] = 5, lower[0] = -5, lower[1] = -5;

  double var = params(0), nu = params(2), beta = params(1);
  double mu1 = params(3), mu2 = params(4), sigma11 = params(5), sigma22 = params(6), sigma12 = params(7);
  
  for (int i = 0; i < m1; ++i) {
    for (int j = 0; j < m2; ++j) {
      
      double x1 = Loc1(i, 0);
      double y1 = Loc2(j, 0);
      double x2 = Loc1(i, 1);
      double y2 = Loc2(j, 1);
      int t1 = Loc1(i, 2);
      int t2 = Loc2(j, 2);

      double val = integrate_velocity_deform(x1, x2, y1, y2, t1, t2, mu1, mu2, sigma11, sigma22, sigma12, var, nu, beta,  lower,  upper);
      cor(i, j) = val;
    }
  }
  return cor;
}
//source: https://stackoverflow.com/questions/53648355/double-integral-in-rcppnumerical
//http://blog.sarantop.com/notes/mvn
