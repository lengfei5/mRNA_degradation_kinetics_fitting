#include <Rcpp.h>
# include <math.h>
using namespace Rcpp;

/*
 * function to calculate s
 */
// [[Rcpp::export]]
NumericVector compute_s_beta(NumericVector t, double Min, double Amp, double phase, double beta){
  
  const double pi = 3.14159265358979323846;
  double w = 2.0*pi/24.0;
  
  int n = t.size();
  NumericVector s(n);
  
  for(int i=0; i<n; i++){
    s(i) = Min + Amp * pow(((1.0+cos(w*(t[i]-phase)))/2.0), beta);
  }
  
  return s;
}
/*compute.s.beta = function(t = zt, Min = 1, Amp = 2, phase = 12, beta = 1)
 *{
 * w = 2*pi/24;
 * s = Min+Amp*((1+cos(w*(t-phase)))/2)^beta
 * return(s)
 *}
 */


/*
 * function to calculate m
 */
// [[Rcpp::export]]
NumericVector compute_m_beta(NumericVector t, 
                             double gamma, double eps_gamma, double phase_gamma,double splicing_k,  
                             double Min, double Amp, double phase, double beta){
  const double pi = 3.14159265358979323846;
  double w = 2.0*pi/24.0;
  
  int n = t.size();
  NumericVector m(n);
  
  for(int i = 0; i < n; i++){
    m(i) = 0.2;
  }
  
  return m;
}




/*
 * compute.m.beta = function(t=seq(0, 94, by=2), gamma=log(2)/5, eps.gamma=0.2, phase.gamma=12, 
 * splicing.k=log(2)/(5/60), Min = 0.5, Amp=5, phase=12, beta=1, simulation.only=FALSE)
 * 
 * 
 */


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// A test function 
// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}
