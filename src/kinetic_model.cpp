//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::depends(RcppNumerical)]]
#include <Rcpp.h>
#include <math.h>
#include <RcppNumerical.h>

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
    s(i) = Min + Amp * pow(((1.0+cos(w*(t(i)-phase)))/2.0), beta);
  }
  
  return s;
  
}

double compute_s_t(double t, double Min, double Amp, double phase, double beta){

  const double pi = 3.14159265358979323846;
  double w = 2.0*pi/24.0;
  double s = Min + Amp * pow(((1.0+cos(w*(t-phase)))/2.0), beta);
  return s;

}

/*compute.s.beta = function(t = zt, Min = 1, Amp = 2, phase = 12, beta = 1)
 *{
 * w = 2*pi/24;
 * s = Min+Amp*((1+cos(w*(t-phase)))/2)^beta
 * return(s)
 *}
 */
class Fun_to_Integrate: public Numer::Func{
private:
  double gamma;
  double eps_gamma; 
  double phase_gamma;
  double Min;
  double Amp;
  double phase;
  double beta;
  
public:
  Fun_to_Integrate(double gamma_, double eps_gamma_, double phase_gamma_,
                   double Min_, double Amp_, double phase_, double beta_) : 
  gamma(gamma_), eps_gamma(eps_gamma_), phase_gamma(phase_gamma_), Min(Min_), Amp(Amp_), phase(phase_), beta(beta_){}
  
  double operator()(const double& t) const
  {
    const double pi = 3.14159265358979323846;
    double w = 2.0*pi/24.0;
    double Gamma = gamma * (t + eps_gamma/w * sin(w * (t - phase_gamma)));
    double s = Min + Amp * pow(((1.0+cos(w*(t-phase)))/2.0), beta);
    return exp(Gamma) * s;  
  }
  
}; 

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
  int i;
  NumericVector m(n);
  
  double lower = 0.0;
  double upper;
  double Tstable = 24*(ceil(log(2000)/gamma/24) + 1.0);
  
  Fun_to_Integrate f(gamma, eps_gamma, phase_gamma, Min, Amp, phase, beta);
  double err_est;
  int err_code;
  double res_integrate;
  double res_Gamma;
  
  for(i = 0; i < n; i++){
    upper = Tstable + t[i];
    
    res_integrate = Numer::integrate(f, lower, upper, err_est, err_code);
    res_Gamma =  gamma * (upper + eps_gamma/w * sin(w * (upper - phase_gamma)));
    m[i] = splicing_k * exp(-res_Gamma) * res_integrate ;
    //m[i] = splicing.k*exp(-Gamma(t = Tstable+time, gamma = gamma, eps.gamma=eps.gamma, phase.gamma=phase.gamma)) *
    //integrate(f2integrate, lower = 0, upper = Tstable+time, par = c(gamma, eps.gamma, phase.gamma, Min, Amp, phase, beta))$value
    
  }
  
  return m;
}

/*
 *  Gamma = function(t = 0, gamma = log(2)/3, eps.gamma = 0.2, phase.gamma = 0)
 {
w = 2*pi/24;
Gamma = gamma*(t + eps.gamma/w * sin(w*(t-phase.gamma)))
#Gamma = gamma*(t+ eps.gamma/w * sin(w*(t-phase.gamma)))
#Gamma = gamma*(t+ eps.gamma/w * sin(w*(t-phase.gamma)))
#Gamma = gamma*t + amp.gamma/2*(t+ sin(w*(t-phase.gamma))/w)
return(Gamma)

*/

double Gamma_copy(double t, double gamma, double eps_gamma, double phase_gamma){
  
  const double pi = 3.14159265358979323846;
  double w = 2.0*pi/24.0;
  double res_Gamma = gamma * (t + eps_gamma/w * sin(w * (t - phase_gamma)));
  
  return res_Gamma;
}

/*
 * f2integrate = function(t, par)
 {
gamma = par[1]; 
eps.gamma= par[2];
phase.gamma= par[3]; 
Min = par[4];
Amp = par[5];
phase = par[6];
beta = par[7];
return(exp(Gamma(t, gamma, eps.gamma, phase.gamma)) * compute.s.beta(t, Min, Amp, phase, beta))
}
 */


/*
 * 
 * integrate.m = function(t=seq(0, 46, by=2), 
 gamma=log(2)/5, eps.gamma=0.2, phase.gamma=12, splicing.k=log(2)/(5/60), 
  Min = 0.5, Amp=2.0, phase=12, beta=1)
  {
  w = 2*pi/24;
  m = rep(0, length(t))
  
  Tstable = 24*(ceiling(log(2000)/gamma/24) +1)
  for(i in 1:length(t))
  {
  time = t[i]; 
  m[i] = splicing.k*exp(-Gamma(t = Tstable+time, gamma = gamma, eps.gamma=eps.gamma, phase.gamma=phase.gamma)) *
 integrate(f2integrate, lower = 0, upper = Tstable+time, par = c(gamma, eps.gamma, phase.gamma, Min, Amp, phase, beta))$value
  }
  return(m)
  }
  }
 * 
 */

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
