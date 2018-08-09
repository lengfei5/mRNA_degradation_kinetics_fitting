//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::depends(RcppNumerical)]]
//#include <Rcpp.h>
#include <math.h>
#include <RcppNumerical.h>

using namespace Rcpp;
class Func
{
public:
  virtual double operator()(const double& x) const = 0;
  virtual void eval(double* x, const int n) const
  {
    for(int i = 0; i < n; i++)
      x[i] = this->operator()(x[i]);
  }
  
  virtual ~Func() {}
};

inline double integrate(
    const Func& f, const double& lower, const double& upper,
    double& err_est, int& err_code,
    const int subdiv = 100, const double& eps_abs = 1e-8, const double& eps_rel = 1e-6,
    const Rcpp Integrator<double>::QuadratureRule rule = Integrator<double>::GaussKronrod41
)
  

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
  int i;
  NumericVector m(n);
  
  for(i = 0; i < n; i++){
    
    m(i) = 0.2 *  Min + Amp * pow(((1.0+cos(w*(t(i)-phase)))/2.0), beta);
    
  }
  
  return m;
  
}

/*
 * integration
 */
class func_to_integrate: public RcppNumerical::Func


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
  m[i] = splicing.k*exp(-Gamma(t = Tstable+time, gamma = gamma, eps.gamma=eps.gamma, phase.gamma=phase.gamma)) *integrate(f2integrate, lower = 0, upper = Tstable+time, par = c(gamma, eps.gamma, phase.gamma, Min, Amp, phase, beta))$value
  }
  return(m)
  }
  f2integrate = function(t, par)
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
#### This is the integral of degradation function
  Gamma = function(t = 0, gamma = log(2)/3, eps.gamma = 0.2, phase.gamma = 0)
  {
  w = 2*pi/24;
  Gamma = gamma*(t + eps.gamma/w * sin(w*(t-phase.gamma)))
#Gamma = gamma*(t+ eps.gamma/w * sin(w*(t-phase.gamma)))
#Gamma = gamma*(t+ eps.gamma/w * sin(w*(t-phase.gamma)))
#Gamma = gamma*t + amp.gamma/2*(t+ sin(w*(t-phase.gamma))/w)
  return(Gamma)
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
