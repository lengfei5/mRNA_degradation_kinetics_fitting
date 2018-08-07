#include <assert.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>

using namespace std;

// Class Model
// =================================
class model {

 public:

  // Functions
  model(char* mynamedata, int mtype, double* myparams, double* myptypes, double* mysigma, double** mypriorparams, vector <int> &dtype, vector <int> &t0id, vector <vector <double> > &mytimes, vector <vector <double> > &mydata);
  void get_solutions();
  void set_postsigm();
  void set_sigmaupdate(double factor);
  void fit_postprop(int numburn);
  void set_params(double* myparams);
  void set_mcmcparams(double logprob);
  void print_fit(string fileout);
  void print_mcmc(string fileout);
  void print_hist(string fileout);
  void print_path(string fileout);
  void print_postprop(string fileout);
  void print_results(int nummodels, double modelprop[]);
  void print_results2(int nummodels, double modelprob[]);
  void print_results3(int nummodels, double modelprob[]);
  void get_mcmc_median(int numburn);  
  void get_mcmc_mean(int numburn);
  void get_mcmc_histograms(int numburn);
  void get_mcmc_functions(int numburn);
  void get_mcmc_maxpost(); 
  void get_mcmc_path();
  void get_covmat(int numburn);
  double* get_newparams(int controlprop);
  double get_postproppdf(int controlprop);
  double get_propratio(int controlprop);
  double get_postproppdf();
  double get_logprob();
  double get_logprior();
  ~model();

 private:
  // Variables model
  int seed;
  int numbins;
  int numdata;
  int modeltype;
  int numparams;
  int numparsamp;
  char *namedata;
  double meanclock;
  double pi;
  double* sigma;
  double* ptypes;
  double* params;
  double* newparams;
  double* mcmc_mean;
  double* mcmc_sigm;
  double* mcmc_median;
  double* mcmc_quan05;
  double* mcmc_quan95;
  double* mcmc_maxpost;
  double* mcmc_maxpara;
  double* hist_dx;
  double* hist_x0;
  double** mcmc_hist;
  double** priorparams;
  double** postpropparams;
  vector <int> datatype;
  vector <int> t0index;
  vector < vector <double> > data;
  vector < vector <double> > times;
  vector < vector <double> > solution;
  vector < vector <double> > mcmc_path;
  vector < vector <double> > mcmc_path2;
  vector < vector <double> > mcmc_params;
  vector <double> mcmc_logprob;

  vector <double> theo_times;
  vector <double> theo_synth;
  vector <double> theo_mrna;
  vector <double> mcmc_func;
  
  static double dt;
  static vector <double> funtimes;
  static vector <double> syntab;
  static vector <double> degtab;

  const gsl_rng_type *gsl_R;
  gsl_rng *gsl_r;

  // Variables ode solver
  int dimension;
  double eps_abs;
  double eps_rel;

  gsl_odeiv_step *step_ptr; 
  gsl_odeiv_control *control_ptr;
  gsl_odeiv_evolve *evolve_ptr;
  gsl_odeiv_system my_system; 
  const gsl_odeiv_step_type *type_ptr;  

  // Functions
  static int diffequation(double t, const double y[], double f[], void *params_ptr);
  static double laurafun(double t, vector <double> &table);
  static double mod(double a, double b);
  void lauratab(vector <double> &table, double a, double b, double T, double f, double p, int control);
  void laurafill();

  void get_accumulation(vector <double> &mytimes, vector <double> &mysolution);
  void get_synthesis(vector <double> &mytimes, vector <double> &mysolution);
  void get_decay(vector <double> &mytimes, vector <double> &mysolutions, int i0);
  void get_sum(int n, double &x2, double &m2, double &xm);
  void get_initvalue(double &y, int &numperiods);
  void set_sigma(double* mysigma); 
  void set_ptypes(double* myptypes); 
  void set_newparams(double* myparams);
  void fit_lognormal(int n, int numburn, double &par1, double &par2);
  void fit_beta(int n, int numburn, double &par1, double &par2);
  void fit_vonmises(int n, int numburn, double &par1, double &par2);
  void get_betaparams(int n, double p, double &s, double &a, double &b);
  string itos(int i);
};
