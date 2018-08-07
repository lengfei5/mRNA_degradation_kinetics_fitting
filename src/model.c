#include <assert.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <time.h>
#include <algorithm>
#include "model.h"
#include "probdist.h"

using namespace std;

double model::dt = 0.5;
vector <double> model::syntab = (*new vector <double>);
vector <double> model::degtab = (*new vector <double>);
vector <double> model::funtimes = (*new vector <double>);

// Model functions
//==================================================================

// Constructor
//**************
model::model(char *mynamedata, int mtype, double* myparams, double* myptypes, double* mysigma, double** mypriorparams, vector <int> &dtype, vector <int> &t0id, vector <vector <double> > &mytimes, vector <vector <double> > &mydata) {
  
  // Variables
  pi = 3.14159265359;
  seed = time(NULL); 
  namedata  = mynamedata;
  numparams = 10;
  dimension = 1;
  modeltype = mtype;
  datatype  = dtype;
  t0index   = t0id;
  data      = mydata;
  times     = mytimes;
  numdata   = data.size();
  priorparams = mypriorparams;
  sigma     = new double[numparams];
  ptypes    = new double[numparams];
  params    = new double[numparams];
  newparams = new double[numparams];
  mcmc_mean = new double[numparams];
  mcmc_sigm = new double[numparams];
  mcmc_median = new double[numparams];
  mcmc_quan05 = new double[numparams];
  mcmc_quan95 = new double[numparams];
  mcmc_maxpost = new double[numparams];
  mcmc_maxpara = new double[numparams];
  hist_dx        = new double[numparams];
  hist_x0        = new double[numparams];
  postpropparams = new double*[numparams];
  mcmc_hist      = new double*[numparams]; 
  mcmc_logprob   = *(new vector <double>);
  mcmc_params    = *(new vector <vector <double> >(numparams));
  solution       = *(new vector <vector <double> >(numdata));
  mcmc_path      = *(new vector <vector <double> >(numdata));
  mcmc_path2     = *(new vector <vector <double> >(numdata));
  numparsamp = 0;
  meanclock = 0;
  numbins = 50;

  // NOT GENERAL! Wrong!
  theo_times = *(new vector <double>);
  theo_synth = *(new vector <double>);
  theo_mrna  = *(new vector <double>);
  mcmc_func  =  *(new vector <double>(6));
  for (double i = 0; i < 24; i += dt) {
    theo_times.push_back(i);
    theo_synth.push_back(0.);
    theo_mrna.push_back(0.);
  }
  
  if (funtimes.empty()) {
    for (double i = 0; i < 24; i += dt) {
      funtimes.push_back(i);
      syntab.push_back(0.);
      degtab.push_back(0.);
    }
  }
  
  set_sigma(mysigma);
  set_ptypes(myptypes);
  set_params(myparams);
  set_newparams(myparams);
  
  for (int n = 0; n < numparams; n++) {
    if (ptypes[n] != 0) {numparsamp++;}
    mcmc_params[n] = (*(new vector <double>));
    postpropparams[n] = new double[2];
    mcmc_hist[n] = new double[numbins]; 
  }

  for (int n = 0; n < numdata; n++) {
    solution[n] = (*(new vector <double>));
    mcmc_path[n] = (*(new vector <double>));
    mcmc_path2[n] = (*(new vector <double>));
    for (int i = 0; i < data[n].size(); i++) {
      solution[n].push_back(0.);
      mcmc_path[n].push_back(0.);
      mcmc_path2[n].push_back(0.);
    }
  }
   
  // gsl rnd
  gsl_rng_env_setup();
  gsl_R = gsl_rng_default;
  gsl_r = gsl_rng_alloc(gsl_R);
  gsl_rng_set(gsl_r,seed);

  // ode solver varibales 
  eps_abs = 1.e-8;       // absolute error requested 
  eps_rel = 1.e-10;      // relative error requested 
  
  // allocate/initialize the stepper, the control function, and the evolution function.
  type_ptr    = gsl_odeiv_step_rkf45;
  step_ptr    = gsl_odeiv_step_alloc(type_ptr,dimension);
  control_ptr = gsl_odeiv_control_y_new(eps_abs,eps_rel);
  evolve_ptr  = gsl_odeiv_evolve_alloc(dimension);
  
  // load values into the my_system structure   
  my_system.function  = diffequation;	// the right-hand-side functions dy[i]/dt
  my_system.dimension = dimension;	// number of diffeq's
  my_system.params    = newparams;	// parameters to pass to rhs and jacobian 
}


// Destructor
//**************
model::~model() {
  
  // Variables
  delete &data;
  delete &times;
  delete &params;
  delete &solution;

  // all done; free up the gsl_odeiv stuff
  if (modeltype == 3 || modeltype == 4) { 
    gsl_odeiv_evolve_free(evolve_ptr);
    gsl_odeiv_control_free(control_ptr);
    gsl_odeiv_step_free(step_ptr);
  }
}

// Initialization of the class
//**************

// set sigma
void model::set_sigma(double* mysigma) 
{
  for (int n = 0; n < numparams; n++) 
  {
    sigma[n] = mysigma[n];
  }
}

// setptypes
void model::set_ptypes(double* myptypes) 
{
  for (int n = 0; n < numparams; n++) {
    ptypes[n] = myptypes[n];
  }
}

// setparams
void model::set_params(double* myparams) 
{
  for (int n = 0; n < numparams; n++) 
  {
    params[n] = myparams[n];
  }
}

// fitting posterior distribution to update the distribution parameters, such as mean and variance for log-normal, alpha and beta for beta distribution, mu and kappa for vonmises  
// and updating the scaling parameter sigma
// transition or preparation for the main code (model jump and mcmc)
//*****************

void model::fit_postprop(int numburn) 
{
  double par1,par2;
  
  for (int n = 0; n < numparams; n++) {
    par1 = 0;
    par2 = 0;
    
    if      (ptypes[n] == 1) {
      fit_lognormal(n,numburn,par1,par2);
    } 
    else if (ptypes[n] == 2) {
      fit_beta(n,numburn,par1,par2);
    }
    else if (ptypes[n] == 3) {
      fit_vonmises(n,numburn,par1,par2);
    }
    postpropparams[n][0] = par1;
    postpropparams[n][1] = par2;
    
    //cerr << n << " " << ptypes[n] << " " << par1 << " " << par2 << endl;
  }
}

void model::fit_lognormal(int n, int numburn, double &par1, double &par2) 
{
  int numiter = mcmc_logprob.size(); 
  int numnorm = numiter-numburn;
  double l = 0;
  
  for (int i = numburn; i < numiter; i++) {
    l = log(mcmc_params[n][i]);
    par1 += l;
    par2 += l*l;
  }
  par1 /= numnorm;
  par2 /= numnorm;
  par2 = sqrt(par2-par1*par1);
}

void model::fit_beta(int n, int numburn, double &par1, double &par2) 
{
  int numiter = mcmc_logprob.size();
  int numnorm = numiter-numburn;
  double p = 0;
  double m = 0;
  double v = 0;
  
  if (numburn != -1) {
    for (int i = numburn; i < numiter; i++) {
      p =  (mcmc_params[n][i]-priorparams[n][0])/(priorparams[n][1]-priorparams[n][0]);
      m += p;
      v += p*p; 
    }
    m /= numnorm;
    v /= numnorm;
    v = v-m*m;
    par1 = (m*m - m*m*m - m*v)/v;
    par2 = (m*m + v - m)*(m-1)/v;
  }
  else {
    par1 = 1;
    par2 = 1;
  }
}

void model::fit_vonmises(int n, int numburn, double &par1, double &par2) 
{
  int numiter = mcmc_logprob.size(); 
  int numnorm = numiter-numburn;
  double c = 0;
  double s = 0;
  double a = 0;
  double r = 0;
  
  if (numburn != -1) {
    for (int i = numburn; i < numiter; i++) 
    {
      a =  mcmc_params[n][i];
      c += cos(a);
      s += sin(a);
    }
    c /= numnorm;
    s /= numnorm;
    r = sqrt(c*c+s*s);
    
    par1 = atan2(s,c);
    if      (par1 < 0)    {par1 += 2*pi;}
    else if (par1 > 2*pi) {par1 -= 2*pi;}
    
    if (r < 0.53) {
      par2 = 2*r + r*r*r + 5*r*r*r*r*r/6;
    }
    else if (r >= 0.53 && r < 0.85) {
      par2 = -0.4 + 1.39*r + 0.43/(1-r);
    }
    else {
      par2 = 1/(r*r*r - 4*r*r + 3*r);
    }
  }
  else {
    par1 = 0;
    par2 = 0;
  }
}

//resetting scaling parameter sigma
void model::set_postsigm() 
{
  double factor = numparsamp;
  double sum = 0;
  double pro = 0;
  
  for (int n = 0; n < numparams; n++) {
    if (ptypes[n] == 1) {
      sigma[n] = postpropparams[n][1]/sqrt(factor);
    }
    else if (ptypes[n] == 2) {
      sigma[n] = factor*(postpropparams[n][0]+postpropparams[n][1]);
    }
    else if (ptypes[n] == 3) {
      sigma[n] = factor*(postpropparams[n][1]);
    }
    //cerr << n << " " << ptypes[n] << " " << sigma[n] << endl; 
  }
}
// if controlmax=1 which is not in the current case
void model::set_sigmaupdate(double factor) 
{
  for (int n = 0; n < numparams; n++) {
    sigma[n] = factor*sigma[n];
  }  
}

// MCMC shared in the transition and the main code
//*******************************

// LIKELIHOOD PATR
///////////

//likelihood
double model::get_logprob() 
{
  double x2,m2,xm;
  double logprob = 0;
  double num = 0;
  double b = 0;
  //double b = 1;
	double Sigma_m = 0.18874;
	double Sigma_s = 0.25514;
  
  get_solutions();// get the solution with the actual parameters
  //cout<<"<<<<<<<<<<"<<endl;
  for (int n = 0; n < numdata; n++) 
  {
    get_sum(n,x2,m2,xm);// calculate the terms of data and solution
    
    num += data[n].size()/2;
    if(n==0)
    {
    	b += (x2+m2-2*xm)/(Sigma_m*Sigma_m); // error function for exon
    }
    else
    {
    	b += (x2+m2-2*xm)/(Sigma_s*Sigma_s); // error function for intron
    }
    
    
    //num = data[n].size()/2;  
    //b *= x2+m2-2*xm;
  }
  logprob = -1.0/2.0*b;
  //cout<<numdata<<endl;
  //cout<<data[1].size()<<" "<<num<<endl;
  logprob += get_logprior();// prior probability
	//cout<<logprob<<" "<<get_logprior()<<" "<<-0.5*b<<endl;
  return(logprob);
}

void model::get_sum(int n, double &x2, double &m2, double &xm) 
{
  x2 = 0;
  m2 = 0;
  xm = 0;
  for (int i = 0; i < data[n].size(); i++) 
  {
    x2 += data[n][i]*data[n][i];
    //cout<<exp(data[n][i])<<" ";
    m2 += solution[n][i]*solution[n][i];
    xm += data[n][i]*solution[n][i];
  }
  //cout<<endl;
}

// prior probability  
double model::get_logprior() 
{
  double logprior = 0;
  for (int n = 0; n < numparams; n++) {
    if  (ptypes[n] == 1) 
    {
      logprior += log(general_log_normal_pdf(newparams[n],priorparams[n][2],priorparams[n][3]));
    }
    else if (ptypes[n] == 2) 
    {
      logprior += log(general_beta_pdf(newparams[n],priorparams[n][0],priorparams[n][1],priorparams[n][2],priorparams[n][3]));
    }
    else if (ptypes[n] == 3) 
    {
      logprior += log(general_von_mises_pdf(newparams[n],priorparams[n][2],priorparams[n][3]));
    }
  }
  return(logprior);
}

// getsolution
void model::get_solutions() 
{
  
  // Filling synthesis and degradation tables
  laurafill();
  
  //for (int i = 0; i < funtimes.size(); i++) {
  //cerr << funtimes[i] << " " << syntab[i] << " " << degtab[i] << endl;
  //}

  //for (double t = 0; t < 24; t+=0.01) {
  //cout << t << " " << laurafun(t,syntab) << " " << laurafun(t,degtab) << endl;
  //}
  //exit(1);
  
  // Getting solutions
  for (int n = 0; n < numdata; n++) 
  {
    if (datatype[n] == 1) 
    {
      get_accumulation(times[n],solution[n]);
    }
    else if (datatype[n] == 2) 
    {
      get_synthesis(times[n],solution[n]);
    }
    else if (datatype[n] == 3) 
    {
      get_decay(times[n],solution[n],t0index[n]);
    }
  }
}

void model::laurafill()
{
  double k0 = newparams[0];
  double ka = newparams[1];
  double kt = newparams[2];
  double kf = newparams[3];
  double kp = newparams[4];
  double g0 = newparams[5];
  double ga = newparams[6];
  double gt = newparams[7];
  double gf = newparams[8];
  double gp = newparams[9];
  
  lauratab(syntab,k0,ka,kt,kf,kp,0);
  lauratab(degtab,g0,ga,gt,gf,gp,1);
}

// rescaling parameters and compute the synthesis and degradation with them 
void model::lauratab(vector <double> &table, double a, double b, double T, double f, double p, int control) 
{
  int numfuntimes = funtimes.size();
  double fun = 0;
  double fun2 = 0.0;
  double up = 0;
  double dn = 0;
  double t  = 0;
  double l = (24*T-4)*f+2;
  double m = (24*T-4)*(1-f)+2;
  double maxF = 0;
  if      (control == 0) {maxF=230;}
  else if (control == 1) {maxF=70;}
                
  p = 3.819718634205488*p;//p = 24*p/(2*pi);
  a = (4.14-0.014)*a+0.014;  
  b = 3.111904508399025*maxF*b;
  //b = -2.07460300559935+1.037301502799675*maxF*b;
  //b = 1000*b;

  for (int n = 0; n < numfuntimes; n++) 
  {
    fun = 0;
    fun2 = 0.0;
    t = funtimes[n];
    for (int i = -2; i <= 2; i++) 
    {
      up = 1/(1+exp(-(8/l)*(l/2-p+24*i+t)));
      dn = 1/(1+exp((8/m)*(-m/2-p+24*i+t)));
      fun += a*(1 + b*(up+dn-1));
      fun2 += (up+dn-1);
    }
    //table[n] = fun/5;
    table[n] = a*(1+b*fun2);
    //cerr<<"diff of two funs "<<fun/5<<" "<<a*(1+b*fun2)<<endl;
    
  }
}


// getsolution including the integration diffequqtion to get mrna accumulation profile
void model::get_accumulation(vector <double> &mytimes, vector <double> &mysolution)
{
  int numtimes = mytimes.size();
  int numperiods = 1;
  double h = 1e-6;
  double t_next;
  double t = mytimes[0];
  double y[1];
  meanclock = 0;

  // step to tmax from tmin 
  get_initvalue(y[0],numperiods);
  for (int n = 0; n < numperiods; n++) 
  {
    meanclock = 0;
    for (int i = 0; i < numtimes; i++) 
    {
      t_next = mytimes[i]+n*48;//WRONG!!
      //t_next = mytimes[i]+n*24;//Wrong!

      while (t < t_next) 
      {
		gsl_odeiv_evolve_apply(evolve_ptr, control_ptr, step_ptr, &my_system, &t, t_next, &h, y);
      }
      mysolution[i] = y[0];
      meanclock += y[0];
    }
  }
  
  meanclock /= numtimes;
  for (int i = 0; i < numtimes; i++) {
    //if (mysolution[i] < 0) {cerr << "HEY!\n";}
    mysolution[i] = log(mysolution[i])-log(meanclock);
    //if (modeltype == 2) {cout << mytimes[i] << " " << mysolution[i] << endl;}
  }
  //if (modeltype == 2) {exit(1);}
}

// get initial conditions
void model::get_initvalue(double &y, int &numperiods) 
{
  numperiods = 2;
  /*
  if (ptypes[9] == 0) {
    y = newparams[2]/newparams[6];
    numperiods = 2;
  }
  else {
    y = newparams[9];
    numperiods = 1;
  }
  */
}

void model::get_synthesis(vector <double> &mytimes, vector <double> &mysolution) 
{
  int numtimes = mytimes.size();
  double mean = 0; 

  for (int i = 0; i < numtimes; i++) 
  {
    mysolution[i] = laurafun(mytimes[i],syntab);
    mean += mysolution[i];
  }
  mean /= numtimes;

  for (int i = 0; i < numtimes; i++) 
  {
    //if (mysolution[i] < 0) {cerr << "UMM!\n";}
    mysolution[i] = log(mysolution[i])-log(mean);
  }
}

void model::get_decay(vector <double> &mytimes, vector <double> &mysolution, int i0) 
{
  // to fill 
}

double model::laurafun(double t, vector <double> &table) 
{
  //cerr << "HEY!" << endl;
  int numfuntimes = funtimes.size();
  double fun = 0;
  int n,n1;
  t = mod(t,24); 
  n = int(t/dt);
  n1 = n+1;
  
  if (n1 == numfuntimes) {n1=0;}
  fun = table[n]+(t-funtimes[n])*(table[n1]-table[n])/(funtimes[n1]-funtimes[n]); 
  
  return(fun);
}

// diffequation
int model::diffequation(double t, const double y[], double f[], void *params_ptr) 
{
  double *myparams = (double *) params_ptr;
    
  /* evaluate the right-hand-side functions at t */
  f[0] = laurafun(t,syntab) - laurafun(t,degtab)*y[0]; 

  return GSL_SUCCESS;		/* GSL_SUCCESS defined in gsl/errno.h as 0 */
}

// integer to string
string model::itos(int i) 
{
  stringstream ss;
  ss << i;
  return(ss.str());
}

double model::mod(double a, double b) 
{
  double m = a-b*int(a/b);
  return(m);
}


// NEW PARAMETER PART 
////////////////////
double* model::get_newparams(int controlprop) 
{
  double pnew = 0;
  double pold = 0;
  double sigm = 0;
  double aold = 0;
  double bold = 0;

  if (controlprop == 0) 
  {
    for (int n = 0; n < numparams; n++) 
    {
      pold = params[n];
      sigm = sigma[n];
      
      // lognormal  
      if      (ptypes[n] == 1) 
      {
		pnew = log_normal_sample(log(pold),sigm,&seed);
      }
    
      // beta 
      else if (ptypes[n] == 2) 
      {
		get_betaparams(n,pold,sigm,aold,bold);
		pnew = general_beta_sample2(priorparams[n][0],priorparams[n][1],aold,bold,&seed);
      }
      
      // van mises 
      else if (ptypes[n] == 3) 
      {
		pnew = general_von_mises_sample(params[n],sigma[n],&seed);
      }
      
      else 
      {
		pnew = params[n];
      }
      // seting newparams
      newparams[n] = pnew;
    }
  }
  else 
  {
    for (int n = 0; n < numparams; n++) {

      // lognormal  
      if      (ptypes[n] == 1) {
	pnew = log_normal_sample(postpropparams[n][0],postpropparams[n][1],&seed);
      }
      
      // beta 
      else if (ptypes[n] == 2) {
	//cerr << modeltype << " " << priorparams[n][0] << " " << priorparams[n][1] << " " << postpropparams[n][0] << " " << postpropparams[n][1] << endl;
	pnew = general_beta_sample2(priorparams[n][0],priorparams[n][1],postpropparams[n][0],postpropparams[n][1],&seed);
      }
      
      // van mises 
      else if (ptypes[n] == 3) {
	pnew = general_von_mises_sample(postpropparams[n][0],postpropparams[n][1],&seed);
      }
      
      else {
	pnew = params[n];
      }
      //cerr << n << " " << ptypes[n] << " " << pnew << endl;

      // seting newparams
      newparams[n] = pnew;
    }
  }
  
  return(newparams);
}

void model::get_betaparams(int n, double p, double &s, double &a, double &b) 
{
  double len = (priorparams[n][1]-priorparams[n][0]);
  p = (p-priorparams[n][0])/len;
  a = s*p+1;
  b = s*(1-p)+1;
}

//ACCEPTANCE RATIO PART
////////////////////

// propratio no jump
double model::get_propratio(int controlprop) 
{
  double proppdf = 1;
  double pold = 0;
  double pnew = 0;
  double sigm = 0;
  double aold = 0;
  double bold = 0;
  double anew = 0;
  double bnew = 0;

  if (controlprop == 0) {
    for (int n = 0; n < numparams; n++) {
      pold = params[n]; 
      pnew = newparams[n];
      sigm = sigma[n];
      
      // lognormal  
      if      (ptypes[n] == 1) {
	proppdf *= (log_normal_pdf(pold,log(pnew),sigm)/log_normal_pdf(pnew,log(pold),sigm));
      }
      
      // beta 
      else if (ptypes[n] == 2) {
	get_betaparams(n,pold,sigm,aold,bold);
	get_betaparams(n,pnew,sigm,anew,bnew);
	proppdf *= (general_beta_pdf(pold,priorparams[n][0],priorparams[n][1],anew,bnew)/general_beta_pdf(pnew,priorparams[n][0],priorparams[n][1],aold,bold));
      }
      
      // van mises is simetric 
      //else if (ptypes[n] == 3) {
      //proppdf *= (von_mises_pdf(pold,pnew,sigm)/von_mises_pdf(pnew,pold,sigm));
      //}
    }
  }
  // WRONG? NOT USED HERE
  else { 
    for (int n = 0; n < numparams; n++) 
    {
      pold = params[n]; 
      pnew = newparams[n];
      
      // lognormal  
      if      (ptypes[n] == 1) {
	proppdf *= (log_normal_pdf(pold,postpropparams[n][0],postpropparams[n][1])/log_normal_pdf(pnew,postpropparams[n][0],postpropparams[n][1]));
      }
      
      // beta 
      else if (ptypes[n] == 2) {
	proppdf *= (general_beta_pdf(pold,priorparams[n][0],priorparams[n][1],postpropparams[n][0],postpropparams[n][1])/general_beta_pdf(pnew,priorparams[n][0],priorparams[n][1],postpropparams[n][0],postpropparams[n][1]));
      }
      
      // van mises is simetric 
      else if (ptypes[n] == 3) {
	proppdf *= (von_mises_pdf(pold,postpropparams[n][0],postpropparams[n][1])/von_mises_pdf(pnew,postpropparams[n][0],postpropparams[n][1]));
      }
      //cerr << "RATIO: " <<  n << " " << ptypes[n] << " " <<  proppdf << endl;
    }
  }
  
  return(proppdf);
}

// model jump
double model::get_postproppdf(int controlprop) 
{
  double proppdf = 1;
  double para = 0;
  
  for (int n = 0; n < numparams; n++) {
    if (controlprop == 0) {para = params[n];}
    else                  {para = newparams[n];}
    
    // lognormal  
    if        (ptypes[n] == 1) {
      proppdf *= log_normal_pdf(para,postpropparams[n][0],postpropparams[n][1]);
    }
      
    // beta 
    else if (ptypes[n] == 2) {
      proppdf *= general_beta_pdf(para,priorparams[n][0],priorparams[n][1],postpropparams[n][0],postpropparams[n][1]);
    }
    
    // van mises is simetric 
    else if (ptypes[n] == 3) {
      proppdf *= general_von_mises_pdf(para,postpropparams[n][0],postpropparams[n][1]);
    }
    //cerr <<  modeltype << " " << n << " " << ptypes[n] << " " << para << " " << proppdf << endl;
  }
  
  return(proppdf);
}
  
// set new params if accepted
void model::set_newparams(double* myparams) 
{
  for (int n = 0; n < numparams; n++) {
    newparams[n] = myparams[n];
  }
}

// save current parameters 
void model::set_mcmcparams(double logprob) 
{
  mcmc_logprob.push_back(logprob);
  for (int n = 0; n < numparams; n++) 
  {
    mcmc_params[n].push_back(params[n]);
  }
}

// compute average trajectory
void model::get_mcmc_path() 
{
  int numtimes = 0;
  double val = 0;
  
  // Wrong! It should not count if iteration < numburn
  for (int n = 0; n < numdata; n++) 
  {
    for (int i = 0; i < times[n].size(); i++) 
    {
      val = solution[n][i];
      mcmc_path[n][i] += val;
      mcmc_path2[n][i] += val*val;
    }
  }
}


//OUTPUT PATRT
/////////////

// computing mean and std 
void model::get_mcmc_mean(int numburn) 
{
  int numiter = mcmc_logprob.size(); 
  double val = 0;
  double m  = 0;
  double m2 = 0;
  double s  = 0;
  
  for (int n = 0; n < numparams; n++) 
  {
    m  = 0;
    m2 = 0;
       
    for (int i = numburn; i < numiter; i++) 
    {
      val = mcmc_params[n][i];
      m  += val;
      m2 += val*val;
    }
    m /= (numiter-numburn);
    m2 /= (numiter-numburn);
    s = sqrt(m2-m*m);
    
    mcmc_mean[n] = m;
    mcmc_sigm[n] = s;
  }
}

// comparing function
bool compangles(double a,double b) 
{
  return(cos(a) > cos(b));
}

// computing mean and std 
void model::get_mcmc_median(int numburn) 
{
  double numiter = (double) (mcmc_logprob.size()-numburn);
  int p05 = int(0.025*numiter);
  int p50 = int(0.500*numiter);
  int p95 = int(0.975*numiter);
  
  for (int n = 0; n < numparams; n++) {
    
    // Sorting 
    if (ptypes[n] == 3) {
      sort(mcmc_params[n].begin()+numburn,mcmc_params[n].end(),compangles);
      //cerr<<" compangles "<<compangles<<endl;
    }
    else {
      sort(mcmc_params[n].begin()+numburn,mcmc_params[n].end());
    }
    
    // Storing median and 95% interval
    mcmc_median[n] = mcmc_params[n][p50+numburn];
    mcmc_quan05[n] = mcmc_params[n][p05+numburn];
    mcmc_quan95[n] = mcmc_params[n][p95+numburn];
    //cerr << n << " " << ptypes[n] << " " << p05 << " " << p50 << " " << p95 << " " << mcmc_median[n] << " " << mcmc_quan05[n] << " " << mcmc_quan95[n] << endl;
  }
}

// histogram
void model::get_mcmc_histograms(int numburn) 
{
  int numiter = mcmc_logprob.size();
  int k = 0;
  
  for (int n = 0; n < numparams; n++) 
  {
    // Initialitation
    for (k = 0; k < numbins; k++) {mcmc_hist[n][k] = 0;}
    k = 0;
    
    // Sorting 
    sort(mcmc_params[n].begin()+numburn,mcmc_params[n].end());
    
    // Computing counts in each bin
    hist_dx[n] = (mcmc_params[n][numiter-1]-mcmc_params[n][numburn])/numbins;
    hist_x0[n] = mcmc_params[n][numburn];
    for (int i = numburn; i < numiter; i++) {
      while (mcmc_params[n][i] > (k+1)*hist_dx[n]+hist_x0[n]) {k++;}
      mcmc_hist[n][k]++;
    }
  }
}

// get maximum
void model::get_mcmc_maxpost() {
  for (int n = 0; n < numparams; n++) {
    mcmc_maxpost[n] = 0;
    mcmc_maxpara[n] = 0;
    for (int k = 0; k < numbins; k++) {
      if (mcmc_maxpost[n] < mcmc_hist[n][k]) {
	mcmc_maxpost[n] = mcmc_hist[n][k];
	mcmc_maxpara[n] = (k+1/2)*hist_dx[n]+hist_x0[n];
      }
    }
  }
}

// get covariance matrix of transformed variables
void model::get_covmat(int numburn) 
{
  int numiter = mcmc_logprob.size();
  int n,m,i;
  double *mcmc_transparams[numparsamp];
  double c = 0;
  
  // Transforming variables
  m = 0;
  for (n = 0; n < numparams; n++) {

    if (ptypes[n] == 1) {
      mcmc_transparams[m] = new double[numiter];
      for (i = 0; i < numiter; i++) {
	//mcmc_transparams[m][i] = log(mcmc_params[n][i]);
	mcmc_transparams[m][i] = mcmc_params[n][i];
      }
      m++;
    }
    else if (ptypes[n] == 2) {
      mcmc_transparams[m] = new double[numiter];
      for (i = 0; i < numiter; i++) {
	//mcmc_transparams[m][i] = log(mcmc_params[n][i]/(1-mcmc_params[n][i]));
	mcmc_transparams[m][i] = mcmc_params[n][i];
      }
      m++;
    }
    else if (ptypes[n] == 3) {
      mcmc_transparams[m] = new double[numiter];
      for (i = 0; i < numiter; i++) {
	c = cos(mcmc_params[n][i]);
	//mcmc_transparams[m][i] = log((1+c)/(1-c));
	mcmc_transparams[m][i] = c;
      }
      m++;
    }
  }
  
  // computing cov matrix
  double covtmp = 0;
  double cov[numparsamp][numparsamp];
  double mu[numparsamp];
  for (n = 0; n < numparsamp; n++) {
    mu[n] = 0;
  }
  for (n = 0; n < numparsamp; n++) {
    mu[n] += mcmc_transparams[n][i];
    for (m = 0; m < numparsamp; m++) {
      covtmp = 0;
      for (i = numburn; i < numiter; i++) {
	covtmp += mcmc_transparams[n][i]*mcmc_transparams[m][i]; 
      }
      mu[n] /= numiter;
      cov[n][m] = covtmp/numiter;
    }
  }
  
  for (n = 0; n < numparsamp; n++) {
    for (m = 0; m < numparsamp; m++) {
      cov[n][m] = cov[n][m]-mu[n]*mu[m];
    }
  }
  for (n = 0; n < numparsamp; n++) {
    for (m = 0; m < numparsamp; m++) {
      fprintf(stderr,"%.7f\t",(cov[n][m]*cov[n][m])/(cov[n][n]*cov[m][m]));
    }
    fprintf(stderr,"\n");
  }
}

// get functions of parameters NOT GENERAL! Wrong!
void model::get_mcmc_functions(int numburn) 
{
  // Setting params
  set_params(mcmc_mean);

  // Getting theoretical curves 
  get_accumulation(theo_times,theo_mrna);
  get_synthesis(theo_times,theo_synth);
  double maxmrna = theo_mrna[0];
  double minmrna = theo_mrna[0];
  double maxsynth = theo_synth[0];
  double minsynth = theo_synth[0];
  double phasemrna = 0;
  double phasesynth = 0;

  // Phase expression
  for (int n = 0; n < theo_times.size(); n++) 
  {
    if (maxmrna  < theo_mrna[n]) 
    {
      phasemrna = theo_times[n];
      maxmrna = theo_mrna[n];
    }
    if (minmrna > theo_mrna[n]) 
    {
      minmrna = theo_mrna[n];
    }
  }
  mcmc_func[0] = phasemrna;
  mcmc_func[1] = exp(maxmrna)-exp(minmrna);

  // Phase synthesis
  for (int n = 0; n < theo_times.size(); n++) 
  {
    if (maxsynth  < theo_synth[n]) 
    {
      phasesynth = theo_times[n];
      maxsynth = theo_synth[n];
    }
    if (minsynth > theo_synth[n]) 
    {
      minsynth = theo_synth[n];
    }
  }
  mcmc_func[2] = phasesynth;
  mcmc_func[3] = exp(maxsynth)-exp(minsynth);

  // Mean log
  double m = 0;
  double m2 = 0;
  double val = 0;
  int numiter = mcmc_logprob.size();
  for (int i = numburn; i < numiter; i++) 
  {
    val = log(mcmc_params[5][i]);
    //cerr<<mcmc_params[0][i]<<" "<<mcmc_params[1][i]<<" "<<mcmc_params[2][i]<<" "<<mcmc_params[3][i]<<" "<<mcmc_params[4][i]<<" "<<mcmc_params[5][i]<<endl;
    //cerr<<log(mcmc_params[5][i])<<endl;
    m  += val;
    m2 += val*val;
  }	
	m /= (numiter-numburn);
  	m2 /= (numiter-numburn);
	//cerr<<m<< "  "<<m2<<endl;
  	mcmc_func[4] = m;
  	mcmc_func[5] = m2;
}

// print resutls3
void model::print_results3(int nummodels, double modelprob[]) 
{

  // printing phases
  fprintf(stdout,"%.1f\t%.3f\t%.1f\t%.3f\t%.7f\t%.7f\t",mcmc_func[0],mcmc_func[1],mcmc_func[2],mcmc_func[3],mcmc_func[4],mcmc_func[5]);
  
  // printing model probabilities
  fprintf(stdout,"%.5f\t%.5f\t%.5f\t%.5f\t%i\t",modelprob[0],modelprob[1],modelprob[2],modelprob[3],modeltype);

  // printing synthesis and degradation paramters
  //fprintf(stdout,"%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t",mcmc_mean[1],mcmc_mean[2],mcmc_mean[3],mcmc_mean[4],mcmc_mean[5],mcmc_mean[6],mcmc_mean[9]);
  fprintf(stdout,"%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t",mcmc_maxpara[1],mcmc_maxpara[2],mcmc_maxpara[3],mcmc_maxpara[4],mcmc_maxpara[5],mcmc_maxpara[6],mcmc_maxpara[9]);

  // printing errors
  fprintf(stdout,"%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t",mcmc_sigm[1],mcmc_sigm[2],mcmc_sigm[3],mcmc_sigm[4],mcmc_sigm[5],mcmc_sigm[6],mcmc_sigm[9]);

  // Name data  
  fprintf(stdout,"%s\n",namedata);
}

// print resutls2
void model::print_results2(int nummodels, double modelprob[]) 
{

  // printing phases
  fprintf(stdout,"%.0f\t%.1f\t%.0f\t%.1f\t",1.0,mcmc_func[2],1.0,mcmc_func[3]);
  
  // printing model probabilities
  fprintf(stdout,"%.5f\t%.5f\t%.5f\t%.5f\t%i\t",modelprob[0],modelprob[1],modelprob[2],modelprob[3],modeltype);

  // printing synthesis and degradation paramters
  fprintf(stdout,"%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.0f\t",mcmc_mean[2],mcmc_mean[3],mcmc_mean[4],mcmc_mean[5],mcmc_mean[6],mcmc_mean[7],mcmc_mean[8],0.0);

  // printing errors
  fprintf(stdout,"%.0f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.0f\t",0.0,mcmc_sigm[2],mcmc_sigm[3],mcmc_sigm[4],mcmc_sigm[5],mcmc_sigm[6],mcmc_sigm[7],mcmc_sigm[8],0.0);

  // Name data  
  fprintf(stdout,"%s\n",namedata);
}

// print results
void model::print_results(int nummodels, double modelprob[]) 
{
  cout << mcmc_func[0] << "\t" << mcmc_func[1] << "\t" << mcmc_func[2] << "\t" << mcmc_func[3] << "\t";

  for (int n = 0; n < nummodels; n++) {
    cout << modelprob[n] << "\t";
  }
  cout << modeltype << "\t";

  for (int n = 0; n < numparams; n++) {
    cout << mcmc_maxpara[n] << "\t"; 
  }
  cout << modeltype << "\t";

  for (int n = 0; n < numparams; n++) {
    cout << mcmc_median[n] << "\t"; 
  }
  cout << modeltype << "\t";

  for (int n = 0; n < numparams; n++) {
    cout << mcmc_mean[n] << "\t"; 
  }
  cout << endl;
}

void model::print_hist(string fileout) 
{
  ofstream out(fileout.c_str());  
  for (int k = 0; k < numbins; k++) {
    for (int n = 0; n < numparams; n++) {
      if (ptypes[n] != 0) {
	out << (k+1/2)*hist_dx[n]+hist_x0[n] << " " << mcmc_hist[n][k] << " ";
      }
    }
    out << endl;
  }
  out.close();
}

// Print average path
void model::print_path(string filetag) 
{
  int numiter = mcmc_logprob.size();
  double val = 0;
  double val2 = 0;
  double sigm = 0;
  string fittag  = "fit.";
  string outtag  = ".out";
  string fileout = "";

  for (int n = 0; n < numdata; n++) {
    fileout = string(filetag+fittag+itos(n)+outtag);
    
    ofstream out(fileout.c_str());
    for (int i = 0; i < times[n].size(); i++) {
      val  = mcmc_path[n][i]/numiter; //Wrong!
      val2 = mcmc_path2[n][i]/numiter;//Wrong!
      sigm = sqrt(val2-val*val);
      out << times[n][i] << "\t" << data[n][i] << "\t" << val << "\t" << sigm << "\n";
    }
    out.close();
  }
}

// print postparams
void model::print_postprop(string fileout) 
{
  ofstream out(fileout.c_str());
  for (int n = 0; n < numparams; n++) {
    out << n << " " << ptypes[n] << " " << postpropparams[n][0] << " " << postpropparams[n][1] << endl;
  }
  out.close();
}

// print solution: Wrong!
void model::print_fit(string fileout) {
  ofstream out(fileout.c_str());
  for (int i = 0; i < times[0].size(); i++) {
    for (int n = 0; n < numdata; n++) {
      out << times[n][i] << " " << data[n][i] << " " << solution[n][i] << " ";
    }
    out << "\n";
  }
  out.close();
}

// print mcmc 
void model::print_mcmc(string fileout) {
  int numiter = mcmc_logprob.size(); 
  
  ofstream out(fileout.c_str());
  for (int i = 0; i < numiter; i++) {
    out << i << " " << mcmc_logprob[i] << " ";
    for (int n = 0; n < numparams; n++) {
      out << mcmc_params[n][i] << " ";
    }
    out << endl;
  }
  out.close();
}


///////Blablabla PART not used for the moment

/*
double model::laurafun(double t, double a, double b, double T, double f, double p) 
{
  double fun = 0;
  double up = 0;
  double dn = 0;
  double l = (24*T-4)*f+2;
  double m = (24*T-4)*(1-f)+2;
  t = t-24*int(t/24);
  p = 3.819718634205488*p;//p = 24*p/(2*pi);
  a = (4.14-0.014)*a+0.014;
  
  for (int i = -1; i <= 1; i++) {
    up = 1/(1+exp(-(8/l)*(l/2-p+24*i+t)));
    dn = 1/(1+exp((8/m)*(-m/2-p+24*i+t)));
    fun += a*(1 + b*(up+dn-1));
  }

  return(fun/3);
}
*/
