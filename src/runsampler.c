#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <algorithm>
#include "model.h"

using namespace std;

#define RNG ((random()+0.1)/2147483648.0)//RAND_MAX
string itos(int i);
void read_models(char *namedata, char *file, vector <model*> &modelset, vector <int> &datatype, vector <int> &t0index, vector <vector <double> > &times, vector <vector <double> > &data);
void read_data(char *file, int &numdata, vector <int> &datatype, vector <int> &t0index, vector <vector <double> > &times, vector <vector <double> > &data);
void read_pmatrix(char *file, vector <vector <double> > &probjump);
void set_beta(double beta[], int numiter, int numburn, int controlmax);
int  modelrnd(int &oldmodel, vector <vector <double> > &probjump, int &controljump);

// ====================================================
//                      MAIN
// ====================================================
int main(int argc, char* argv[11]) 
{ 
  //cout << "Jesus is my Lord"<<endl;
  //abort();
      
  // Read Arguments
  // ==============================================================
  if(argc != 11) {cerr << "Usage: <namedata> <filedata> <filemodels> <fileprobmatrix> <fileout> <numiterpost> <numiteradap> <numiterjmp> <numburn> <controlmax>\n";exit(1);} 
  char* namedata  = argv[1];
  char* filedata  = argv[2];
  char* filemods  = argv[3];
  char* filepmat  = argv[4];
  char* filetag   = argv[5];
  int numiterpost = atoi(argv[6]);
  int numiteradap = atoi(argv[7]);
  int numiterjmp  = atoi(argv[8]);
  int numburn     = atoi(argv[9]);
  int controlmax  = atoi(argv[10]);
  
  // Time and Seed
  time_t start, finish;
  double deltatime;
  start = time(NULL);
  srandom(start);
  
  // Random generator
  const gsl_rng_type * R;
  gsl_rng * r;
  gsl_rng_env_setup();
  R = gsl_rng_default;
  r = gsl_rng_alloc(R);
  gsl_rng_set(r,start);
  
  // Files
  string fileout(filetag);
  string outmcmc("mcmc.out");
  string outpost(".post.out");
  string outfit("fit.out"); 
  string outhist("hist.out");
  fileout = fileout + ".";


  // Initialitation
  // ==============================================================
  vector < vector <double> > probjump;
  vector < vector <double> > times;
  vector < vector <double> > data;
  vector <model*> modelset;
  vector <int> datatype;
  vector <int> t0index;
  int numdata = 0;

  // Reading data file
  read_data(filedata,numdata,datatype,t0index,times,data);    
  //cout<<"You are my Lord"<<endl;
  //abort();
  // Setting models
  read_models(namedata,filemods,modelset,datatype,t0index,times,data);

  // Reading prob matrix
  read_pmatrix(filepmat,probjump);


  // Running MCMC
  // ==============================================================
  int nummodels  = modelset.size();   
  int numiterall = numiterpost*nummodels;
  int numitertot = numiterall+numiterjmp;
  int numparams = 10; //wrong
  int controljump = 0;
  int newmodel = 0;
  int oldmodel = 0;
  int controla = 1;
  int acept = 0;
  int a = 0;
  double facnew = 1;
  double facold = 1; 
  double ratio = 0;
  double ratio1 = 0;
  double ratio2 = 0;
  double logprob = 0;
  double newlogprob = 0;
  double ratioprop = 0;
  double newproppdf = 0;
  double proppdf = 0;
  double *newparams = new double[numparams];
  double modelcounts[nummodels];
  for (int n = 0; n < nummodels; n++) {modelcounts[n]=0;}

  // Setting temperature function
  double beta[numitertot];
  set_beta(beta,numitertot,numburn,controlmax);  	
  //double time11 = time(NULL);
  // Running MCMC
  logprob = modelset[newmodel] -> get_logprob();
  
	for (int n = 0; n < numitertot; n++) {
     
      // transition before model jumping
      // fitting posterior distribution and update scaling parameter sigma
      //====================================
      if (controla == 1) {a++;}
      if (a == numiteradap && n < numiterall) {
          modelset[newmodel] -> fit_postprop(numburn);
          modelset[newmodel] -> set_postsigm();
          ratio1 = (double) acept/numiteradap;
          controla = 0;
          acept = 0;
          a = 0;
          //cout<<"HEY "<<" "<<n<<endl;
          /*
           string modeltag = itos(newmodel+1);
           string filepost = string(fileout + modeltag + outpost);
           modelset[newmodel] -> print_postprop(filepost);
           */
      }
      
      // Updating sigma if controlmax=1
      if (controlmax == 1) {
          if (n <= numburn) {facnew = 1;}
          else              {facnew = n/1000.0;}
          modelset[newmodel] -> set_sigmaupdate(1/facold);
          modelset[newmodel] -> set_sigmaupdate(facnew);
          facold = facnew;
      }
      
      if (newmodel < ((n+1)/numiterpost) && n < numiterall) {
          modelset[newmodel] -> fit_postprop(numiteradap);
          ratio2 = (double) acept/(numiterpost-numiteradap);
          cerr << "ACEPT: " << ratio1 << " " << ratio2 << endl;
          //cout<<"Hey, here !!"<<endl;
          //cout<<newmodel<<" "<<n<<endl;
          controla = 1;
          acept = 0;
          
          /*
           string modeltag = itos(newmodel+1);
           string filefit  = string(fileout + modeltag + outfit);
           string filemcmc = string(fileout + modeltag + outmcmc);
           string filepost = string(fileout + modeltag + outpost);
           
           modelset[newmodel] -> print_postprop(filepost);
           modelset[newmodel] -> print_fit(filefit);
           modelset[newmodel] -> print_mcmc(filemcmc);
           */
          
          if (newmodel < nummodels-1) {
              newmodel++;
              oldmodel = newmodel;
              logprob  = modelset[newmodel] -> get_logprob();
          }
      }
      
	//After transition 
    // Selecting model	
    if (n >= numiterall){
      newmodel = modelrnd(oldmodel,probjump,controljump);
      modelcounts[oldmodel]++;
    }
      
    // New params and prop if no model jump
    if (controljump == 0) 
    {
      newparams = modelset[newmodel] -> get_newparams(0);
      ratioprop = modelset[newmodel] -> get_propratio(0);
    }
      
    // New params and prop if model jump
    else {
      newparams  = modelset[newmodel] -> get_newparams(1);
      newproppdf = modelset[newmodel] -> get_postproppdf(1);
      proppdf    = modelset[oldmodel] -> get_postproppdf(0);
      ratioprop  = proppdf/newproppdf;
    }
    // New prob
    newlogprob = modelset[newmodel] -> get_logprob();

    // Ratio
    ratio = exp(beta[n]*(newlogprob-logprob))*ratioprop;
      
    // Aceptance
    if (ratio > RNG) {
      modelset[newmodel] -> set_params(newparams);
      logprob  = newlogprob;
      oldmodel = newmodel;
      acept += 1;
    }
    
    // Storing current params
    modelset[oldmodel] -> set_mcmcparams(logprob);

    // Computing average path
    modelset[oldmodel] -> get_mcmc_path();

    // Printing status
    //if ((n % (int) (numiterpost/10)) == 0) {cerr << n << " " << oldmodel+1 << " " << newmodel+1 << " " << logprob << endl;}
    //cerr << n << " " << oldmodel+1 << " " << newmodel+1 << " " << logprob << endl;
    
  }
  cerr << "ACEPT: " <<  (double) acept/numiterjmp << endl;
     
  // Processing MCMC results
  // ============================================================== 
  // Getting max model
  int maxmodel = 0;
  double maxprob = 0;
  for (int n = 0; n < nummodels; n++) {
    modelcounts[n] /= numiterjmp; 
    if (maxprob < modelcounts[n]) {
      maxprob = modelcounts[n];
      maxmodel = n;
    }
  }

  // Printing mcmc trajectories
  //modelset[maxmodel] -> print_mcmc(fileout+outmcmc);
  
  // Cmputing median and mean params
  if (numiteradap > numburn) {numburn=numiteradap;}
  modelset[maxmodel] -> get_mcmc_mean(numburn);
  modelset[maxmodel] -> get_mcmc_median(numburn);
  modelset[maxmodel] -> get_mcmc_histograms(numburn);
  modelset[maxmodel] -> get_mcmc_maxpost();
  //modelset[maxmodel] -> get_covmat(numburn);

  // Computing phases
  modelset[maxmodel] -> get_mcmc_functions(numburn);
  
  // Printing out results
  modelset[maxmodel] -> print_results3(nummodels,modelcounts);
  
  // Printing histograms
  modelset[maxmodel] -> print_hist(string(fileout+outhist));

  // Printing average path
  //modelset[maxmodel] -> print_path(string(fileout+outfit));
  modelset[maxmodel] -> print_path(fileout);
  
  
  for (int n = 0; n < nummodels; n++) {
    string dot(".");
    string modeltag = itos(n+1);
    string filemcmc = string(fileout + modeltag + dot + outmcmc);
    //string filefit  = string(fileout + modeltag + dot + outfit);
    //string filepost = string(fileout + modeltag + dot + outpost);
    //modelset[n] -> print_postprop(filepost);
    //modelset[n] -> print_fit(filefit);
    modelset[n] -> print_mcmc(filemcmc);
   }
    
  // Time!
  // ============================================================== 
  finish = time(NULL);
  deltatime = (finish - start)/60.0;
  cerr << "TOTAI TIME: " << deltatime << " min" << endl;
}


// ====================================================
//                      FUNCTIONS
// ====================================================

// Selecting model 
int modelrnd(int &oldmodel, vector <vector <double> > &probjump, int &controljump) {
  double rnum = RNG;
  int n = 0;
  
  while (rnum > probjump[oldmodel][n]) {
    n++;
  }

  if (n == oldmodel) {controljump = 0;}
  else               {controljump = 1;}
  
  return(n);
}


// Reading models 
void read_models(char *namedata, char *file, vector <model*> &modelset, vector <int> &datatype, vector <int> &t0index, vector <vector <double> > &times, vector <vector <double> > &data) {
  string modelclass;
  int  modeltype;
  int  numparams;
  int  tmp;

  ifstream in(file);
  if(!in) {
    cerr << "ERROR: file could not be opened" << endl;
    exit(1);
  }
  
  while (!in.eof()) {
    if (in >> modelclass >> modeltype >> numparams) {
      string *paraname  = new string[numparams];
      double *paratype  = new double[numparams];
      double *paravalue = new double[numparams];
      double *parasigma = new double[numparams];
      double **paraprior = new double*[numparams];
      
      for (int n = 0; n < numparams; n++) {
	paraprior[n] = new double[4];
	in >> paraname[n] >> paravalue[n] >> paratype[n] >> paraprior[n][0] >> paraprior[n][1] >> paraprior[n][2] >> paraprior[n][3] >> parasigma[n];
      }
      model *mymodel = new model(namedata,modeltype,paravalue,paratype,parasigma,paraprior,datatype,t0index,times,data);
      modelset.push_back(mymodel);
    }
  }
  in.close();
}

// Reading data file 
void read_data(char *file, int &numdata, vector <int> &datatype, vector <int> &t0index, vector <vector <double> > &times, vector <vector <double> > &data)
{
    ifstream in(file);
    string dataname;
    double tvalue = 0;
    double dvalue = 0;
    int timeindex = 0;
    int numpoints = 0;
    int dtype = 0;

  if(!in) {
    cerr << "ERROR: file could not be opened" << endl;
    exit(1);
  }

  while (!in.eof())
  {
    if (in >> dataname >> dtype >> numpoints >> timeindex)
    {
        datatype.push_back(dtype);
        t0index.push_back(timeindex);
        times.push_back(*(new vector <double>(numpoints)));
        data.push_back(*(new vector <double>(numpoints)));
        //cout<<"dataname  "<<dataname<<" "<<"dtype "<<dtype<<" "<<"numpoints"<<numpoints<<" "<<endl;
        for (int n = 0; n < numpoints; n++)
        {
            in >> tvalue >> dvalue;
            times[numdata][n] = tvalue;
            data[numdata][n] = log(dvalue);
        }
        numdata++;
    }
  }
in.close();
}

// Reading prob matrix
void read_pmatrix(char *file, vector <vector <double> > &probjump) {
  int nummodels;
  string word;
  double prob;
  double num;
  int n = 0;
  int m = 0;
  
  ifstream in(file);
  if(!in) {
    cerr << "ERROR: file could not be opened" << endl;
    exit(1);
  }

  in >> word >> nummodels;
  for (n = 0; n < nummodels; n++) {
    probjump.push_back(*(new vector <double>));
  }
  
  n = 0;
  while (n < nummodels) {
    prob = 0;
    for (m = 0; m < nummodels; m++) {
      in >> num;
      prob += num;
      probjump[n].push_back(prob);
    }
    n++;
  }
  in.close();
  
  //for (n = 0; n < nummodels; n++) {
  //for (m = 0; m < nummodels; m++) {
  //cerr << probjump[n][m] << " ";
  //}
  //cerr << endl;
  //}
}

// Seting temperature function
void set_beta(double beta[], int numiter, int numburn, int controlmax) {
  for (int n = 0; n < numiter; n++) {
    if (controlmax == 1 && n > numburn) {
      beta[n] = 1+0.02*(n-numburn);
    }
    else {
      beta[n] = 1;
    }
  }
}

// integer to string
string itos(int i) {
  stringstream ss;
  ss << i;
  return(ss.str());
}
