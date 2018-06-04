##########################################################################
##########################################################################
## Project:
## Script purpose: Set the boundary and initialize parameters for each model
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Jun  4 12:23:08 2018
##########################################################################
##########################################################################
set.bounds = function(model = 4, parametrization =c('cosine.beta'), absolute.signal = TRUE)
{
  ### minimum and maximum of gamma
  gamma.min = log(2)/24
  gamma.max = log(2)/(10/60);  
  
  ### relative amplitudes of gamma (normalized)
  eps.gamma.min = 0;
  eps.gamma.max = 1;
  
  ### phases of gamma (normalized)
  phase.gamma.min = 0;
  phase.gamma.max = 24; 
  
  ### ratio between splicing rate and degradation rate   
  splicing.k.max = 10^5 ## 1s of splicing time;
  splicing.k.min = 10^(-1) ## 30 min of splicing time;
  
  Min.int.max = 2^20; ## intron signal mean 2^9;
  Min.int.min = 2^(-20);
  
  Amp.int.min = 0;
  Amp.int.max = 200;
  
  phase.int.min = 0;
  phase.int.max = 24; 
  
  beta.int.min = 1;
  beta.int.max = 5;
  
  param.synthesis.upper = c(Min.int.max, Amp.int.max, phase.int.max, beta.int.max)
  param.synthesis.lower = c(Min.int.min, Amp.int.min, phase.int.min, beta.int.min)
  
  if(absolute.signal)
  {
    if(model==1)
    {
      upper = c(gamma.max, splicing.k.max, mean.int.max)
      lower = c(gamma.min, splicing.k.min, mean.int.min)
    }
    if(model==2)
    {
      upper = c(gamma.max, splicing.k.max, param.synthesis.upper)
      lower = c(gamma.min, splicing.k.min, param.synthesis.lower)
    }
    if(model==3)
    {
      upper = c(gamma.max, eps.gamma.max, phase.gamma.max, splicing.k.max, Min.int.max)
      lower = c(gamma.min, eps.gamma.min, phase.gamma.min, splicing.k.min, Min.int.min)
    }
    if(model==4)
    {
      upper = c(gamma.max, eps.gamma.max, phase.gamma.max, splicing.k.max,  param.synthesis.upper)
      lower = c(gamma.min, eps.gamma.min, phase.gamma.min, splicing.k.min,  param.synthesis.lower)
    }
  }
  return(list(lower = lower, upper = upper))
}


