##########################################################################
##########################################################################
## Project:
## Script purpose: Set the boundary and initialize parameters for each model
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Jun  4 12:23:08 2018
##########################################################################
##########################################################################
Gamma.Initiation = function(eps.gamma.init, min.half.life=0.5, max.half.life=6, w=2*pi/24)
{
  # min.half.life=0.5; max.half.life=6; w=2*pi/24
  eps.gamma.init = as.numeric(eps.gamma.init);
  gamma.init = rep(log(2)/(10/60), length(eps.gamma.init))
  lss = lseq(log(2)/max.half.life, log(2)/min.half.life, length=100)
  for(n in 1:length(gamma.init))
  {
    kk = which(lss>=sqrt(w^2/(1/eps.gamma.init[n]^2-1)));
    if(length(kk)>0) gamma.init[n] = sample(lss[kk], 1);
  }
  return(gamma.init)
  #gamma.init = c(rep(log(2)/lseq(max.half.life, min.half.life, length = Nfit.M), Nfit.M%/%Nfit.M), rep(log(2)/5,Nfit.M%%Nfit.M))
}

sigmoid.bound.contraint = function(eps.gamma)
{
  return(10^4/(1+exp(-100*(eps.gamma-1.0))));
  
  smooth.bound.constraint.function = FALSE
  if(smooth.bound.constraint.function)
  {
    xx = lseq(0.001, 1.2, length.out = 1000)
    x0 = 1.0;y0=10^4
    yy = y0/(1+exp(-100*(xx-x0)))
    plot(xx, yy, type='l', col='blue', log='', ylim=c(0, y0));abline(h=1, col='red');abline(h=0, col='red');
    abline(v=x0, col='black');abline(v=(x0-0.05), col='black');abline(v=1, col='black');
  }
}

####################
## main function for parameter boundaries 
####################
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
    
  upper = c(gamma.max, eps.gamma.max, phase.gamma.max, splicing.k.max,  param.synthesis.upper)
  lower = c(gamma.min, eps.gamma.min, phase.gamma.min, splicing.k.min,  param.synthesis.lower)
  
  names(upper) = c("gamma", "eps.gamma", "phase.gamma", "splicing.k", 
                           "Min.int", "Amp.int", "phase.int", "beta.int")
  names(lower) =  names(upper)
  
  if(absolute.signal){
    if(model==1){ upper = upper[c(1, 4, 5)];  lower = lower[c(1, 4, 5)]; }
    if(model==2){ upper = upper[c(1, 4, 5:8)];  lower = lower[c(1, 4, 5:8)];}
    if(model==3){ upper = upper[c(1:5)]; lower = lower[c(1:5)]; }
  }
  
  return(list(lower = lower, upper = upper))
}


