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
## function for general parameter boundaries 
####################
set.general.bounds.int = function(lower.user = NULL, 
                                  upper.user = NULL,
                                  parametrization =c('cosine.beta'), 
                                  absolute.signal = TRUE)
{
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
  names(param.synthesis.upper) = c("Min.int", "Amp.int", "phase.int", "beta.int")
  names(param.synthesis.lower) = names(param.synthesis.upper);

  return(list(lower = param.synthesis.lower, upper = param.synthesis.upper))
}

set.general.bounds.degr.splicing = function(lower.user = NULL, 
                                  upper.user = NULL,
                                  parametrization =c('cosine.beta'), 
                                  absolute.signal = TRUE)
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
  
  degr.splicing.upper = c(gamma.max, eps.gamma.max, phase.gamma.max, splicing.k.max)
  degr.splicing.lower = c(gamma.min, eps.gamma.min, phase.gamma.min, splicing.k.min)
  
  names(degr.splicing.upper) = c("gamma", "eps.gamma", "phase.gamma", "splicing.k")
  names(degr.splicing.lower) =  names(degr.splicing.upper)
  
  return(list(lower = degr.splicing.lower, upper = degr.splicing.upper))
  
}

set.bounds.general = function(model = 4, lower.user = NULL, upper.user = NULL,
                              parametrization =c('cosine.beta'), absolute.signal = TRUE)
{
  bounds.int = set.general.bounds.int();
  bounds.degr.splicing = set.general.bounds.degr.splicing();
  
  upper = c(bounds.degr.splicing$upper, bounds.int$upper);
  lower = c(bounds.degr.splicing$lower, bounds.int$lower);
  #upper = c(gamma.max, eps.gamma.max, phase.gamma.max, splicing.k.max,  param.synthesis.upper)
  #lower = c(gamma.min, eps.gamma.min, phase.gamma.min, splicing.k.min,  param.synthesis.lower)
  #names(upper) = c("gamma", "eps.gamma", "phase.gamma", "splicing.k", 
  #                         "Min.int", "Amp.int", "phase.int", "beta.int")
  #names(lower) =  names(upper)
  
  if(absolute.signal){
    if(model==1){ upper = upper[c(1, 4, 5)];  lower = lower[c(1, 4, 5)]; }
    if(model==2){ upper = upper[c(1, 4, 5:8)];  lower = lower[c(1, 4, 5:8)];}
    if(model==3){ upper = upper[c(1:5)]; lower = lower[c(1:5)]; }
  }
  
  return(list(lower = lower, upper = upper))
}

####################
## Set gene-specific parameter boundaries and inital values for fitting pre-mRNAs
## Initial values were smapled 
####################
Sampling.Initial.Values.for.fitting.S = function(S, Nfit.S = 4, zt = seq(0,94,by = 2)) 
{
  set.seed(8675309);
  
  bounds.general = set.bounds.general(model = 2);
  upper.general = bounds.general$upper; 
  lower.general = bounds.general$lower;
  
  Min.init = rep(min(S), Nfit.S)
  
  Amp.init = rep((max(S)-min(S)), Nfit.S)
  
  phase.init = zt[which.max(S)]
  phase.init = (rep(phase.init, Nfit.S)+rnorm(Nfit.S,sd = 3))%%24
  
  beta.min = lower.general[which(names(lower.general)=="beta.int")]; 
  beta.max = upper.general[which(names(upper.general)=="beta.int")];
  beta.init = lseq(beta.min, beta.max, length = Nfit.S)
  
  PAR.INIT.S = cbind(Min.init, Amp.init, phase.init, beta.init)
  colnames(PAR.INIT.S) = c('Min.int', 'Amp.int', "phase.int", "beta.int")
  
  return(PAR.INIT.S)
  
}

set.bounds.gene.s = function(S, range_scalingFactor=5)
{
  bounds.int = set.general.bounds.int();
  
  lower.g.s = c(min(S)/range_scalingFactor, (max(S)-min(S))/range_scalingFactor, bounds.int$lower[c(3:4)]);
  upper.g.s = c(max(S), (max(S)-min(S))*range_scalingFactor, bounds.int$upper[c(3:4)])
  
  return(list(upper = upper.g.s, lower = lower.g.s))
  
}

####################
## Set gene-specific parameter boundaries and inital values for fitting mRNAs
## Initial values were smapled 
####################
Sampling.Initial.Values.for.fitting.M = function(M, S, Nfit.M = 6, zt = seq(0,94,by = 2)) 
{
  set.seed(8675309);
  eps.m = min((max(M) - min(M))/mean(M)/2, 1);
  if(Nfit.M%%2==0){
    eps.gamma.init = c(sample(seq(0.1, 0.6, length = 10), Nfit.M/2, replace = TRUE), rep(eps.m/1.25, length=Nfit.M/2));
  }else{
    eps.gamma.init = c(sample(seq(0.1, 0.6, length = 10), (Nfit.M-1)/2, replace = TRUE), rep(eps.m/1.25, length=(Nfit.M+1)/2));
  }
  if(any(eps.gamma.init>0.8)) eps.gamma.init[which(eps.gamma.init>0.8)] = 0.8;
    
  phase.m = zt[which.max(M)]
  phase.gamma.init =  (rep((phase.m+12), Nfit.M)+rnorm(Nfit.M, sd = 3))%%24
    
  gamma.init = Gamma.Initiation(eps.gamma.init, 0.5, 6);
  #if(debug){cat('gamma inital values ', gamma.init, '\n')};
  
  # define the ratio between splicing rate and degratation rate
  a.init = rep(mean(M)/mean(S), Nfit.M)
  
  PAR.INIT.M = cbind(gamma.init, eps.gamma.init, phase.gamma.init, a.init)
  colnames(PAR.INIT.M) = c('gamma', 'eps.gamma', "phase.gamma", "splicing.k")
  
  return(PAR.INIT.M)
  
}

set.bounds.gene.m = function(M, S, range_scalingFactor=5)
{
  bounds.gamma.k = set.general.bounds.degr.splicing();
  a = mean(M)/mean(S);
  
  lower.gamma.k = c(bounds.gamma.k$lower[c(1:3)], a/range_scalingFactor)
  upper.gamma.k = c(bounds.gamma.k$upper[c(1:3)], a*range_scalingFactor)
  
  return(list(upper = upper.gamma.k, lower = lower.gamma.k))
  
}

####################
## Set gene-specific parameter boundaries and inital values for fitting mRNA and pre-mRNA together
## Initial values were smapled 
####################
Sampling.Initial.Values.for.fitting.M.S = function(M, S, model = 4, Nfit = 6, zt = seq(0,94,by = 2)) 
{
  set.seed(8675309);
  
  ### Define the initial values for degradation and splicing parameters
  a = mean(M)/mean(S) # define the ratio between splicing rate and degratation rate
  a.init = lseq(max(0.1, a/5), min(10^5, a*2), length=Nfit)
  max.half.life = 12;
  min.half.life = 30/60;
  
  if(model==2){
    Min.init = rep(res.fit.s[1], Nfit)
    Amp.init = seq(res.fit.s[2]/2, res.fit.s[2]*2, length=Nfit)
    phase.init = (rep(res.fit.s[3],Nfit)+rnorm(Nfit,sd = 2))%%24
    beta.init = lseq(max(1, (res.fit.s[4]-2)), min((res.fit.s[4]+2), 10), length=Nfit)
    
    gamma.init = c(rep(log(2)/lseq(max.half.life, min.half.life, length = Nfit), Nfit%/%Nfit), rep(log(2)/5,Nfit%%Nfit))
    
    PAR.INIT = cbind(gamma.init, a.init, Min.init, Amp.init, phase.init, beta.init); 
    colnames(PAR.INIT) = c('gamma','splicing.k','Min.int','Amp.int','phase.int','beta.int');
    
  }
  if(model==3){	
    ### initialize decay parameters for model 3
    eps.m = min((max(M) - min(M))/mean(M)/2, 1); # close to mRNA amplitude
    eps.gamma.init = c(sample(seq(0.1, 0.6, length = 10), Nfit/2, replace = TRUE), rep(eps.m, length=Nfit/2));
    if(length(which(eps.gamma.init>0.8))>=1) eps.gamma.init[which(eps.gamma.init>0.8)] = 0.8;
    phase.m = zt[which.max(M)]
    phase.gamma.init =  (rep((phase.m+12), Nfit)+rnorm(Nfit, sd = 2))%%24 ## antiphasic to mRNA
    gamma.init = Gamma.Initiation(eps.gamma.init, 0.5, 6);
    
    Min.init = rep(mean(S), Nfit)
    PAR.INIT = cbind(gamma.init, eps.gamma.init, phase.gamma.init, a.init, Min.init)
    colnames(PAR.INIT) = c('gamma','eps.gamma','phase.gamma','splicing.k', 'Min.int');
   
  }
  if(model==4){
    ### initialize decay parameters for model 4
    eps.gamma.init = seq(res.fit.m[2]/4, min(res.fit.m[2]*4, 1), length=Nfit);
    if(length(which(eps.gamma.init>0.8))>=1) eps.gamma.init[which(eps.gamma.init>0.8)] = 0.8;
    gamma.init = Gamma.Initiation(eps.gamma.init, 0.5, 8);
    phase.gamma.init = (rep(res.fit.m[3],Nfit)+rnorm(Nfit,sd = 1.5))%%24;
    a.init = seq(max(0.1, res.fit.m[4]/1.5), min(10^5, res.fit.m[4]*1.5), length=Nfit)
    
    Min.init = rep(res.fit.s[1], Nfit)
    Amp.init = seq(res.fit.s[2]/2, res.fit.s[2]*2, length=Nfit)
    phase.init = (rep(res.fit.s[3],Nfit)+rnorm(Nfit,sd = 2))%%24
    beta.init = lseq(max(1, (res.fit.s[4]-2)), min((res.fit.s[4]+2), 10), length=Nfit)
    
    PAR.INIT = cbind(gamma.init, eps.gamma.init, phase.gamma.init, a.init, Min.init, Amp.init, phase.init, beta.init);
    colnames(PAR.INIT) = c('gamma','eps.gamma','phase.gamma','splicing.k', 'Min.int','Amp.int','phase.int','beta.int')
  }
  
}

set.bounds.gene = function(M, S, model = 4, range_scalingFactor=5)
{
  a = mean(M)/mean(S) # define the ratio between splicing rate and degratation rate
  
  if(model==2){
    upper[2] = a*range_scalingFactor;   lower[2] = a/range_scalingFactor;
    lower[3] = min(S)/range_scalingFactor;     upper[3] = max(S)*range_scalingFactor;
    lower[4] = (max(S)-min(S))/range_scalingFactor;     upper[4] = (max(S)-min(S))*range_scalingFactor;
  }
  if(model==3){	
    upper[4] = a*range_scalingFactor;     lower[4] = a/range_scalingFactor;
    lower[5] = min(S)/range_scalingFactor;     upper[5] = max(S)*range_scalingFactor;
  }
  if(model==4){
    upper[4] = a*range_scalingFactor; lower[4] = a/range_scalingFactor;
    lower[5] = min(S)/range_scalingFactor;  upper[5] = max(S)*range_scalingFactor;
    lower[6] = (max(S)-min(S))/range_scalingFactor; upper[6] = (max(S)-min(S))*range_scalingFactor;
  }
  
  return(list(upper = upper, lower = lower))
  
}


