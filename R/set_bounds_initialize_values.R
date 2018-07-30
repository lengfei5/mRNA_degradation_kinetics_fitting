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
  
}

####################
## function for set gene-specific parameter boundaries 
####################
set.bounds.gene.m = function(GeneDataSet, range_scalingFactor=5)
{
  M = GeneDataSet$Norm.m;
  S = GeneDataSet$Norm.s;
  bounds.gamma.k = set.general.bounds.degr.splicing();
  a = mean(M)/mean(S);
  
  lower.gamma.k = c(bounds.gamma.k$lower[c(1:3)], a/range_scalingFactor)
  upper.gamma.k = c(bounds.gamma.k$upper[c(1:3)], a*range_scalingFactor)
  
  return(list(upper = upper.gamma.k, lower = lower.gamma.k))
  
}

set.bounds.gene.s = function(GeneDataSet, range_scalingFactor=5)
{
  S = GeneDataSet$Norm.s;
  
  bounds.int = set.general.bounds.int();
  
  lower.g.s = c(Min.int = min(S)/range_scalingFactor, Amp.int = (max(S)-min(S))/range_scalingFactor, bounds.int$lower[c(3:4)]);
  upper.g.s = c(Min.int = max(S), Amp.int = (max(S)-min(S))*range_scalingFactor, bounds.int$upper[c(3:4)])
  
  return(list(upper = upper.g.s, lower = lower.g.s))
  
}

set.bounds.gene = function(GeneDataSet, model = 4, range_scalingFactor=5)
{
  M = GeneDataSet$Norm.m;
  S = GeneDataSet$Norm.s;
  
  bounds.general = set.bounds.general(model = model)
  upper = bounds.general$upper;
  lower = bounds.general$lower;
  
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

####################
##  set inital values for fitting pre-mRNAs
## Initial values were smapled 
####################
Sampling.Initial.Values.for.fitting.S = function(GeneDataSet, Nfit.S = 4) 
{
  S = GeneDataSet$Norm.s;
  zt = GeneDataSet$zt;
  
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

####################
## Set gene-specific parameter boundaries and inital values for fitting mRNAs
## Initial values were smapled 
####################
Sampling.Initial.Values.for.fitting.M = function(GeneDataSet, Nfit.M = 6)
{
  M = GeneDataSet$Norm.m;
  S = GeneDataSet$Norm.s;
  zt = GeneDataSet$zt;
  
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


####################
## Set gene-specific parameter boundaries and inital values for fitting mRNA and pre-mRNA together
## Initial values were smapled 
####################
Sampling.Initial.Values.for.fitting.M.S = function(GeneDataSet, model = 4, Nfit = 6, res.fit.s = NULL, res.fit.m = NULL) 
{
  M = GeneDataSet$Norm.m;
  S = GeneDataSet$Norm.s;
  zt = GeneDataSet$zt;
  
  set.seed(8675309);
  ### Define the initial values for degradation and splicing parameters
  a = mean(M)/mean(S) 
  a.init = lseq(max(0.1, a/5), min(10^5, a*5), length=Nfit)
  
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
  
  return(PAR.INIT)
}



