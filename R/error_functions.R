##########################################################################
##########################################################################
## Project:
## Script purpose: all error functions used in the optimization step
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Jun  4 12:27:15 2018
##########################################################################
##########################################################################
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

###########
## error for model 1
###########
calculate.error.for.flat.model = function(GeneDataSet, debug = FALSE, outliers = FALSE, 
                                          parametrization = c('cosine.beta'), absolute.signal = TRUE)
{
  # unwrap GeneDataSet
  zt = unlist(GeneDataSet$zt)
  M = unlist(GeneDataSet$Norm.m);
  S = unlist(GeneDataSet$Norm.s);
  
  outlier.m = unlist(GeneDataSet$outlier.m)
  outlier.s = unlist(GeneDataSet$outlier.s)
  
  if(GeneDataSet$mode == "NB"){
    R.m = unlist(GeneDataSet$R.m) #R.m = unlist(T[gene.index, i.ex]) ## nb of reads for exon
    R.s = unlist(GeneDataSet$R.s) #R.s = unlist(T[gene.index, i.int]) ## nb of reads for intron
    L.m = GeneDataSet$L.m # L.m = T$length.mRNA[gene.index];
    L.s = GeneDataSet$L.s  #L.s = T$length.premRNA[gene.index];
    
    alpha.m = unlist(GeneDataSet$alpha.m)
    alpha.s = unlist(GeneDataSet$alpha.s)
    mu.m = convert.nb.reads(rep(mean(M), length(R.m)), L.m);
    mu.s = convert.nb.reads(rep(mean(S), length(R.s)), L.s);
    
    err = NB.error(R.m = R.m, R.s = R.s, alpha.m = alpha.m, alpha.s = alpha.s, 
                   mu.m = mu.m, mu.s = mu.s, outlier.m = outlier.m, outlier.s=outlier.s, species = 'both')
    
  }else{
    var.m = unlist(GeneDataSet$var.m);  var.s = unlist(GeneDataSet$var.s);
    mu.m = rep(mean(M), length(M)); mu.s = rep(mean(S), length(S));
    
    err = Gaussian.error(M = M, S = S, var.m = var.m, var.s = var.s,
                   mu.m = mu.m, mu.s = mu.s, outlier.m = outlier.m, outlier.s=outlier.s, species = 'both')
  }
  
  names(err) = 'error.m1';
  return(err)
  
}

########
### ERROR function for premRNA profile
########
f2min.int = function(par.int, GeneDataSet,
                     parametrization ='cosine.beta', debug = FALSE)
{
  zt = unlist(GeneDataSet$zt)
  M = unlist(GeneDataSet$Norm.m);
  S = unlist(GeneDataSet$Norm.s);
  
  outlier.m = unlist(GeneDataSet$outlier.m)
  outlier.s = unlist(GeneDataSet$outlier.s)
  
  s = compute.s.beta(t = zt, Min = par.int[1], Amp = par.int[2], phase = par.int[3], beta= par.int[4]);
  
  if(GeneDataSet$mode == "NB"){
    # R.s=R.s, L.s=L.s, alpha.s=alpha.s, outlier.s=outlier.s, zt = zt;
    R.s = as.numeric(unlist(GeneDataSet$R.s));
    L.s = as.numeric(GeneDataSet$L.s);
    alpha.s = as.numeric(unlist(GeneDataSet$alpha.s));
    
    mu.s = convert.nb.reads(s, L.s); ### convert normalized rpkm calculated from model into nb of reads
    err = NB.error(R.s = R.s, alpha.s = alpha.s, mu.s = mu.s, outlier.s = outlier.s, species = 'premRNA');
  }else{
    var.s = unlist(GeneDataSet$var.s)
    mu.s = s;
    err = Gaussian.error(S = S, var.s = var.s, mu.s = mu.s, outlier.s=outlier.s, species = 'premRNA')
    
  }
  
  return(err)
}

###
### ERROR function for mRNA profile
###
f2min.mrna = function(par.init.m, res.fit.s, GeneDataSet, debug = FALSE, norm.params=TRUE, parametrization ='cosine.beta')
{
  w = 2*pi/24;
  # R.m=R.m, L.m=L.m, alpha.m=alpha.m, outlier.m=outlier.m, zt = zt, norm.params=norm.params,
  zt = unlist(GeneDataSet$zt)
  M = unlist(GeneDataSet$Norm.m);
  S = unlist(GeneDataSet$Norm.s);
  
  outlier.m = unlist(GeneDataSet$outlier.m)
  outlier.s = unlist(GeneDataSet$outlier.s)
  
  res.fit.s = as.numeric(unlist(res.fit.s));
  m = compute.m.beta(t = zt, par.init.m[1], 
                     par.init.m[2]*sqrt(1+w^2/par.init.m[1]^2), (par.init.m[3]-atan2(w, par.init.m[1])/w), par.init.m[4]*par.init.m[1],
                     Min = res.fit.s[1], Amp = res.fit.s[2], phase = res.fit.s[3], beta= res.fit.s[4]); 
  
  if(GeneDataSet$mode == "NB"){
    R.m = as.numeric(unlist(GeneDataSet$R.m));
    L.m = as.numeric(GeneDataSet$L.m);
    alpha.m = as.numeric(unlist(GeneDataSet$alpha.m));
    
    mu.m = convert.nb.reads(m, L.m); ### convert normalized rpkm calculated from model into nb of reads
    err = NB.error(R.m = R.m, alpha.m = alpha.m, mu.m = mu.m, outlier.m = outlier.m, species = 'mRNA');
  }else{
    var.m = unlist(GeneDataSet$var.m)
    mu.m = m;
    err = Gaussian.error(M = M, var.m = var.m, mu.m = mu.m, outlier.m = outlier.m, species = 'mRNA')
  }
  
  eps.non.scaled = par.init.m[2]*sqrt(1+w^2/par.init.m[1]^2); 
  err = err + sigmoid.bound.contraint(eps.non.scaled);
  
  return(err)
}

####
#### ERROR function for premRNA and mRNA
####
f2min = function(par, GeneDataSet, model=4, debug = FALSE, parametrization =c('cosine.beta'), norm.params = TRUE)
{
  w = 2*pi/24;
  
  # initialization of the parameters
  gamma = par[1];
  if(model==2){
    eps.gamma = 0.0;
    phase.gamma = 12.0; 
    splicing.k = par[2];
    param.synthesis.1 = par[3];
    param.synthesis.2 = par[4];  
    param.synthesis.3 = par[5]; 
    param.synthesis.4 = par[6]; 
  }
  if(model==3) {
    eps.gamma = par[2];
    phase.gamma = par[3];
    splicing.k = par[4];
    param.synthesis.1 = par[5];
    param.synthesis.2 = 0;  
    param.synthesis.3 = 0; 
    param.synthesis.4 = 1; 
  }
  if(model==4) {
    eps.gamma = par[2];
    phase.gamma = par[3];
    splicing.k = par[4];
    param.synthesis.1 = par[5];
    param.synthesis.2 = par[6];  
    param.synthesis.3 = par[7]; 
    param.synthesis.4 = par[8]; 
  }
  
  # unwrap the GeneDataSet
  zt = unlist(GeneDataSet$zt)
  M = unlist(GeneDataSet$Norm.m);
  S = unlist(GeneDataSet$Norm.s);
  
  outlier.m = unlist(GeneDataSet$outlier.m)
  outlier.s = unlist(GeneDataSet$outlier.s)
  
  s = compute.s.beta(t = zt, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4)
  m = compute.m.beta(t = zt, gamma, eps.gamma*sqrt(1+w^2/gamma^2), (phase.gamma-atan2(w, gamma)/w), splicing.k*gamma,
                     param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4);
  
  if(GeneDataSet$mode == "NB"){
    R.m = as.numeric(unlist(GeneDataSet$R.m));
    L.m = as.numeric(GeneDataSet$L.m);
    alpha.m = as.numeric(unlist(GeneDataSet$alpha.m));
    R.s = as.numeric(unlist(GeneDataSet$R.s));
    L.s = as.numeric(GeneDataSet$L.s);
    alpha.s = as.numeric(unlist(GeneDataSet$alpha.s));
    
    ### convert rpkm into mean of read numbers
    mu.m = convert.nb.reads(m, L.m);
    mu.s = convert.nb.reads(s, L.s);
   
    err.fit = NB.error(R.m = R.m, R.s = R.s, alpha.m = alpha.m, alpha.s = alpha.s, mu.m = mu.m, mu.s = mu.s, 
                       outlier.m = outlier.m, outlier.s = outlier.s, species = 'both');
  }else{
    var.s = unlist(GeneDataSet$var.s)
    var.m = unlist(GeneDataSet$var.m)
    mu.s = s;
    mu.m = m;
    err.fit = Gaussian.error(M = M, S = S, var.m = var.m, var.s = var.s, 
                             mu.m = mu.m, mu.s = mu.s, outlier.m = outlier.m, outlier.s = outlier.s, species = 'both')
  }
  
  eps.non.scaled = eps.gamma*sqrt(1+w^2/gamma^2);
  err.fit = err.fit + sigmoid.bound.contraint(eps.non.scaled);  
  #print(eps.non.scaled)
  #print(err.fit)
  
  return(err.fit)
}

###############################
# ERROR function for profile likelihood (identifiability analysis)
###############################
f2min.profile = function(param.start, GeneDataSet, which.to.fix, value.to.fix, model=4)
{
  ## convert back the list of parameters
  new.starts = rep(NA, length(c(value.to.fix, param.start)))
  new.starts[which.to.fix] = value.to.fix;
  if(which.to.fix==1) {
    new.starts = c(value.to.fix, param.start)
  }else{
    if(which.to.fix==(length(param.start)+1)) {
      new.starts = c(param.start, value.to.fix)
    }else{
      new.starts = c(param.start[c(1:(which.to.fix-1))], value.to.fix, param.start[c(which.to.fix:length(param.start))])
    }
  }
  
  err.profile = f2min(new.starts, GeneDataSet, model=model)
    
  return(err.profile)
  
}

###############################
# function for loglikelihood contribution from each data
# used for outlier detection  
###############################
calculate.loglike.contribution = function(param.fits.results, GeneDataSet, model = 4, w = 2*pi/24)
{
  w = 2*pi/24;
  #loglike.m = matrix(NA, nrow=3, ncol=48)
  #loglike.s = matrix(NA, nrow=3, ncol=48)
  
  ## extract data from GeneDataSet list
  zt = unlist(GeneDataSet$zt)
  M = unlist(GeneDataSet$Norm.m);
  S = unlist(GeneDataSet$Norm.s);
  outlier.m = unlist(GeneDataSet$outlier.m)
  outlier.s = unlist(GeneDataSet$outlier.s)
  
  ## extract fitted parameters for M4
  #model = 4
  gamma = param.fits.results[which(names(param.fits.results)==paste('gamma.m', model, sep=''))];
  eps.gamma = param.fits.results[which(names(param.fits.results)==paste('eps.gamma.m', model, sep=''))];
  phase.gamma = param.fits.results[which(names(param.fits.results)==paste('phase.gamma.m', model, sep=''))];
  splicing.k = param.fits.results[which(names(param.fits.results)==paste('splicing.k.m', model, sep=''))];
  
  Min.int = param.fits.results[which(names(param.fits.results)==paste('Min.int.m', model, sep=''))];
  Amp.int = param.fits.results[which(names(param.fits.results)==paste('Amp.int.m', model, sep=''))];
  phase.int = param.fits.results[which(names(param.fits.results)==paste('phase.int.m', model, sep=''))];
  beta.int = param.fits.results[which(names(param.fits.results)==paste('beta.int.m', model, sep=''))];
  
  #m = compute.m.beta(zt, gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
  m = compute.m.beta(t = zt, gamma, eps.gamma*sqrt(1+w^2/gamma^2), (phase.gamma-atan2(w, gamma)/w), splicing.k*gamma,
                     Min.int, Amp.int, phase.int, beta.int);
  s = compute.s.beta(zt, Min.int, Amp.int, phase.int, beta.int);
  
  if(GeneDataSet$mode == "NB"){
    R.m = unlist(GeneDataSet$R.m) #R.m = unlist(T[gene.index, i.ex]) ## nb of reads for exon
    R.s = unlist(GeneDataSet$R.s) #R.s = unlist(T[gene.index, i.int]) ## nb of reads for intron
    L.m = GeneDataSet$L.m # L.m = T$length.mRNA[gene.index];
    L.s = GeneDataSet$L.s  #L.s = T$length.premRNA[gene.index];
    alpha.m = unlist(GeneDataSet$alpha.m)
    alpha.s = unlist(GeneDataSet$alpha.s)
    
    mu.m = convert.nb.reads(m, L.m);
    mu.s = convert.nb.reads(s, L.s);
    
    loglike.m = -2*dnbinom(as.numeric(R.m), size=1/alpha.m, mu=as.numeric(mu.m), log = TRUE)
    loglike.s = -2*dnbinom(as.numeric(R.s), size=1/alpha.s, mu=as.numeric(mu.s), log = TRUE)
    
  }else{
    var.s = unlist(GeneDataSet$var.s)
    var.m = unlist(GeneDataSet$var.m)
    
    loglike.m = (log(M)-log(m))^2/var.m;
    loglike.s = (log(S)-log(s))^2/var.s
  }
  
  return(list(loglike.m = loglike.m, loglike.s = loglike.s))
  
}



########################################################
# Section: Error functions for NB and Gaussian modes
########################################################
NB.error = function(R.m=rep(1000, 48), R.s=c(100, 48), alpha.m=rep(0.02, 48), alpha.s=rep(0.03, 48), 
                    mu.m=rep(1000, 48), mu.s=rep(100, 48), outlier.m = c(), outlier.s=c(), species = 'both', intense.debug=FALSE)
{
  if(species=='both'|species=='mRNA'){
    R.m = as.numeric(unlist(R.m));
    alpha.m = as.numeric(alpha.m);
    mu.m = as.numeric(unlist(mu.m));
    remain.m = setdiff(c(1:length(R.m)), outlier.m)
    
    if(any(R.m[remain.m]<0)|any(mu.m<=0)){
      error.M = 10^10;
    }else{
      #ptm <- proc.time()
      error.M = -2*sum(dnbinom(R.m[remain.m], size=1/alpha.m[remain.m], mu=mu.m[remain.m], log = TRUE))
      #error.M = -2*sum(log((1-q)*dnbinom(R.m[remain.m], size=1/alpha.m[remain.m], mu=mu.m[remain.m], log = FALSE) + q*exp(-10)))
      #proc.time() - ptm
    }
  }
  if(species=='both'|species=='premRNA') {
    R.s = as.numeric(unlist(R.s));
    alpha.s = as.numeric(alpha.s);
    mu.s = as.numeric(unlist(mu.s));
    remain.s = setdiff(c(1:length(R.s)), outlier.s)
    
    if(any(R.s[remain.s]<0)|any(mu.s<=0)) {
      error.S = 10^10;
    }else{
      #error.S = -2*sum(R.s*log(mu.s/(mu.s+1/alpha.s))-1.0/alpha.s*log(mu.s+1.0/alpha.s))
      #error.S = -2*sum(R.s*log(mu.s*alpha.s/(1.0+mu.s*alpha.s)) - 1.0/alpha.s*log(1.0+alpha.s*mu.s)
      #                 + lgamma(R.s+1.0/alpha.s) -lgamma(1.0/alpha.s) - lfactorial(R.s));
      error.S = -2*sum(dnbinom(R.s[remain.s], size=1/alpha.s[remain.s], mu=mu.s[remain.s], log = TRUE))
      #error.S = -2*sum(log((1-q)*dnbinom(R.s[remain.s], size=1/alpha.s[remain.s], mu=mu.s[remain.s], log = FALSE) + q*exp(-10)))
      #cat(error.S, error.S2, '\n')
    }
  }
  
  ## -2loglike
  if(species == 'mRNA') { error = error.M; }
  if(species == 'premRNA') { error = error.S;}
  if(species == 'both') {error = error.M + error.S; }
    
  return(error)
  
}

### error for the Gaussian mode
Gaussian.error = function(M = re(100, 48), S = re(10, 48), var.m = rep(0.05, 48), var.s = rep(0.1, 48), 
               mu.m = rep(90, 48), mu.s = rep(12, 48), outlier.m = c(), outlier.s=c(), species = 'both', intense.debug=FALSE)
{
  if(species=='both'|species=='mRNA'){
    M = as.numeric(unlist(M));
    s2.m = as.numeric(unlist(var.m));
    mu.m = as.numeric(unlist(mu.m));
    
    remain.m = setdiff(c(1:length(M)), outlier.m)
    if(any(M[remain.m]<0)|any(mu.m<=0)){
      error.M = 10^10;
    }else{
      error.M = sum((log(M[remain.m])-log(mu.m[remain.m]))^2/s2.m[remain.m]);
    }
  }
  if(species=='both'|species=='premRNA') {
    S = as.numeric(unlist(S));
    s2.s = as.numeric(var.s);
    mu.s = as.numeric(unlist(mu.s));
    
    remain.s = setdiff(c(1:length(S)), outlier.s)
    
    if(any(S[remain.s]<0)|any(mu.s<=0)){
      error.S = 10^10;
    }else{
      error.S = sum((log(S[remain.s])-log(mu.s[remain.s]))^2/s2.s[remain.s]);
    }
  }
  
  ## -2loglike in logNormal distributed noise
  if(species == 'mRNA') { error = error.M; }
  if(species == 'premRNA') { error = error.S;}
  if(species == 'both') { error = error.M + error.S; }
  
  return(error)
  
}

# old gaussian error function used in function_RPKM_v3.R
# not used here
error.from.functions_RPKM_v3.R = function(S,s,M,m, sigma.ss=rep(0.2, length(S)), sigma.mm=rep(0.1, length(M)), log = TRUE, intense.debug=FALSE)
{
  S = as.numeric(unlist(S));
  s = as.numeric(unlist(s));
  M =  as.numeric(unlist(M));
  m =  as.numeric(unlist(m));
  
  if(log) ## calculate separately the errors of premRNAs and mRNAs
  {
    if(any(S<=0)|any(s<=0)){error.S = 10^10;
    }else{
      error.S = sum((log(S)-log(s))^2/sigma.ss^2);
    }
    if(any(M<=0)|any(m<=0)){
      error.M = 10^10
    }else{
      error.M = sum((log(M)-log(m))^2/sigma.mm^2);
    }
  }
  
  ### Here we estimate gene-specific variances with replicates
  error = error.M + error.S ## -2loglike
  
  #if(intense.debug){cat('______________________________-------- error.S = ',error.S,'\n')}
  #if(intense.debug){cat('______________________________-------- error.M = ',error.M,'\n')}
  return(error)  
}

