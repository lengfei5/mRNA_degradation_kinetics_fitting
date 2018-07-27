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
  R.m = unlist(GeneDataSet$R.m) #R.m = unlist(T[gene.index, i.ex]) ## nb of reads for exon
  R.s = unlist(GeneDataSet$R.s) #R.s = unlist(T[gene.index, i.int]) ## nb of reads for intron
  L.m = GeneDataSet$L.m # L.m = T$length.mRNA[gene.index];
  L.s = GeneDataSet$L.s  #L.s = T$length.premRNA[gene.index];
  
  alpha.m = unlist(GeneDataSet$alpha.m)
  alpha.s = unlist(GeneDataSet$alpha.s)
  
  if(outliers) {
    outlier.m = unlist(GeneDataSet$outlier.m)
    outlier.s = unlist(GeneDataSet$outlier.s)
  }else{
    outlier.m = c();
    outlier.s = c();
  }
  
  M = norm.RPKM(R.m, L.m)
  S = norm.RPKM(R.s, L.s)
  
  mu.m = convert.nb.reads(rep(mean(M), length(R.m)), L.m);
  mu.s = convert.nb.reads(rep(mean(S), length(R.s)), L.s);
  
  err = NB.error(R.m = R.m, R.s = R.s, alpha.m = alpha.m, alpha.s = alpha.s, 
                 mu.m = mu.m, mu.s = mu.s, outlier.m = outlier.m, outlier.s=outlier.s, specie = 'both')
  names(err) = 'error.m1';
  
  return(err)
  
}

########
### ERROR function for premRNA profile
########
f2min.int = function(par.int, R.s, L.s=1000, alpha.s=c(0.02, 48),  outlier.s = c(), zt = seq(0,94,by = 2), 
                     parametrization ='cosine.beta', debug = FALSE)
{
  R.s = as.numeric(unlist(R.s));
  alpha.s = as.numeric(unlist(alpha.s));
  #S = as.numeric(unlist(S)); 
  s = compute.s.beta(t = zt, Min = par.int[1], Amp = par.int[2], phase = par.int[3], beta= par.int[4]);
  mu.s = convert.nb.reads(s, L.s); ### convert normalized rpkm calculated from model into nb of reads
  
  #err = sum((log(S)-log(s))^2/sigma.ss^2);
  err = NB.error(R.s = R.s, alpha.s = alpha.s, mu.s = mu.s, outlier.s = outlier.s, specie = 'premRNA');
  #remain.s = setdiff(c(1:48), outlier.s)
  #error.S = -2*sum(dnbinom(R.s[remain.s], size=1/alpha.s[remain.s], mu=mu.s[remain.s], log = TRUE))
  #cat(par.int, '\n')
  #cat(error.S, '\n')
  #cat(err, '\n')
  return(err)
}

###
### ERROR function for mRNA profile
###
f2min.mrna = function(par.init.m, res.fit.s, R.m, L.m=1000, alpha.m=c(0.02, 48), outlier.m = c(), zt = seq(0,94,by = 2), norm.params=TRUE,
                      parametrization ='cosine.beta', debug = FALSE)
{
  ### par.int.m = res.fit.m
  R.m = as.numeric(unlist(R.m));
  alpha.m = as.numeric(unlist(alpha.m));
  res.fit.s = as.numeric(unlist(res.fit.s));
  #S = as.numeric(unlist(S));
  
  w = 2*pi/24;
  m = compute.m.beta(t = zt, par.init.m[1], 
                     par.init.m[2]*sqrt(1+w^2/par.init.m[1]^2), (par.init.m[3]-atan2(w, par.init.m[1])/w), par.init.m[4]*par.init.m[1],
                     Min = res.fit.s[1], Amp = res.fit.s[2], phase = res.fit.s[3], beta= res.fit.s[4]); 
  mu.m = convert.nb.reads(m, L.m); ### convert normalized rpkm calculated from model into nb of reads
  
  #err = sum((log(S)-log(s))^2/sigma.ss^2);
  err = NB.error(R.m = R.m, alpha.m = alpha.m, mu.m = mu.m, outlier.m = outlier.m, specie = 'mRNA'); 
  eps.non.scaled = par.init.m[2]*sqrt(1+w^2/par.init.m[1]^2); 
  err = err + sigmoid.bound.contraint(eps.non.scaled);  
  #if(eps.non.scaled>1){err = 10^10;}
  #cat(par.init.m, '\n');
  #cat(err, '\n');
  return(err)
}

####
#### ERROR function for premRNA and mRNA
####
f2min = function(par, R.m, R.s, L.m, L.s, alpha.m, alpha.s, outlier.m = c(), outlier.s = c(), model=4, zt = seq(0,94,by = 2), 
                 parametrization =c('cosine.beta'), debug = FALSE, norm.params = TRUE)
{
  # initialization of the parameters
  gamma = par[1];
  if(model==2)	
  {
    eps.gamma = 0.0;
    phase.gamma = 12.0; 
    splicing.k = par[2];
    param.synthesis.1 = par[3];
    param.synthesis.2 = par[4];  
    param.synthesis.3 = par[5]; 
    param.synthesis.4 = par[6]; 
  }
  if(model==3)
  {
    eps.gamma = par[2];
    phase.gamma = par[3];
    splicing.k = par[4];
    param.synthesis.1 = par[5];
    param.synthesis.2 = 0;  
    param.synthesis.3 = 0; 
    param.synthesis.4 = 1; 
  }
  if(model==4)
  {
    eps.gamma = par[2];
    phase.gamma = par[3];
    splicing.k = par[4];
    param.synthesis.1 = par[5];
    param.synthesis.2 = par[6];  
    param.synthesis.3 = par[7]; 
    param.synthesis.4 = par[8]; 
  }
  
  s = compute.s.beta(t = zt, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4)
  #print(par);
  #print(compute.m.beta(t = zt, gamma, eps.gamma, phase.gamma, splicing.k, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4, simulation.only=TRUE))
  w = 2*pi/24;
  m = compute.m.beta(t = zt, gamma, eps.gamma*sqrt(1+w^2/gamma^2), (phase.gamma-atan2(w, gamma)/w), splicing.k*gamma,
                     param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4);
  ### convert rpkm into mean of read numbers
  mu.m = convert.nb.reads(m, L.m);
  mu.s = convert.nb.reads(s, L.s);
  #print(m);
  err.fit = NB.error(R.m = R.m, R.s = R.s, alpha.m = alpha.m, alpha.s = alpha.s, mu.m = mu.m, mu.s = mu.s, 
                     outlier.m = outlier.m, outlier.s = outlier.s, specie = 'both');
  
  eps.non.scaled = eps.gamma*sqrt(1+w^2/gamma^2);
  err.fit = err.fit + sigmoid.bound.contraint(eps.non.scaled);  
  #print(eps.non.scaled)
  #print(err.fit)
  
  return(err.fit)
}

NB.error = function(R.m=rep(1000, 48), R.s=c(100, 48), alpha.m=rep(0.02, 48), alpha.s=rep(0.03, 48), 
                    mu.m=rep(1000, 48), mu.s=rep(100, 48), outlier.m = c(), outlier.s=c(), specie = 'both', intense.debug=FALSE)
{
  if(specie=='both'|specie=='mRNA')
  {
    R.m = as.numeric(unlist(R.m));
    alpha.m = as.numeric(alpha.m);
    mu.m = as.numeric(unlist(mu.m));
    remain.m = setdiff(c(1:48), outlier.m)
    
    if(any(R.m[remain.m]<0)|any(mu.m<=0))
    {
      error.M = 10^10;
    }else{
      #error.M = -2*sum(R.m*log(mu.m/(mu.m+1/alpha.m))-1.0/alpha.m*log(mu.m+1.0/alpha.m))
      #ptm <- proc.time()
      #error.M = -2*sum(R.m*log(mu.m*alpha.m/(1.0+mu.m*alpha.m)) - 1.0/alpha.m*log(1.0+alpha.m*mu.m)
      #                 + lgamma(R.m+1.0/alpha.m) -lgamma(1.0/alpha.m) - lfactorial(R.m));
      #proc.time() - ptm
      #kk = 3
      #ptm <- proc.time()
      error.M = -2*sum(dnbinom(R.m[remain.m], size=1/alpha.m[remain.m], mu=mu.m[remain.m], log = TRUE))
      #error.M = -2*sum(log((1-q)*dnbinom(R.m[remain.m], size=1/alpha.m[remain.m], mu=mu.m[remain.m], log = FALSE) + q*exp(-10)))
      #proc.time() - ptm
      #cat(error.M, error.M2, '\n')
    }
  }
  if(specie=='both'|specie=='premRNA')
  {
    R.s = as.numeric(unlist(R.s));
    alpha.s = as.numeric(alpha.s);
    mu.s = as.numeric(unlist(mu.s));
    remain.s = setdiff(c(1:48), outlier.s)
    
    if(any(R.s[remain.s]<0)|any(mu.s<=0))
    {
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
  if(specie == 'mRNA') 
  {
    error = error.M;
  }
  if(specie == 'premRNA') 
  {
    error = error.S;
  }
  if(specie == 'both')
  {
    error = error.M + error.S; 
  }
  #error = error.M
  #if(intense.debug){cat('-------- error.S = ',error.S)}
  #if(intense.debug){cat('-------- error.M = ',error.M)}
  #if(intense.debug){cat('-------- error.M = ',error,'\n')}
  return(error)
  
}
