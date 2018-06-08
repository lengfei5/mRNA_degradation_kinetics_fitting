##########################################################################
##########################################################################
## Project:
## Script purpose: to detect outliers after mode fitting
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Tue Jun  5 12:04:24 2018
##########################################################################
##########################################################################

####################
## utility functions for outlier detection 
####################
index.outliers.loglike = function(data.xx, c=1.5)
{
  #c = 3
  #data.xx = c(2, 3, 6, 9, 13, 18, 21, 106)
  #data.xx = data.xx[which(!is.na(data.xx)==TRUE)]
  Q1 = quantile(data.xx, 0.25,type=5)
  Q3 = quantile(data.xx, 0.75, type=5)
  IQD = Q3 - Q1
  lower = Q1 - c*IQD
  upper = Q3 + c*IQD
  index = which(data.xx>upper)
  #boxplot(data.xx);abline(h=Q1);abline(h=Q3);
  return(index)
}

index.outliers = function(data.xx)
{
  c = 1.5
  #data.xx = c(2, 3, 6, 9, 13, 18, 21, 106)
  Q1 = quantile(data.xx, 0.25,type=5)
  Q3 = quantile(data.xx, 0.75, type=5)
  IQD = Q3 - Q1
  lower = Q1 - c*IQD
  upper = Q3 + c*IQD
  index = which(data.xx<lower|data.xx>upper)
  #boxplot(data.xx);abline(h=Q1);abline(h=Q3);
}

## testing the scope of variables in R
test.funciton = function()
{
  cat("test ----------\n")
  cat("intially defined L.m -- ", L.m , "\n")
  L.m = 10;
  cat("newly defefined L.m -- ", L.m, "\n")
  cat("L.s -- ", L.s /2, "\n")
  #g <- function(y) { Lm <- 100; f(y); }
  
}

####################
## main function for outlier detection 
####################
detect.ouliters.loglike = function(param.fits.results,
                                   R.m,
                                   R.s,
                                   a.m, 
                                   a.s,
                                   L.m,
                                   L.s,
                                   outlier.m = c(), 
                                   outlier.s = c(),
                                   nb.additonal.m = 1,
                                   nb.additonal.s = 1
                                   )
{
  loglike.m = matrix(NA, nrow=3, ncol=48)
  loglike.s = matrix(NA, nrow=3, ncol=48)
  
  for(model in c(2:4))
  {
    #cat(' Model ', model, '\n');
    gamma = param.fits.results[which(names(param.fits.results)==paste('gamma.m', model, sep=''))];
    splicing.k = param.fits.results[which(names(param.fits.results)==paste('splicing.k.m', model, sep=''))];
    Min = param.fits.results[which(names(param.fits.results)==paste('Min.int.m', model, sep=''))];
    if(model == 2){
      eps.gamma = 0;
      phase.gamma = 12;
      Amp=param.fits.results[which(names(param.fits.results)==paste('Amp.int.m', model, sep=''))];
      phase=param.fits.results[which(names(param.fits.results)==paste('phase.int.m', model, sep=''))];
      beta=param.fits.results[which(names(param.fits.results)==paste('beta.int.m', model, sep=''))];
    }
    if(model == 3)
    {
      eps.gamma = param.fits.results[which(names(param.fits.results)==paste('eps.gamma.m', model, sep=''))];
      phase.gamma = param.fits.results[which(names(param.fits.results)==paste('phase.gamma.m', model, sep=''))];
      Amp=0; phase=12; beta=1;
    }
    if(model == 4){
      eps.gamma = param.fits.results[which(names(param.fits.results)==paste('eps.gamma.m', model, sep=''))];
      phase.gamma = param.fits.results[which(names(param.fits.results)==paste('phase.gamma.m', model, sep=''))];
      Amp=param.fits.results[which(names(param.fits.results)==paste('Amp.int.m', model, sep=''))];
      phase=param.fits.results[which(names(param.fits.results)==paste('phase.int.m', model, sep=''))];
      beta=param.fits.results[which(names(param.fits.results)==paste('beta.int.m', model, sep=''))];
    }
    
    m = compute.m.beta(zt, gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
    s = compute.s.beta(zt, Min, Amp, phase, beta);
    read.m = convert.nb.reads(m, L.m)
    read.s = convert.nb.reads(s, L.s)
    
    loglike = -2*dnbinom(as.numeric(R.m), size=1/a.m, mu=as.numeric(read.m), log = TRUE)
    #cat(index.outliers.loglike(loglike), '\n')
    loglike.m[(model-1), ] = loglike
    loglike = -2*dnbinom(as.numeric(R.s), size=1/a.s, mu=as.numeric(read.s), log = TRUE)
    #cat(index.outliers.loglike(loglike), '\n')
    loglike.s[(model-1), ] = loglike
  }
  
  ## identified outlier index
  #additional.m = intersect(index.outliers.loglike(loglike.m[3,]), intersect(index.outliers.loglike(loglike.m[1,]), index.outliers.loglike(loglike.m[2,])))
  #additional.s = intersect(index.outliers.loglike(loglike.s[3,]), intersect(index.outliers.loglike(loglike.s[1,]), index.outliers.loglike(loglike.s[2,])))
   
  index.outliers.m = index.outliers.loglike(loglike.m[3,]); 
  index.outliers.s = index.outliers.loglike(loglike.s[3,]);
    
  if(length(outlier.m)==0){
    outlier.m = index.outliers.m;
  }else{
    index.outliers.m = setdiff(index.outliers.m, outlier.m);
    outlier.m = c(outlier.m, index.outliers.m);
  }
  
  if(length(outlier.s)==0){
    outlier.s = index.outliers.s;
  }else{
    index.outliers.s = setdiff(index.outliers.s, outlier.s);
    outlier.s = c(outlier.s, index.outliers.s);
  }
  
  outlier.m = outlier.m[order(outlier.m)];
  outlier.s = outlier.s[order(outlier.s)];
  
  res.outliers.detection = list(nb.newOutliers.m = length(index.outliers.m), 
                                nb.newOutliers.s = length(index.outliers.s), 
                                outlier.m = outlier.m,
                                outlier.s = outlier.s);
  
  return(res.outliers.detection)
  
}
