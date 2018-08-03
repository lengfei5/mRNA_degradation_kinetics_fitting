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
#source("R/error_functions.R", local = TRUE)
index.outliers.loglike = function(data.xx, c=1.5)
{
  #c = 3;
  #data.xx = c(2, 3, 6, 9, 13, 18, 21, NA, NA, 106)
  #data.xx = data.xx[which(!is.na(data.xx)==TRUE)]
  Q1 = quantile(data.xx, 0.25, type=5, na.rm = TRUE)
  Q3 = quantile(data.xx, 0.75, type=5, na.rm = TRUE)
  IQD = Q3 - Q1
  lower = Q1 - c*IQD
  upper = Q3 + c*IQD
  # boxplot(data.xx);abline(h=c(Q3 + 1.5*IQD), col='blue'); abline(h=c(Q3 + 3*IQD), col='red');
  return(which(data.xx>upper))
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

####################
# main function for outlier detection based on the contribution of -2loglikelihood  
# this iterative outlier detection based on the loglikelihood, because the contribution to the loglikehood is relevant to
# model selection (BIC) used afterward, since some outliers will dominant the BIC score if they contribute too much
# in addition, the outlier is detected based the most complicated model, i.e. M4; the idea is that if data are identified as outliers 
# in M4, and they are also outliers in M1-M3
####################
detect.ouliters.loglike = function(param.fits.results, GeneDataSet)
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
  model = 4
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
  
  ## not consider the outliers identified already
  if(length(outlier.m)>0) loglike.m[outlier.m] = NA;
  if(length(outlier.s)>0) loglike.s[outlier.s] = NA;
  ## identified outlier index
  newOutliers.m = index.outliers.loglike(loglike.m); 
  newOutliers.s = index.outliers.loglike(loglike.s);
  
  ## if there are more than 3 outliers for P or M; then pick the first three 
  if(length(newOutliers.m)>3){
   newOutliers.m = newOutliers.m[order(-loglike.m[newOutliers.m])];
   newOutliers.m = newOutliers.m[c(1:3)]
  }
  if(length(newOutliers.s)>3){
    newOutliers.s = newOutliers.s[order(-loglike.s[newOutliers.s])];
    newOutliers.s = newOutliers.s[c(1:3)]
  }
  
  res.outliers.detection = list(newOutliers.m = newOutliers.m, 
                                newOutliers.s = newOutliers.s)
                                
  
  return(res.outliers.detection)
  
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

