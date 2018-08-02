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

####################
## main function for outlier detection 
####################
detect.ouliters.loglike = function(param.fits.results, GeneDataSet)
{
  # extract fitted parameters
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
  }else{
    var.s = unlist(GeneDataSet$var.s)
    var.m = unlist(GeneDataSet$var.m)
  }
  
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
    if(model == 3){
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
    
    loglike = -2*dnbinom(as.numeric(R.m), size=1/alpha.m, mu=as.numeric(read.m), log = TRUE)
    #cat(index.outliers.loglike(loglike), '\n')
    loglike.m[(model-1), ] = loglike
    loglike = -2*dnbinom(as.numeric(R.s), size=1/alpha.s, mu=as.numeric(read.s), log = TRUE)
    #cat(index.outliers.loglike(loglike), '\n')
    loglike.s[(model-1), ] = loglike
  }
  
  ## identified outlier index
  new.outliers.m = index.outliers.loglike(loglike.m[3,]); 
  new.outliers.s = index.outliers.loglike(loglike.s[3,]);
    
  if(length(outlier.m)==0){
    outlier.m = new.outliers.m;
  }else{
    new.outliers.m = setdiff(new.outliers.m, outlier.m);
    outlier.m = c(outlier.m, new.outliers.m);
  }
  
  if(length(outlier.s)==0){
    outlier.s = new.outliers.s;
  }else{
    new.outliers.s = setdiff(new.outliers.s, outlier.s);
    outlier.s = c(outlier.s, new.outliers.s);
  }
  
  outlier.m = outlier.m[order(outlier.m)];
  outlier.s = outlier.s[order(outlier.s)];
  
  res.outliers.detection = list(nb.newOutliers.m = length(new.outliers.m), 
                                nb.newOutliers.s = length(new.outliers.s), 
                                outlier.m = outlier.m,
                                outlier.s = outlier.s);
  
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

