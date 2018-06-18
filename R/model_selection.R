##########################################################################
##########################################################################
## Project:
## Script purpose:
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Jun  4 13:48:08 2018
##########################################################################
##########################################################################
######################
# Functions for model selecion
######################
prob.best.model = function(BIC)
{
  if(length(which(is.na(BIC)==TRUE))>0){
    return(rep(NA, 5));
  }else{
    bic = BIC-min(BIC)
    prob.model = exp(-0.5*bic)
    prob.model = prob.model/sum(prob.model)
    
    #prob.model = exp(-0.5*BIC)
    #prob.model = prob.model/sum(prob.model)
    return(c(prob.model, which.min(BIC)));
  }
}

nb.outliers = function(param.fit, index.outliers = c(182, 183))
{
  if(length(index.outliers)==2){
    return(length(unlist(strsplit(param.fit[index.outliers[1]], ';'))) + length(unlist(strsplit(param.fit[index.outliers[2]], ';'))))
  }else{
    return('error')
  }
}

####################
## Model Selection for one specific gene
####################
my.model.selection.one.gene.loglike = function(param.fits.results, method = 'BIC', 
                                               outliers.m = c(), outliers.s = c(), absolute.signal=TRUE)
{
  ## method = c('BIC', 'AIC', 'AICc');
  set.nb.data.param(absolute.signal=TRUE);
  
  index = match(c('error.m1', 'error.m2', 'error.m3', 'error.m4'), names(param.fits.results))
  error.m1 = param.fits.results[index[1]]
  error.m2 = param.fits.results[index[2]]
  error.m3 = param.fits.results[index[3]]
  error.m4 = param.fits.results[index[4]]
  
  nb.data = nb.data - length(outliers.m[which(!is.na(outliers.m)==TRUE)]) - length(outliers.s[which(!is.na(outliers.s)==TRUE)]);
  
  ## the formula of BIC used here is chi-square+k*ln(n)==error/sigma^2+k(ln(n)) in which sigma of noise is supported to be known.
  if(method == 'BIC') {
    score.m1 = (error.m1) + log(nb.data)*n.param[1];
    score.m2 = (error.m2) + log(nb.data)*n.param[2];
    score.m3 = (error.m3) + log(nb.data)*n.param[3];
    score.m4 = (error.m4) + log(nb.data)*n.param[4];
  }
  
  if(method == 'AIC'|method == 'AICc') {
    score.m1 = (error.m1) + 2*n.param[1] 
    score.m2 = (error.m2) + 2*n.param[2]
    score.m3 = (error.m3) + 2*n.param[3] 
    score.m4 = (error.m4) + 2*n.param[4] 
    if(method== 'AICc') {
      score.m1 = score.m1 + 2*n.param[1]*(n.param[1]+1)/(nb.data-n.param[1]-1)
      score.m2 = score.m2 + 2*n.param[2]*(n.param[2]+1)/(nb.data-n.param[2]-1)
      score.m3 = score.m3 + 2*n.param[3]*(n.param[3]+1)/(nb.data-n.param[3]-1)
      score.m4 = score.m4 + 2*n.param[4]*(n.param[4]+1)/(nb.data-n.param[4]-1)
    }
  }
  scores = c(score.m1, score.m2, score.m3, score.m4)
  scores.relavtive = scores-min(scores)
  prob.model = exp(-0.5*scores.relavtive)
  prob.model = prob.model/sum(prob.model)
  
  res.MS = c(scores, which.min(scores), prob.model[which.min(scores)])
  names(res.MS) = paste(method, c('.m1', '.m2', '.m3', '.m4', '.best.model', '.prob.best.model'), sep='')
  
  return(res.MS)
}

####################
## Model Selection for all genes with the matrix of parameters as input
####################
my.model.selection.loglike = function(T, nb.data=96, method = 'BIC',  outliers = FALSE, absolute.signal=TRUE)
{
  ## method = c('BIC', 'AIC', 'AICc'); outliers = TRUE
  set.nb.param(absolute.signal)
  
  index = match(c('error.m1', 'error.m2', 'error.m3', 'error.m4'), colnames(T))
  error.m1 = T[,index[1]]
  error.m2 = T[,index[2]]
  error.m3 = T[,index[3]]
  error.m4 = T[,index[4]]
  
  if(outliers)
  {
    index.outliers = match(c('outlier.m', 'outlier.s'), colnames(T));
    nb.data = 96 - apply(T, 1, nb.outliers, index.outliers);
  }
  
  ## the formula of BIC used here is chi-square+k*ln(n)==error/sigma^2+k(ln(n)) in which sigma of noise is supported to be known.
  if(method == 'BIC'){
    score.m1 = (error.m1) + log(nb.data)*n.param[1];
    score.m2 = (error.m2) + log(nb.data)*n.param[2];
    score.m3 = (error.m3) + log(nb.data)*n.param[3];
    score.m4 = (error.m4) + log(nb.data)*n.param[4];
  }
  if(method == 'AIC'|method == 'AICc'){
    score.m1 = (error.m1) + 2*n.param[1] 
    score.m2 = (error.m2) + 2*n.param[2]
    score.m3 = (error.m3) + 2*n.param[3] 
    score.m4 = (error.m4) + 2*n.param[4] 
    if(method== 'AICc'){
      score.m1 = score.m1 + 2*n.param[1]*(n.param[1]+1)/(nb.data-n.param[1]-1)
      score.m2 = score.m2 + 2*n.param[2]*(n.param[2]+1)/(nb.data-n.param[2]-1)
      score.m3 = score.m3 + 2*n.param[3]*(n.param[3]+1)/(nb.data-n.param[3]-1)
      score.m4 = score.m4 + 2*n.param[4]*(n.param[4]+1)/(nb.data-n.param[4]-1)
    }
  }
  
  scores = cbind(score.m1, score.m2, score.m3, score.m4)
  prob.model = t(apply(as.matrix(scores), 1, prob.best.model));
  res.MS = cbind(scores, prob.model)
  colnames(res.MS) = paste(method, c('.m1', '.m2', '.m3', '.m4', '.prob.m1', '.prob.m2', '.prob.m3', '.prob.m4', '.best.model'), sep='')
  
  return(res.MS)
}


