##########################################################################
##########################################################################
## Project:
## Script purpose:
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Jun  4 13:48:08 2018
##########################################################################
##########################################################################
####################
## Model Selection for one specific gene
####################
my.model.selection.one.gene.loglike = function(param.fits.results, method = 'BIC', 
                                               outlier.m = c(), outlier.s = c(), absolute.signal=TRUE)
{
  ## method = c('BIC', 'AIC', 'AICc');
  index = match(c('error.m1', 'error.m2', 'error.m3', 'error.m4'), names(param.fits.results))
  error.m1 = param.fits.results[index[1]]
  error.m2 = param.fits.results[index[2]]
  error.m3 = param.fits.results[index[3]]
  error.m4 = param.fits.results[index[4]]
  
  cat('\t nb of data --', nb.data, "\n");
  nb.data.without.ouliers = nb.data - length(outlier.m[which(!is.na(outlier.m)==TRUE)]) - length(outlier.s[which(!is.na(outlier.s)==TRUE)]);
  
  ## the formula of BIC used here is chi-square+k*ln(n)==error/sigma^2+k(ln(n)) in which sigma of noise is supported to be known.
  if(method == 'BIC') {
    cat('\t nb of data to use for model selection --', nb.data.without.ouliers, "\n");
    score.m1 = (error.m1) + log(nb.data.without.ouliers)*nb.param[1];
    score.m2 = (error.m2) + log(nb.data.without.ouliers)*nb.param[2];
    score.m3 = (error.m3) + log(nb.data.without.ouliers)*nb.param[3];
    score.m4 = (error.m4) + log(nb.data.without.ouliers)*nb.param[4];
  }
  
  if(method == 'AIC'|method == 'AICc') {
    score.m1 = (error.m1) + 2*nb.param[1] 
    score.m2 = (error.m2) + 2*nb.param[2]
    score.m3 = (error.m3) + 2*nb.param[3] 
    score.m4 = (error.m4) + 2*nb.param[4] 
    if(method== 'AICc') {
      score.m1 = score.m1 + 2*nb.param[1]*(nb.param[1]+1)/(nb.data.without.outliers-nb.param[1]-1)
      score.m2 = score.m2 + 2*nb.param[2]*(nb.param[2]+1)/(nb.data.without.outliers-nb.param[2]-1)
      score.m3 = score.m3 + 2*nb.param[3]*(nb.param[3]+1)/(nb.data.without.outliers-nb.param[3]-1)
      score.m4 = score.m4 + 2*nb.param[4]*(nb.param[4]+1)/(nb.data.without.outliers-nb.param[4]-1)
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

