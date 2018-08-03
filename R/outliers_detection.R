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
source("R/error_functions.R", local = TRUE)

detect.ouliters.loglike = function(param.fits.results, GeneDataSet)
{
  outlier.m = unlist(GeneDataSet$outlier.m)
  outlier.s = unlist(GeneDataSet$outlier.s)
  
  loglike.contribution = calculate.loglike.contribution(param.fits.results, GeneDataSet)
  loglike.m = loglike.contribution$loglike.m;
  loglike.s = loglike.contribution$loglike.s;
  
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
