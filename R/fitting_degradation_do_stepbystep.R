##########################################################################
##########################################################################
## Project:
## Script purpose:
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Tue Jun  5 11:03:16 2018
##########################################################################
##########################################################################
make.fits.with.all.models.for.one.gene.remove.outliers = function(mds, 
                                                                  gene.index = 1, 
                                                                  debug = FALSE,
                                                                  outliers.removal = FALSE, 
                                                                  identifiablity.analysis.gamma = FALSE, 
                                                                  parametrization = c('cosine.beta'), 
                                                                  absolute.signal = TRUE)
{
  ####################
  ## define here global variables and global functions that are accessed for all steps;
  ####################
  # set global functions
  source("R/utilities_generalFunctions.R");
  set.nb.param();
  set.time.points(mds$zt) #actually this can be also a global parameter, because it is gene-independent
  
  # set scaling factors here as global variables for NB mode
  if(mds$mode == "NB"){set.scaling.factors(mds$scaling.factors);}
  
  ## extract the data relevant for one specific gene and wrap it in a list
  GeneDataSet = ExtractGeneData(mds, gene.index)
   
  ####################
  ## fitting the data for each model
  ####################
  source("R/optimization_params.R", local = TRUE)
  
  if(!outliers.removal){ 
    # ignore outliers 
    # without outlier detection and removal
    if(debug){cat('starting optimization without outliers \n')}
    param.fits.results = make.fits.with.all.models.for.one.gene(GeneDataSet = GeneDataSet, debug = debug);
    #param.fits.results = make.fits.with.all.models.for.one.gene();
    outlier.m = NA;
    outlier.s = NA;
    
  }else{ 
    ## outlier detection and removal
    source("R/outliers_detection.R", local = TRUE)
    if(debug){cat('starting optimization with outlier detection ----------\n ');}
    
    newOutliers.m = c(1);
    newOutliers.s = c(1);
    
    ## Todo: threshold in which the outlier detection stops should be revised ...
    while((length(newOutliers.m) > 0 | length(newOutliers.s) > 0) & length(GeneDataSet$outlier.m) <= length(GeneDataSet$zt)/8 &
          length(GeneDataSet$outlier.s) <= length(GeneDataSet$zt)/8)
    {
      if(debug){cat('-- outlier index of mRNA : ', paste0(GeneDataSet$outlier.m, collapse = ","), "\n");  
                cat('-- outlier index of premRNA : ', paste0(GeneDataSet$outlier.s, collapse = ","),  "\n");}
            
      param.fits.results = make.fits.with.all.models.for.one.gene(GeneDataSet = GeneDataSet, debug = debug); 
      res.outliers.detection = detect.ouliters.loglike(param.fits.results, GeneDataSet);
      
      newOutliers.m = res.outliers.detection$newOutliers.m
      newOutliers.s = res.outliers.detection$newOutliers.s
      
      # Important Note: change outlier records in the matrix T which is used to pass the outlier index for optimization module
      if(length(newOutliers.m)>0){GeneDataSet$outlier.m = unique(c(GeneDataSet$outlier.m, unlist(newOutliers.m))); } 
      if(length(newOutliers.s)>0) {GeneDataSet$outlier.s = unique(c(GeneDataSet$outlier.s, unlist(newOutliers.s))); }
      
    }
    
    if(length(GeneDataSet$outlier.m)==0) {outlier.m = NA; 
    }else{outlier.m = GeneDataSet$outlier.m; outlier.m = outlier.m[order(outlier.m)]}
    if(length(GeneDataSet$outlier.s)==0) {outlier.s = NA;
    }else{outlier.s = GeneDataSet$outlier.s; outlier.s = outlier.s[order(outlier.s)]}
    
  }
  
  ####################
  ## Analyze non-identifiability using Profile Likelihood approache
  ####################
  if(identifiablity.analysis.gamma){
    if(debug){cat('starting non-identifiability analysis for gamma \n')}
    
    source("R/identifiability_analysis.R", local = TRUE)
    res.nonident.gamma = Identifiablity.analysis.gamma.all.models(param.fits.results, GeneDataSet);
                                                                                      
  }else{
    res.nonident.gamma = rep(NA, 8);
    names(res.nonident.gamma) = paste0(c('non.identifiability.gamma.L.m', 'non.identifiability.gamma.R.m'), 
                                                          rep(c(1:4), each=2))
  }
  
  ####################
  ## model selection  
  ####################
  if(debug){cat('starting model selection \n')}
  set.nb.data(GeneDataSet);
  source("R/model_selection.R", local = TRUE);
  res.model.sel = my.model.selection.one.gene.loglike(param.fits.results, method = 'BIC', 
                                                                 outlier.m = outlier.m, outlier.s = outlier.s) 
  
  ####################
  ## parameter transformation 
  ####################
  if(debug){cat('starting parameter transformation \n')}
  source("R/params_transformation_cleaning.R", local = TRUE);
  
  # param.fits.results.save4test = param.fits.results;
  param.transformed.cleaned = transform.parameter.combinations.cleaning(param.fits.results,
                                                                        res.model.sel,
                                                                        res.nonident.gamma)
  
  ####################
  ## output  
  ####################
  if(debug){cat("final result is being wrapped --\n"); cat("---------- Done !!! ----------\n")}
  
  return(list(norm.RPKM = list(norm.S = GeneDataSet$Norm.s, norm.M = GeneDataSet$Norm.s),
              param.combinations = list(m1 = param.fits.results[grep('.m1', names(param.fits.results))],
                                m2 = param.fits.results[grep('.m2', names(param.fits.results))],
                                m3 = param.fits.results[grep('.m3', names(param.fits.results))],
                                m4 = param.fits.results[grep('.m4', names(param.fits.results))]),
              nonident.analysis.for.gamma = res.nonident.gamma,
              outliers = list(outlier.m = paste(outlier.m, sep='', collapse = ','), outlier.s = paste(outlier.s, sep='', collapse = ',')),
              model.selection = res.model.sel,
              param.fit.cleaned = param.transformed.cleaned)
         )
  
}