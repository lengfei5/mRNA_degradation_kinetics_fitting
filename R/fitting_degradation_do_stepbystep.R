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
  # import and set global functions
  source("R/utilities_generalFunctions.R");
  
  set.time.points(mds$zt) #actually this can be also a global parameter, because it is gene-independent
  if(mds$mode == "NB"){
    # set scaling factors here as global variables
    set.scaling.factors(mds$scaling.factors)
  }
  
  ## extract the data relevant for one specific gene and wrap it in a list
  GeneDataSet = ExtractGeneData(mds, gene.index)
   
  ####################
  ## fitting the data for each model
  ####################
  source("R/optimization_params.R", local = TRUE)
  
  if(!outliers.removal){ 
    # ignore outliers 
    # without outlier detection and removal
    
    if(debug){cat('starting optimization without outliers \n ')}
    param.fits.results = make.fits.with.all.models.for.one.gene(GeneDataSet = GeneDataSet, debug = debug);
    #param.fits.results = make.fits.with.all.models.for.one.gene();
    outlier.m = NA;
    outlier.s = NA;
    
  }else{ ## outlier detection and removal
    source("R/outliers_detection.R", local = TRUE)
    
    if(debug){cat('starting optimization with outlier detection ----------\n ');}
    
    outlier.m = GeneDataSet$outlier.m;
    outlier.s = GeneDataSet$outlier.s;
    nb.newOutliers.m = 1;
    nb.newOutliers.s = 1;
    while((nb.newOutliers.m > 0 | nb.newOutliers.s > 0) & length(c(outlier.m, outlier.s)) <= 12)
    {
      if(debug){ cat('--outlier index of mRNA :', paste0(outlier.m, collapse = ","), "\n" );  
        cat('-- outlier index of premRNA : ', paste0(outlier.s, collapse = ","),  '\n'); }
      
      param.fits.results = make.fits.with.all.models.for.one.gene(GeneDataSet = GeneDataSet, outliers = TRUE, debug = debug); 
      
      res.outliers.detection = detect.ouliters.loglike(param.fits.results, GeneDataSet);
      
      nb.newOutliers.m = res.outliers.detection$nb.newOutliers.m
      nb.newOutliers.s = res.outliers.detection$nb.newOutliers.s
      outlier.m = res.outliers.detection$outlier.m;
      outlier.s = res.outliers.detection$outlier.s;
      
      # Important Note: change outlier records in the matrix T which is used to pass the outlier index for optimization module
      GeneDataSet$outlier.m = outlier.m;
      GeneDataSet$outlier.s = outlier.s;
    }
    
    if(length(outlier.m)==0) outlier.m = NA; 
    if(length(outlier.s)==0) outlier.s = NA;
  }
  
  ####################
  ## Analyze non-identifiability using Profile Likelihood approache
  ####################
  if(identifiablity.analysis.gamma){
    source("R/identifiability_analysis.R", local = TRUE)
    if(debug){cat('starting non-identifiability analysis for gamma \n')}
    res.nonident.analysis.gamma.all.models = Identifiablity.analysis.gamma.all.models(param.fits.results, GeneDataSet);
                                                                                      
  }else{
    res.nonident.analysis.gamma.all.models = rep(NA, 8);
    names(res.nonident.analysis.gamma.all.models) = paste0(c('non.identifiability.gamma.L.m', 'non.identifiability.gamma.R.m'), 
                                                          rep(c(1:4), each=2))
  }
  
  ####################
  ## model selection  
  ####################
  if(debug){cat('starting model selection \n')}
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
                                                                        res.nonident.analysis.gamma.all.models)
  
  ####################
  ## output  
  ####################
  if(debug){cat("final result is being wrapped --\n"); cat("----------Done !!!----------\n")}
  
  return(list(norm.RPKM = list(norm.S = GeneDataSet$Norm.s, norm.M = GeneDataSet$Norm.s),
              param.combinations = list(m1 = param.fits.results[grep('.m1', names(param.fits.results))],
                                m2 = param.fits.results[grep('.m2', names(param.fits.results))],
                                m3 = param.fits.results[grep('.m3', names(param.fits.results))],
                                m4 = param.fits.results[grep('.m4', names(param.fits.results))]),
              nonident.analysis.for.gamma = res.nonident.analysis.gamma.all.models,
              outliers = list(outlier.m = paste(outlier.m, sep='', collapse = ','), outlier.s = paste(outlier.s, sep='', collapse = ',')),
              model.selection = res.model.sel,
              param.fit.cleaned = param.transformed.cleaned)
         )
  
}