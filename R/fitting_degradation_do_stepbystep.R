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
  # import global functions
  source("R/utilities_generalFunctions.R"); 
  
  # set scaling factors here as global variables
  set.scaling.factors(mds$scaling.factors)
  
  zt  = mds$zt
  # gene.name = T$gene[gene.index];
  R.m = mds$M[gene.index, ]
  R.s = mds$P[gene.index, ]
  L.m = T$length.mRNA[gene.index];
  L.s = T$length.premRNA[gene.index];
  
  alpha.m = mds$dispersions.M[gene.index, ]
  alpha.s = mds$dispersions.P[gene.index, ]
  
  M = norm.RPKM(R.m, L.m)
  S = norm.RPKM(R.s, L.s)
  
  #identifiablity.analysis.gamma = Identifiablity.Analysis
  
  ####################
  ## fitting the data for each model
  ####################
  source("R/optimization_params.R", local = TRUE)
  
  if(!outliers.removal){
    ## without outlier detection and removal
    if(debug){cat('starting optimization without outliers \n ')}
    
    param.fits.results = make.fits.with.all.models.for.one.gene(R.m = R.m, R.s = R.s, 
                                                                alpha.m = alpha.m, alpha.s = alpha.s, 
                                                                L.m = L.m, L.s = L.s, 
                                                                zt = zt,
                                                                debug = debug); 
    #param.fits.results = make.fits.with.all.models.for.one.gene();
    outlier.m = NA;
    outlier.s = NA;
    
  }else{
    ## outlier detection and removal
    source("R/outliers_detection.R", local = TRUE)
    
    outlier.m = c(); 
    outlier.s = c();
    nb.newOutliers.m = 1; 
    nb.newOutliers.s = 1;
    T$mRNA.outlier[gene.index] = '';  T$premRNA.outlier[gene.index] = '';
    
    while((nb.newOutliers.m > 0 | nb.newOutliers.s > 0) & length(c(outlier.m, outlier.s)) <= 12)
    {
      if(debug){cat('starting optimization with outlier detection ----------\n ');
        cat('-- outlier index of mRNA :', paste0(outlier.m, collapse = ",") );  
        cat('-- outlier index of premRNA : ', paste0(outlier.s, collapse = ","),  '\n');
      }
      
      param.fits.results = make.fits.with.all.models.for.one.gene(T = T, gene.index = gene.index, debug = debug, zt = zt, 
                                                                  i.ex = ZT.ex, i.int = ZT.int, outliers = TRUE);
      
      res.outliers.detection = detect.ouliters.loglike(param.fits.results, R.m, R.s, alpha.m, alpha.s, L.m, L.s,  
                                                                  outlier.m = outlier.m, outlier.s = outlier.s, 
                                                                  nb.additonal.m = nb.additonal.m, 
                                                                  nb.additonal.s = nb.additonal.s);
      
      nb.newOutliers.m = res.outliers.detection$nb.newOutliers.m
      nb.newOutliers.s = res.outliers.detection$nb.newOutliers.s
      outlier.m = res.outliers.detection$outlier.m;
      outlier.s = res.outliers.detection$outlier.s;
      
      # Important Note: change outlier records in the matrix T which is used to pass the outlier index for optimization module
      T$mRNA.outlier[gene.index] = paste(outlier.m, sep='', collapse = ';')
      T$premRNA.outlier[gene.index] = paste(outlier.s, sep='', collapse = ';')
    
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
    
    res.nonident.analysis.gamma.all.models = Identifiablity.analysis.gamma.all.models(param.fits.results,
                                                                                      R.m, R.s, 
                                                                                      L.m, L.s,
                                                                                      alpha.m = alpha.m, 
                                                                                      alpha.s = alpha.s, 
                                                                                      outlier.m = outlier.m, 
                                                                                      outlier.s = outlier.s,
                                                                                      zt = zt);
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
  if(debug){cat('final result is ----------\n')}
  
  return(list(gene.length = list(L.s = L.s, L.m = L.m),
              readCounts = list(R.s = R.s, R.m = R.m),
              norm.RPKM = list(norm.S = S, norm.M = M),
              dispersions.param = list(dispersion.s = alpha.s, dispersion.m = alpha.m),
              param.combinations = list(m1 = param.fits.results[grep('.m1', names(param.fits.results))],
                                m2 = param.fits.results[grep('.m2', names(param.fits.results))],
                                m3 = param.fits.results[grep('.m3', names(param.fits.results))],
                                m4 = param.fits.results[grep('.m4', names(param.fits.results))]),
              nonident.analysis.for.gamma = res.nonident.analysis.gamma.all.models,
              outliers = list(outlier.m = paste(outlier.m, sep='', collapse = ';'), outlier.s = paste(outlier.s, sep='', collapse = ';')),
              model.sel = res.model.sel,
              param.fit.cleaned = param.transformed.cleaned
              )
         )
  
}