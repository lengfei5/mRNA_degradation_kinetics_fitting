##########################################################################
##########################################################################
## Project:
## Script purpose:
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Tue Jun  5 11:03:16 2018
##########################################################################
##########################################################################
make.fits.with.all.models.for.one.gene.remove.outliers = function(T = T, 
                                                                  gene.index = 1, 
                                                                  zt = seq(0,94,by = 2), 
                                                                  i.ex = ZT.ex, 
                                                                  i.int = ZT.int, 
                                                                  debug = FALSE,
                                                                  outliers.removal = FALSE, 
                                                                  Identifiablity.Analysis.by.Profile.Likelihood.gamma = FALSE, 
                                                                  parametrization = c('cosine.beta'), 
                                                                  absolute.signal = TRUE)
{
  ####################
  ## define here global variables for all steps; 
  ####################
  gene.name = T$gene[gene.index];
  R.m = T[gene.index, ZT.ex]
  R.s = T[gene.index, ZT.int]
  a.m = rep(as.numeric(T[gene.index, grep('alpha.mRNA.ZT', colnames(T))]), 4);
  a.s = rep(as.numeric(T[gene.index, grep('alpha.premRNA.ZT', colnames(T))]), 4);
  L.m = T$length.mRNA[gene.index];
  L.s = T$length.premRNA[gene.index];
  
  ####################
  ## fitting the data for each model
  ####################
  source("R/optimization_params.R", local = TRUE)
  
  if(!outliers.removal){
    ## without outlier detection and removal
    if(debug){cat('starting optimization without outliers \n ')}
    
    param.fits.results = make.fits.with.all.models.for.one.gene(T = T, gene.index = gene.index, debug = debug, zt = zt, 
                                                                i.ex = ZT.ex, i.int = ZT.int, outliers = outliers.removal); 
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
        cat('-- outlier index of mRNA :', paste0(outlier.m, collapse = ",") );  cat('-- outlier index of premRNA : ', paste0(outlier.s, collapse = ","),  '\n');
      }
      
      param.fits.results = make.fits.with.all.models.for.one.gene(T = T, gene.index = gene.index, debug = debug, zt = zt, 
                                                                  i.ex = ZT.ex, i.int = ZT.int, outliers = outliers.removal);
      
      res.outliers.detection = detect.ouliters.loglike(param.fits.results, R.m, R.s, a.m, a.s, L.m, L.s,  
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
  if(Identifiablity.Analysis.by.Profile.Likelihood.gamma){
    source("R/identifiability_analysis.R", local = TRUE)
    if(debug){cat('starting non-identifiability analysis for gamma \n')}
    
    res.nonident.analysis.gamma.all.models = Identifiablity.analysis.gamma.all.models(param.fits.results,
                                                                                      R.m, R.s, 
                                                                                      L.m, L.s,
                                                                                      alpha.m = a.m, 
                                                                                      alpha.s = a.s, 
                                                                                      outlier.m = outlier.m, 
                                                                                      outlier.s = outlier.s,
                                                                                      zt = zt);
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
  
  
  ####################
  ## output  
  ####################
  if(debug){cat('final result is \n')}
  
  return(list(param.fits = list(m1 = param.fits.results[grep('.m1', names(param.fits.results))],
                                m2 = param.fits.results[grep('.m2', names(param.fits.results))],
                                m3 = param.fits.results[grep('.m3', names(param.fits.results))],
                                m4 = param.fits.results[grep('.m4', names(param.fits.results))]),
              nonident.analysis.for.gamma = res.nonident.analysis.gamma.all.models,
              outliers = list(outlier.m = paste(outlier.m, sep='', collapse = ';'), outlier.s = paste(outlier.s, sep='', collapse = ';')),
              model.sel = res.model.sel
              
              )
         )
  
}