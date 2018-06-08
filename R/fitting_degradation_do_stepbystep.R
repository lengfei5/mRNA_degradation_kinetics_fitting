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
  ## check the parameters
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
  if(Identifiablity.Analysis.by.Profile.Likelihood.gamma)
  {
    source("R/identifiability_analysis.R", local = TRUE)
    if(debug){cat('starting non-identifiability analysis for gamma \n')}
    
    M = norm.RPKM(R.m, L.m)
    S = norm.RPKM(R.s, L.s)
    a = mean(M)/mean(S) # ratio between splicing rate and degratation rate
    
    for(model in c(2:4))
    {
      cat('Model ', model, '\n');
      bounds = set.bounds(model = model);
      upper = bounds$upper; 
      lower = bounds$lower;
      error.fit = param.fits.results[which(names(param.fits.results)==paste('error.m', model, sep=''))];
      res.fit = param.fits.results[grep(paste('.m', model, sep=''), names(param.fits.results))];
      res.fit = res.fit[-1];
      if(model==2){
        res.fit = res.fit[1:6];
        lower[2] = a/5;upper[2] = a*5;
        lower[3] = min(S)/5; upper[3] = max(S)*5;
        lower[4] = (max(S)-min(S))/5; upper[4] = (max(S)-min(S))*5;
      }
      if(model==3){
        res.fit = res.fit[1:5];
        lower[4] = a/5;upper[4] = a*5; 
        lower[5] = min(S)/5; upper[5] = max(S)*5;
      } 
      if(model==4){
        res.fit = res.fit[1:8];
        lower[4] = a/5;upper[4] = a*5;
        lower[5] = min(S)/5;upper[5] = max(S)*5;
        lower[6] = (max(S)-min(S))/5;upper[6] = (max(S)-min(S))*5;
      } 
      res.nonident.analysis.gamma = Identifiablity.analysis.gamma.each.model(error.opt=error.fit, params.opt =res.fit, 
                                                                             lower = lower, upper = upper, 
                                                                             R.m = R.m, R.s = R.s, L.m=L.m, L.s = L.s, alpha.m=a.m, alpha.s=a.s, 
                                                                             outlier.m = outlier.m, outlier.s = outlier.s, model = model, zt = zt, 
                                                                             PLOT.PL = PLOT.Ident.analysis, gene2plot=gene.name);
      
      names(res.nonident.analysis.gamma) = paste(names(res.nonident.analysis.gamma), '.m', model, sep='');
      param.fits.results = c(param.fits.results, res.nonident.analysis.gamma);
      #names(param.fits.results)[length(param.fits.results)] = paste('non.identifiability.gamma.m', model, sep = '') 
    }
    
  }
  
  ####################
  ## model selection  
  ####################
  if(debug){cat('starting model selection \n')}
  
  
  ####################
  ## parameter transformation 
  ####################
  if(debug){cat('starting parameter transformation \n')}
  
  
  ####################
  ## output  
  ####################
  if(debug){cat('final result is \n')}
  
  return(c(param.fits.results, outlier.m = paste(outlier.m, sep='', collapse = ';'), outlier.s = paste(outlier.s, sep='', collapse = ';')));
  
}