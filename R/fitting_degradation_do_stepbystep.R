##########################################################################
##########################################################################
## Project:
## Script purpose:
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Tue Jun  5 11:03:16 2018
##########################################################################
##########################################################################
source("R/configure.R")
source("R/preparation_before_modelfitting.R")
make.fits.with.all.models.for.one.gene.remove.outliers = function(T = T, 
                                                                  gene.index = 1, 
                                                                  zt = seq(0,94,by = 2), 
                                                                  i.ex = ZT.ex, 
                                                                  i.int = ZT.int, 
                                                                  debug = FALSE,
                                                                  outliers = FALSE, 
                                                                  Identifiablity.Analysis.by.Profile.Likelihood.gamma = TRUE, 
                                                                  PLOT.Ident.analysis = FALSE, 
                                                                  parametrization = c('cosine.beta'), 
                                                                  absolute.signal = TRUE)
{
  ####################
  ## check the parameters
  ####################
  #gg = 'Per2'; gene.index = which(T$gene==gg); source('functions.R');zt = seq(0,94,by = 2); debug=TRUE; outliers = TRUE;
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
  source("R/optimization_params.R")
  if(!outliers)
  {
    param.fits.results = make.fits.with.all.models.for.one.gene(T = T, gene.index = gene.index, debug = debug, zt = zt, 
                                                                i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE); 
    outlier.m = NA; outlier.s = NA;
    
  }else{
    ####################
    ## outlier detection and loop over the mode fitting if there are outliers 
    ####################
    ## source()
    outlier.m = c(); outlier.s = c();
    nb.additonal.m = 1; nb.additonal.s = 1;
    T$mRNA.outlier[gene.index] = '';  T$premRNA.outlier[gene.index] = '';
    
    while((nb.additonal.m>0 | nb.additonal.s>0) & sum(c(nb.additonal.m, nb.additonal.s))<=12)
    {
      if(debug){
        cat('Outlier of mRNA : ', T$mRNA.outlier[gene.index],  '\n');  cat('Outlier of premRNA : ', T$premRNA.outlier[gene.index],  '\n');
        #missing.data = length(unlist(strsplit(T$mRNA.outlier[j], ';'))) + length(unlist(strsplit(T$premRNA.outlier[j], ';')))
        #print(c(param.fits.results, my.model.selection.one.gene.loglike(param.fits.results, nb.data = (96-missing.data), method = 'BIC'),
        #        my.model.selection.one.gene.loglike(param.fits.results, nb.data = (96-missing.data), method = 'AICc'), 
        #        my.model.selection.one.gene.loglike(param.fits.results, nb.data = (96-missing.data), method = 'AIC')));
        #print(param.fits.results);
      }
      
      param.fits.results = make.fits.with.all.models.for.one.gene(T = T, gene.index = gene.index, debug = debug, zt = zt, 
                                                                  i.ex = ZT.ex, i.int = ZT.int, outliers = TRUE);
      ### chekc if there are outliers
      loglike.m = matrix(NA, nrow=3, ncol=48)
      loglike.s = matrix(NA, nrow=3, ncol=48)
      for(model in c(2:4))
      {
        #cat(' Model ', model, '\n');
        gamma = param.fits.results[which(names(param.fits.results)==paste('gamma.m', model, sep=''))];
        splicing.k = param.fits.results[which(names(param.fits.results)==paste('splicing.k.m', model, sep=''))];
        Min = param.fits.results[which(names(param.fits.results)==paste('Min.int.m', model, sep=''))];
        if(model == 2){
          eps.gamma = 0;phase.gamma = 12;
          Amp=param.fits.results[which(names(param.fits.results)==paste('Amp.int.m', model, sep=''))];
          phase=param.fits.results[which(names(param.fits.results)==paste('phase.int.m', model, sep=''))];
          beta=param.fits.results[which(names(param.fits.results)==paste('beta.int.m', model, sep=''))];
        }
        if(model == 3)
        {
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
        
        loglike = -2*dnbinom(as.numeric(R.m), size=1/a.m, mu=as.numeric(read.m), log = TRUE)
        #cat(index.outliers.loglike(loglike), '\n')
        loglike.m[(model-1), ] = loglike
        loglike = -2*dnbinom(as.numeric(R.s), size=1/a.s, mu=as.numeric(read.s), log = TRUE)
        #cat(index.outliers.loglike(loglike), '\n')
        loglike.s[(model-1), ] = loglike
      }
      #additional.m = intersect(index.outliers.loglike(loglike.m[3,]), intersect(index.outliers.loglike(loglike.m[1,]), index.outliers.loglike(loglike.m[2,])))
      #additional.s = intersect(index.outliers.loglike(loglike.s[3,]), intersect(index.outliers.loglike(loglike.s[1,]), index.outliers.loglike(loglike.s[2,])))
      additional.m = index.outliers.loglike(loglike.m[3,]); 
      additional.s = index.outliers.loglike(loglike.s[3,]);
      #apply(cbind(loglike.s, loglike.m), 1, sum)
      #apply(cbind(loglike.s[, setdiff(c(1:48), outlier.s)], loglike.m[, setdiff(c(1:48), outlier.m)]), 1, sum)
      if(length(outlier.m)==0)
      {
        outlier.m = additional.m;
      }else{
        additional.m = setdiff(additional.m, outlier.m);
        outlier.m = c(outlier.m, additional.m);
      }
      nb.additonal.m = length(additional.m);
      
      if(length(outlier.s)==0)
      {
        outlier.s = additional.s;
      }else{
        additional.s = setdiff(additional.s, outlier.s);
        outlier.s = c(outlier.s, additional.s);
      }
      nb.additonal.s = length(additional.s);
      
      outlier.m = outlier.m[order(outlier.m)];
      outlier.s = outlier.s[order(outlier.s)];
      
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
    source("R/identifiability_analysis.R")
    if(debug){cat('\t\t starting non-identifiability analysis for gamma \n')}
    
    M = norm.RPKM(R.m, L.m)
    S = norm.RPKM(R.s, L.s)
    a = mean(M)/mean(S) # ratio between splicing rate and degratation rate
    
    if(PLOT.Ident.analysis){
      pdf.name = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Total/Total_counts/identifiability_test_Real_Data/gamma_nonidentif_',
                       gene.name, '_all_models.pdf', sep=''); 
      pdf(pdf.name, width = 10 , height = 3);
      par(cex = 0.7, las = 1, mgp = c(2.0,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3);
      par(mfrow=c(1,3),cex=1.0)
      #### FAST analysis
      #ptm <- proc.time()
    }
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
    if(PLOT.Ident.analysis)  dev.off()
    #proc.time() - ptm;
    
  }
  
  ####################
  ## model selection  
  ####################
  
  ####################
  ## parameter transformation 
  ####################
  
  ####################
  ## output  
  ####################
  
  return(c(param.fits.results, outlier.m = paste(outlier.m, sep='', collapse = ';'), outlier.s = paste(outlier.s, sep='', collapse = ';')));
  
}