##########################################################################
##########################################################################
## Project:
## Script purpose:
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Jun  4 13:46:18 2018
##########################################################################
##########################################################################
Identifiablity.analysis.gamma.each.model = function (error.opt, params.opt, lower, upper, R.m, R.s, L.m, L.s, alpha.m, alpha.s, model = 3, zt,
                                                     outlier.m = c(), outlier.s = c(), PLOT.PL = FALSE, gene2plot='Per1') 
{
  # error.opt=errors.fit[imin]; params.opt =res.fit;PLOT.PL = TRUE; gene2plot=gene.name;
  #cat('FAST identifiability analysis for gamma \n');
  if(length(outlier.m)>0) outlier.m = outlier.m[which(!is.na(outlier.m)==TRUE)];
  if(length(outlier.s)>0) outlier.s = outlier.s[which(!is.na(outlier.s)==TRUE)];
  params = params.opt;
  gamma.opt = params[1];
  #gamma.opt = log(2)/(24)
  nb.test.iden = 20;
  list.par = lseq(log(2)/24, log(2)/(10/60), length=nb.test.iden);
  #list.par = unique(c(list.par, gamma.opt))
  #list.par = list.par[order(list.par)];
  #ii = which(list.par==gamma.opt);
  
  #res.pl = rep(NA, Nfit.pl)
  profiles = matrix(NA, nrow =length(list.par) , ncol=2);
  colnames(profiles) = c('gamma', 'pls');
  
  profiles[,1] = list.par;
  ptm <- proc.time()
  for(n in 1:nrow(profiles))
  {
    opt.pl = optim(params, par.fixed = profiles[n, 1], index.fixed = 1, f2min.profile, 
                   R.m = R.m, R.s = R.s, L.m=L.m, L.s = L.s, alpha.m=alpha.m, alpha.s=alpha.s, outlier.m = outlier.m, outlier.s = outlier.s, model = model, 
                   zt = zt, method = 'L-BFGS-B', lower = lower, upper = upper)
    profiles[n, 2] = opt.pl$value;
  }
  proc.time() - ptm;
  
  #### plot the identifiablity analysis (cf. Jens Timmer 2009)
  PL = data.frame(profiles, stringsAsFactors = FALSE)
  #qchisq(p, df, ncp = 0, lower.tail = TRUE)
  if(model==2) df = 6; if(model==3) df = 5; if(model==4) df = 8;
  threshold1 = qchisq(p=0.95, df=df, lower.tail = TRUE) + error.opt;
  threshold2 = qchisq(p=0.95, df=1, lower.tail = TRUE) + error.opt;
  
  ## the following threshod is used for the current code
  threshold3 = qchisq(p=0.68, df=1, lower.tail = TRUE) + error.opt;
  
  if(PLOT.PL){ ## NOT Used here
    #pdf.name = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Total/Total_counts/identifiability_test_Real_Data/gamma_nonidentif_',
    #                   gene2plot, '_m', model, '.pdf', sep=''); 
    #pdf(pdf.name, width = 5 , height = 3);
    #par(cex = 0.7, las = 1, mgp = c(2.0,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3);
    lims = range(c(error.opt, threshold1, threshold2, (threshold1+5)));
    plot(PL$gamma, PL$pls, type='l', col='darkblue', log='', lwd=1.5, xlab='gamma', ylab='-2log(PL)', ylim = lims, main = paste('m', model, sep='')); 
    points(PL$gamma, PL$pls, type='p', cex=0.8, col='darkblue'); 
    points(params[1], error.opt, col='darkred', cex=1.5, pch=8); 
    abline(h=threshold1, col='red', lwd=1.5, lty=2);abline(h=threshold2, col='red', lwd=1.5, lty=1);abline(h=threshold3, col='red', lwd=1.5, lty=1)  
    #dev.off()
  }
  
  bound.left = c(); bound.right = c();
  #bound.left = length(which(PL$gamma<gamma.opt & PL$pls>threshold2));
  #bound.right = length(which(PL$gamma>gamma.opt & PL$pls>threshold2));
  bound.left = max(PL$pls[which(PL$gamma<=gamma.opt)], na.rm = TRUE) - error.opt;
  bound.right = max(PL$pls[which(PL$gamma>=gamma.opt)], na.rm = TRUE) - error.opt;
  res.bounds = c(bound.left, bound.right);
  names(res.bounds) = c('non.identifiability.gamma.L', 'non.identifiability.gamma.R');
  return(res.bounds);
  
}

Identifiablity.analysis.gamma.all.models = function(param.fits.results, GeneDataSet)
{
  zt = unlist(GeneDataSet$zt)
  R.m = unlist(GeneDataSet$R.m) #R.m = unlist(T[gene.index, i.ex]) ## nb of reads for exon
  R.s = unlist(GeneDataSet$R.s) #R.s = unlist(T[gene.index, i.int]) ## nb of reads for intron
  L.m = GeneDataSet$L.m # L.m = T$length.mRNA[gene.index];
  L.s = GeneDataSet$L.s  #L.s = T$length.premRNA[gene.index];
  
  alpha.m = unlist(GeneDataSet$alpha.m)
  alpha.s = unlist(GeneDataSet$alpha.s)
  
  outlier.m = unlist(GeneDataSet$outlier.m)
  outlier.s = unlist(GeneDataSet$outlier.s)
  
  M = norm.RPKM(R.m, L.m)
  S = norm.RPKM(R.s, L.s)
  a = mean(M)/mean(S) # ratio between splicing rate and degratation rate
  
  res.iden.all = c();
  for(model in c(2:4))
  {
    cat('\t model ', model, '\n');
    error.fit = param.fits.results[which(names(param.fits.results)==paste('error.m', model, sep=''))];
    res.fit = param.fits.results[grep(paste('.m', model, sep=''), names(param.fits.results))];
    res.fit = res.fit[-1];
    
    bounds.gene = set.bounds.gene(GeneDataSet = GeneDataSet, model);
    upper = bounds.gene$upper; 
    lower = bounds.gene$lower;
    
    if(model==2) res.fit = res.fit[1:6];
    if(model==3) res.fit = res.fit[1:5];
    if(model==4) res.fit = res.fit[1:8];
    
    res.nonident.analysis.gamma.spec.model = Identifiablity.analysis.gamma.each.model(error.opt=error.fit, 
                                                                           params.opt =res.fit, 
                                                                           lower = lower, 
                                                                           upper = upper, 
                                                                           R.m = R.m, 
                                                                           R.s = R.s, 
                                                                           L.m=L.m, 
                                                                           L.s = L.s, 
                                                                           alpha.m=alpha.m, 
                                                                           alpha.s=alpha.s, 
                                                                           outlier.m = outlier.m, 
                                                                           outlier.s = outlier.s, 
                                                                           model = model, 
                                                                           zt = zt);
    
    names(res.nonident.analysis.gamma.spec.model) = paste(names(res.nonident.analysis.gamma.spec.model), '.m', model, sep='');
    res.iden.all = c(res.iden.all, res.nonident.analysis.gamma.spec.model);
  }
  
  return(res.iden.all)
  
}

