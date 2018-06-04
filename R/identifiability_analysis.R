##########################################################################
##########################################################################
## Project:
## Script purpose:
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Jun  4 13:46:18 2018
##########################################################################
##########################################################################
###
### ERROR function for profile likelihood (identifiability analysis)
###
f2min.proflie = function(par, par.fixed, index.fixed, R.m, R.s, L.m, L.s, alpha.m, alpha.s, outlier.m = c(), outlier.s = c(), model=3,
                         zt = seq(0,94,by = 2), norm.params=TRUE)
{
  par[index.fixed] = par.fixed; ### fixe one parameter
  gamma = par[1];
  if(model==2)	
  {
    eps.gamma = 0.0;
    phase.gamma = 12.0; 
    splicing.k = par[2];
    param.synthesis.1 = par[3];
    param.synthesis.2 = par[4];  
    param.synthesis.3 = par[5]; 
    param.synthesis.4 = par[6]; 
  }
  if(model==3)
  {
    eps.gamma = par[2];
    phase.gamma = par[3];
    splicing.k = par[4];
    param.synthesis.1 = par[5];
    param.synthesis.2 = 0;  
    param.synthesis.3 = 0; 
    param.synthesis.4 = 1; 
  }
  if(model==4)
  {
    eps.gamma = par[2];
    phase.gamma = par[3];
    splicing.k = par[4];
    param.synthesis.1 = par[5];
    param.synthesis.2 = par[6];  
    param.synthesis.3 = par[7]; 
    param.synthesis.4 = par[8]; 
  }
  if(model==3)
  {
    param.synthesis.2 = 0;  
    param.synthesis.3 = 0; 
    param.synthesis.4 = 1; 
  }
  if(model==4)
  {
    param.synthesis.2 = par[6];  
    param.synthesis.3 = par[7]; 
    param.synthesis.4 = par[8]; 
  }
  
  s = compute.s.beta(t = zt, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4)
  #print(par);
  #print(compute.m.beta(t = zt, gamma, eps.gamma, phase.gamma, splicing.k, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4, simulation.only=TRUE))
  w = 2*pi/24;
  m = compute.m.beta(t = zt, gamma, eps.gamma*sqrt(1+w^2/gamma^2), (phase.gamma-atan2(w, gamma)/w), splicing.k*gamma,
                     param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4);
  ### convert rpkm into mean of read numbers
  mu.m = convert.nb.reads(m, L.m);
  mu.s = convert.nb.reads(s, L.s);
  
  #print(m);
  err.profile = NB.error(R.m = R.m, R.s = R.s, alpha.m = alpha.m, alpha.s = alpha.s, mu.m = mu.m, mu.s = mu.s, 
                         outlier.m = outlier.m, outlier.s = outlier.s, specie = 'both');
  eps.non.scaled = eps.gamma*sqrt(1+w^2/gamma^2); 
  err.profile = err.profile + sigmoid.bound.contraint(eps.non.scaled);  
  #if(eps.non.scaled>1){err.profile = 10^10;}
  #print(par);
  
  return(err.profile)
}


Identifiablity.analysis.gamma.each.model = function (error.opt, params.opt, lower, upper, R.m, R.s, L.m, L.s, alpha.m, alpha.s, model = 3, zt,
                                                     outlier.m = c(), outlier.s = c(), PLOT.PL = FALSE, gene2plot='gene_test') 
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
    opt.pl = optim(params, par.fixed = profiles[n, 1], index.fixed = 1, f2min.proflie, 
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
  threshold3 = qchisq(p=0.68, df=1, lower.tail = TRUE) + error.opt;
  #qchisq(p=0.68, df=1, lower.tail = TRUE)
  if(PLOT.PL)
  {
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

Identifiablity.analysis.M3.M4.Example = function(error.opt, params, lower, upper, ## initial conditions and boundaries for maximum profile likelihood
                                                 R.m, R.s, L.m, L.s, alpha.m, alpha.s, outlier.m, outlier.s, model, zt, PLOT = FALSE, gene.name = 'gene2test') ## all arguments for optimization function
{
  # error.opt = error.fit; params = res.fit; 
  #ptm <- proc.time()
  ### 3 parameters to test identifiability
  nb.test.iden = c(50, 50, 50)
  #res.pl = rep(NA, Nfit.pl)
  keep = matrix(NA, nrow = sum(nb.test.iden), ncol=5);
  
  colnames(keep) = c('index.par', 'pls', 'par1', 'par2', 'par3')
  for(index.par in c(1:3))
  {
    cat(index.par, '\n');
    ## index.par = 1;
    nb.test = nb.test.iden[index.par];
    if(index.par==1)  
    {
      list.par = lseq(lower[index.par], upper[index.par], length=nb.test);
    }else{
      list.par = seq(lower[index.par], upper[index.par], length=nb.test);
    }
    
    #keep = matrix(NA, ncol=6, nrow=nb.test)
    for(nn in 1:length(list.par))
    {
      # nn = 1;
      # par.test = list.par[1]
      par.test = list.par[nn];
      #ptm <- proc.time()
      opt.pl = optim(params, par.fixed = par.test, index.fixed = index.par, f2min.proflie, 
                     R.m = R.m, R.s = R.s, L.m=L.m, L.s = L.s, alpha.m=alpha.m, alpha.s=alpha.s, outlier.m = outlier.m, outlier.s = outlier.s, model = model, 
                     zt = zt, method = 'L-BFGS-B', lower = lower, upper = upper)
      #proc.time() - ptm
      
      if(index.par==1) keep[nn, ] = c(index.par, opt.pl$value, par.test, opt.pl$par[c(2, 3)])
      if(index.par==2) keep[(nn+sum(nb.test.iden[1])), ] = c(index.par, opt.pl$value, opt.pl$par[1], par.test, opt.pl$par[3])
      if(index.par==3) keep[(nn+sum(nb.test.iden[c(1:2)])), ] = c(index.par, opt.pl$value, opt.pl$par[c(1, 2)], par.test)
    }
  }
  
  #### plot the identifiablity analysis (cf. Jens Timmer 2009)
  keep = data.frame(keep, stringsAsFactors = FALSE)
  
  #qchisq(p, df, ncp = 0, lower.tail = TRUE)
  if(model==3) df = 5;
  if(model==4) df = 8;
  threshold1 = qchisq(p=0.95, df=df, lower.tail = TRUE) + error.opt;
  threshold2 = qchisq(p=0.95, df=1, lower.tail = TRUE) + error.opt;
  #qchisq(p=0.68, df=1, lower.tail = TRUE)
  if(PLOT)
  {
    pdf.name = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Total/Total_counts/identifiability_test/Iden_Analysis_norm_parames_Model3_', 
                     gene.name, '.pdf', sep=''); 
    pdf(pdf.name, width = 12 , height = 3);
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3);
  }
  par(mfrow=c(1,3),cex=1.0)
  index.p = which(keep[,1]==1); lims = range(c(error.opt, threshold1, threshold2, (error.opt+15)));
  plot(keep$par1[index.p], keep$pls[index.p], type='l', log='', cex=0.5, xlab='gamma', ylab='-2log(PL)', ylim = lims); 
  points(params[1], error.opt, col='darkred', cex=1.5, pch=8);abline(h=threshold1, col='red', lwd=1.5, lty=2);abline(h=threshold2, col='red', lwd=1.5, lty=2)
  
  index.p = which(keep[,1]==2); 
  plot(keep$par2[index.p], keep$pls[index.p], type='l', log='', cex=0.5, xlab='epsilon.norm', ylab='-2log(PL)',ylim=lims);
  points(params[2], error.opt, col='darkred', cex=1.5, pch=8);abline(h=threshold1, col='red', lwd=1.5, lty=2);abline(h=threshold2, col='red', lwd=1.5, lty=2)
  index.p = which(keep[,1]==3);
  plot(keep$par3[index.p], keep$pls[index.p], type='l', log='', cex=0.5, xlab='phase.norm', ylab='-2log(PL)',ylim=lims); 
  points(params[3], error.opt, col='darkred', cex=1.5, pch=8);abline(h=threshold1, col='red', lwd=1.5, lty=2);abline(h=threshold2, col='red', lwd=1.5, lty=2)
  #proc.time() - ptm
  if(PLOT) dev.off()
  
  Test.function.relation = FALSE
  if(Test.function.relation)
  {
    plot(keep$list.par[kk], keep$par2[kk], type='b', log=log, cex=0.5)
    points(keep$list.par[kk], 0.75*sqrt(w^2+keep$list.par[kk]^2)/keep$list.par[kk], type='l', col='red', lwd=2.0)
    plot(keep$list.par[kk], keep$par3[kk], type='b', log=log, cex=0.5)
    points(keep$list.par[kk], 1.65-atan2(w, keep$list.par[kk])/2/pi*24, type='l', col='red', lwd=2.0)
    plot(keep$list.par[kk], keep$par4[kk], type='b', log='', cex=0.5)
    
    plot(keep$list.par[kk], keep$par1[kk], type='b', log=log, cex=0.5)
    #points(keep$list.par[kk], sqrt(w^2+keep$list.par[kk]^2)/keep$list.par[kk]*0.2, type='l', col='red', lwd=2.0)
    plot(keep$list.par[kk], keep$par3[kk], type='b', log=log, cex=0.5)
    #points(keep$list.par[kk], 16-atan2(w, keep$list.par[kk])/2/pi*24, type='l', col='red', lwd=2.0)
    plot(keep$list.par[kk], keep$par4[kk], type='b', log='', cex=0.5)
    
    
    plot(keep$list.par[kk], keep$par1[kk], type='b', log=log, cex=0.5)
    #points(keep$list.par[kk], sqrt(w^2+keep$list.par[kk]^2)/keep$list.par[kk]*0.2, type='l', col='red', lwd=2.0)
    plot(keep$list.par[kk], keep$par2[kk], type='b', log=log, cex=0.5)
    #points(keep$list.par[kk], 16-atan2(w, keep$list.par[kk])/2/pi*24, type='l', col='red', lwd=2.0)
    plot(keep$list.par[kk], keep$par4[kk], type='b', log='', cex=0.5) 
  }
  
}

Identifiablity.analysis.M3.with.plot.result = function(T, gene.index = 389, model=3, debug = FALSE, zt = seq(0,94,by = 2), i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE)
{
  i = gene.index; w = 2*pi/24;
  # i = j; gene.index=j;zt =  seq(0,94,by = 2); i.ex = ZT.ex; i.int = ZT.int;absolute.signal = TRUE; Nfit=NA; debug = TRUE; model = 3;outliers = TRUE; 
  #norm.params = TRUE;
  #model  = T$BIC.best.model[i];
  if(model!=3)
  {
    cat('ONLY for Model 3');
  }else{
    R.m = unlist(T[i, i.ex]) ## nb of reads for exon
    R.s = unlist(T[i, i.int]) ## nb of reads for intron
    L.m = T$length.mRNA[i];
    L.s = T$length.premRNA[i];
    M = norm.RPKM(R.m, L.m)
    S = norm.RPKM(R.s, L.s)
    #alpha.m = T$alpha.mRNA[i];
    #alpha.s = T$alpha.premRNA[i];
    alpha.m = rep(as.numeric(T[i, grep('alpha.mRNA.ZT', colnames(T))]), 4);
    alpha.s = rep(as.numeric(T[i, grep('alpha.premRNA.ZT', colnames(T))]), 4);
    
    if(outliers)
    {
      outlier.m = as.numeric(unlist(strsplit(as.character(T$outlier.m[i]), ';')))
      outlier.s = as.numeric(unlist(strsplit(as.character(T$outlier.s[i]), ';'))) 
    }else{
      outlier.m = c();
      outlier.s = c();
    }
    
    bounds = set.bounds(model = model);
    upper = bounds$upper; 
    lower = bounds$lower;
    
    param.fit = T[gene.index, grep('.m3', colnames(T))];
    error.fit = param.fit[1];
    res.fit = param.fit[c(2:6)];
    
    ptm <- proc.time()
    Identifiablity.analysis.M3.M4.Example(error.opt = error.fit, params = res.fit, lower = lower, upper = upper, 
                                          R.m = R.m, R.s = R.s, L.m=L.m, L.s = L.s, alpha.m=alpha.m, alpha.s=alpha.s, 
                                          outlier.m = outlier.m, outlier.s = outlier.s, model = model, zt = zt, PLOT = TRUE, gene.name=T$gene[gene.index]); 
    proc.time() - ptm
    
  }
}


Identifiablity.analysis.4gamma.M3.with.plot.result = function(T, gene.index = 389, model=3, i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE, nb.test.iden = 20,
                                                              folder='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/', PLOT.PL=TRUE )
{
  #gene.index = kk[2];outliers = FALSE; nb.test.iden = 50; 
  #i.ex = ZT.ex; i.int = ZT.int;folder='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/analysis_non_identifiability/';
  i = gene.index; w = 2*pi/24;
  # i = j; gene.index=j;zt =  seq(0,94,by = 2); i.ex = ZT.ex; i.int = ZT.int;absolute.signal = TRUE; Nfit=NA; debug = TRUE; model = 3;outliers = TRUE; 
  #norm.params = TRUE;
  #model  = T$BIC.best.model[i];
  gene.name = T$gene[gene.index];
  if(model!=3)
  {
    cat('ONLY for Model 3');
  }else{
    cat('Identifiability analysis for Model ', model, '\n');
    R.m = unlist(T[i, i.ex]) ## nb of reads for exon
    R.s = unlist(T[i, i.int]) ## nb of reads for intron
    L.m = T$length.mRNA[i];
    L.s = T$length.premRNA[i];
    M = norm.RPKM(R.m, L.m)
    S = norm.RPKM(R.s, L.s)
    a = mean(M)/mean(S) # ratio between splicing rate and degratation rate
    #alpha.m = T$alpha.mRNA[i];
    #alpha.s = T$alpha.premRNA[i];
    alpha.m = rep(as.numeric(T[i, grep('alpha.mRNA.ZT', colnames(T))]), 4);
    alpha.s = rep(as.numeric(T[i, grep('alpha.premRNA.ZT', colnames(T))]), 4);
    if(outliers){
      outlier.m = as.numeric(unlist(strsplit(as.character(T$outlier.m[i]), ';')))
      outlier.s = as.numeric(unlist(strsplit(as.character(T$outlier.s[i]), ';'))) 
    }else{
      outlier.m = c();
      outlier.s = c();
    }
    bounds = set.bounds(model = model);
    upper = bounds$upper; 
    lower = bounds$lower;
    #error.fit = param.fits.results[which(names(param.fits.results)==paste('error.m', model, sep=''))];
    #res.fit = param.fits.results[grep(paste('.m', model, sep=''), names(param.fits.results))];
    #res.fit = res.fit[-1];
    param.fit = T[gene.index, grep('.m3', colnames(T))];
    error.fit = param.fit[1];
    res.fit = param.fit[c(2:6)];
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
    
    ##################
    ## computation of PL
    ##################
    error.opt=error.fit;params = res.fit;
    gamma.opt = params[1];
    gamma.true = T$gamma[gene.index]
    #gamma.opt = log(2)/(24)
    hls = lseq(10/60, 16, length.out = nb.test.iden)
    list.par = log(2)/hls;
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
      opt.pl = optim(params, par.fixed = profiles[n, 1], index.fixed = 1, f2min.proflie, 
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
    threshold3 = qchisq(p=0.68, df=1, lower.tail = TRUE) + error.opt;
    #qchisq(p=0.68, df=1, lower.tail = TRUE)
    if(PLOT.PL)
    {
      pdf.name = paste(folder, 'non_identifiablilty_analysis', gene.name, '_m', model, '.pdf', sep=''); 
      pdf(pdf.name, width = 1.8, height = 1.5);
      par(cex = 0.7, las = 1, mgp = c(2.0,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3);
      
      lims = range(c(error.opt, threshold2, (threshold2+1)));
      cex.axis = 0.7;
      plot(log(2)/PL$gamma, PL$pls, type='l', col='darkblue',lwd=1., xlab=NA, ylab=NA, ylim = lims, main = NA, log='x', axes = FALSE); 
      #points(log(2)/PL$gamma, PL$pls, type='p', cex=0.7, col='darkblue'); 
      points(log(2)/params[1], error.opt, col='red', cex=1.2, pch=8); 
      abline(h=threshold1, col='red', lwd=1., lty=3);
      abline(h=threshold2, col='darkgray', lwd=0.8, lty=2);abline(h=threshold3, col='darkgray', lwd=0.8, lty=2)
      abline(v=log(2)/gamma.true, lwd=1.2, col='green')
      #text(3, threshold3+1.5, paste('eps : ', signif(T$eps.gamma[i], d=2), ' -> ', signif(T$eps.gamma.m3[i], d=2), sep=''), cex=0.6)
      #text(3, threshold3+1., paste('phase : ', signif(T$phase.gamma[i], d=2), 'h -> ', signif(T$phase.gamma.m3[i], d=2), 'h', sep=''), cex=0.6)
      axis(1, at=c(0.2, 0.5, 2, 10), cex.axis = cex.axis)
      axis(2, las=1,cex.axis = cex.axis)
      box()
      dev.off()
    }
  }
}

