##########################################################################
##########################################################################
## Project:
## Script purpose:
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Jun  4 13:49:17 2018
##########################################################################
##########################################################################
#################################################################################################################
#################################################################################################################
###########  Main Function for model selection and fitting
#################################################################################################################
#################################################################################################################
### Model Selection for one gene
my.model.selection.one.gene.loglike = function(param.fits.results, nb.data=96, method = 'BIC', absolute.signal=TRUE)
{
  ## method = c('BIC', 'AIC', 'AICc');
  set.nb.param(absolute.signal=TRUE);
  index = match(c('error.m1', 'error.m2', 'error.m3', 'error.m4'), names(param.fits.results))
  error.m1 = param.fits.results[index[1]]
  error.m2 = param.fits.results[index[2]]
  error.m3 = param.fits.results[index[3]]
  error.m4 = param.fits.results[index[4]]
  
  ## the formula of BIC used here is chi-square+k*ln(n)==error/sigma^2+k(ln(n)) in which sigma of noise is supported to be known.
  if(method == 'BIC')
  {
    score.m1 = (error.m1) + log(nb.data)*n.param[1];
    score.m2 = (error.m2) + log(nb.data)*n.param[2];
    score.m3 = (error.m3) + log(nb.data)*n.param[3];
    score.m4 = (error.m4) + log(nb.data)*n.param[4];
  }
  if(method == 'AIC'|method == 'AICc')
  {
    score.m1 = (error.m1) + 2*n.param[1] 
    score.m2 = (error.m2) + 2*n.param[2]
    score.m3 = (error.m3) + 2*n.param[3] 
    score.m4 = (error.m4) + 2*n.param[4] 
    if(method== 'AICc')
    {
      score.m1 = score.m1 + 2*n.param[1]*(n.param[1]+1)/(nb.data-n.param[1]-1)
      score.m2 = score.m2 + 2*n.param[2]*(n.param[2]+1)/(nb.data-n.param[2]-1)
      score.m3 = score.m3 + 2*n.param[3]*(n.param[3]+1)/(nb.data-n.param[3]-1)
      score.m4 = score.m4 + 2*n.param[4]*(n.param[4]+1)/(nb.data-n.param[4]-1)
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

#########
#### fitting all models for one genes and remove outliers with iterations
#########
make.fits.with.all.models.for.one.gene.remove.outliers.4specific.model = function(T = T, gene.index = j, model=3, debug = TRUE,
                                                                                  zt = zt, i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE)
{
  #make.fit.spec.model = function(T = T, gene.index = 1, model = 1, debug = FALSE, zt = seq(0,46,by = 2), 
  #                               i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE)
  param.fits.results = make.fit.spec.model(T = T, gene.index = gene.index, model = model,debug = debug, zt = zt, 
                                           i.ex = ZT.ex, i.int = ZT.int, outliers = outliers); 
  return(param.fits.results);
}

make.fits.with.all.models.for.one.gene.remove.outliers = function(T = T, gene.index = 1, debug = FALSE, zt = seq(0,94,by = 2), i.ex = ZT.ex, i.int = ZT.int, 
                                                                  outliers = FALSE, Identifiablity.Analysis.by.Profile.Likelihood.gamma = TRUE, PLOT.Ident.analysis = FALSE, 
                                                                  parametrization = c('cosine.beta'), absolute.signal = TRUE)
{
  #gg = 'Per2'; gene.index = which(T$gene==gg); source('functions.R');zt = seq(0,94,by = 2); debug=TRUE; outliers = TRUE;
  gene.name = T$gene[gene.index];
  R.m = T[gene.index, ZT.ex]
  R.s = T[gene.index, ZT.int]
  a.m = rep(as.numeric(T[gene.index, grep('alpha.mRNA.ZT', colnames(T))]), 4);
  a.s = rep(as.numeric(T[gene.index, grep('alpha.premRNA.ZT', colnames(T))]), 4);
  L.m = T$length.mRNA[gene.index];
  L.s = T$length.premRNA[gene.index];
  
  #######
  #### fitting and outlier removal
  #######
  if(!outliers)
  {
    param.fits.results = make.fits.with.all.models.for.one.gene(T = T, gene.index = gene.index, debug = debug, zt = zt, 
                                                                i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE); 
    outlier.m = NA; outlier.s = NA;
  }else{
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
  
  ##### Analyze non-identifiability using Profile Likelihood approache
  if(Identifiablity.Analysis.by.Profile.Likelihood.gamma)
  {
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
  
  return(c(param.fits.results, outlier.m = paste(outlier.m, sep='', collapse = ';'), outlier.s = paste(outlier.s, sep='', collapse = ';')));
}

make.fits.with.all.models.for.one.gene = function(T = T, gene.index = 1, debug = FALSE, zt = seq(0,46,by = 2),
                                                  i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE, parametrization = c('cosine.beta'), absolute.signal = TRUE)
{
  #T = T; gene.index = j; debug = TRUE; parametrization = 'cosine.beta';  zt = zt; i.ex = ZT.ex; i.int = ZT.int; absolute.signal = TRUE
  for(model in 1:4)
  {
    if(debug){cat('\t starting model ',model,'\n');}
    
    param.fit = make.fit.spec.model(T = T, gene.index = gene.index, model = model, debug = debug, zt = zt, i.ex = i.ex, i.int = i.int, outliers=outliers);
    
    if(model == 1)
    {
      Param.fit.for.gene = param.fit
    }else{
      Param.fit.for.gene = c(Param.fit.for.gene, param.fit)
    }
    if(debug){cat('\t model ',model,' finished \n')};
  }
  
  return(Param.fit.for.gene)
}

make.fit.spec.model = function(T = T, gene.index = 1, model = 1, debug = FALSE, zt = seq(0,46,by = 2), 
                               i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE, parametrization = c('cosine.beta'), absolute.signal = TRUE)
{
  param.fit = NA
  ## Model 1 premrna and mran are both flat
  if(model == 1)
  {
    ## gene.index = j; model=1; i.ex = ZT.ex; i.int=ZT.int;
    R.m = unlist(T[gene.index, i.ex]) ## nb of reads for exon
    R.s = unlist(T[gene.index, i.int]) ## nb of reads for intron
    L.m = T$length.mRNA[gene.index];
    L.s = T$length.premRNA[gene.index];
    M = norm.RPKM(R.m, L.m)
    S = norm.RPKM(R.s, L.s)
    alpha.m = rep(as.numeric(T[gene.index, grep('alpha.mRNA.ZT', colnames(T))]), 4);
    alpha.s = rep(as.numeric(T[gene.index, grep('alpha.premRNA.ZT', colnames(T))]), 4);
    
    if(outliers)
    {
      outlier.m = as.numeric(unlist(strsplit(as.character(T$mRNA.outlier[gene.index]), ';')))
      outlier.s = as.numeric(unlist(strsplit(as.character(T$premRNA.outlier[gene.index]), ';'))) 
    }else{
      outlier.m = c();
      outlier.s = c();
    }
    
    mu.m = convert.nb.reads(rep(mean(M), length(R.m)), L.m);
    mu.s = convert.nb.reads(rep(mean(S), length(R.s)), L.s);
    
    err = NB.error(R.m = R.m, R.s = R.s, alpha.m = alpha.m, alpha.s = alpha.s, 
                   mu.m = mu.m, mu.s = mu.s, outlier.m = outlier.m, outlier.s=outlier.s, specie = 'both')
    param.fit = err;
    names(param.fit) = 'error.m1';
  }
  ## parameter estimations for model 2,3,4
  if(model > 1)
  {
    param.fit = make.optimization(T = T, i = gene.index, model = model, Nfit = NA, debug = debug, zt = zt, i.ex = i.ex, i.int = i.int, outliers = outliers)
  }
  return(param.fit)
}

###############
#### Function to estimat parameters for model 2, 3, 4
###############
norm.RPKM = function(nb.reads, length)
{
  set.scaling.factors();
  return(nb.reads/length/scaling.factors*10^9);
}
convert.nb.reads = function(rpkm, length)
{
  set.scaling.factors();
  return(rpkm*length*scaling.factors/10^9);
}

norm.RPKM.libary.size = function(nb.reads, length)
{
  load(file='Libary_size_48_samples.Rdata')
  return(nb.reads/length/ss*10^9);
}
convert.nb.reads.libary.size = function(rpkm, length)
{
  load(file='Libary_size_48_samples.Rdata')
  return(rpkm*length*ss/10^9);
}

Gamma.Initiation = function(eps.gamma.init, min.half.life=0.5, max.half.life=6, w=2*pi/24)
{
  # min.half.life=0.5; max.half.life=6; w=2*pi/24
  eps.gamma.init = as.numeric(eps.gamma.init);
  gamma.init = rep(log(2)/(10/60), length(eps.gamma.init))
  lss = lseq(log(2)/max.half.life, log(2)/min.half.life, length=100)
  for(n in 1:length(gamma.init))
  {
    kk = which(lss>=sqrt(w^2/(1/eps.gamma.init[n]^2-1)));
    if(length(kk)>0) gamma.init[n] = sample(lss[kk], 1);
  }
  return(gamma.init)
  #gamma.init = c(rep(log(2)/lseq(max.half.life, min.half.life, length = Nfit.M), Nfit.M%/%Nfit.M), rep(log(2)/5,Nfit.M%%Nfit.M))
}
sigmoid.bound.contraint = function(eps.gamma)
{
  return(10^4/(1+exp(-100*(eps.gamma-1.0))));
  
  smooth.bound.constraint.function = FALSE
  if(smooth.bound.constraint.function)
  {
    xx = lseq(0.001, 1.2, length.out = 1000)
    x0 = 1.0;y0=10^4
    yy = y0/(1+exp(-100*(xx-x0)))
    plot(xx, yy, type='l', col='blue', log='', ylim=c(0, y0));abline(h=1, col='red');abline(h=0, col='red');
    abline(v=x0, col='black');abline(v=(x0-0.05), col='black');abline(v=1, col='black');
  }
}

make.optimization = function(T = T, i = 1, model = 4, Nfit = NA, debug = FALSE, zt =  seq(0,46,by = 2), 
                             i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE,parametrization =c('cosine.beta'), norm.params = TRUE, absolute.signal = TRUE)
{
  # i = j; zt =  seq(0,94,by = 2); i.ex = ZT.ex; i.int = ZT.int;absolute.signal = TRUE; Nfit=NA; debug = TRUE; model = 3;outliers = TRUE; 
  #norm.params = TRUE;
  w = 2*pi/24;
  gene2opt = T$gene[i];
  param.fit = NA
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
    outlier.m = as.numeric(unlist(strsplit(as.character(T$mRNA.outlier[i]), ';')))
    outlier.s = as.numeric(unlist(strsplit(as.character(T$premRNA.outlier[i]), ';'))) 
  }else{
    outlier.m = c();
    outlier.s = c();
  }
  
  fitting.factor = 1; ### define the nb of initial values in optimization
  bounds = set.bounds(model = model);
  upper = bounds$upper; 
  lower = bounds$lower;
  
  #######
  ####### FITTING the absolute signal
  #######
  ### Here we want first to fit the pre-mRNA profile to identify parameters mean.int, fold.change.int, phase.int, beta.int
  ## Prefit premRNA for model2 and model4 in order to have the good initial values for premRNA parameters
  #ptm = proc.time();
  prefit.S = TRUE 
  if((model==2|model==4) & prefit.S) 
  {
    Nfit.S = 4*fitting.factor;
    #index = i
    
    ### Choose different initial values of parameters for S fitting
    Min.init = rep(min(S), Nfit.S)
    Amp.init = rep((max(S)-min(S)), Nfit.S)
    phase.init = zt[which.max(S)]
    phase.init = (rep(phase.init, Nfit.S)+rnorm(Nfit.S,sd = 3))%%24
    beta.min = 1;
    beta.max = 5;
    beta.init = lseq(beta.min, beta.max, length = Nfit.S)
    PAR.INIT.S = cbind(Min.init, Amp.init, phase.init, beta.init)
    
    errors.fit.s = rep(NA, Nfit.S)
    
    limit.factor = 5
    for(fit.nb.s in 1:Nfit.S)
    {
      # fit.nb.s = 1;
      par.init.s = PAR.INIT.S[fit.nb.s,]
      opt.s = optim(par.init.s, f2min.int, R.s=R.s, L.s=L.s, alpha.s=alpha.s, outlier.s=outlier.s, zt = zt, method = 'L-BFGS-B', 
                    lower = c(min(S)/limit.factor, (max(S)-min(S))/limit.factor, 0, 1), upper = c(max(S), (max(S)-min(S))*limit.factor, 24, 5))
      #print(opt.s$convergence)
      res.fit.s = opt.s$par
      errors.fit.s[fit.nb.s] = opt.s$value
      eval(parse(text = paste('res.fit.s.', fit.nb.s, ' = res.fit.s', sep = '')))
    }
    ## choose the best-fitting parameters for S
    imin.s = which.min(errors.fit.s); 
    eval(parse(text = paste('res.fit.s = res.fit.s.', imin.s, sep = '')))			
  }
  #proc.time() - ptm
  
  #ptm = proc.time();
  prefit.M = TRUE
  if(model==4 & prefit.M)
  {
    Nfit.M = 6*fitting.factor;
    
    ### Choose different initial values of parameters for M fitting
    eps.m = min((max(M) - min(M))/mean(M)/2, 1);
    eps.gamma.init = c(sample(seq(0.1, 0.6, length = 10), Nfit.M/2, replace = TRUE), rep(eps.m/1.25, length=Nfit.M/2));
    #if(debug){cat('eps.gamma inital values ... ', eps.gamma.init, '\n')};
    if(length(which(eps.gamma.init>0.8))>=1) eps.gamma.init[which(eps.gamma.init>0.8)] = 0.8;
    phase.m = zt[which.max(M)]
    phase.gamma.init =  (rep((phase.m+12), Nfit.M)+rnorm(Nfit.M, sd = 3))%%24
    #if(debug){cat('eps.gamma inital values ', eps.gamma.init, '\n')};
    gamma.init = Gamma.Initiation(eps.gamma.init, 0.5, 6);
    #if(debug){cat('gamma inital values ', gamma.init, '\n')};
    a = mean(M)/mean(S) # define the ratio between splicing rate and degratation rate
    a.init = rep(a, Nfit.M)
    PAR.INIT.M = cbind(gamma.init, eps.gamma.init, phase.gamma.init, a.init)
    
    errors.fit.m = rep(NA, Nfit.M)
    for(fit.nb.m in 1:Nfit.M)
    {
      # fit.nb.m = 1;
      par.init.m = PAR.INIT.M[fit.nb.m,]
      opt.m = optim(par.init.m, f2min.mrna, res.fit.s=res.fit.s, R.m=R.m, L.m=L.m, alpha.m=alpha.m, outlier.m=outlier.m, zt = zt, 
                    norm.params=norm.params, method = 'L-BFGS-B', lower = c(lower[c(1:3)], a/10), upper = c(upper[c(1:3)], a*5))
      res.fit.m = opt.m$par
      errors.fit.m[fit.nb.m] = opt.m$value
      eval(parse(text = paste('res.fit.m.', fit.nb.m, ' = res.fit.m', sep = '')))
    }
    
    ## choose the best-fitting parameters for S
    imin.m = which.min(errors.fit.m); 
    eval(parse(text = paste('res.fit.m = res.fit.m.', imin.m, sep = '')))
    
    if(debug){cat('prefix parameter for M4', res.fit.m, res.fit.s, '\n')};
    #cat('true parameters for M4', unlist(T[j, c(142, 146, 147, 141, 140, 143, 144, 145)]), '\n')
    #cat('-2loglike = ', errors.fit.m[imin.m], errors.fit.s[imin.s], '\n sum of -2loglike = ', sum(errors.fit.m[imin.m], errors.fit.s[imin.s]), '\n')
  }
  #proc.time() - ptm
  ######
  ### Now fit mRNA and premRNA all together for model2, model3 and model4
  ######
  if(is.na(Nfit))
  {
    if(model==2) Nfit = fitting.factor*10;
    if(model==3) Nfit = fitting.factor*12;
    if(model==4) Nfit = fitting.factor*12;
  }
  
  ### Define the initial values for degradation and splicing parameters
  a = mean(M)/mean(S) # define the ratio between splicing rate and degratation rate
  a.init = lseq(max(0.1, a/5), min(10^5, a*2), length=Nfit)
  max.half.life = 12;
  min.half.life = 30/60;
  if(model==2)
  {
    Min.init = rep(res.fit.s[1], Nfit)
    Amp.init = seq(res.fit.s[2]/2, res.fit.s[2]*2, length=Nfit)
    phase.init = (rep(res.fit.s[3],Nfit)+rnorm(Nfit,sd = 2))%%24
    beta.init = lseq(max(1, (res.fit.s[4]-2)), min((res.fit.s[4]+2), 10), length=Nfit)
    
    gamma.init = c(rep(log(2)/lseq(max.half.life, min.half.life, length = Nfit), Nfit%/%Nfit), rep(log(2)/5,Nfit%%Nfit))
    
    PAR.INIT = cbind(gamma.init, a.init, Min.init, Amp.init, phase.init, beta.init); 
    colnames(PAR.INIT) = c('gamma','splicing.k','Min.int','Amp.int','phase.int','beta.int');
    upper[2] = a*5 
    lower[2] = a/5;
    lower[3] = min(S)/5;
    upper[3] = max(S)*5;
    lower[4] = (max(S)-min(S))/5;
    upper[4] = (max(S)-min(S))*5;
  }
  if(model==3)
  {	
    ### initialize decay parameters for model 3
    #gamma.init = c(rep(log(2)/lseq(6, 0.5, length = Nfit), Nfit%/%Nfit), rep(log(2)/5,Nfit%%Nfit))
    eps.m = min((max(M) - min(M))/mean(M)/2, 1); # close to mRNA amplitude
    eps.gamma.init = c(sample(seq(0.1, 0.6, length = 10), Nfit/2, replace = TRUE), rep(eps.m, length=Nfit/2));
    if(length(which(eps.gamma.init>0.8))>=1) eps.gamma.init[which(eps.gamma.init>0.8)] = 0.8;
    phase.m = zt[which.max(M)]
    phase.gamma.init =  (rep((phase.m+12), Nfit)+rnorm(Nfit, sd = 2))%%24 ## antiphasic to mRNA
    gamma.init = Gamma.Initiation(eps.gamma.init, 0.5, 6);
    
    Min.init = rep(mean(S), Nfit)
    PAR.INIT = cbind(gamma.init, eps.gamma.init, phase.gamma.init, a.init, Min.init)
    colnames(PAR.INIT) = c('gamma','eps.gamma','phase.gamma','splicing.k', 'Min.int');
    upper[4] = a*5;
    lower[4] = a/5;
    lower[5] = min(S)/5;
    upper[5] = max(S)*5;
  }
  if(model==4)
  {
    ### initialize decay parameters for model 4
    #gamma.init = c(seq(max(log(2)/max.half.life, res.fit.m[1]/2), min(log(2)/min.half.life, res.fit.m[1]*2), length=(Nfit-Nfit%/%3)),
    #               log(2)/lseq(10, 1, length = Nfit%/%3));
    eps.gamma.init = seq(res.fit.m[2]/4, min(res.fit.m[2]*4, 1), length=Nfit);
    if(length(which(eps.gamma.init>0.8))>=1) eps.gamma.init[which(eps.gamma.init>0.8)] = 0.8;
    gamma.init = Gamma.Initiation(eps.gamma.init, 0.5, 8);
    phase.gamma.init = (rep(res.fit.m[3],Nfit)+rnorm(Nfit,sd = 1.5))%%24;
    a.init = seq(max(0.1, res.fit.m[4]/1.5), min(10^5, res.fit.m[4]*1.5), length=Nfit)
    
    Min.init = rep(res.fit.s[1], Nfit)
    Amp.init = seq(res.fit.s[2]/2, res.fit.s[2]*2, length=Nfit)
    phase.init = (rep(res.fit.s[3],Nfit)+rnorm(Nfit,sd = 2))%%24
    beta.init = lseq(max(1, (res.fit.s[4]-2)), min((res.fit.s[4]+2), 10), length=Nfit)
    
    PAR.INIT = cbind(gamma.init, eps.gamma.init, phase.gamma.init, a.init, Min.init, Amp.init, phase.init, beta.init);
    colnames(PAR.INIT) = c('gamma','eps.gamma','phase.gamma','splicing.k', 'Min.int','Amp.int','phase.int','beta.int')
    upper[4] = a*5;
    lower[4] = a/5;
    lower[5] = min(S)/5;
    upper[5] = max(S)*5;
    lower[6] = (max(S)-min(S))/5;
    upper[6] = (max(S)-min(S))*5;
  }
  
  #################################
  #### Fit the S and M all together 
  #################################
  #ptm = proc.time();
  errors.fit = rep(NA, Nfit)
  #rank.hess = rep(NA, Nfit)
  #aic.fit = rep(NA, Nfit)
  # fit.number = 1; 
  #lower[1] = log(2)/24;upper[1]=log(2)/(10/60);
  #while(fit.number<(Nfit+1))
  ptm <- proc.time()
  for(fit.number in 1:Nfit)
  {
    ## cat(fit.number, '\n'); alpha.m = 0.002; alpha.s = 0.2;
    par.init = PAR.INIT[fit.number,]
    if(debug){cat('\t \t optimization # ',fit.number,' started : \t',par.init,'\n')}
    
    #ptm <- proc.time()
    opt = optim(par.init, f2min, R.m = R.m, R.s = R.s, L.m=L.m, L.s = L.s, alpha.m=alpha.m, alpha.s=alpha.s, 
                outlier.m = outlier.m, outlier.s = outlier.s, model = model, debug = debug, zt = zt, norm.params=norm.params,
                method = 'L-BFGS-B', lower = lower, upper = upper, control = list(maxit=500,trace=0), hessian = FALSE)
    #cat(proc.time() - ptm);
    ## extract the parameter of optimization results and errors
    errors.fit[fit.number] = opt$value;
    eval(parse(text = paste('res.fit.', fit.number, ' = opt$par', sep = '')))
    #ttry = try(sqrt(diag(solve(0.5*opt$hessian))), silent = TRUE); 
    #if(!inherits(ttry, "try-error")) { opt$stderr = sqrt(diag(solve(0.5*opt$hessian)));}else{ opt$stderr = rep(NA, length(opt$par));}
    #eval(parse(text = paste('res.fit.stderr.', fit.number, ' = opt$stderr', sep = '')))
    #res.fit.stderr = opt$stderr
    #aic.fit[fit.number] = errors.fit[fit.number] + rank.hess[fit.number]*2
    if(debug){cat('\t\t optimization # ',fit.number,' finished : \t', opt$par, '\t',opt$value, '\n')}
  }
  proc.time() - ptm;
  #### First choose the best estimated parameters
  imin = which.min(errors.fit); 
  #imin.test = which.min(aic.fit)
  eval(parse(text = paste('res.fit = res.fit.', imin, sep = '')))
  #eval(parse(text = paste('res.fit.stderr = res.fit.stderr.', imin, sep = '')))
  #if(debug){cat('\t\t select optimal fitting ', res.fit, '\n')};
  
  ### compute the Standard Error of estimates using hessian function
  hess = hessian(func=f2min, R.m = R.m, R.s = R.s, L.m=L.m, L.s = L.s, alpha.m=alpha.m, alpha.s=alpha.s, 
                 outlier.m = outlier.m, outlier.s = outlier.s, model = model, zt = zt, x=res.fit);
  #if(debug){cat('\t\t calculate hessian matrix..... \n')};
  ttry = try(sqrt(diag(solve(0.5*hess))), silent = FALSE);
  if(!inherits(ttry, "try-error")){res.fit.stderr = sqrt(diag(solve(0.5*hess)));}else{res.fit.stderr = rep(NA, length(res.fit))}
  #if(debug){cat('\t\t calculate standard error of estimates..... \n')};
  param.fit = c(errors.fit[imin], res.fit, res.fit.stderr);
  names(param.fit) = paste(c('error', colnames(PAR.INIT), paste(colnames(PAR.INIT), '.stderr', sep='')),'.m',model,sep = ''); 
  
  #####
  ##### identifiability analysis (for M2, M3 AND M4) using rank of hessian matrix or profile likelihood
  #Identifiablity.Analysis.by.Profile.Likelihood = FALSE
  #if(Identifiablity.Analysis.by.Profile.Likelihood)
  #{
  #  if(debug){cat('\t\t start non-identifiability analysis for gamma \n')}
  #### FAST analysis 
  #  ptm <- proc.time()
  #  res.nonident.analysis.gamma = Identifiablity.analysis.gamma.each.model(error.opt=errors.fit[imin], params.opt =res.fit, 
  #                                                                         lower = lower, upper = upper, 
  #                                     R.m = R.m, R.s = R.s, L.m=L.m, L.s = L.s, alpha.m=alpha.m, alpha.s=alpha.s, 
  #                                   outlier.m = outlier.m, outlier.s = outlier.s, model = model, zt = zt, PLOT.PL = TRUE, gene2plot=gene2opt); 
  
  # proc.time() - ptm;
  #param.fit = c(param.fit, res.nonident.analysis.gamma);
  #names(param.fit)[length(param.fit)] = paste('non.identifiability.gamma.m', model, sep = '')
  
  #### Detail analysis with individual examples
  #Detail.analysis = FALSE
  #if(Detail.analysis){
  #ptm <- proc.time()
  #Identifiablity.analysis.M3.M4.Example(error.opt = errors.fit[imin], params = res.fit, lower = lower, upper = upper, 
  #                                     R.m = R.m, R.s = R.s, L.m=L.m, L.s = L.s, alpha.m=alpha.m, alpha.s=alpha.s, 
  #                                     outlier.m = outlier.m, outlier.s = outlier.s, model = model, zt = zt); 
  #proc.time() - ptm
  #}
  #}
  #if(debug){cat('\t\t\t optimization finished\n')};
  #if(debug){cat('\t\t fitting result ', param.fit, '\n')};
  return(param.fit);
}
