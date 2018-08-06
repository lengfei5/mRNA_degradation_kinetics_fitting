##########################################################################
##########################################################################
## Project:
## Script purpose:
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Jun  4 13:46:18 2018
##########################################################################
##########################################################################
source("R/error_functions.R", local = TRUE)
source("R/set_bounds_initialize_values.R", local = TRUE)

Identifiablity.analysis.gamma.all.models = function(param.fits.results, GeneDataSet)
{
  res.ident.all = c();
  
  for(model in c(2:4))
  {
    cat('\t model ', model, '\n');
    
    id.gamma.spec.model = Identifiablity.analysis.gamma.each.model(param.opt = param.fits.results, GeneDataSet=GeneDataSet, model = model)
    names(id.gamma.spec.model) = paste(names(id.gamma.spec.model), '.m', model, sep='');
    
    res.iden.all = c(res.iden.all, id.gamma.spec.model);
  }
  
  return(res.iden.all)
  
}


###############################
# model-wise in case checking gamma identifiability for particular model 
###############################
Identifiablity.analysis.gamma.each.model = function(param.opt, GeneDataSet, model = 4,
                                                    quick.check.boundary = TRUE, profile.plotting = FALSE) 
{
  # param.opt = param.fits.results; model = 4; quick.check.boundary = TRUE;
  if(quick.check.boundary){ nb.profile = 2;}
  if(profile.plotting) {nb.profile = 10; }
  
  ## specify the boundaries for parameters in the same way as optimization_params.R
  bounds.gene = set.bounds.gene(GeneDataSet = GeneDataSet, model);
  upper = bounds.gene$upper; 
  lower = bounds.gene$lower;
  
  list.gamma = lseq(lower[which(names(lower)=="gamma")], 
                    upper[which(names(upper)=="gamma")], length=nb.profile);
  #list.par = unique(c(list.par, gamma.opt))
  #list.par = list.par[order(list.par)];
  #ii = which(list.par==gamma.opt);
  
  profiles = matrix(NA, nrow =length(list.gamma) , ncol=2);
  colnames(profiles) = c('gamma', 'pls');
  profiles[,1] = list.gamma;
  
  ## extracting the optimal parameters and error for specific model
  index.error = which(names(param.opt)==paste('error.m', model, sep=''));
  index.param = setdiff(grep(paste('.m', model, sep=''), names(param.opt)), c(index.error, grep("stderr", names(param.opt)))) 
  error.fit = param.opt[index.error];
  param.fit = param.opt[index.param];
  
  
  ptm <- proc.time()
  for(n in 1:nrow(profiles))
  {
    if(profiles[n, 1] == param.fit[1]){
      profiles[n, 2] = error.fit;
    }else{
      param.start = param.fit; 
      param.start[1] = profiles[n, 1];
      #opt = optim(par.init, f2min, GeneDataSet = GeneDataSet, model = model, debug = debug, norm.params=norm.params,
      #            method = 'L-BFGS-B', lower = bounds.g$lower, upper = bounds.g$upper, 
      #            control = list(maxit=500,trace=0), hessian = FALSE)
      opt.pl = optim(param.start, f2min.profile, GeneDataSet = GeneDataSet, which.to.fix =1,  model = model, 
                    method = 'L-BFGS-B', lower = lower, upper = upper)
      profiles[n, 2] = opt.pl$value;
    }
    
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
  
  if(profile.plotting){ ## NOT Used here
    #par(cex = 0.7, las = 1, mgp = c(2.0,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3);
    lims = range(c(error.opt, threshold1, threshold2, (threshold1+5)));
    plot(PL$gamma, PL$pls, type='l', col='darkblue', log='', lwd=1.5, xlab='gamma', ylab='-2log(PL)', ylim = lims, main = paste('m', model, sep='')); 
    points(PL$gamma, PL$pls, type='p', cex=0.8, col='darkblue'); 
    points(params[1], error.opt, col='darkred', cex=1.5, pch=8); 
    abline(h=threshold1, col='red', lwd=1.5, lty=2);abline(h=threshold2, col='red', lwd=1.5, lty=1);abline(h=threshold3, col='red', lwd=1.5, lty=1)  
    
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
