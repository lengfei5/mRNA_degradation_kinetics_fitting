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
    
    res.ident.all = c(res.ident.all, id.gamma.spec.model);
  }
  
  return(res.ident.all)
  
}


###############################
# model-wise in case checking gamma identifiability for particular model 
###############################
Identifiablity.analysis.gamma.each.model = function(param.opt, GeneDataSet, model = 4,
                                                    quick.check.boundary = TRUE, profile.plotting = FALSE) 
{
  # param.opt = param.fits.results; model = 4; quick.check.boundary = TRUE; profile.plotting = FALSE
  if(quick.check.boundary){ nb.profile = 2;}
  if(profile.plotting) {nb.profile = 10; }
  
  ## specify the boundaries for parameters in the same way as optimization_params.R
  bounds.gene = set.bounds.gene(GeneDataSet = GeneDataSet, model);
  upper = bounds.gene$upper; 
  lower = bounds.gene$lower;
  
  ###############################
  # gamma sampling :
  # for the quick.check.boundary option, we just check the left and right boundary
  # for the profile.plotting option, we sampled 10 values for gamma
  # To-do: maybe improvement can be done 
  ###############################
  list.gamma = lseq(lower[which(names(lower)=="gamma")], 
                    upper[which(names(upper)=="gamma")], length=nb.profile);
    
  profiles = matrix(NA, nrow =length(list.gamma) , ncol=2);
  colnames(profiles) = c('gamma', 'pls');
  profiles[,1] = list.gamma;
  
  ## extracting the optimal parameters and error for specific model
  index.error = which(names(param.opt)==paste('error.m', model, sep=''));
  index.param = setdiff(grep(paste('.m', model, sep=''), names(param.opt)), c(index.error, grep("stderr", names(param.opt)))) 
  error.fit = param.opt[index.error];
  param.fit = param.opt[index.param];
  
  which.to.fix = 1
  
  ptm <- proc.time()
  for(n in 1:nrow(profiles))
  {
    if(profiles[n, 1] == param.fit[1]){
      profiles[n, 2] = error.fit;
    }else{
      param.start = param.fit[-which.to.fix];
      #param.start[1] = profiles[n, 1];
      #opt = optim(par.init, f2min, GeneDataSet = GeneDataSet, model = model, debug = debug, norm.params=norm.params,
      #            method = 'L-BFGS-B', lower = bounds.g$lower, upper = bounds.g$upper, 
      #            control = list(maxit=500,trace=0), hessian = FALSE)
      opt.pl = optim(param.start, f2min.profile, GeneDataSet = GeneDataSet, which.to.fix = 1, value.to.fix = profiles[n, 1], model = model, 
                    method = 'L-BFGS-B', lower = lower, upper = upper)
      profiles[n, 2] = opt.pl$value;
    }
    
  }
  proc.time() - ptm;
  
  #### plot the identifiablity analysis (cf. Jens Timmer 2009)
  PL = data.frame(profiles, stringsAsFactors = FALSE)
  
  if(profile.plotting){ ## NOT Used here
    thresholds = calculate.zscore.thresholds() + error.fit;
    #par(cex = 0.7, las = 1, mgp = c(2.0,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3);
    lims = range(c(error.fit, threshold1, threshold2, (threshold1+5)));
    plot(PL$gamma, PL$pls, type='l', col='darkblue', log='', lwd=1.5, xlab='gamma', ylab='-2log(PL)', ylim = lims, main = paste('m', model, sep='')); 
    points(PL$gamma, PL$pls, type='p', cex=0.8, col='darkblue'); 
    points(params[1], error.fit, col='darkred', cex=1.5, pch=8); 
    abline(h=thresholds, col='red', lwd=1.5, lty=2)
    
  }
  
  #z.left = c(); z.right = c();
  z.left = max(PL$pls[which(PL$gamma<=param.fit[which.to.fix])], na.rm = TRUE) - error.fit;
  z.right = max(PL$pls[which(PL$gamma>=param.fit[which.to.fix])], na.rm = TRUE) - error.fit;
  zscores = c(z.left, z.right);
  names(zscores) = c('gamma.profile.L', 'gamma.profile.R');
  z.test = c(evaluate.id.zscore(z.left), evaluate.id.zscore(z.right))
  names(z.test) = c('gamma.identifiability.L', 'gamma.identifiability.R');
  return(c(zscores, z.test));
}

calculate.zscore.thresholds = function()
{
  #if(model==2) df = 6; if(model==3) df = 5; if(model==4) df = 8;
  #df = nb.param[model]
  threshold1 = qchisq(p=0.68, df=1, lower.tail = TRUE)  ## the following threshod is used for the current code
  threshold2 = qchisq(p=0.95, df=1, lower.tail = TRUE) 
  #threshold3 = qchisq(p=0.95, df=df, lower.tail = TRUE)
  
  return(c(threshold1, threshold2))
  
}

evaluate.id.zscore = function(z)
{
  if(z<=0) { return("optimization.errror");
  }else{
    cutoffs = calculate.zscore.thresholds();
    if(z<cutoffs[1]) return("nonidentifiable")
    if(z>=cutoffs[1] & z<cutoffs[2]) return("pass")
    if(z>cutoffs[2]) return("good")
  }
}

