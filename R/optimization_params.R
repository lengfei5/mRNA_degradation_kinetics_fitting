##########################################################################
##########################################################################
## Project:
## Script purpose: Main Function for model fitting
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Jun  4 13:49:17 2018
##########################################################################
##########################################################################

## import function dependencies
source("R/error_functions.R", local = TRUE)
source("R/kinetic_model.R", local = TRUE)
source("R/set_bounds_initialize_values.R", local = TRUE)

#########
#### fitting all models for one genes and remove outliers with iterations
#########
make.fits.with.all.models.for.one.gene = function(R.m, R.s, alpha.m, alpha.s, L.m, L.s, 
                                                  zt = seq(0,94,by = 2), 
                                                  debug = FALSE, outliers = FALSE, 
                                                  parametrization = c('cosine.beta'), absolute.signal = TRUE)
{
  #T = T; gene.index = j; debug = TRUE; parametrization = 'cosine.beta';  zt = zt; i.ex = ZT.ex; i.int = ZT.int; absolute.signal = TRUE
  
  Param.fit.for.gene = c();
  for(model in 1:4)
  {
    if(debug){cat('\t starting model ',model,'\n');}
    
    param.fit = make.fit.spec.model(T = T, gene.index = gene.index, model = model, debug = debug, 
                                    zt = zt, i.ex = i.ex, i.int = i.int, outliers=outliers);
    Param.fit.for.gene = c(Param.fit.for.gene, param.fit)
    
    if(debug){cat('\t model ',model,' finished \n')};
  }
  
  return(Param.fit.for.gene)

}

make.fit.spec.model = function(T = T, gene.index = 1, model = 1, debug = FALSE, zt = seq(0,46,by = 2), 
                               i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE, parametrization = c('cosine.beta'), absolute.signal = TRUE)
{
  param.fit = NA;
  
  ## error for model 1
  if(model == 1){
    param.fit = calculate.error.for.flat.model(T = T, 
                                               gene.index = gene.index, 
                                               debug = debug,  
                                               zt = zt, i.ex = i.ex, i.int = i.int, 
                                               outliers = outliers, 
                                               parametrization = parametrization, 
                                               absolute.signal = absolute.signal)
  }
  
  ## parameter estimations for model 2,3,4
  if(model > 1){
    param.fit = make.optimization(T = T, i = gene.index, model = model, Nfit = NA, debug = debug, zt = zt, i.ex = i.ex, i.int = i.int, 
                                  outliers = outliers)
  }
  
  return(param.fit)
  
}

###########
## error for model 1
###########
calculate.error.for.flat.model = function(T = T, 
                                          gene.index = 1, 
                                          debug = FALSE, 
                                          zt = seq(0,46,by = 2), 
                                          i.ex = ZT.ex, 
                                          i.int = ZT.int, 
                                          outliers = FALSE, 
                                          parametrization = c('cosine.beta'), 
                                          absolute.signal = TRUE)
{
  R.m = unlist(T[gene.index, i.ex]) ## nb of reads for exon
  R.s = unlist(T[gene.index, i.int]) ## nb of reads for intron
  L.m = T$length.mRNA[gene.index];
  L.s = T$length.premRNA[gene.index];
 
  alpha.m = rep(as.numeric(T[gene.index, grep('alpha.mRNA.ZT', colnames(T))]), 4);
  alpha.s = rep(as.numeric(T[gene.index, grep('alpha.premRNA.ZT', colnames(T))]), 4);
  
  M = norm.RPKM(R.m, L.m)
  S = norm.RPKM(R.s, L.s)
  
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
  names(err) = 'error.m1';
  
  return(err)
  
}

####################
## main optimization function for Model 2, 3 and 4
####################
make.optimization = function(T = T, 
                             i = 1, 
                             model = 4, 
                             Nfit = NA,
                             zt =  seq(0,94,by = 2), 
                             i.ex = ZT.ex, 
                             i.int = ZT.int, 
                             outliers = FALSE,
                             parametrization =c('cosine.beta'), 
                             debug = FALSE,
                             norm.params = TRUE, 
                             absolute.signal = TRUE)
{
  # i = gene.index; zt =  seq(0,94,by = 2); i.ex = ZT.ex; i.int = ZT.int;absolute.signal = TRUE; Nfit=NA; debug = TRUE; model = 4;outliers = TRUE; 
  
  ####################
  ## prepare parameters for the optimization 
  ####################
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
  alpha.m = rep(as.numeric(T[i, grep('alpha.mRNA.ZT', colnames(T))]), 4);  # dispersion parameter alpha for each time points
  alpha.s = rep(as.numeric(T[i, grep('alpha.premRNA.ZT', colnames(T))]), 4); 
  
  if(outliers){
    outlier.m = as.numeric(unlist(strsplit(as.character(T$mRNA.outlier[i]), ';')))
    outlier.s = as.numeric(unlist(strsplit(as.character(T$premRNA.outlier[i]), ';'))) 
  }else{
    outlier.m = c();
    outlier.s = c();
  }
  
  ####################
  ## Prefix parameters in optimization process
  ####################
  prefit.S = TRUE # optimization warm-up by fitting pre-mRNA individually to have good initial values
  prefit.M = TRUE # optimization warm-up by fitting mRNA individually 
  
  fitting.factor = 1;  
  Nfit.S = 4*fitting.factor; # nb of inital values for pre-mRNA
  Nfit.M = 6*fitting.factor; # nb of initial values for mRNA
  if(is.na(Nfit)) # nb of intial valeus for pre-mRNA and mRNA together
  {
    if(model==2) Nfit = fitting.factor*4;
    if(model==3) Nfit = fitting.factor*6;
    if(model==4) Nfit = fitting.factor*8;
  }
  
  ####################
  ## fit premRNA individually for M2 and M4 with the aim of having good initial values for its parameters 
  ####################
  #ptm = proc.time();
  if((model==2|model==4) & prefit.S) 
  {
    ## set initial values of parameters for S fitting
    PAR.INIT.S = Sampling.Initial.Values.for.fitting.S(S , Nfit.S, zt)
    bounds.g.s = set.bounds.gene.s(S, range_scalingFactor=5)

    errors.fit.s = rep(NA, Nfit.S)
    
    for(fit.nb.s in 1:Nfit.S)
    {
      par.init.s = PAR.INIT.S[fit.nb.s,]
      opt.s = optim(par.init.s, f2min.int, R.s=R.s, L.s=L.s, alpha.s=alpha.s, outlier.s=outlier.s, zt = zt, method = 'L-BFGS-B', 
                    lower = bounds.g.s$lower, upper = bounds.g.s$upper)
      
      res.fit.s = opt.s$par
      errors.fit.s[fit.nb.s] = opt.s$value
      eval(parse(text = paste('res.fit.s.', fit.nb.s, ' = res.fit.s', sep = '')))
    }
    
    ## choose the best-fitting parameters for S
    imin.s = which.min(errors.fit.s); 
    eval(parse(text = paste('res.fit.s = res.fit.s.', imin.s, sep = '')))
  }
  #proc.time() - ptm
  
  ####################
  ## fit mRNA individually for M4 with the aim of having good initial values for its parameters     
  ####################
  #ptm = proc.time();
  if(model==4 & prefit.M)
  {
    ### sampling initial values of parameters for M fitting
    PAR.INIT.M = Sampling.Initial.Values.for.fitting.M(M, S, Nfit, zt = zt)
    bounds.g.m = set.bounds.gene.m(M, S, range_scalingFactor=5)
    
    errors.fit.m = rep(NA, Nfit.M)
    for(fit.nb.m in 1:Nfit.M)
    {
      par.init.m = PAR.INIT.M[fit.nb.m,]
      opt.m = optim(par.init.m, f2min.mrna, res.fit.s=res.fit.s, R.m=R.m, L.m=L.m, alpha.m=alpha.m, outlier.m=outlier.m, zt = zt, 
                    norm.params=norm.params, method = 'L-BFGS-B', lower = bounds.g.m$lower, upper = bounds.g.m$upper)
      res.fit.m = opt.m$par
      errors.fit.m[fit.nb.m] = opt.m$value
      eval(parse(text = paste('res.fit.m.', fit.nb.m, ' = res.fit.m', sep = '')))
    }
    
    ## choose the best-fitting parameters for S
    imin.m = which.min(errors.fit.m); 
    eval(parse(text = paste('res.fit.m = res.fit.m.', imin.m, sep = '')))
    
    if(debug){cat('\t fitting warm-up for M4 ', res.fit.m, res.fit.s, '\n')};
 
  }
  #proc.time() - ptm
  
  #################################
  #### Fit both S (pre-mRNA) and M (mRNA) together for Model 2, 3, 4 
  #################################
  ### sampling initial values for parameters and set boundaries
  if(model == 2) {res.fit.m = NULL;}
  if(model == 3) {res.fit.s = NULL; res.fit.m = NULL}
  
  PAR.INIT = Sampling.Initial.Values.for.fitting.M.S(M, S, model, Nfit, zt, res.fit.s, res.fit.m)
  bounds.g = set.bounds.gene(M, S, model = model)
  
  errors.fit = rep(NA, Nfit)
  
  for(fit.number in 1:Nfit)
  {
    # cat(fit.number, '\n'); alpha.m = 0.002; alpha.s = 0.2;
    par.init = PAR.INIT[fit.number,]
    if(debug){cat('\t \t optimization # ',fit.number,' started : \t',par.init,'\n')}
    
    opt = optim(par.init, f2min, R.m = R.m, R.s = R.s, L.m=L.m, L.s = L.s, alpha.m=alpha.m, alpha.s=alpha.s, 
                outlier.m = outlier.m, outlier.s = outlier.s, model = model, debug = debug, zt = zt, norm.params=norm.params,
                method = 'L-BFGS-B', lower = bounds.g$lower, upper = bounds.g$upper, 
                control = list(maxit=500,trace=0), hessian = FALSE)
    
    # extract the parameter of optimization results and errors
    errors.fit[fit.number] = opt$value;
    eval(parse(text = paste('res.fit.', fit.number, ' = opt$par', sep = '')))
        
    if(debug){cat('\t\t optimization # ',fit.number,' finished : \t', opt$par, '\t',opt$value, '\n')}
  }
    
  ## choose the best estimated parameters
  imin = which.min(errors.fit); 
  eval(parse(text = paste('res.fit = res.fit.', imin, sep = '')))
    
  ## compute the Standard Error of estimates using hessian function
  hess = hessian(func=f2min, R.m = R.m, R.s = R.s, L.m=L.m, L.s = L.s, alpha.m=alpha.m, alpha.s=alpha.s, 
                 outlier.m = outlier.m, outlier.s = outlier.s, model = model, zt = zt, x=res.fit);
  ttry = try(sqrt(diag(solve(0.5*hess))), silent = FALSE);
  
  if(!inherits(ttry, "try-error")){res.fit.stderr = sqrt(diag(solve(0.5*hess)));}else{res.fit.stderr = rep(NA, length(res.fit))}
  
  param.fit = c(errors.fit[imin], res.fit, res.fit.stderr);
  names(param.fit) = paste(c('error', colnames(PAR.INIT), paste(colnames(PAR.INIT), '.stderr', sep='')),'.m',model,sep = ''); 
  
  return(param.fit);
}
