##########################################################################
##########################################################################
## Project:
## Script purpose: transform parameters back 
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Jun 18 15:37:17 2018
##########################################################################
##########################################################################
transform.parameter.combinations.cleaning = function(param.fits.results, res.model.sel, res.nonident.analysis.gamma.all.models)
{
  param.t = transform.parameter.combinations(param.fits.results, res.model.sel, res.nonident.analysis.gamma.all.models);
  param.tc = parameter.cleaning(param.t);
  
  return(param.tc)
  
}

transform.parameter.combinations = function(param.fits.results, res.model.sel, res.nonident.analysis.gamma.all.models)
{
  w = 2*pi/24;
  cat("\t transforming back parameter combination\n")
  
  param.t = rep(NA, 20)
  names(param.t) = c('best.model', 'prob.best.model', 'non.identifiability.gamma.L', 'non.identifiability.gamma.R',
                     'half.life', 'cv.half.life', 'splicing.time', 'cv.splicing.time',
                     'eps.gamma', 'cv.eps.gamma', 'phase.gamma', 'sem.phase.gamma', 
                     'Min.int', 'cv.Min.int','Amp.int', 'cv.Amp.int', 'phase.int', 'sem.phase.int', 'beta.int', 'cv.beta.int')
  
  parameter.list = c('Min.int','Amp.int','phase.int','beta.int', 'gamma','eps.gamma','phase.gamma','splicing.k')
  models.p = cbind(c(0,2,3,4), c(0,2,0,4), c(0,2,0,4),c(0,2,0,4), c(0,2,3,4),c(0,0,3,4),c(0,0,3,4),c(0,2,3,4))
  
  ### extract parameter estimation from the fitting result
  param.t["best.model"] = res.model.sel["BIC.best.model"]
  param.t["prob.best.model"] = res.model.sel["BIC.prob.best.model"]
  
  model = param.t["best.model"]

  ## analysis of non-identifiability
  eval(parse(text = paste('param.t["non.identifiability.gamma.L"] = res.nonident.analysis.gamma.all.models["non.identifiability.gamma.L.m', 
                          model, '"]', sep='')));
  eval(parse(text = paste('param.t["non.identifiability.gamma.R"] = res.nonident.analysis.gamma.all.models["non.identifiability.gamma.R.m', 
                          model, '"]', sep=''))); 
  
  param.t = data.frame(t(param.t))
  param.T = data.frame(t(param.fits.results))
  
  if(model>1){
    ### 8 parameters
    for(n.p in 1:length(parameter.list))
    {
      if(models.p[model, n.p]>0)
      {
        # n.p = 8;
        par = parameter.list[n.p];
        
        if(n.p<=4){
          eval(parse(text = paste('param.t$', par, ' = param.T$', par, '.m', model, sep='')));
          if(n.p!=3){
            eval(parse(text = paste('param.t$cv.', par, ' = as.numeric(param.T$', par, '.stderr.m', model, ')/param.T$', par, '.m', 
                                    model, sep='')));
          }else{
            ## standard error for estimated pre-mRNA phase
            eval(parse(text = paste('param.t$sem.', par, ' = param.T$', par, '.stderr.m', model, sep = '')));
          }
        }else{
          gamma = c(); delta.gamma=c(); 
          eval(parse(text = paste('gamma = param.T$gamma.m', model, sep='')));
          eval(parse(text = paste('delta.gamma = param.T$gamma.stderr.m', model, sep='')));
          if(n.p==5){param.t$half.life = log(2)/gamma; param.t$cv.half.life = delta.gamma/gamma;}
          if(n.p==6){
            eval(parse(text = paste('epsilon = param.T$', par, '.m', model, sep='')));
            eval(parse(text = paste('delta.epsilon = param.T$', par, '.stderr.m', model, sep='')));
            Y = epsilon*sqrt(w^2/gamma^2+1); 
            delta.Y = sqrt((sqrt(w^2/gamma^2+1)*delta.epsilon)^2 + (epsilon*w^2/gamma^3*delta.gamma/sqrt(w^2/gamma^2+1))^2);
            param.t$eps.gamma = Y; param.t$cv.eps.gamma = delta.Y/Y;
          }
          if(n.p==7){
            eval(parse(text = paste('phi = param.T$', par, '.m', model, sep='')));
            eval(parse(text = paste('delta.phi = param.T$', par, '.stderr.m', model, '', sep='')));
            ## standard error for estimation of phase.gamma
            Y = (phi - atan2(w, gamma)/w)%%24; 
            delta.Y = sqrt((delta.phi)^2 + (delta.gamma/(gamma^2+w^2))^2);
            param.t$phase.gamma = Y; 
            param.t$sem.phase.gamma = delta.Y;
          }
          if(n.p==8){
            eval(parse(text = paste('a = param.T$', par, '.m', model, sep='')));
            eval(parse(text = paste('delta.a = param.T$', par, '.stderr.m', model, sep='')));
            Y = log(2)/a/gamma*60; delta.Y = log(2)*60*sqrt((delta.gamma/a/gamma^2)^2 + (delta.a/a^2/gamma)^2)
            param.t$splicing.time = Y; param.t$cv.splicing.time = delta.Y/Y;
          }
        }
      }
    } 
  }
  
  return(unlist(param.t))
  
}

####################
## original code for transforming estimated parameters and making clean table for figures  
####################
parameter.cleaning = function(param.t)
{
  ## Filter and clean the estimated parameters
  ## not sure to use it or not 
  param.tc = param.t;
  Filter.keep = FALSE
  if(Filter.keep)
  {
    keep.old = keep;
    
    source('functions.R')
    keep = T;
    
    ###
    ### select genes with probability of best model large than 0.5
    kk = which(keep$prob.best.model>0.5)
    sels =  keep[kk,]
    
    kk = which(sels$best.model==2 | sels$best.model==4);length(kk)
    cutoff = 0.4
    jj1 = kk[which(sels$cv.Min.int[kk]<cutoff & sels$cv.Amp.int[kk]<cutoff & sels$sem.phase.int[kk]<1. & sels$cv.beta.int[kk]<1)]; length(jj1);
    jj2 = kk[which((sels$cv.Min.int[kk]<cutoff|is.na(sels$cv.Min.int[kk])==TRUE) 
                   & (sels$cv.Amp.int[kk]<cutoff |is.na(sels$cv.Amp.int[kk])==TRUE)
                   & (sels$sem.phase.int[kk]<1 | is.na(sels$sem.phase.int[kk])==TRUE)
                   & (sels$cv.beta.int[kk]<1 | is.na(sels$cv.beta.int[kk])==TRUE))]; 
    length(jj2);
    length(jj2)/length(kk)
    ii = which(sels$best.model==3);length(ii);
    ii1 = ii[which(sels$cv.Min.int[ii]<cutoff)];length(ii1);
    ii2 = ii[which(sels$cv.Min.int[ii]<cutoff | is.na(sels$cv.Min.int[ii])==TRUE)];length(ii2);
    
    sels = sels[c(jj2, ii2), ]
    
    ## select genes with identifiable parameter gamma
    threshold2 = qchisq(p=0.95, df=1, lower.tail = TRUE)
    threshold3 = qchisq(p=0.68, df=1, lower.tail = TRUE)
    threshold = threshold3;
    kk = which(sels$best.model==4);length(kk);
    
    length(which(sels$non.identifiability.gamma.L[kk]>threshold & sels$non.identifiability.gamma.R[kk]>threshold));
    kk = which(sels$best.model==3);length(kk);
    length(which(sels$non.identifiability.gamma.L[kk]>threshold & sels$non.identifiability.gamma.R[kk]>threshold));
    kk = which(sels$best.model==2);length(kk);
    length(which(sels$non.identifiability.gamma.L[kk]>threshold & sels$non.identifiability.gamma.R[kk]>threshold));
    
    xx = sels[which(sels$gamma.identifiability.L & sels$gamma.identifiability.R), ]
    
    cutoff = 1
    kk = which((xx$cv.half.life<cutoff | is.na(xx$cv.half.life))
               & (xx$cv.eps.gamma<cutoff | is.na(xx$cv.eps.gamma))
               & (xx$cv.splicing.time<cutoff | is.na(xx$cv.splicing.time))
               & (xx$sem.phase.gamma<2 | is.na(xx$sem.phase.gamma)))
    xx = xx[kk, ]
    
    w = 2*pi/24;
    yy = sels[which( sels$gamma.identifiability.L & !sels$gamma.identifiability.R
                     & sels$half.life<(log(2)/(5*w))), ]
    kk = which(yy$best.model==3)
    test = data.frame(T$eps.gamma.m3[yy$index[kk]], yy$eps.gamma[kk], T$phase.gamma.m3[yy$index[kk]], yy$phase.gamma[kk])
    yy$cv.half.life = 100;
    yy$cv.splicing.time = 100;
    yy$cv.eps.gamma = NA;
    yy$sem.phase.gamma = NA;
    
    keep = rbind(xx, yy)
    
    o1 = order(keep$best.model)
    keep = keep[o1,]
  }
  
  cat("\t parameter cleaning or filtering (but NOTHING done here)\n")
  return(param.tc)
  
  
}

