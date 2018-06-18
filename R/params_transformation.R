##########################################################################
##########################################################################
## Project:
## Script purpose: transform parameters back 
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Jun 18 15:37:17 2018
##########################################################################
##########################################################################
transform.parameter.combinations = function(T = T)
{
  kk = which(!is.na(T$BIC.best.model>1)==TRUE)
  keep = matrix(data=NA, nrow=length(kk), ncol=26)
  colnames(keep) = c('gene', 'index', 'best.model', 'prob.best.model', 'non.identifiability.gamma.L', 'non.identifiability.gamma.R',
                     'half.life', 'cv.half.life', 'splicing.time', 'cv.splicing.time',
                     'eps.gamma', 'cv.eps.gamma', 'phase.gamma', 'sem.phase.gamma', 
                     'Min.int', 'cv.Min.int','Amp.int', 'cv.Amp.int', 'phase.int', 'sem.phase.int', 'beta.int', 'cv.beta.int',
                     'Fried', 'Shar', 'Schwa', 'score.params')
  #keep[,4] = T$BIC.prob.[kk]
  keep = data.frame(T$gene[kk], keep[,-1], stringsAsFactors=FALSE);colnames(keep)[1] = 'gene'
  keep$index = kk; keep$best.model = T$BIC.best.model[kk];
  
  load(file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/mRNAs_half_lives_databases.Rdata')
  parameter.list = c('Min.int','Amp.int','phase.int','beta.int', 'gamma','eps.gamma','phase.gamma','splicing.k')
  models.p = cbind(c(0,2,3,4), c(0,2,0,4), c(0,2,0,4),c(0,2,0,4), c(0,2,3,4),c(0,0,3,4),c(0,0,3,4),c(0,2,3,4))
  
  ### extract parameter estimation from the fitting result
  for(n in 1:nrow(keep))
    #for(n in c(1:10))
  {
    #n = 1
    cat(n, '\n')
    model = keep$best.model[n]; index = keep$index[n];
    eval(parse(text = paste('keep$prob.best.model[n] = T$BIC.prob.m', model, '[index]', sep=''))); ## prob of best model
    
    if(model>1){
      ## analysis of non-identifiability
      eval(parse(text = paste('keep$non.identifiability.gamma.L[n] = T$non.identifiability.gamma.L.m', model, '[index]', sep=''))); 
      eval(parse(text = paste('keep$non.identifiability.gamma.R[n] = T$non.identifiability.gamma.R.m', model, '[index]', sep='')));
      ### 8 parameters
      w = 2*pi/24;
      for(n.p in 1:length(parameter.list))
      {
        if(models.p[model, n.p]>0)
        {
          # n.p = 8;
          par = parameter.list[n.p];
          
          if(n.p<=4){
            eval(parse(text = paste('keep$', par, '[n] = T$', par, '.m', model, '[index]', sep='')));
            if(n.p!=3){
              eval(parse(text = paste('keep$cv.', par, '[n] = as.numeric(T$', par, '.stderr.m', model, '[index])/T$', par, '.m', model, '[index]', sep='')));
            }else{
              ## standard error for estimated pre-mRNA phase
              eval(parse(text = paste('keep$sem.', par, '[n] = T$', par, '.stderr.m', model, '[index]', sep = '')));
            }
          }else{
            gamma = c(); delta.gamma=c(); eval(parse(text = paste('gamma = T$gamma.m', model, '[index]', sep='')));
            eval(parse(text = paste('delta.gamma = T$gamma.stderr.m', model, '[index]', sep='')));
            if(n.p==5){keep$half.life[n] = log(2)/gamma; keep$cv.half.life[n] = delta.gamma/gamma;}
            if(n.p==6){
              eval(parse(text = paste('epsilon = T$', par, '.m', model, '[index]', sep='')));
              eval(parse(text = paste('delta.epsilon = T$', par, '.stderr.m', model, '[index]', sep='')));
              Y = epsilon*sqrt(w^2/gamma^2+1); 
              delta.Y = sqrt((sqrt(w^2/gamma^2+1)*delta.epsilon)^2 + (epsilon*w^2/gamma^3*delta.gamma/sqrt(w^2/gamma^2+1))^2);
              keep$eps.gamma[n] = Y; keep$cv.eps.gamma[n] = delta.Y/Y;
            }
            if(n.p==7){
              eval(parse(text = paste('phi = T$', par, '.m', model, '[index]', sep='')));
              eval(parse(text = paste('delta.phi = T$', par, '.stderr.m', model, '[index]', sep='')));
              ## standard error for estimation of phase.gamma
              Y = (phi - atan2(w, gamma)/w)%%24; 
              delta.Y = sqrt((delta.phi)^2 + (delta.gamma/(gamma^2+w^2))^2);
              keep$phase.gamma[n] = Y; 
              keep$sem.phase.gamma[n] = delta.Y;
            }
            if(n.p==8){
              eval(parse(text = paste('a = T$', par, '.m', model, '[index]', sep='')));
              eval(parse(text = paste('delta.a = T$', par, '.stderr.m', model, '[index]', sep='')));
              Y = log(2)/a/gamma*60; delta.Y = log(2)*60*sqrt((delta.gamma/a/gamma^2)^2 + (delta.a/a^2/gamma)^2)
              keep$splicing.time[n] = Y; keep$cv.splicing.time[n] = delta.Y/Y;
            }
          }
        }
      } 
    }
    ## scoring the estimation of parameters
    #if(length(test)==length(which(!is.na(test)))) keep$score.params[n] = max(test)
  }
  ### add half-lived measured from other data sets
  mm = match(keep$gene, fried[,1]); keep$Fried = as.numeric(fried[mm, 2]);
  mm = match(keep$gene, shar[,1]); keep$Shar = as.numeric(shar[mm, 2]);
  mm = match(keep$gene, schwa[,1]); keep$Schwa = as.numeric(schwa[mm, 2]);
  
  return(keep);
}