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

####################
## original code for transforming estimated parameters and making clean table for figures  
####################
Make.Table.Main.Results = FALSE
if(Make.Table.Main.Results)
{
  source('functions.R')
  
  Make.Summary.Results = FALSE
  if(Make.Summary.Results)
  {
    res.version.old = '_total_counts_all_genes_norm_params_v6';
    res.version = '_total_counts_all_genes_norm_params_final_v7';
        
    load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version.old, '_OLD.Rdata',sep = ''))
    load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
    
    load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_RNA_seq_total_counts_table_sx', res.version.old, '.Rdata',sep = ''))
    load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_RNA_seq_total_counts_table_sx_All', res.version.old, '.Rdata',sep = ''))
    
    names = colnames(keep)[1:25]
    index = setdiff(c(1:nrow(T)), keep$index)
    
    kk = index
    source('functions.R')
    xx = Transform.Scoring.parameters(T)
    xx = data.frame(T, xx, stringsAsFactors = FALSE)
    
    kk = c(138:187, 190:200, 224)
    yy = xx[, -kk]
    mm = match(yy$gene, R.norm$gene)
    
    ii = grep('.rpkm.', colnames(yy)); print(colnames(yy)[ii])
    jj = grep('.log2rpkm.', colnames(R.norm)); print(colnames(R.norm)[jj])
    yy[, ii] = R.norm[mm, jj];colnames(yy)[ii] = colnames(R.norm)[jj]
    
    colnames(yy)[grep('outlier', colnames(yy))] = c('outlier.mRNA', 'outlier.premRNA')
    colnames(yy)[c(160:162)] = paste0(colnames(yy)[c(160:162)], '.half.life')
    threshold.identifiability = qchisq(p=0.68, df=1, lower.tail = TRUE)
    yy$non.identifiability.gamma.L = yy$non.identifiability.gamma.L>threshold.identifiability;
    yy$non.identifiability.gamma.R = yy$non.identifiability.gamma.R>threshold.identifiability;
    colnames(yy)[grep('identifiability', colnames(yy))] = c('gamma.identifiability.L','gamma.identifiability.R')
    
    yy$best.model[which(yy$prob.best.model<=0.5)] = NA
    
    T = yy;
    res.version = '_total_counts_all_genes_norm_params_final_v7_All_summary';
    res.version.old = '_total_counts_all_genes_norm_params_v6_OLD';
    save(T, file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
    save(T.old, keepAll, keep, file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', 
                                          res.version.old, '.Rdata',sep = ''))
  }
  
    
  res.version = '_total_counts_all_genes_norm_params_final_v7_All_summary';
  res.version.old = '_total_counts_all_genes_norm_params_v6_OLD';
  load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
  load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', 
                  res.version.old, '.Rdata',sep = ''))
  T$BIC.best.model = T$best.model;
  
  #jj = setdiff(which(xx$best.model>1), keep$index)
  xx = T
  
  jj = grep('.count.', colnames(xx))
  xx = xx[, -jj]
  ii = grep('.log2rpkm.', colnames(xx))
  xx = xx[, -ii]
  xx = xx[, -ncol(xx)]
  
  colnames(xx)[which(colnames(xx)=='best.model')] = 'BIC.best.model'
  
  jj = which(xx$sem.phase.gamma>12)
  kk = grep("gamma.identifiability", colnames(xx))
  
  xx[jj, kk]
  
  xx$sem.phase.gamma[which(xx$sem.phase.gamma>12)] =12 
  xx$sem.phase.int[which(xx$sem.phase.int>12)]=12
  
  #xx = xx[, -c(2, 26)]
  supDir = "/Users/jiwang/Dropbox/mRNA_decay/Revision/Tables_supp/"
  write.table(xx, file=paste0(supDir, "Table_S1_2.txt"), sep='\t', col.names = TRUE, row.names = FALSE, quote=FALSE)
  
  #### Filter the estimated parameters
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
    
    ###
    ### filter parameters estimated for premRNAs
    TEST.plot = FALSE
    if(TEST.plot)
    {
      par(mfrow=c(2, 2))
      hist(sels$cv.Min.int, breaks=100, xlim = c(0, 2));abline(v=c(0.2, 0.3, 0.4, 0.5), col='red')
      hist(sels$cv.Amp.int, breaks=100, xlim = c(0, 2));abline(v=c(0.2, 0.3, 0.4, 0.5), col='red')
      hist(sels$sem.phase.int, breaks=500, xlim = c(0, 3));abline(v=c(0.5, 1, 1.5, 2), col='red')
      hist(sels$cv.beta.int, breaks=1000, xlim = c(0, 2));abline(v=c(0.5, 1, 1.5, 2), col='red')
      
      par(mfrow=c(1, 1))
      plot(sels$beta.int, sels$cv.beta.int, cex=0.2, ylim = c(0, 5))
      kk = which(sels$best.model==2 | sels$best.model==4);length(kk)
      jj = kk[which(sels$cv.Min.int[kk]<0.3)]; length(jj);
      jj = kk[which(sels$cv.Min.int[kk]<0.3 & sels$cv.Amp.int[kk]<0.3)]; length(jj);
      jj = kk[which(sels$cv.Min.int[kk]<0.3 & sels$cv.Amp.int[kk]<0.3 & sels$sem.phase.int[kk]<1.5)]; length(jj);
    }
    
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
    
    examples = c('Per2', 'Cry1', 'Per3','Per1', 'Npas2', 'Arntl', 'Clock', 'Cry2', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Rorc','Rora',    
                 'Nfil3', 'Bhlhe40', 'Bhlhe41')
    mm = match(examples, sels$gene);
    examples[which(is.na(mm))]
    
    ###
    ### select genes with identifiable parameter gamma
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
    
    par(mfrow=c(2, 2))
    hist(xx$cv.half.life, breaks=100, xlim = c(0, 2));abline(v=c(0.2, 0.3, 0.4, 0.5, 1), col='red')
    hist(xx$cv.splicing.time, breaks=100, xlim = c(0, 2));abline(v=c(0.2, 0.3, 0.4, 0.5, 1), col='red')
    hist(xx$cv.eps.gamma, breaks=100, xlim = c(0, 3));abline(v=c(0.2, 0.3, 0.5, 1), col='red')
    hist(xx$sem.phase.gamma, breaks=200, xlim = c(0, 2));abline(v=c(0.5, 1, 1.5, 2), col='red')
    
    keep = rbind(xx, yy)
    
    o1 = order(keep$best.model)
    keep = keep[o1,]
    
    res.version = '_total_counts_all_genes_norm_params_final_v7_All_summary';
    save(T, keep, file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
    
  }
  
  #save(keep, file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_RNA_seq_total_counts_table_sx', res.version, '.Rdata',sep = ''))
  #keep$score.params = apply(keep, 1, scoring.estimation.params)
  
  Add.phase.amp.4mrna = FALSE
  if(Add.phase.amp.4mrna)
  {
    load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_RNA_seq_total_counts_table_sx', res.version, '.Rdata',sep = ''))
    index = keep$index;
    keep$qv.premRNA = T$qv.rpkm.premRNA[index]; keep$qv.mRNA = T$qv.rpkm.mRNA[index];
    
    ptm <- proc.time()
    
    #library(parallel)
    #no_cores <- detectCores()-1; cl <- makeCluster(no_cores);  # Initiate cluster
    source('functions.R')
    #parLapply(cl, index[1:10],  compute.phase.amp.4premrna.mrna.with.fitting.res, T=T)
    res = t(apply(as.matrix(index), 1, compute.phase.amp.4premrna.mrna.with.fitting.res, T=T));
    proc.time() - ptm
    #aa$relamp.premRNA = T$rel.amp.rpkm.premRNA[index];aa$relamp.mRNA = T$rel.amp.rpkm.mRNA[index];
    #aa$phase.premRNA = T$phase.rpkm.premRNA[index];aa$phase.mRNA = T$phase.rpkm.mRNA[index];
    xx = data.frame(keep, res, stringsAsFactors = FALSE)
    
    ### quick check
    plot(xx$relamp.mRNA, T$rel.amp.rpkm.mRNA[index], cex=0.2);abline(0, 1, col='red')
    plot(xx$relamp.premRNA, T$rel.amp.rpkm.premRNA[index], cex=0.2);abline(0, 1, col='red')
    plot(xx$phase.mRNA, T$phase.rpkm.mRNA[index], cex=0.2);abline(0, 1, col='red')
    plot(xx$phase.premRNA, T$phase.rpkm.premRNA[index], cex=0.2);abline(0, 1, col='red')
    
    #keep = xx;
    #save(keep, file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_RNA_seq_total_counts_table_sx', res.version, '.Rdata',sep = ''))
  }
  
  #xx = keep[, -c(2)]
  #write.table(xx, file="/Users/jiwang/Dropbox/mRNA_decay/Manuscript/Tables_supp/Table_S2_3.txt", sep='\t', col.names = TRUE, row.names = FALSE, quote=FALSE)
  
}
