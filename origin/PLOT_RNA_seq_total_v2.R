################################################################################
################################################################################
######## Fit Results (model selection and parameter estimation) with Read data of read counts
################################################################################
################################################################################
Fitting.Results = FALSE
if(Fitting.Results)
{
  #res.version = '_total_counts_all_genes_norm_params_v3';
  #res.version = '_total_counts_all_genes_norm_params_v4';
  #res.version = '_total_counts_all_genes_norm_params_v5';
  load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
  source('functions.R')
  BIC = my.model.selection.loglike(T=T, method='BIC', outliers = TRUE)
  T = data.frame(T, BIC, stringsAsFactors = FALSE);
  T = T[which(!is.na(T$BIC.best.model)==TRUE), ]
  Tt = T
  
  write.table(keep, file = '/Users/jiwang/Dropbox/mRNA_decay/analysis/Fitting_results_total_RNA_seq_read_counts_NB_estimated_parameters_v2.txt', 
              quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  write.table(R.norm, file = '/Users/jiwang/Dropbox/mRNA_decay/analysis/Table_RPKM_premRNA_mRNA.txt', 
              quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  #write.table(T, file = '/Users/jiwang/Dropbox/mRNA_degradatation/analysis/Fitting_results_total_RNA_seq_read_counts_NB.txt', 
  #quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  
  #### Result-check part
  Quick.Check.Results = FALSE
  if(Quick.Check.Results)
  {
    length(which(Tt$BIC.best.model==1));length(which(Tt$BIC.best.model==2));length(which(Tt$BIC.best.model==3));length(which(Tt$BIC.best.model==4))
    length(which(Tt$BIC.best.model==2))/length(which(Tt$BIC.best.model!=1))
    length(which(Tt$BIC.best.model==3 | Tt$BIC.best.model==4))/length(which(Tt$BIC.best.model!=1))
    length(which(Tt$BIC.best.model==3))/length(which(Tt$BIC.best.model!=1))
    length(which(Tt$BIC.best.model==4))/length(which(Tt$BIC.best.model!=1))
    
    examples = c('Npas2', 'Arntl', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry1', 'Cry2', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Nfil3', 'Bhlhe40', 'Bhlhe41', 'Nampt', 'Parp1', 'Prkaa2', 'Prkag1', 'Gsk3a', 'Gsk3b', 'Csnk1d')
    mm = match(examples, T[,1])
    cbind(T$gene[mm], T$BIC.best.model[mm])
    
    #### Check Half-lives, half-life amplitudes, splicing.k
    source('model_RNA_seq_total/functions.R')
    bounds = set.bounds(model=4)
    upper = bounds$upper
    lower = bounds$lower
    
    gamma.max = log(2)/(1/60);
    gamma.min = log(2)/24;
    splicing.k.max = upper[4]
    splicing.k.min = lower[4]
    amp.gamma.max = upper[2]
    
    ### 10 min half-life boundary
    length(which(Tt$BIC.best.model==2 & (abs(Tt$gamma.m2-gamma.max)<0.00001)))/length(which(Tt$BIC.best.model==2))
    length(which(Tt$BIC.best.model==3 & (abs(Tt$gamma.m3-gamma.max)<0.00001)))/length(which(Tt$BIC.best.model==3))
    length(which(Tt$BIC.best.model==4 & (abs(Tt$gamma.m4-gamma.max)<0.00001)))/length(which(Tt$BIC.best.model==4))
    #### 24 h half-life boundary
    length(which(Tt$BIC.best.model==2 & (abs(Tt$gamma.m2-gamma.min)<0.00001)))/length(which(Tt$BIC.best.model==2))
    length(which(Tt$BIC.best.model==3 & (abs(Tt$gamma.m3-gamma.min)<0.00001)))/length(which(Tt$BIC.best.model==3))
    length(which(Tt$BIC.best.model==4 & (abs(Tt$gamma.m4-gamma.min)<0.00001)))/length(which(Tt$BIC.best.model==4))
    #### amplitude of gamma boundary
    length(which(Tt$BIC.best.model==3 & (abs(Tt$amp.gamma.m3-amp.gamma.max)<0.00001)))/length(which(Tt$BIC.best.model==3))
    length(which(Tt$BIC.best.model==4 & (abs(Tt$amp.gamma.m4-amp.gamma.max)<0.00001)))/length(which(Tt$BIC.best.model==4))
    #### 1s splicing.k boundary
    length(which(Tt$BIC.best.model==2 & (abs(Tt$splicing.k.m2-splicing.k.max)<0.00001)))/length(which(Tt$BIC.best.model==2))
    length(which(Tt$BIC.best.model==3 & (abs(Tt$splicing.k.m3-splicing.k.max)<0.00001)))/length(which(Tt$BIC.best.model==3))
    length(which(Tt$BIC.best.model==4 & (abs(Tt$splicing.k.m4-splicing.k.max)<0.00001)))/length(which(Tt$BIC.best.model==4))
    #### 30min splicing.k boundary
    length(which(Tt$BIC.best.model==2 & (abs(Tt$splicing.k.m2-splicing.k.min)<0.00001)))/length(which(Tt$BIC.best.model==2))
    length(which(Tt$BIC.best.model==3 & (abs(Tt$splicing.k.m3-splicing.k.min)<0.00001)))/length(which(Tt$BIC.best.model==3))
    length(which(Tt$BIC.best.model==4 & (abs(Tt$splicing.k.m4-splicing.k.min)<0.00001)))/length(which(Tt$BIC.best.model==4))
    
    T[which(T$gene=='Arhgef18'), ]
    
    ##### Check individual examples
    if(!file.exists('myplots/Total/gene_examples/')){dir.create('myplots/Total/gene_examples/')}
    
    source('model_RNA_seq_total/functions.R')
    
    for(m in c(2:4))
    {
      index = which(T$BIC.best.model==m)
      index = index[order(T$pval.mRNA[index])]
      #index = kk[order(T$BIC.best.model[kk])]
      #index = c(1:10)
      pdfname = paste('myplots/Total/gene_examples/fitting_results_all_model_', m, '.pdf', sep='')
      plot.genes.examples(index, pdfname, T)
    }
    
    #### Check half-lives of mRNAs
    source('model_RNA_seq_total/functions.R')
    pdfname = 'myplots/Total/gene_examples/Comparison_half_lives_v2.pdf'
    plot.half.life.comparison(index, pdfname, T)
    
    ############
    ### check the noise for the genes with bad fitting
    ############
    noise.check.4fitting = FALSE
    if(noise.check.4fitting)
    {
      exon = log(as.matrix(T[,c(2:49)]))
      intron = log(as.matrix(T[,c(50:97)]))
      test = cbind(apply(exon, 1, variance.estimate.replicates, nb.replicates=4), apply(intron, 1, variance.estimate.replicates, nb.replicates=4))
      
      plot(test[,1], test[,2], xlab='Exon sigma', ylab='Intron sigma', log='', cex=0.5)
      mm = match(c('Hlf', 'Ppara'), T$gene)
      points(test[mm, 1], test[mm, 2], cex=0.7, col='red')
      abline(0, 1, col='red', lwd=2.0)
    }
    
  }
}

#######################################
#######################################
####### Plot Figures 
#######################################
#######################################
Figure_1 = FALSE
if(Figure_1)
{
  ### total RNA-seq
  data.version = '_total_counts_v2';
  load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_sel_alphas', data.version, '.Rdata', sep=''))
  #source("/Users/jiwang/Degradation_Liver/Main_Code_Degradation/f24_modified_1.0.r");
  attribute.global.variable.colors();
  load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_sel_alphas_log2_RPKM', 
                  data.version, '.Rdata', sep=''))
  T = R.norm;
  source('functions.R')

  ZT.ex = grep('.count.mRNA', colnames(T))
  ZT.int = grep('.count.premRNA', colnames(T))
  zt = seq(0,94,by = 2)
  
  ##############################
  ##### Overview of the Method
  ##############################  
  
  Overview.Method.with.different.scenarios = FALSE
  if(Overview.Method.with.different.scenarios)
  {
    folder = "/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/"
    ### distribution of mRNA/pre-mRNA ratios for cyclic and rhythmic genes
    data.version = '_total_counts_v2'
    #load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_sel_alphas', data.version, '.Rdata', sep=''))
    #source("/Users/jiwang/Degradation_Liver/Main_Code_Degradation/f24_modified_1.0.r")
    #T = R;
    
    ### representative example
    zt = seq(0, 94, by=2); Min = 5;
    gamma = log(2)/2; eps.gamma = 0.3; phase.gamma = 18; splicing.k = 5*gamma;
    param.synthesis.1 = Min; param.synthesis.2 = 10; param.synthesis.3 = 6; param.synthesis.4 = 1;
    par = c(gamma, eps.gamma, phase.gamma, splicing.k, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4)
    s = compute.s.beta(t = zt, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4)
    m = compute.m.beta(t = zt, gamma, eps.gamma*sqrt(1+w^2/gamma^2), (phase.gamma-atan2(w, gamma)/w), splicing.k*gamma,
                       param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4);
    set.seed(42); m = m + rnorm(n=length(m), mean = 0, sd=1); s = s + rnorm(n=length(s), mean = 0, sd=0.4);
    lims = range(s, m)
    
    pdfname = paste(folder, 'temporal_profiles_representative_example.pdf', sep='')
    pdf(pdfname, width=2.2, height=1.8)
    #pdf(pdf.name, width=2.4, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    plot(zt, s, type='l', lwd=1., lty=1, col=col.int, ylim=lims, xlab=NA, ylab=NA, main=NA, axes=FALSE, log='y')
    points(zt, s, col=col.int, type='p', cex=0.5, pch=5)
    points(zt, m, type='l', lwd=1., col=col.ex, lty=1)
    points(zt, m, type='p', cex=0.7, col=col.ex, pch=16)
    #box()
    axis(2, at=seq(5, 30, by=5), cex.axis=1.0, labels = FALSE, tick = TRUE, tck=-0.02)
    axis(1, at=seq(0, 96, by=6), las=1,cex.axis = 1.0, labels = FALSE, tck=-0.02)
    par(mgp = c(0, 0.2,0))
    axis(1, at=seq(0, 96, by=24), las=1,cex.axis = 1.0, tick = FALSE)
    #abline(v=6, col='gray', lwd=2.0)
    #abline(v=12, col='gray', lwd=2.0)
    dev.off()
    
    zt.p = seq(0, 96, by=0.1);
    for(model in c(1:4))
    {
      pdfname = paste(folder, 'temporal_profiles_representative_example_model', model, '.pdf', sep='')
      pdf(pdfname, width=0.7, height=0.5)
      #pdf(pdf.name, width=2.4, height=2.0)
      par(cex = 0.7, las = 1, mgp = c(0.2,0.2,0), mar = c(0.2,0.2,0.2,0.2), tcl = -0.3)
      
      if(model==1){
        s = rep(2, length(zt.p)); deg = rep(1, length(zt.p));
      }else{
        if(model==2){
          eps.gamma = 0.0; phase.gamma = 0.0; splicing.k = par[2];
          param.synthesis.1 = par[5]; param.synthesis.2 = par[6]; param.synthesis.3 = par[7]; param.synthesis.4 = par[8];
        }
        if(model==3){
          eps.gamma = par[2]; phase.gamma = par[3]; splicing.k = par[4]; param.synthesis.1 = par[5];
          param.synthesis.2 = 0; param.synthesis.3 = 0; param.synthesis.4 = 1;
        }
        if(model==4){
          eps.gamma = par[2]; phase.gamma = par[3]; splicing.k = par[4];
          param.synthesis.1 = par[5]; param.synthesis.2 = par[6]; param.synthesis.3 = par[7]; param.synthesis.4 = par[8];
        }
        #print(c(param.synthesis.3, phase.gamma))
        w = 2*pi/24;
        s = compute.s.beta(t = zt.p, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4)
        s = s/mean(s)+1
        deg = 1*(1+eps.gamma*cos(w*(zt.p-phase.gamma)));
      }

      plot(zt.p, s, type='l', lwd=1, lty=1, col=col.int, ylim=c(0.6, 2.6), xlab=NA, ylab=NA, main=NA, axes=FALSE)
      points(zt.p, deg, type='l', lwd=1., col=col.deg, lty=1)
      abline(h=1, col=col.deg, lwd=1.0, lty=2);abline(h=2, col=col.int, lwd=1.0, lty=2);
      axis(2, cex.axis=1.0, labels = FALSE, tick = TRUE, lwd.ticks = 0.0)
      axis(1, cex.axis = 1.0, labels = FALSE, tick = TRUE, lwd.ticks = 0.0)
      box()
      #abline(v=12, col='gray', lwd=2.0)
      dev.off()
    }
    
    ### Illustration of model selection
    pdfname = paste(folder, 'Illustration_model_selection.pdf', sep='')
    pdf(pdfname, width=1., height=0.8)
    #pdf(pdf.name, width=2.4, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.0,0.4,0), mar = c(0.4,1.6,0.4,0.2), tcl = -0.3)
    #par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(1.0,1.0,1.0,3)+0.1, tcl = -0.3)
    pie.sales= c(0.05, 0.1, 0.05, (1-0.05-0.1-0.05))
    #names(pie.sales) <- NA;
    cols= c(1:4)
    barplot(pie.sales, col = cols, axes = FALSE, ylim=c(0, 1))
    axis(2, las=1, at= seq(0, 1, by=0.5), cex.axis=0.7, labels = TRUE, tick = TRUE, tcl=-0.2)
    #axis(1, cex.axis = 1.0, labels = FALSE, tick = TRUE, lwd.ticks = 0.0)
    #pie(pie.sales,col=cols, cex=0.7)
    #legend(0.65, 1.25, legend = c(' ',' ', ' ', ' ') , cex=0.9, pt.cex=1.5, pt.lwd=1.0, fill= cols, border = NA, bty = 'n')
    #title(main="Distribution of proteins (4016 proteins of 5827) for each models", cex.main=1.8, font.main=2)
    #text(0, -0.98, paste("M.1: Constant mRNA and Constant Protein \n M.2: Rhythmic mRNA and Constant Protein \n M.3: Constant mRNA and Rhythmic Protein \n M.4: Rhythmic mRNA and Rhythmic Protein"))
    dev.off()
    
  }
  
  ##########
  ### Overview of Data
  ##########
  Check.rhythmicy.total.RNA_seq = FALSE
  if(Check.rhythmicy.total.RNA_seq)
  {
    data.version = '_total_counts_v2'
    attribute.global.variable.colors();
    load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_sel_alphas_log2_RPKM', 
                    data.version, '.Rdata', sep=''))
    
    folder = "/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/"
    
    jj = which(R.norm$qv.log2rpkm.mRNA<0.05 & R.norm$rel.amp.log2rpkm.mRNA>0.1)
    phases = R.norm$phase.log2rpkm.mRNA[jj]
    
    pdfname = paste(folder, 'Phases_Rhythmic_mRNAs_all.pdf', sep='')
    pdf(pdfname, width=1.8, height=1.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    circular_phase24H_histogram(phases, col=col.ex, cex.axis=0.5, cex.lab=0.01, lwd=0.7)
    dev.off()
    
    amps = R.norm$amp.log2rpkm.mRNA[jj]
    pdfname = paste(folder, 'Amplitudes_Rhythmic_mRNAs_all.pdf', sep='')
    pdf(pdfname, width=1.8, height=1.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    breaks=seq(0, 8, by=0.5);
    h = hist(amps, breaks=breaks, plot=FALSE)
    y = h$counts
    y[which(y<=0)] = 0.5;
    lwd = 1.
    plot(h$breaks, c(NA,y), type='S', ylim=c(1, max(y)), main=NA, col=col.ex, xlab=NA, ylab=NA, axes=FALSE, log='y', lwd=lwd)
    axis(1,at=seq(0, 8, by=1),cex.axis =1.0)
    axis(2, at= c(1, 10, 100, 1000), las=1,cex.axis = 1.0)
    lines(h$breaks, c(h$counts,NA), type='s', lwd=lwd, col=col.ex)
    lines(h$breaks, c(NA,h$counts), type='h', lwd=lwd, col=col.ex)
    lines(h$breaks, c(h$counts,NA), type='h',lwd=lwd, col=col.ex)
    lines(h$breaks, rep(0,length(h$breaks)), type='S', col=col.ex)
    invisible(h)
    dev.off()
    
    jj = which(R.norm$qv.log2rpkm.premRNA<0.05 & R.norm$rel.amp.log2rpkm.premRNA>0.1)
    phases = R.norm$phase.log2rpkm.premRNA[jj]
    pdfname = paste(folder, 'Phases_Rhythmic_premRNAs_all.pdf', sep='')
    pdf(pdfname, width=1.8, height=1.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    circular_phase24H_histogram(phases, col=col.int, cex.axis=0.5, cex.lab=0.01, lwd=0.7)
    dev.off()
    
    amps = R.norm$amp.log2rpkm.premRNA[jj]
    pdfname = paste(folder, 'Amplitudes_Rhythmic_premRNAs_all.pdf', sep='')
    pdf(pdfname, width=1.8, height=1.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    breaks=seq(0, 8, by=0.5);
    h = hist(amps, breaks=breaks, plot=FALSE)
    y = h$counts
    y[which(y<=0)] = 0.5;
    lwd = 1.
    plot(h$breaks, c(NA,y), type='S', ylim=c(1, max(y)), main=NA, col=col.int, xlab=NA, ylab=NA, axes=FALSE, log='y', lwd=lwd)
    axis(1,at=seq(0, 8, by=1),cex.axis =1.0)
    axis(2, at= c(1, 10, 100, 1000), las=1,cex.axis = 1.0)
    lines(h$breaks, c(h$counts,NA), type='s', lwd=lwd, col=col.int)
    lines(h$breaks, c(NA,h$counts), type='h', lwd=lwd, col=col.int)
    lines(h$breaks, c(h$counts,NA), type='h',lwd=lwd, col=col.int)
    lines(h$breaks, rep(0,length(h$breaks)), type='S', col=col.int)
    invisible(h)
    dev.off()
    
    jj = which(R.norm$qv.log2rpkm.premRNA<0.05 & R.norm$qv.log2rpkm.mRNA<0.05 & R.norm$rel.amp.log2rpkm.premRNA>0.1 & R.norm$rel.amp.log2rpkm.mRNA>0.1)
    phases1 = R.norm$phase.log2rpkm.premRNA[jj]
    phases2 = R.norm$phase.log2rpkm.mRNA[jj]
    dd = (phases2 - phases1)%%24;
    
    amp1 = R.norm$rel.amp.log2rpkm.premRNA[jj]
    amp2 = R.norm$rel.amp.log2rpkm.mRNA[jj]
    
    kk = which(dd>12);dd[kk] = dd[kk]-24;
    circular.error=function(p,x,y)
    {
      sum(1-cos(2*pi/24*(y-x-p)))
    }
    fit=nlm(f=circular.error, p=2, x=phases1, y=phases2)
    
    pdfname = paste(folder, 'Phases_Rhythmic_mRNAs_vs_premRNAs.pdf', sep='')
    pdf(pdfname, width=2.0, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3, pty='s')
    #circular_phase24H_histogram(phases, col=col.int, cex.axis=0.5, cex.lab=0.25, lwd=0.5)
    plot(phases1, phases2, cex=0.03, col='black', xlab=NA, ylab=NA, axes=FALSE)
    abline(fit$estimate%%24, 1, col='blue', lwd=1.5, lty=1);
    #abline((fit$estimate%%24-24), 1, col='red', lwd=1.0);abline((fit$estimate%%24+24), 1, col='red', lwd=1.0)
    abline(0, 1, col='red', lty=1, lwd=1.0)
    abline(-24, 1, col='red', lty=1, lwd=1.0)
    abline(24, 1, col='red', lty=1, lwd=1.0)
    
    axis(2, at=seq(0, 24, by=6), las=1,cex.axis = 1.0, tck=-.04)
    axis(1, at=seq(0, 24, by=6), las=1,cex.axis = 1.0, tck=-.04)
    box()
    
    dev.off()
    
    pdfname = paste(folder, 'RelAmp_Rhythmic_mRNAs_vs_premRNAs.pdf', sep='')
    pdf(pdfname, width=2.0, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3, pty='s')
    #circular_phase24H_histogram(phases, col=col.int, cex.axis=0.5, cex.lab=0.25, lwd=0.5)
    plot(amp1, amp2, cex=0.03, col='black', xlab=NA, ylab=NA, axes=FALSE, log='xy')
    #abline(fit$estimate%%24, 1, col='blue', lwd=1.5, lty=1);
    #abline((fit$estimate%%24-24), 1, col='red', lwd=1.0);abline((fit$estimate%%24+24), 1, col='red', lwd=1.0)
    abline(0, 1, col='red', lty=1, lwd=1.0)
    #abline(-24, 1, col='red', lty=1, lwd=1.0)
    #abline(24, 1, col='red', lty=1, lwd=1.0)
    axis(2, at=c(0.1, 0.2, 0.5, 1), las=1,cex.axis = 1.0, tck=-.04)
    axis(1, at=c(0.1, 0.2, 0.5, 1), las=1,cex.axis = 1.0, tck=-.04)
    box()
    
    dev.off()
    
  }
  
  ####### Compare total RNA-seq with nascent RNA-seq
  Compare.nascent.RNA.seq.micro_exons = FALSE
  if(Compare.nascent.RNA.seq.micro_exons)
  {
    ### import nascent RNA-seq and microarray data
    nascent = read.table('/Users/jiwang/RNA_seq_Data/Total-RNA-seq/Menet_Nascent_RNAs.txt',sep='\t', header=TRUE, na.strings = "-Inf")
    source("/Users/jiwang/Degradation_Liver/Main_Code_Degradation/f24_modified_1.0.r");
    
    ### process nascent rna-seq data
    zt.nascent = grep('_Nascent', colnames(nascent))[1:12]
    zt.mrna = grep('mRNA', colnames(nascent))[1:12]
    xx = data.frame(nascent[,1], nascent[,zt.nascent], nascent[, zt.mrna])
    colnames(xx)[1] = 'gene'
    nascent = xx
    ii = grep('_Nascent', colnames(nascent))
    data = as.matrix(nascent[,ii])
    time = c(0:11)*4
    res = t(apply(data, 1, f24_R2_alt2, t=time))
    res = cbind(res, qv=qvals(res[,6]))
    #res[,4] = t(apply(2^data, 1, f24_R2_alt2, t=time))
    colnames(res) = paste(colnames(res), '.nascent', sep='')
    nascent = data.frame(nascent, res, stringsAsFactors=FALSE)
    ii = grep('_mRNA', colnames(nascent))
    data = as.matrix(nascent[,ii])
    time = c(0:11)*4+2
    res = t(apply(data, 1, f24_R2_alt2, t=time))
    res = cbind(res, qv=qvals(res[,6]))
    #res[,4] = t(apply(2^data, 1, f24_R2_alt2, t=time))
    colnames(res) = paste(colnames(res), '.mrna', sep='')
    nascent = data.frame(nascent, res, stringsAsFactors=FALSE)
    o1 = order(nascent$qv.mrna)
    nascent = nascent[o1,]
    
    ## Here we use the linear signals normalized by their means
    zt.nascent = grep('Nascent', colnames(nascent))
    zt.mrna = grep('mRNA', colnames(nascent))
    ZT.mRNA = grep('.rpkm.mRNA', colnames(R.norm))[1:48]
    ZT.premRNA = grep('.rpkm.premRNA', colnames(R.norm))[1:48]
    
    #examples = c('Per2', 'Cry1', 'Per3','Per1', 'Npas2', 'Arntl', 'Clock', 'Cry2', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Rorc','Rora',    
    #             'Nfil3', 'Bhlhe40', 'Bhlhe41')
    examples = c('Arntl', 'Per1', 'Per2', 'Cry1', 'Cry2', 'Nr1d1', 'Nr1d2', 'Rorc')
    list = examples[which(!is.na(match(examples, R.norm$gene)))];
    folder = "/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/Examples/Compare_profiles"
    source('functions.R')
    #for(n in 1:nrow(nascent))
    for(name in list)
    {
      pdfname = paste(folder, '/comparison_total_nascent_exonarray_', name, '.pdf', sep='')
      pdf(pdfname, width = 1.8, height = 1.2)
      par(cex = 0.7, las = 1, tcl = -0.3)
      #par(mfcol = c(2,1))
      #layout(matrix(c(1,2),nrow = 2, ncol = 1), heights=c(1.1,1.))
      ### pre-mRNAs
      kk = which(R.norm$gene==name);yy0=as.numeric(R.norm[kk, ZT.premRNA]);yy0 = yy0/mean(yy0,na.rm=TRUE);
      averg0 = mean.err(yy0, interval = 2)[1,]; err0 = mean.err(yy0, interval = 2)[2,];
      kk1 = which(nascent$gene==name)
      if(length(kk1)>0){
        yy1 = as.numeric(nascent[kk1, zt.nascent]); yy1 = yy1/mean(yy1, na.rm = TRUE);
        averg1 = mean.err(yy1, interval = 4)[1,]; err1 = mean.err(yy1, interval = 4)[2,];
      }else{
        averg1 = rep(NA, 6); err1 = rep(NA, 6);
      }
      
      lims = range(c(averg0-err0, averg0+err0, averg1-err1, averg1+err1), na.rm=TRUE)
      time = c(0:11)*2
      #par(mgp = c(0.5,0.3,0.), mar = c(1.2,1.0,0.5,0.5));
      par(mgp = c(0.5,0.3,0), mar = c(1.4,1.4,0.7,0.5));
      plot(time, averg0, type='l', ylim=lims, xlim=c(0, 24), col=col.int, lwd=1.2, pch=16, main=name, cex.main=0.7, axes=FALSE, xlab=NA, ylab=NA)
      points(time, averg0, type='p', cex=0.6, col=col.int, pch=1)
      arrows(time, averg0-err0, time, averg0+err0, length=0.02, angle=90, code=3, col=col.int, lwd=1.0, cex=0.6)
      #points(time, averg2, type='l', ylim=lims, col='green', lwd=1.0, lty=2)
      #points(time, averg2, type='p', cex=0.3, col='green', pch=2)
      #arrows(time, averg2-err2, time, averg2+err2, length=0.02, angle=90, code=3, col='green', lwd=0.7, cex=0.3)
      time = seq(0, 20, by=4)
      points(time, averg1, type='l', ylim=lims, col='darkorange', lwd=1.2, lty=2)
      points(time, averg1, type='p', cex=0.6, col='darkorange', pch=1)
      arrows(time, averg1-err1, time, averg1+err1, length=0.02, angle=90, code=3, col='darkorange', lwd=1.0, cex=0.6)
      #abline(h=1, col='darkgray', lwd=1.0)
      axis(2, las=1,cex.axis = 0.8, tck=-.04)
      axis(1, at=seq(0, 24, by=6), las=1,cex.axis = 0.8, tck=-.04, labels = FALSE)
      axis(1, at=seq(0, 24, by=6), las=1,cex.axis = 0.8, tck=-.0, line = -0.2, lwd = 0)
      #box();
      abline(h=1, col='gray', lwd=1.0)
      dev.off()
      
    }
    
    #### Comparison of pre-mRNA phases measured by total and nascent RNA-seq 
    jj = which(R.norm$qv.log2rpkm.premRNA<0.05 & R.norm$rel.amp.log2rpkm.premRNA>0.1)
    mm = match(R.norm$gene[jj], nascent$gene)
    #kk = which(nascent$pval.nascent[mm]<0.1);jj=jj[kk];mm=mm[kk];
    
    library(circular)
    x=R.norm$phase.log2rpkm.premRNA[jj];
    y=nascent$phase.nascent[mm];
    kk = which(!is.na(y)==TRUE);x=x[kk];y=y[kk];
    x=x/24*2*pi
    y=y/24*2*pi
    
    x = circular(x, modulo="asis", template="clock24", type = 'angles', units = 'radians', rotation = 'clock')
    y = circular(y, modulo="asis", template="clock24", type = 'angles', units = 'radians', rotation = 'clock')
    
    cor.circular(x, y, test=TRUE)
    
    pdfname = paste(folder, '_phase_comparison_Total_nascent.pdf', sep='')
    pdf(pdfname, width=2.0, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3, pty='s')
    
    plot(R.norm$phase.log2rpkm.premRNA[jj], nascent$phase.nascent[mm], cex=0.1, col=col.int, xlab=NA, ylab=NA, axes=FALSE)
    abline(0, 1, col='red', lty=1, lwd=1.0)
    abline(-24, 1, col='red', lty=1, lwd=1.0)
    abline(24, 1, col='red', lty=1, lwd=1.0)
    axis(2, at=seq(0, 24, by=6), las=1,cex.axis = 1.0, tck=-.04)
    axis(1, at=seq(0, 24, by=6), las=1,cex.axis = 1.0, tck=-.04)
    box()
    dev.off()
    
    #### phase delay between mRNA and pre-mRNA in Menet et al.
    jj = which(nascent$pval.nascent<0.05 & nascent$pval.mrna<0.05)
    phases1 = nascent$phase.nascent[jj]
    phases2 = nascent$phase.mrna[jj]
    dd = (phases2 - phases1)%%24;
    kk = which(dd>12);dd[kk] = dd[kk]-24;
    circular.error=function(p,x,y)
    {
      sum(1-cos(2*pi/24*(y-x-p)))
    }
    fit=nlm(f=circular.error, p=2, x=phases1, y=phases2)
    
    pdfname = paste(folder, 'Phases_Rhythmic_mRNAs_vs_premRNAs.pdf', sep='')
    pdf(pdfname, width=2.0, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3, pty='s')
    #circular_phase24H_histogram(phases, col=col.int, cex.axis=0.5, cex.lab=0.25, lwd=0.5)
    plot(phases1, phases2, cex=0.03, col='black', xlab=NA, ylab=NA, axes=FALSE)
    abline(fit$estimate%%24, 1, col='blue', lwd=1.5, lty=1);
    #abline((fit$estimate%%24-24), 1, col='red', lwd=1.0);abline((fit$estimate%%24+24), 1, col='red', lwd=1.0)
    abline(0, 1, col='red', lty=1, lwd=1.0)
    abline(-24, 1, col='red', lty=1, lwd=1.0)
    abline(24, 1, col='red', lty=1, lwd=1.0)
    
    axis(2, at=seq(0, 24, by=6), las=1,cex.axis = 1.0, tck=-.04)
    axis(1, at=seq(0, 24, by=6), las=1,cex.axis = 1.0, tck=-.04)
    box()
    
    dev.off()
    
    #### amplitude comparison between total and nascent RNA-seq
    pdfname = paste(folder, '_relative_amps_comparison_Total_nascent.pdf', sep='')
    pdf(pdfname, width=1.5, height=1.5)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3, pty='s')
    
    plot(R.norm$rel.amp.log2rpkm.mRNA[jj], nascent$relamp.nascent[mm], cex=0.05, xlim=c(0, 1.2), ylim=c(0, 1.2), 
         col=col.int, log='', xlab=NA, ylab=NA, axes=FALSE)
    abline(0, 1, col='red', lty=1, lwd=1)
    axis(2, at=seq(0., 1.2, by=0.4), las=1,cex.axis = 0.7, tck=-.04)
    axis(1, at=seq(0., 1.2, by=0.4), las=1,cex.axis = 0.7, tck=-.04)
    box()
    dev.off()
    
  }
  
  #############
  #### Apply method for total RNA-seq data
  #############
  Illustrative.examples.for.method = FALSE
  if(Illustrative.examples.for.method)
  {
    ##### Import fitting result
    res.version = '_total_counts_all_genes_norm_params_v6';
    load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
    source('functions.R')
    BIC = my.model.selection.loglike(T=T, method='BIC', outliers = TRUE)
    T = data.frame(T, BIC, stringsAsFactors = FALSE);
    T = T[which(!is.na(T$BIC.best.model)==TRUE), ]
    Tt = T
    
    #### caculate parameters of mRNA from fitting reulsts instead of harmony repression
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
    
    res.version = '_total_counts_all_genes_norm_params_v6';
    #load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
    #res.version = '_total_counts_all_genes_norm_params_v6';
    load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
    source('functions.R')
    BIC = my.model.selection.loglike(T=T, method='BIC', outliers = TRUE)
    T = data.frame(T, BIC, stringsAsFactors = FALSE);
    T = T[which(!is.na(T$BIC.best.model)==TRUE), ]
    Tt = T
    
    #source('functions.R')
    load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_RNA_seq_total_counts_table_sx', res.version, '.Rdata',sep = ''))
    attribute.global.variable.colors();
    
    examples = c('Npas2', 'Arntl', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry1', 'Cry2', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef')
    mm = match(examples, keep$gene); 
    cbind(examples, keep$best.model[mm], keep$half.life[mm], keep$cv.half.life[mm]);
    
    examples = c('Acly', 'Actb', 'Ass1', 'Fasn', 'G6pc', 'Pck1', 'Srebf1', 'Insr', 'Gcgr', 'Nlrp6', 'Mlxipl')
    mm = match(examples, keep$gene); #mm = mm[which(!is.na(mm)==TRUE)]
    cbind(examples, keep$best.model[mm], log(2)/keep$half.life[mm], keep$cv.half.life[mm]);

    ##### Examples (positive controls; core clock and clock-related mRNA, previously reproted mRNAs)
    Examples4Figure = FALSE
    if(Examples4Figure)
    {
      ### representive examples of M2
      folder = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/Examples/';
      #threshold2 = qchisq(p=0.95, df=1, lower.tail = TRUE)
      #threshold3 = qchisq(p=0.68, df=1, lower.tail = TRUE)
      ### examples clock gene (short half-life)
      kk = match(c('Rorc', 'Cp','Fus', 'Per3', 'Cbs', 'Clock', 'Srebf1', 'Hnrnpc', 'Nr1d1', 'Arntl', 'Cry1'), keep$gene)
      index = keep$index[kk]
      #diff = T$rel.amp.rpkm.premRNA[index]-T$rel.amp.rpkm.mRNA[index]
      #index = index[order(-diff)]
      
      source('functions.R')
      plot.genes.examples(index=index, T=T, Figure=TRUE, folder=folder, summary.plot=TRUE);
      
      ### long half-life
      #kk = which(keep$best.model==2 & keep$qv.mRNA<0.05 & keep$qv.premRNA<0.01 & keep$half.life>=6 & keep$half.life<10
      #           & keep$cv.half.life<0.4 & keep$relamp.premRNA>0.5);
      #index = index[which(T$qv.rpkm.mRNA[index]<10^(-5) & T$rel.amp.rpkm.mRNA[index]>0.1 
      #                    & T$qv.rpkm.premRNA[index]<10^(-5) & T$rel.amp.rpkm.premRNA[index]>0.3
      #                   & (T$rel.amp.rpkm.premRNA[index]-T$rel.amp.rpkm.mRNA[index])>0.2)]
      #mm = match(examples, keep$gene); 
      #index = keep$index[kk]
      #index = which(T$gene=='Srebf1')
      #source('functions.R')
      #plot.genes.examples(index=index, T=T, Figure=TRUE, folder=folder, summary.plot=TRUE);
      
      #kk = which(keep$best.model==2 & keep$half.life>8 & keep$cv.half.life<0.3 & keep$qv.premRNA<0.001 & keep$relamp.premRNA>0.3)
      
      
      ### representive examples of M3
      folder = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/Examples/model_3/';
      #kk = which(Tt$BIC.best.model==3)
      mm = match(c('Fus', 'Hdac3', 'Loxl1', 'Exog', 'Hnrnpc'), keep$gene)
      kk = c(mm[which(!is.na(mm)==TRUE)])
      #kk = which(keep$best.model==3 & (keep$cv.half.life<0.4 | is.na(keep$cv.half.life)))
      index = c(keep$index[kk], which(T$gene == 'Hdac3'))
      which(keep$gene=='Hnrnpc')
      
      source('functions.R')
      index = match(c('Fus', 'Hdac3', 'Loxl1', 'Exog', 'Hnrnpc'), T$gene) 
      #index = which(T$gene=='Fus');
      plot.genes.examples(index=index, T=T, Figure=TRUE, folder=folder, summary.plot = TRUE);
      #index = which(T$gene=='Hnrnpc'); plot.genes.examples(index, T, Figure=TRUE, folder=folder, summary.plot = TRUE, shape.test.M3=FALSE);
      
      kk = which(keep$best.model==3 & keep$qv.mRNA<0.05 & keep$qv.premRNA>0.4 & keep$relamp.mRNA>0.25);
      index = keep$index[kk]
      source('functions.R')
      plot.genes.examples(index, T, Figure=TRUE, folder=folder, summary.plot=TRUE);
      
      ### representive examples of M4
      ### positive controls
      source('functions.R')
      folder = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/Examples/model_4/'
      mm = match(c('Arntl', 'Npas2','Per1', 'Per3', 'Cry1','Nr1d1', 'Gys2', 'Ppard', 'Eps8l2', 'Rbm3', 'Insig2', 'Ddo', 'Trfc'), keep$gene)
      kk = c(mm[which(!is.na(mm)==TRUE)])
      index = keep$index[kk]
      source('functions.R');index = which(T$gene=='Per3')
      
      source('functions.R');
      #index = match(c('Per1', 'Per3', 'Cry1', 'Npas2', 'Arntl', 'Nr1d1'), T.cc$gene) 
      index = match(c('Per3', 'Per1', 'Cry1', 'Npas2', 'Arntl', 'Nr1d1', 'Tef'), T$gene) 
      plot.genes.examples(index=index, T=T, Figure=TRUE, folder=folder, summary.plot = TRUE);
      
      dd = (keep$phase.mRNA-keep$phase.premRNA)%%24; rr = keep$relamp.mRNA/keep$relamp.premRNA;
      jj = which(keep$best.model==4 & keep$relamp.premRNA>0.2 & ((dd>=6 & dd<22) | rr>=4) & keep$qv.mRNA<0.01 & keep$qv.premRNA<0.01)
      kk = unique(c(kk, jj))
      index = keep$index[kk]
      
      jj = which(keep$best.model==4 & keep$relamp.premRNA>0.1 & (rr>=1.5) & keep$qv.mRNA<0.05 & keep$qv.premRNA<0.05)
      jj = jj[order(-rr[jj])]
      index = keep$index[jj]
      
      source('functions.R');
      index = match(c('Fas', 'Ddo'), T$gene)
      #index = index[which(T$qv.rpkm.mRNA[index]<10^(-4) & T$rel.amp.rpkm.mRNA[index]>0.2
      #                    & T$qv.rpkm.premRNA[index]<10^(-4) & T$rel.amp.rpkm.premRNA[index]>0.1)]
      #relamps = T$rel.amp.rpkm.mRNA[index]
      #diff.amps = T$rel.amp.rpkm.mRNA[index] /T$rel.amp.rpkm.premRNA[index]
      #diff.phase = (T$phase.rpkm.mRNA[index] - T$phase.rpkm.premRNA[index])%%24
      #index = index[which(diff.phase>5)]
      plot.genes.examples(index, T, Figure=TRUE, folder=folder, summary.plot = TRUE);
      
      yy = (T$phase.rpkm.mRNA - T$phase.rpkm.premRNA)%%24 
      yy[which(T$gene=='Cbs')]
      
    }
  }
  
}

##########################
##########################
#### Figure S2 : Valid method with simulation and parameter identifiability analysis 
##########################
##########################
Figure_S2 = FALSE
if(Figure_S2)
{
  #### Valid model selection and parameter estimation with simulated data
  Method.validation.with.simulation = FALSE
  if(Method.validation.with.simulation)
  {
    folder = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/';
    
    res.version = '_total_counts_s1_fitting_norm_parameters_constrained_gamma_v4';
    #res.version = '_total_counts_s5_all_v1';
    load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
    source('functions.R')
    
    BIC = my.model.selection.loglike(T=T, method='BIC', outliers = FALSE)
    #AIC = my.model.selection.loglike(T=T, method='AIC', outliers = TRUE)
    #AICc = my.model.selection.loglike(T=T, method='AICc', outliers = TRUE)
    T = data.frame(T, BIC, stringsAsFactors = FALSE);
    T = T[which(!is.na(T$BIC.best.model)==TRUE), ]
    #Tt = T
    
    ###### Calculate the Parameter Combination of true values which is to be compared afterwards 
    normalize.parameteres = TRUE
    if(normalize.parameteres)
    {
      xx = T$splicing.k/T$gamma
      T$splicing.k = xx;
      #head(T$gamma)
      w = 2*pi/24;
      xx = T$eps.gamma*T$gamma/sqrt(w^2+T$gamma^2);
      T$eps.gamma = xx;
      xx = T$phase.gamma+atan2(w, T$gamma)/w;
      kk = which(T$true.model==1 |T$true.model==2)
      xx[kk] = 0;
      T$phase.gamma = xx;
    }
    
    for(model in c(1:4))
    {
      print(length(which(T$BIC.best.model==model & T$true.model==model))/length(which(T$true.model==model)))
      #kk = which(T$true.model==model & T$BIC.best.model!=2)
      print(length(which(T$BIC.best.model==model & T$pval.rpkm.mRNA<0.01 & T$true.model==model))/length(which(T$true.model==model  & T$pval.rpkm.mRNA<0.01)))
      #hist(log(2)/T$gamma[kk])
    }
    
    #### PLOT for model selection reulst and compared true value of parameters with estimated 
    Plot.Results = TRUE
    if(Plot.Results)
    {
      source('functions.R')
      library(stringr)
      #selecting.approach = c('BIC.best.model','AICc.best.model','AIC.best.model')
      selecting.approach = c('BIC.best.model')
      palette(c('gray','green3','tomato','black'))
      #palette(c('gray','green3','tomato','black'))
      #source('functions.R')
      
      pdf.name = paste(folder, 'my_simulated_data_with_fits_results_BIC_', res.version, '.pdf', sep='')
      pdf(pdf.name, width = 2., height = 2.)
      par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
      
      Prediction.Quality = matrix(NA, nrow = 4, ncol =4)
      for(model in 1:4){
        #h = hist(T$LRT.best.model[((model-1)*X+1):(model*X)], breaks = seq(0.5,4.5,by = 1), plot = FALSE)
        h = hist(T$BIC.best.model[which(T$true.model==model)], breaks = seq(0.5,4.5,by = 1), plot = FALSE)
        Prediction.Quality[model,] = h$counts/sum(h$counts)
      }
      barplot(t(Prediction.Quality), col = 1:4, names.arg = c('', '', '', ''), border = NA, main = NA)
      dev.off()
      
      
      ###### wrongly calssied genes in M3 and M4
      pdf.name = paste(folder, 'my_simulated_data_with_fits_results_BIC_model_3.pdf', sep='')
      pdf(pdf.name, width = 1.8, height = 1.8)
      par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3, pty='s')
      
      jj = which(T$true.model==3)
      plot(T$rel.amp.rpkm.mRNA[jj], -log10(T$pval.rpkm.mRNA[jj]), type='n', xlab=NA, ylab=NA, main=NA, axes=FALSE)
      for(model in c(1:4))
      {
        kk = which(T$true.model==3 & T$BIC.best.model==model);
        points(T$rel.amp.rpkm.mRNA[kk], -log10(T$pval.rpkm.mRNA[kk]),  col=model, cex=0.5)
      }
      abline(h=(-log10(0.05)), col='darkblue', lty=1, lwd=1.0)
      abline(v=0.05, col='darkblue', lty=1, lwd=1.0)
      axis(1,cex.axis =0.7, at=seq(0, 1, by=0.2))
      axis(2, at=seq(0, 25, by=5), las=1,cex.axis = 0.7)
      box()
      dev.off()
      
      pdf.name = paste(folder, 'my_simulated_data_with_fits_results_BIC_model_4.pdf', sep='')
      pdf(pdf.name, width = 1.8, height = 1.8)
      par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3, pty='s')
      cex.axis = 0.7
      jj = which(T$true.model==4)
      plot((T$phase.rpkm.mRNA[jj]-T$phase.rpkm.premRNA[jj])%%24, T$rel.amp.rpkm.mRNA[jj]/T$rel.amp.rpkm.premRNA[jj], type='n', ylim = c(0, 2.0),
           xlab=NA, ylab=NA, main=NA, axes=FALSE, log='')
      #plot(T$rel.amp.rpkm.mRNA[kk], -log10(T$pval.rpkm.mRNA[kk]), type='p', xlab=NA, ylab=NA, main=NA, axes=TRUE)
      for(model in c(1:4)){
        kk = which(T$true.model==4 & T$BIC.best.model==model);
        ## transform eps combination back to real relative amplitude 
        yy = T$eps.gamma[kk]*sqrt(1+(w/T$gamma[kk])^2)/T$rel.amp.rpkm.premRNA[kk]
        #points((T$phase.gamma[kk]-T$phase.int[kk])%%24, yy, col=model, cex=0.5)
        points((T$phase.rpkm.mRNA[kk]-T$phase.rpkm.premRNA[kk])%%24, T$rel.amp.rpkm.mRNA[kk]/T$rel.amp.rpkm.premRNA[kk], col=model, cex=0.5)
      }
      xx = seq(0, 6, by=0.1);yy=cos(xx*w);points(xx, yy, col='blue', lwd=1.5, type='l')
      #abline(h=1, col='blue', lty=1, lwd=1.0)
      #abline(v=6, col='blue', lty=1, lwd=1.0)
      axis(1, at=seq(0, 24, by=6), cex.axis = cex.axis)
      axis(2, at=seq(0, 4, by=0.5), las=1,cex.axis = cex.axis)
      box()
      dev.off()
      
      ##################
      ### Parameter Comparison including degradation rate : gamma
      ##################
      #sd.noise = c(0.1, 0.25, 0.5)
      #set.global.sigma();
      #filter = TRUE
      source('functions.R')
      pdf.name = paste(folder, 'my_simulated_data_with_fits_results_BIC_AIC_LRT', res.version, '.pdf', sep='')
      pdf(pdf.name, width = 4.5, height = 3.5)
      
      compare.parameters.RNA.seq.total.fake(T=T, plot.parameters = TRUE)
      
      dev.off()
    }
    
  }
  
  #### Search.parameter.regions.non.identifiability for M3
  Search.parameter.regions.non.identifiability = FALSE
  if(Search.parameter.regions.non.identifiability)
  {
    #res.version = '_total_counts_s3_nonident_model3_v1';
    res.version = '_total_counts_s4_nonident_model3_v1';
    load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', 
                      res.version, '.Rdata',sep = ''))
    w = 2*pi/24;
    Transform.parameters = TRUE
    if(Transform.parameters){
      #eps.gamma.r = eps.gamma*sqrt(1+w^2/gamma^2); phase.gamma.r = (phase.gamma-atan2(w, gamma)/w)%%24;
      T$eps.gamma.m3 = T$eps.gamma.m3*sqrt(1+w^2/T$gamma.m3^2); T$phase.gamma.m3 = (T$phase.gamma.m3-atan2(w, T$gamma.m3)/w)%%24;
    }
    
    ### verify transformed phase and rel.amp for M3
    normalize.parameteres = FALSE
    if(normalize.parameteres)
    {
      xx = T$eps.gamma*T$gamma/sqrt(w^2+T$gamma^2); T$eps.gamma = xx;
      xx = (T$phase.gamma+atan2(w, T$gamma)/w)%%24; T$phase.gamma = xx;
      
      par(mfcol=c(1, 2))
      plot(T$eps.gamma, T$eps.gamma.m3, cex=0.2)
      abline(0, 1, col='red')
      
      kk = grep('alpha_', T$gene)
      plot(T$phase.gamma[kk], T$phase.gamma.m3[kk], cex=0.2)
      abline(0, 1, col='red')
      abline(-24, 1, col='red')
      abline(24, 1, col='red')
    }
    
    #### Plot selected mdoels and estimated parameters
    Plot.Non.Identifiability.Results = FALSE
    if(Plot.Non.Identifiability.Results)
    {
      #a.dispersion = 
      alphas = c(10^-4, 10^-3, 10^-2);
      source('functions.R')
     
      #### Here we just show one case with alpha=0.0001 for M3
      m0 = 1; model=3;
      cat(m0, '\n');
      
      index0 = grep(paste('alpha_', m0, sep=''), T$gene);
      index.par = unique(unlist(strsplit(as.character(T$gene[index0]), '_'))[5+6*c(0:(length(index0)-1))])
      grids = matrix(NA, nrow = length(index.par), ncol = 6);
      for(m1 in 1:length(index.par))
      {
        #m1 = 185
        ii = index.par[m1];
        kk = grep(paste('alpha_', m0, '_fake_m', model,'_', ii, '_r', sep=''), T$gene)
        eps.true = T$eps.gamma[kk[1]];
        gamma.true = log(2)/T$gamma[kk[1]];
        #gamma.true = T$gamma[kk[1]];
        ## error meric
        errs = (log(2)/T$gamma.m3[kk]/gamma.true);
        #errs = T$gamma.m3[kk]/gamma.true;
        #errs = errs[which(!is.na(errs)==TRUE)];
        #errs = errs[which(!is.na(errs)==TRUE)];
        #err.estimate = mean(abs(T$gamma.m3[kk]-gamma.true)/gamma.true)
        err.estimate = mean(log2(errs), na.rm = TRUE);
        grids[m1, ] = c(gamma.true, eps.true, err.estimate, sd(errs), min(T$gamma.m3[kk], na.rm = TRUE), max(T$gamma.m3[kk], na.rm = TRUE))
      }
      
      ##tranform half-life in log2 scale 
      yy = data.frame(grids)
      yy[, 1] = log2(yy[, 1]) 
      colnames(yy) = c('gamma.true', 'eps.true', 'err.gamma.estimate', 'sd.err', 'min.gamma', 'max.gamma')
      
      grids = yy;
      library(lattice)
      
      folder = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/analysis_non_identifiability/';
      pdf.name = paste(folder, 'non_ident_simulated_data_', res.version, '_alpha_', alphas[m0], '_half_life.pdf', sep='')
      pdf(pdf.name, width = 5, height = 2.4)
      par(cex =0.7, mar = c(1.0,1,1,0.8)+0.1, mgp = c(1.,0.5,0),las = 0, tcl = -0.3)
      
      x.scale <- list(at=seq(-2, 4, by=1))
      y.scale <- list(at=seq(0.1, 1, by=0.1))
      w = 2*pi/24;
      rr1 = seq(log(2)/(5/60), log(2)/24 , length.out=200); zz1 = 0.1*rr1; 
      rr2 = seq(log(2)/(5/60), log(2)/24 , length.out=200); zz2 = 0.05/rr2*sqrt((w^2+rr2^2)+1) #zz2 = 0.05*rr2*(1+0.5/(rr2*w/log(2)));  
      #levelplot(z ~ x * y, scales=list(x=x.scale, y=y.scale))
      levelplot(err.gamma.estimate~gamma.true*eps.true, grids, xlab='', ylab='', main=NA, sub=NA, region = TRUE, contour=FALSE,  
                scales=list(x=x.scale, y=y.scale, xlab=list(cex=.01), cex.axis=0.01, cex.label=0.0),
                col.regions=colorRampPalette(c("blue","green", "red")), 
                #at=seq(range(grids[,3])[1], range(grids[,3])[2], by=0.5)
                colorkey = list(at=seq(3, -4, by=-0.5), width=1, tck=0.5, labels=list(cex=0.7)),
                at = c(range(grids[,3])[1], -2, -1, -0.5,0, 0.5, 1, 2, range(grids[,3])[2]),
                panel = function(...){
                  panel.levelplot(...)
                  #panel.abline(v = log10(2))
                  #panel.abline(v = 0.4)
                  #panel.xyplot((rr1), zz1, col='orange', type='l', lwd=2.0)
                  #panel.xyplot((rr2), zz2, col='orange', type='l', lwd=2.0)
                }
      )
      
      dev.off()
      
      ### Profilelikelihood plot for 3 different examples
      ratios = log2(T$gamma.m3/T$gamma)
      #kk = which(T$gamma>3 & T$eps.gamma<0.2)
      #jj = match(c('gamma', 'eps.gamma', 'phase.gamma', 'gamma.m3', 'eps.gamma.m3', 'phase.gamma.m3'), colnames(T))
      #head(T[kk, c(1, jj)])
      ZT.int = grep('.count.premRNA', colnames(T))
      ZT.ex = grep('.count.mRNA', colnames(T))
      examples = paste('alpha_1_fake_m3_', c('78_r29', '82_r41', '17_r10', '51_r10'), sep='')
      index = match(examples, T$gene)
      
      folder = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/analysis_non_identifiability/'
      source('functions.R')
      
      for(kk in index) 
      {
        Identifiablity.analysis.4gamma.M3.with.plot.result(T=T, gene.index = kk, model=3, folder=folder, i.ex = ZT.ex, i.int = ZT.int, nb.test.iden = 40)
      }
      
      ###########
      ### Find the differences between identifiable and non-identifiable regions
      ###########
      Try.shape.character.from.non.identifiable.examples = FALSE
      if(Try.shape.character.from.non.identifiable.examples)
      {
        model = 3;
        zt.p = seq(0, 24, by=0.1)
        params = c(0.5, 0.5, 20, 10, 0.2, 0, 13, 1)
        eps.gamma = 0.6; phase.gamma=12;
        
        gamma = 
          eps.gamma.r = eps.gamma*sqrt(1+w^2/gamma^2); phase.gamma.r = (phase.gamma-atan2(w, gamma)/w)%%24;
        #s = compute.s.beta(t = zt, Min = params[1], Amp = params[2], phase = params[3], beta = params[4]); ## always the same
        rr = log(2)/2;m1 = compute.m.beta(t=zt.p, rr, params[2],  params[3], params[4]*rr,
                                          params[5], params[6], params[7], params[8])
        rr = log(2)/(0.25);m2 = compute.m.beta(t=zt.p, rr, params[2],  params[3], params[4]*rr,
                                               params[5], params[6], params[7], params[8])
        rr = log(2)/(8);m3 = compute.m.beta(t=zt.p, rr, params[2],  params[3], params[4]*rr,
                                            params[5], params[6], params[7], params[8])
        generate.cosin = function(m, zt.p){
          m.mean = mean(m);
          m = m/m.mean;
          m.rel = (max(m)-min(m))/2; m.phi = zt.p[which(m==max(m))][1];
          return(m.mean*(1+m.rel*cos(2*pi/24*(zt.p-m.phi))));
        }
        lims = range(c(m1, m2, m3));
        par(mfrow=c(1, 3))
        plot(zt.p, m1, type='l', col='green', lwd=2.0);points(zt.p, generate.cosin(m=m1, zt.p), type='l', col='green', lty=2, lwd=2.0);
        plot(zt.p, m2, type='l', col='blue', lwd=2.0);points(zt.p, generate.cosin(m=m2, zt.p), type='l', col='blue', lty=2, lwd=2.0);
        plot(zt.p, m3, type='l', col='red', lwd=2.0);points(zt.p, generate.cosin(m=m3, zt.p), type='l', col='red', lty=2, lwd=2.0);
      }
     
      
    }
    
    
  }
  
}

#############
#### Figure 2: Model selection results and comparison with PA-test (Sara)
#############
Figure_2 = FALSE
if(Figure_2)
{
  #res.version = '_total_counts_all_genes_norm_params_v6';
  #load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
  #res.version = '_total_counts_all_genes_norm_params_v6';
  #load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
  source('functions.R')
  #BIC = my.model.selection.loglike(T=T, method='BIC', outliers = TRUE)
  #T = data.frame(T, BIC, stringsAsFactors = FALSE);
  #T = T[which(!is.na(T$BIC.best.model)==TRUE), ]
  #Tt = T
  res.version = '_total_counts_all_genes_norm_params_final_v7_All_summary';
  res.version.old = '_total_counts_all_genes_norm_params_v6_OLD';
  load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
  load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', 
                                        res.version.old, '.Rdata',sep = ''))
  T$BIC.best.model = T$best.model;
  #load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_RNA_seq_total_counts_table_sx', res.version, '.Rdata',sep = ''))
  attribute.global.variable.colors();
  
  ##########
  #### Partition of models and percentages of mRNA post-transcriptionally regulated
  ##########
  Tt = T.old;
  source('functions.R')
  length(which(Tt$BIC.best.model==1))/nrow(Tt);
  length(which(Tt$BIC.best.model==2))/nrow(Tt);
  length(which(Tt$BIC.best.model==3))/nrow(Tt);
  length(which(Tt$BIC.best.model==4))/nrow(Tt);
  table(T$best.model)/sum(table(T$best.model))
  
  #cols=c('bisque','olivedrab','palevioletred3','purple4')
  cols = c(1:4)
  counts = matrix(NA, nrow=3, ncol=3)
  counts[1,1] = length(which(Tt$BIC.best.model==2))
  counts[2,1] = length(which(Tt$BIC.best.model==3))
  counts[3,1] = length(which(Tt$BIC.best.model==4))
  
  counts[1,2] = length(which(Tt$BIC.best.model==2 & Tt$qv.rpkm.mRNA<0.05))
  counts[2,2] = length(which(Tt$BIC.best.model==3 & Tt$qv.rpkm.mRNA<0.05))
  counts[3,2] = length(which(Tt$BIC.best.model==4 & Tt$qv.rpkm.mRNA<0.05))
  
  counts[1,3] = length(which(Tt$BIC.best.model==2 & Tt$qv.rpkm.mRNA<0.05 & Tt$rel.amp.rpkm.mRNA>0.1))
  counts[2,3] = length(which(Tt$BIC.best.model==3 & Tt$qv.rpkm.mRNA<0.05 & Tt$rel.amp.rpkm.mRNA>0.1))
  counts[3,3] = length(which(Tt$BIC.best.model==4 & Tt$qv.rpkm.mRNA<0.05 & Tt$rel.amp.rpkm.mRNA>0.1))
  #counts <- table(mtcars$vs, mtcars$gear)
  counts.old = counts[, 2]
  
  counts = c(length(which(T$best.model==2 & T$qv.log2rpkm.mRNA <0.05)), 
             length(which(T$best.model==3 & T$qv.log2rpkm.mRNA <0.05)),
             length(which(T$best.model==4 & T$qv.log2rpkm.mRNA <0.05)))
  counts = cbind(counts, counts/sum(counts))
  
  
  pdf.name = paste("/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/Partition_models.pdf", sep='')
  pdf(pdf.name, width=2.6, height=2.2)
  par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,1.5,3)+0.1, tcl = -0.3)
  
  palette(c('gray','green3','tomato','black'))
  #colnames(counts) = c('M2', 'M3', 'M4')
  barplot(counts[,2], main=NA, ylab= NA, names.arg=c("", "", ""), col='white', space = 0.7, ylim = c(0, 0.7),
          axes = FALSE)
  axis(side = 4)
  mtext(side = 4, line = 2, '', las=3, cex=1.0)
  par(new=TRUE)
  barplot(counts[,1], main=NA, ylab= NA, names.arg=c("", "", ""), col=cols[c(2:4)], space = 0.7, ylim = c(0, 4000),
          legend.text = c('RS-CD', 'CS-RD', 'RS-RD'), args.legend = list(bty='n'), axes = FALSE)
  axis(2,at = c(0, 1000, 2000, 3000, 4000),las=2,cex.axis = 1.)
  
  dev.off()
  
  ### Comparison with Westermark results
  Compare.Our.results.vs.Westermark = FALSE
  if(Compare.Our.results.vs.Westermark)
  {
    #### Westermark method
    #ggs = T$gene[which((T$BIC.best.model==3|T$BIC.best.model==4) & T$qv.rpkm.mRNA<0.05)]
    com1 = read.table('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/Tables/PAtest_result_v2.txt', header = TRUE, sep='\t', as.is = c(1))
    com2 = read.table('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/Tables/PAtest_result_withEstimated_HL.txt', header = TRUE, sep='\t', as.is = c(1))
    com1 = data.frame(com1, stringsAsFactors = FALSE);  
    com2 = data.frame(com2, stringsAsFactors = FALSE);
    overlap = intersect(com1$gene, com2$gene);
    com1 = com1[match(overlap, com1$gene), ]
    com2 = com2[match(overlap, com2$gene), ]
    #com1 = com1[which(!is.na(match(com1$MGI_symbol, T$gene[which(T$qv.rpkm.mRNA<0.05)]))==TRUE), ]
    mm1 = match(com1$gene, T$gene); 
    mm2 = match(com2$gene, T$gene);  
    yy1 = data.frame(T$best.model[mm1], log10(com1$PA.test.p.val), com1$halflife);
    yy2 = data.frame(T$best.model[mm2], log10(com2$PA.test.p.val), com2$est.HL);
    #kk = which(T$qv.rpkm.mRNA[mm1]<0.05); yy1 = yy1[kk,]; 
    #kk = which(T$qv.rpkm.mRNA[mm2]<0.05); yy2 = yy2[kk,]; 
    #yy = yy[which(yy[,3]<15), ]
    #selecting.genes.with.half.lives = FALSE
    #if(selecting.genes.with.half.lives) {}
    #mm = match(com1$gene, T$gene)
    #yy = data.frame(T$BIC.best.model[mm], log10(com1$PA.test.p.val))
    pvals = seq(-8, 0, by=0.4)
    percent = c()
    for(pval in pvals)
    {
      #gg1 = com1$MGI_symbol[which(com1$Deg_Q_Val<qv & com1$Abun_Q_Val<0.05)]
      #gg1 = com1$MGI_symbol[which(com1$Deg_Q_Val<qv & T$qv.rpkm.mRNA)]
      #index = which(yy[,2]<=pval);
      #percent = c(percent, length(which(yy[index, 1]>2))/length(index));
      #percent = rbind(percent, 
      #                c(pval, length(which(yy1[,2]<=pval & yy1[,1]>2))/(length(which(yy1[,2]<=pval))*length(which(yy1[,1]>2))/nrow(yy1)), 
      #                  length(which(yy2[,2]<=pval & yy2[,1]>2))/(length(which(yy2[,2]<=pval))*length(which(yy2[,1]>2))/nrow(yy2))))
      percent = rbind(percent, 
                      c(pval, 100*length(which(yy1[,2]<=pval & yy1[,1]>2))/(length(which(yy1[,2]<=pval))), 
                        100*length(which(yy2[,2]<=pval & yy2[,1]>2))/(length(which(yy2[,2]<=pval)))))
      #percent2 = c(percent2, length(which(yy2[,2]<=pval & yy2[,1]>2))/(length(which(yy2[,2]<=pval))*length(which(yy2[,1]>2))/nrow(yy2)))
    }
    
    folder = "/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/"
    pdfname = paste(folder, 'Sarah_Comparison_enrichment_v5.0.pdf', sep='')
    pdf(pdfname, width=2.2, height=2.2)
    par(cex = 0.7, las = 0, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    plot(percent[, c(1:2)], type='p', lty=2, main=NA, log='', xlab=NA, ylab=NA, ylim = c(0, 100),
         col=c('black'), axes=FALSE,  cex=0.5)
    points(percent[, c(1:2)], type='l', lty=1, lwd=1.5, col='black')
    points(percent[, c(1,3)], type='p', lty=1, cex=0.5, col='darkblue', pch=1)
    points(percent[, c(1,3)], type='l', lty=1, lwd=1.5, col='darkblue')
    abline(h=100*length(which(yy1[,2]<=1 & yy1[,1]>2))/length(which(yy1[,2]<=1)),col='gray',lwd=1.,lty=2)
    axis(1,at=seq(-8, 0, by=2),cex.axis =1.0)
    axis(2, at=seq(0, 100, by=20), las=1,cex.axis = 1.0)
    box()
    
    legend(-8, 90, legend = c('', ''), cex=1.0, lwd=1.5, lty=1, col = c('black', 'darkblue'), border = NA, bty = 'n')
    
    dev.off()
    
  }
  
  
  #### Compare phases, amplitudes for mRNAs and pre-mRNAs in M2-M4
  Compare.Phase.amplitudes.mRNA.premRNA.M2.M3.M4 = FALSE
  if(Compare.M2.M3.M4)
  {
    #### Here we used phases, amplitudes and relative amplitudes computed with log2 scale of mRNA RPKM
    #load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_sel_alphas_log2_RPKM', 
    #                data.version, '.Rdata', sep=''))
    folder = "/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/"
  
    #sel1 = match(T$gene[which(T$BIC.best.model==1)], R.norm$gene);
    #sel2 = match(T$gene[which(T$BIC.best.model==2 & T$BIC.prob.m2>0.5 &T$qv.rpkm.mRNA<0.05)], R.norm$gene);
    #sel3 = match(T$gene[which(T$BIC.best.model==3 & T$BIC.prob.m3>0.5 &T$qv.rpkm.mRNA<0.05)], R.norm$gene);
    #sel4 = match(T$gene[which(T$BIC.best.model==4 & T$BIC.prob.m4>0.5 &T$qv.rpkm.mRNA<0.05)], R.norm$gene);
    
    #dd2 = (R.norm$phase.log2rpkm.mRNA[sel2] -  R.norm$phase.log2rpkm.premRNA[sel2])%%24
    #dd4 = (R.norm$phase.log2rpkm.mRNA[sel4] -  R.norm$phase.log2rpkm.premRNA[sel4])%%24
    #cols = c(1:4)
    
    ### Phases of mRNA and pre-mRNA
    for(model in c(2:4))
    {
      #phases = eval(parse(text=paste('R.norm$phase.log2rpkm.mRNA[sel', model, ']',  sep='')));
      #pdfname = paste(folder, 'Phases_rhythmic_mRNAs_M', model, '.pdf', sep='')
      #pdf(pdfname, width=1.5, height=1.5)
      #par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
      #circular_phase24H_histogram(phases, col=cols[model], cex.axis=0.4, cex.lab=0.25, lwd=0.5)
      #dev.off()
      kk = which(T$BIC.best.model==model & T$qv.log2rpkm.mRNA<0.05)
      if(model!=3)
      {
        phases = as.numeric(T$phase.log2rpkm.premRNA[kk])
        pdfname = paste(folder, 'phases_rhythmic_premRNAs_model_', model, '.pdf', sep='')
        pdf(pdfname, width=1.8, height=1.8)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        circular_phase24H_histogram(phases, col= col.int, cex.axis=0.4, cex.lab=0.25, lwd=0.7)
        
        dev.off()
      }
      
      phases = T$phase.log2rpkm.mRNA[kk]
      pdfname = paste(folder, 'phases_rhythmic_mRNAs_model_', model, '.pdf', sep='')
      pdf(pdfname, width=1.8, height=1.8)
      par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
      circular_phase24H_histogram(phases, col= col.ex, cex.axis=0.5, cex.lab=0.01, lwd=0.7)
      dev.off()
      
      Density.amplitudes = FALSE
      if(Density.amplitudes){
        amps.m = R.norm$amp.log2rpkm.mRNA[kk];amps.p = R.norm$amp.log2rpkm.premRNA[kk];
        
        pdfname = paste(folder, 'Ampls_density_mRNAs_premRNA_model_', model, '.pdf', sep='')
        pdf(pdfname, width=2.0, height=1.4)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        #circular_phase24H_histogram(phases, col= model, cex.axis=0.25, cex.lab=0.01, lwd=0.4)
        
        plot(density(amps.m), type='l', main=NA, log='', xlab=NA, ylab=NA, xlim = c(0, 4),
             col=col.ex, axes=FALSE, lwd=1.5)
        abline(v=median(amps.m), col=col.ex, lwd=1.5)
        if(model!=3) {points(density(amps.p), type='l', lty=1, lwd=1.5, col=col.int);abline(v=median(amps.p), col=col.int,lwd=1.5)}
        
        axis(2, at=seq(0, 3, by=1), las=1,cex.axis = 1.0)
        axis(1,cex.axis =1.0, tick = TRUE, labels = FALSE, lwd.ticks = 0)
        #box()
        if(model==4)  axis(1,at=seq(0, 5, by=1),cex.axis =1.0)
        dev.off() 
      }
    }
    
    ### Amplitudes
    sel1 = which(T$BIC.best.model==1);
    sel2 = which(T$BIC.best.model==2 & T$qv.log2rpkm.mRNA<0.05);
    sel3 = which(T$BIC.best.model==3 & T$qv.log2rpkm.mRNA<0.05);
    sel4 = which(T$BIC.best.model==4 & T$qv.log2rpkm.mRNA<0.05);
    
    yy = data.frame(c(rep(2, length(sel2)), rep(3, length(sel3)), rep(4, length(sel4))), T$amp.log2rpkm.mRNA[c(sel2, sel3, sel4)]);
    
    VioPLOT = FALSE
    pdfname = paste(folder, 'Amplitudes_log2_Rhythmic_mRNA_M2_3_4.pdf', sep='')
    pdf(pdfname, width=2., height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    if(VioPLOT)
    {
      require('vioplot')
      plot(yy[, 1], yy[,2], type='n', xlab = NA, ylab=NA, main=NA, axes = FALSE, xlim=c(1.5, 4.5), ylim = c(0, 3))
      vioplot(yy[which(yy[,1]==2),2],  col = 2+4, at=2, rectCol = 2, cexMed=2, boxwex = 2,
              names=NA, lwd=0.7, colMed = 2, pchMed = '-', horizontal=FALSE, wex=0.8, lty=1, add=TRUE)
      vioplot(yy[which(yy[,1]==3),2],  col = 3+4, at = 3, rectCol = 3, 
              names=NA, lwd=0.7, colMed = 3, pchMed = '-', horizontal=FALSE, wex=0.8, lty=1, add=TRUE)
      vioplot(yy[which(yy[,1]==4),2], col = 4+4, at = 4, rectCol = 4, 
             names=NA, lwd=0.7, colMed = 4, pchMed = '-', horizontal=FALSE, wex=0.8, lty=1, add=TRUE)
      axis(2, las=1,cex.axis = 1.)
      
    }else{
      boxplot(yy[,2] ~ yy[,1], main=NA, col=col.ex, names=c('', '', ''), ylim=c(0, 4), lwd=0.5, cex = 0.8)
      axis(2, las=1,cex.axis = 1.)
    }
    
    dev.off()
    
    ### mean expression (in log2) for pre-mRNAs and mRNAs
    yy = data.frame(c(rep(1, length(sel1)), rep(2, length(sel2)), rep(3, length(sel3)), rep(4, length(sel4))), 
                    T$mean.log2rpkm.premRNA[c(sel1, sel2, sel3, sel4)]);
    pdfname = paste(folder, 'mean_log2_premRNA.pdf', sep='')
    pdf(pdfname, width=2.0, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    boxplot(yy[,2] ~ yy[,1], main=NA, cex.main=0.7, lwd=0.7, cex=0.5, col=cols[1:4], names=c('M1', 'M2', 'M3', 'M4'))
    
    dev.off()      
    
    yy = data.frame(c(rep(1, length(sel1)), rep(2, length(sel2)), rep(3, length(sel3)), rep(4, length(sel4))), 
                    T$mean.log2rpkm.mRNA[c(sel1, sel2, sel3, sel4)]);
    pdfname = paste(folder, 'mean_log2_mRNA.pdf', sep='')
    pdf(pdfname, width=2.0, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    boxplot(yy[,2] ~ yy[,1], main=NA, cex.main=0.7, lwd=0.7, cex=0.5, col=cols[1:4], names=c('M1', 'M2', 'M3', 'M4'))
    
    dev.off()      
    
    ### relative amplitudes of premRNA 
    yy = data.frame(c(rep(2, length(sel2)), rep(3, length(sel3)), rep(4, length(sel4))), R.norm$rel.amp.log2rpkm.premRNA[c(sel2, sel3, sel4)]);
    
    pdfname = paste(folder, 'rel_amplitudes_log2_premRNA.pdf', sep='')
    pdf(pdfname, width=1.8, height=1.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    boxplot(yy[,2] ~ yy[,1], main=' rel.ampl of premRNA', cex.main=0.7, col=cols[2:4], names=c('M2', 'M3', 'M4'), ylim=c(0, 1), lwd=0.5)
    
    dev.off()
    
    ### rel.amplitude ratios
    yy = data.frame(c(rep(2, length(sel2)), rep(4, length(sel4))), 
                    R.norm$rel.amp.log2rpkm.mRNA[c(sel2, sel4)]/R.norm$rel.amp.log2rpkm.premRNA[c(sel2, sel4)]);
    
    pdfname = paste(folder, 'relampl_ratio_mRNA_vs_premRNA.pdf', sep='')
    pdf(pdfname, width=1.7, height=1.7)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    boxplot(yy[,2] ~ yy[,1], main=' rel.ampl ratio (mRNA/pre-mRNA)', cex.main=0.7, col=cols[c(2, 4)], names=c('M2', 'M4'), ylim=c(0, 4), lwd=0.5)
    abline(h=1, lwd=1.0, col='darkgray')
    
    dev.off()
    
  }
 

  Compare.Our.results.vs.PolyA.dynamics = FALSE
  if(Compare.Our.results.vs.PolyA.dynamics)
  {
    #T[which(T$gene=='Cyp2a4'), ]
    
    ### polyA length dynamics
    #com2 = read.csv('/Users/jiwang/Degradation_Liver/Database_2compare/polyA_length_dynamics.csv', header = FALSE)
    #com2 = com2[-c(1:4), ]
    #colnames(com2) = c('Probe.set', 'gene', paste('ZT', seq(0, 20, by=4), sep=''), 'pval', 'amp', 'phase', 'R2')
    #write.table(com2, file='/Users/jiwang/Degradation_Liver/Database_2compare/polyA_length_dynamics.txt', quote = FALSE, col.names = TRUE, row.names = FALSE)
    
    com2 = read.table(file='/Users/jiwang/Degradation_Liver/Database_2compare/polyA_length_dynamics.txt', header = TRUE, as.is = c(1, 2))
    mm = match(com2$gene, T$gene);
    yy = data.frame(T$BIC.best.model[mm], -log10(as.numeric(com2$pval)))
    yy = yy[which(!is.na(yy[,1])==TRUE),]
    #kk = which(T$qv.log2rpkm.mRNA[mm]<0.05)
    #yy = yy[kk,]
    
    cols = c(1:4)
    pdfname = paste(folder, 'Kojima_Comparison_v2.pdf', sep='')
    pdf(pdfname, width=2.0, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    boxplot(yy[,2] ~ yy[,1], main=NA, ylab=NA, cex.main=0.7, 
            col=cols[1:4], names=c('M1', 'M2', 'M3', 'M4'),lwd=0.7)
    abline(h=2, lwd=1., lty=1, col='blue')
    
    dev.off()
    
    #com2 = com2[which(!is.na(match(com2$gene, T$gene))==TRUE), ]
    index = match(com2$gene, T$gene)
    com2 = com2[which(T$qv.rpkm.mRNA[index]<0.05), ]
    mm = match(com2$gene, T$gene)
    counts = c(length(which(T$BIC.best.model[mm]==2))/length(mm), 
               length(which(T$BIC.best.model[mm]==3))/length(mm),
               length(which(T$BIC.best.model[mm]==4))/length(mm))
    counts = data.frame(counts)
    
    pdf.name = paste("/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Total/Total_counts/Overlap_polyA_dynamics.pdf", sep='')
    pdf(pdf.name, width=2.5, height=2.5)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,1.5,0.8)+0.1, tcl = -0.3)
    
    palette(c('gray','green3','tomato','black'))
    #colnames(counts) = c('M2', 'M3', 'M4')
    barplot(counts[,1], main=NA, ylab= 'nb of rhythmic mRNAs', names.arg=c("M2", "M3", "M4"), col=c(2:4), space = 0.7, ylim = c(0, 0.6),
            axes = FALSE, legend.text = c('RS-CD', 'CS-RD', 'RS-RD'), args.legend = list(bty='n'))
    axis(2, at = seq(0, 0.5, by=0.1),las=3,cex.axis = 1.)
    #box()
    
    dev.off()
    
  }
  
}

######################################
######################################
### Figure 3 estimated mRNA half-life and processing time and predicted mRBPs involved in RD 
######################################
######################################
FIGURE_3 = FALSE
if(FIGURE_3)
{
  #res.version = '_total_counts_all_genes_norm_params_v6';
  #### import tables selected
  #load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_RNA_seq_total_counts_table_sx', res.version, '.Rdata',sep = ''))
  
  #res.version.old = '_total_counts_all_genes_norm_params_v6_OLD';
  #load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', 
  #                res.version.old, '.Rdata',sep = ''))
  res.version = '_total_counts_all_genes_norm_params_final_v7_All_summary';
  load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
  
  #save(keep, file="Table_clean_to_Laura_usedforFigure.Rdata")
  
  T$BIC.best.model = T$best.model;
  
  source('functions.R')
  attribute.global.variable.colors();
  folder = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/'
  
  
  ###############
  ### plot half-life and processing time distribution and compare with datasets from other papers
  ###############
  Distribution.half.life.processing.time = FALSE
  if(Distribution.half.life.processing.time)
  {
    examples = c('Npas2', 'Arntl', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry1', 'Cry2', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Nfil3', 
                 'Bhlhe40', 'Bhlhe41')
    mm = match(examples, keep$gene); 
    cbind(examples, keep$half.life[mm], keep$cv.half.life[mm]);
    cutoff = 0.4;
    sel2 = which(keep$best.model==2 & (keep$cv.half.life<cutoff|is.na(keep$cv.half.life)) & keep$gamma.identifiability.L & keep$gamma.identifiability.R) #sel2 = which(keep$best.model==2 & keep$cv.half.life<cutoff & is.na(keep$non.identifiability.gamma))
    sel3 = which(keep$best.model==3 & keep$cv.half.life<cutoff & keep$gamma.identifiability.L & keep$gamma.identifiability.R)
    sel4 = which(keep$best.model==4 & (keep$cv.half.life<cutoff|is.na(keep$cv.half.life)) & keep$gamma.identifiability.L & keep$gamma.identifiability.R) #sel4 = which(keep$best.model==4 & keep$cv.half.life<cutoff & is.na(keep$non.identifiability.gamma))
    index = c(intersect(mm, sel2),  mm4 = intersect(mm, sel4));
    exs = data.frame(keep$gene[index], keep$best.model[index], keep$half.life[index], keep$cv.half.life[index], stringsAsFactors = FALSE)
    o1 = order(exs[,3])
    exs = exs[o1, ]
    
    cat(length(sel2), length(sel4), '\n');
    dens2 = density(keep$half.life[sel2], from = 0);
    dens3 = density(keep$half.life[sel3], from = 0);
    dens4 = density(keep$half.life[sel4], from = 0);
    sel1 = which(T$BIC.best.model==1 & T$BIC.prob.m1>0.5);
    #dens1 = density(T$mean.rpkm.mRNA[sel1]/T$mean.rpkm.premRNA[sel1]*10/60, from = 0)
    
    pdfname = paste(folder, 'Half_life_distribution_model_2_3_4_', res.version, '.pdf', sep='')
    pdf(pdfname, width=2.8, height=2.4)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    plot(dens2, type='l', main=NA, log='', xlab=NA, ylab=NA, xlim = c(-2, 20),
         col=2, axes=FALSE, lwd=2.0)
    points(dens4, type='l', col=4, lwd=2.0)
    #points(dens3, type='l', col=3, lwd=2.0)
    abline(v=median(keep$half.life[sel2]), lwd=1.5, col=2, lty=1)
    abline(v=median(keep$half.life[sel4]), lwd=1.5, col=4, lty=1)
    for(n in 1:nrow(exs))
    {
      points(exs[n, 3], 0.08+0.017*(n-1), type='p', pch=21, bg=exs[n, 2], col='white', cex=1.1);
      pos=2;offset=0.4;
      gg = exs[n,1]
      #if(n<=5){pos=4;}
      #if(gg=='Per2') {offset=0.5}
      #if(gg=='Bhlhe40') {offset=0.2; pos=4}
      text(exs[n, 3], 0.08+0.017*(n-1), exs[n, 1], col='black', cex=0.7, offset = offset, pos=pos);
    }
    axis(1, at=seq(0, 24, by=4), las=1,cex.axis = 1.0)
    axis(2,at=seq(0, 0.3, by=0.1), cex.axis=1.0)
    
    box()
    dev.off()
    
    cutoff = 0.4;
    sel2 = which(keep$best.model==2 & (keep$cv.splicing.time<cutoff|is.na(keep$cv.splicing.time & keep$gamma.identifiability.L & keep$gamma.identifiability.R))) #sel2 = which(keep$best.model==2 & keep$cv.half.life<cutoff & is.na(keep$non.identifiability.gamma))
    sel4 = which(keep$best.model==4 & (keep$cv.splicing.time<cutoff|is.na(keep$cv.splicing.time & keep$gamma.identifiability.L & keep$gamma.identifiability.R)))
    cat(length(sel2), length(sel4), '\n');
    dens2 = density(keep$splicing.time[sel2], from = 0);
    dens4 = density(keep$splicing.time[sel4], from = 0);
    
    pdfname = paste(folder, 'Processing_time_distribution_model_2_4_', res.version, '.pdf', sep='')
    pdf(pdfname, width=2.8, height=2.4)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    plot(dens2, type='l', main=NA, log='', xlab=NA, ylab=NA, xlim = c(0, 60),
         col=2, axes=FALSE, lwd=2.0)
    points(dens4, type='l', col=4, lwd=2.0)
    #points(dens1, type='l', col=1, lwd=2.0)
    abline(v=median(keep$splicing.time[sel2]), lwd=1.5, col=2, lty=1)
    abline(v=median(keep$splicing.time[sel4]), lwd=1.5, col=4, lty=1)
    axis(1, at=seq(0, 60, by=10), las=1,cex.axis = 1.0)
    axis(2,at=seq(0, 0.1, by=0.05), las=3, cex.axis=1.0)
    box()
    dev.off()
  }
  
  ###########################
  ### Compare estimated half-lives with other database
  ###########################
  Compare.half.life.database = FALSE
  if(Compare.half.life.database)
  {
    ###### comapre our estimated with other database
    databases = c('Fried', 'Shar', 'Schwa')
    cex = 0.2;
    for(db in databases)
    {
      #db = databases[1];
      pdfname = paste(folder, 'half-life_comparison_', db, res.version, '.pdf', sep='')
      cutoff = 0.4
      jj = eval(parse(text=paste('which(keep$best.model!=3 & keep$gamma.identifiability.L & keep$gamma.identifiability.R & keep$cv.half.life<', cutoff, '& !is.na(keep$', db, ')==TRUE)', 
                                 sep='')))
      jj1 = eval(parse(text=paste('which(keep$best.model!=3 &  keep$cv.half.life<', cutoff, '& !is.na(keep$', db, ')==TRUE)', 
                                 sep='')))
      print(length(jj))
      print(length(jj1))
      #jj = which(!is.na(keep[,kk])==TRUE)
      method = 'pearson'; R = eval(parse(text=paste('signif(cor(keep$half.life[jj], keep$', db, '[jj], method = method), d=2)', sep='')));
      cat(R, '\n');
      pdf(pdfname, width=1.8, height=1.8)
      par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3, pty='s')
      for(m in c(2:4))
      {
        kk = jj[which(keep$best.model[jj]==m)]
        yy1 = keep$half.life[kk]
        xx1 = eval(parse(text=paste('keep$', db, '[kk]', sep='')));
        if(m==2){
          plot(xx1, yy1, cex=cex, col=m, cex.lab=0.7, cex.main=0.7, xlim=c(0.5,24), ylim=c(0.5, 24),
               xlab=paste(db, sep=''), ylab=NA, main=paste('R = ', R, sep=''), log='xy', axes=FALSE)
        }else{
          points(xx1, yy1, cex=cex, col=m)
        }
      }
      box()
      axis(2,at=c(1, 2, 5, 10, 24), cex.axis=0.8)
      axis(1, at=c(1, 2, 5, 10, 24, 50), las=1,cex.axis = 0.8)
      abline(0, 1, lwd=1.5, col='darkgray')
      dev.off()
    }
    
    ### Compare different half-life databases
    load(file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/mRNAs_half_lives_databases.Rdata')
    method = 'pearson'; 
    mm = match(fried$gene, shar$gene); cor(fried$half.life, shar$half.life[mm], use = "na.or.complete")
    mm = match(fried$gene, schwa$gene); cor(fried$half.life, schwa$half.life[mm], use = "na.or.complete")
    mm = match(shar$gene, schwa$gene); cor(shar$half.life, schwa$half.life[mm], use = "na.or.complete")
    
    pdfname = paste(folder, 'half-life_comparison_Fried_vs_Shar.pdf', sep='')
    pdf(pdfname, width=1.8, height=1.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3, pty='s')
    mm = match(fried$gene, shar$gene); R = signif(cor(fried$half.life, shar$half.life[mm], use = "na.or.complete"), d=2)
    plot(fried$half.life, shar$half.life[mm], cex=0.03,  cex.lab=0.7, cex.main=0.7, xlim=c(0.5,24), ylim=c(0.5, 24),
         xlab=NA, ylab=NA, col='black',
         main=paste('R = ', R, sep=''), log='xy', axes=FALSE)
    axis(2,at=c(1, 2, 5, 10, 24), cex.axis=0.8)
    axis(1, at=c(1, 2, 5, 10, 24, 50), las=1,cex.axis = 0.8)
    abline(0, 1, lwd=1.2, col='darkgray', lty=1)
    box()
    
    dev.off()
    
  }
  
  #####################
  ### Compare estimated half-life and processing time with features of exon and intron
  #####################
  Correlation.half.life.processing.time.with.exon.intron.features = FALSE
  if(Correlation.half.life.processing.time.with.exon.intron.features)
  {
    ### import table containing half-life, processing time and all relavent features
    load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/half_life_splicing_time_genes_introns_exons_UTRs.Rdata',sep = ''))
    neweff = read.table('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/Tables/TE_nobiais_genes.txt', sep='\t', header = TRUE)
    #mm = match(neweff[,1], yy);
    #neweff$gene = xx[mm]
    #write.table(neweff, file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/Tables/TE_nobiais_genes.txt', sep='\t', col.names = TRUE, row.names = FALSE, quote=FALSE)
    mm = match(aa$gene, neweff$gene)
    aa$translation.effs = neweff$Translation_Efficiency[mm]
    aa$translation.effs = 2^aa$translation.effs
    
    cutoff = 0.4;
    sel2 = which(aa$best.model==2 & (aa$cv.half.life<cutoff | is.na(aa$cv.half.life))) #sel2 = which(keep$best.model==2 & keep$cv.half.life<cutoff & is.na(keep$non.identifiability.gamma))
    sel4 = which(aa$best.model==4 & (aa$cv.half.life<cutoff | is.na(aa$cv.half.life))) #sel4 = which(keep$best.model==4 & keep$cv.hal
    kk = c(sel2, sel4);
    sels = c(sel2, sel4)
    cor(log(aa$translation.effs[kk]), log(aa$half.life[kk]), use = "na.or.complete")
    cor((aa$translation.effs[kk]), (aa$half.life[kk]), use = "na.or.complete")
    cor((aa$translation.effs[kk]), (log(2)/aa$half.life[kk]), use = "na.or.complete")
    #cor((aa$translation.effs[kk]), aa$half.life[kk], use = "na.or.complete")
    cor(log(aa$translation.effs[kk]), log(1/(aa$half.life[kk])), use = "na.or.complete")
    #cor(log(aa$translation.effs[kk]), log(1/(aa$half.life[kk])), use = "na.or.complete")
    
    pdfname = paste(folder, 'Half_life_vs_translation_efficiency.pdf', sep='')
    pdf(pdfname, width=2.2, height=1.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    xx = log2(aa$translation.effs[sels]); yy = aa$half.life[sels];
    jj = which(xx<6 & !is.na(xx)==TRUE & !is.na(yy)==TRUE);xx=xx[jj];yy=yy[jj];
    cor(xx, yy, use = 'na.or.complete')
    plot(xx, yy, log='y', xlab=NA, ylab=NA, cex=0.1, axes = FALSE)
    #boxplot(yy[,2] ~ yy[,1], main=' rel.ampl of premRNA', cex.main=0.7, col=cols[2:4], names=c('M2', 'M3', 'M4'), ylim=c(0, 1), lwd=0.5
    #boxplot(log2(aa$translation.effs[kk]) ~ groups, ylim=c(-3.6, 4.5), axes=FALSE)
    #require('MASS');library(sfsmisc)
    fit = lm(log2(yy) ~ xx)
    tt = seq(-4, 5, length.out = 100); 
    points(tt, 2^(fit$coefficients[1]+fit$coefficients[2]*tt), type='l', lwd=1.0, col='blue')
    axis(1, at=seq(-4, 4, by=2), las=1,cex.axis = 1.0)
    axis(2,at=c(0.5, 1, 2, 5, 10, 20), cex.axis=1.0)
    box()
    
    dev.off()
    
    
    EXPLORATION_half.life.splicing.time.correlation.intron.exons = FALSE
    if(EXPLORATION.half.life.splicing.time.correlation.intron.exons)
    {
      Exon.Intron.table = FALSE
      if(Exon.Intron.table)
      {
        transcripts = read.table('/Users/jiwang/RNA_seq_Data/usful_mapping_files/transcript_mapping.txt', sep='\t', header = FALSE, as.is=c(1:3))
        exons = read.table('/Users/jiwang/RNA_seq_Data/usful_mapping_files/mm10_ens_exons_ucsc.bed', sep = '\t', header = FALSE, as.is = c(1,4,5, 6))
        introns = read.table('/Users/jiwang/RNA_seq_Data/usful_mapping_files/mm10_ens_introns_ucsc.bed', sep = '\t', header = FALSE, as.is = c(1,4,5, 6))
        utr = read.table('/Users/jiwang/RNA_seq_Data/mm10_annotation/mm10_3UTR.txt', sep='\t', header = FALSE, as.is = c(1, 4:6))
        extract_transcript = function(x){
          return(unlist(strsplit(as.character(x), '[_]'))[1]);
        }
        collapse.names = function(x){
          return(paste(x, sep='', collapse = '_'));
        }
        
        ggs = unlist(sapply(exons[, 4], extract_transcript)); exons$transcripts = ggs; exons$genes = transcripts[match(exons$transcripts, transcripts[,1]),3]
        ggs = unlist(sapply(introns[, 4], extract_transcript)); introns$transcripts = ggs; introns$genes = transcripts[match(introns$transcripts, transcripts[,1]),3]
        ggs = unlist(sapply(utr[, 4], extract_transcript)); utr$transcripts = ggs; utr$genes = transcripts[match(utr$transcripts, transcripts[,1]),3]
        
        
        aa = matrix(NA, nrow=nrow(keep), ncol = 17); aa = data.frame(aa);
        colnames(aa) = c('gene', 'index', 'best.model','half.life', 'cv.half.life', 'splicing.time', 'cv.splicing.time', 
                         'length.gene','median.length.transcript', 'nb.exon', 'median.length.exon',  'length.3UTR',
                         'nb.intron', 'median.length.intron', 'total.length.intron', 'length.last.intron', 'length.longest.intron');
        aa[, c(1:7)] = keep[, match(colnames(aa)[1:7], colnames(keep))];
        
        ## gene length
        xx = read.table('/Users/jiwang/RNA_seq_Data/usful_mapping_files/mm10_ens_genes_all.txt', sep='\t', header = FALSE, as.is = c(1, 4, 5))
        ggs = unlist(strsplit(as.character(xx[,4]), '[|]'))[2*c(1:nrow(xx))]; lls = abs(xx[,2] - xx[,3]);
        mm = match(aa$gene, ggs); aa$length.gene = lls[mm];
        
        for(n in 1:nrow(aa))
          #for(n in 1:1000)
        {
          #n = 3
          cat(n, '\n');
          gg = aa$gene[n];
          jj = which(transcripts[,3]==aa$gene[n]);
          if(length(jj)>0)
          {
            sign = transcripts[jj[1], 7];
            aa$median.length.transcript[n] = median(transcripts[jj, 6]);
            kk = which(exons$genes==gg);
            names.exons = apply(exons[kk, c(1:3)], 1, collapse.names);
            kk = kk[match(unique(names.exons), names.exons)]
            aa$nb.exon[n] = length(kk); 
            aa$median.length.exon[n]=median(abs(as.numeric(exons[kk, 3]) - as.numeric(exons[kk, 2])));
            ii = which(utr$genes==gg); aa$length.3UTR[n] = median(abs(utr[ii, 3]- utr[ii, 2]));
            mm = which(introns$genes ==gg);
            if(length(mm)>0){
              names.introns = apply(introns[mm, c(1:3)], 1, collapse.names);
              mm = mm[match(unique(names.introns), names.introns)]
              aa$nb.intron[n] = length(mm); aa$median.length.intron[n] = median(abs(introns[mm, 3] - introns[mm, 2]));
              aa$total.length.intron[n] = sum(abs(introns[mm, 3] - introns[mm, 2])); 
              aa$length.longest.intron[n] = max(abs(introns[mm, 3] - introns[mm, 2]));
              if(sign>0) {o1 = order(-introns[mm, 3]);
              }else {o1 = order(introns[mm, 2]);}
              aa$length.last.intron[n] = (abs(introns[mm, 3] - introns[mm, 2]))[1];
            }
          }
        }
        
        xx = read.csv('/Users/jiwang/Degradation_Liver/Database_2compare/Atger_translation_efficiency.csv', header = TRUE)
        effs = xx$mean_RFP_WT - xx$mean_Exon_WT;
        aa$translation.effs = effs[match(aa$gene, xx$Gene_Symbol)]
        aa$translation.effs = 2^aa$translation.effs;
        
        #save(aa, file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/half_life_splicing_time_genes_introns_exons_UTRs.Rdata',sep = ''))
        
        #### 
        #library(GenomicFeatures)
        #library(dplyr)
        #refSeq  <- makeTxDbFromUCSC(genom="mm10",tablename="refGene")   
        #aa = keep[c(sel2, sel4), ]
      }
      
      #####################
      ### explore correlations 
      #####################
      load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/half_life_splicing_time_genes_introns_exons_UTRs.Rdata',sep = ''))
      
      #boxplot(log10(aa$length.3UTR) ~ aa$best.model)
      cutoff = 0.4;
      sel2 = which(aa$best.model==2 & (aa$cv.half.life<cutoff|is.na(aa$cv.half.life))) #sel2 = which(keep$best.model==2 & keep$cv.half.life<cutoff & is.na(keep$non.identifiability.gamma))
      sel4 = which(aa$best.model==4 & (aa$cv.half.life<cutoff | is.na(aa$cv.half.life))) #sel4 = which(keep$best.model==4 & keep$cv.hal
      #sel2 = which(keep$best.model==2 & (keep$cv.half.life<cutoff|is.na(keep$cv.half.life))) #sel2 = which(keep$best.model==2 & keep$cv.half.life<cutoff & is.na(keep$non.identifiability.gamma))
      #sel4 = which(keep$best.model==4 & (keep$cv.half.life<cutoff | is.na(keep$cv.half.life))) #sel4 = which(keep$best.model==4 & keep$cv.half.life<cutoff & is.na(keep$non.identifiability.gamma))
      #sel2 = which(keep$best.model==2 & (keep$cv.half.life<cutoff)) #sel2 = which(keep$best.model==2 & keep$cv.half.life<cutoff & is.na(keep$non.identifiability.gamma))
      #sel4 = which(keep$best.model==4 & (keep$cv.half.life<cutoff))
      kk = c(sel2, sel4);
      
      #########
      #### Half-life 
      #########
      boxplot(log10(aa$length.3UTR) ~ aa$best.model)
      cor((aa$half.life[kk]), log(aa$length.gene[kk]), use = "na.or.complete")
      cor((aa$half.life[kk]), (aa$median.length.transcript[kk]), use = "na.or.complete")
      plot((aa$half.life[kk]), log(aa$median.length.transcript[kk]), cex=0.5)
      
      hist(log10(aa$median.length.transcript[kk]))
      hist(aa$nb.exon, breaks=100, xlim=c(0, 25))
      hist(aa$nb.intron, breaks=100, xlim=c(0, 25))
      
      cor((aa$half.life[kk]), (aa$nb.exon[kk]), use = "na.or.complete")
      cor((aa$half.life[kk]), aa$median.length.exon[kk], use = "na.or.complete")
      cor((aa$half.life[kk]), log10(aa$length.3UTR[kk]), use = "na.or.complete")
      plot(log10(aa$half.life[kk]), log10(aa$length.3UTR[kk]), cex=0.6)
      
      cor(log10(aa$length.3UTR[kk]), log(aa$median.length.transcript[kk]), use = "na.or.complete")
      cor(log10(aa$length.3UTR[kk]), log(aa$length.gene[kk]), use = "na.or.complete")
      cor(log(aa$length.gene[kk]),log(aa$median.length.transcript[kk]), use = "na.or.complete")
      
      hist(aa$half.life[kk], breaks = 50)
      plot(density(aa$half.life[kk]))
      
      groups = rep(NA, length(kk)); xx = aa$half.life[kk];
      cuts = c(2, 4, 6, 12);
      for(n in 1:(length(cuts)+1)){
        if(n==1) groups[which(xx<=cuts[n])] = 1;
        if(n<=length(cuts) & n>1) groups[which(xx<=cuts[n] & xx>cuts[(n-1)])] = n;
        if(n>length(cuts)) groups[which(xx>cuts[n-1])] = n;
      }
      
      pdfname = paste(folder, 'translation_efficiency_vs_half_life.pdf', sep='')
      pdf(pdfname, width=2., height=2.0)
      par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
      
      #boxplot(yy[,2] ~ yy[,1], main=' rel.ampl of premRNA', cex.main=0.7, col=cols[2:4], names=c('M2', 'M3', 'M4'), ylim=c(0, 1), lwd=0.5
      boxplot(log2(aa$translation.effs[kk]) ~ groups, ylim=c(-3.6, 4.5), axes=FALSE)
      
      axis(1, at=seq(1, 5, by=1), las=1, cex.axis = 1.0, labels = FALSE)
      axis(2, cex.axis=1.0)
      box()
      dev.off()
      
      #require('vioplot')
      #vioplot(log2(aa$translation.effs[kk[which(groups==1 & !is.na(aa$translation.effs[kk])==TRUE)]]),  
      #        log2(aa$translation.effs[kk[which(groups==2 & !is.na(aa$translation.effs[kk])==TRUE)]]), 
      #       log2(aa$translation.effs[kk[which(groups==3 & !is.na(aa$translation.effs[kk])==TRUE)]]), 
      #        log2(aa$translation.effs[kk[which(groups==4 & !is.na(aa$translation.effs[kk])==TRUE)]]),
      #        log2(aa$translation.effs[kk[which(groups==5 & !is.na(aa$translation.effs[kk])==TRUE)]]), 
      #        col = 'lightblue2',
      #        names=c(), lwd=0.5, colMed = 'white', pchMed = '.', horizontal=FALSE, wex=1, lty=2)
      #abline(h=-log10(0.05), col='darkgray', lwd=0.7)
      
      ### test for nb of exons
      boxplot(log10(aa$length.3UTR[kk]) ~ groups)
      
      test = c()
      cuts.nb.exons = c(5, 10, 15, 20, 30)
      yy = aa$nb.exon[kk];
      for(n in 1:(length(cuts.nb.exons)+1))
      {
        #nb.total = length(groups)
        if(n==1) controls = (which(yy<cuts.nb.exons[n]))
        if(n<=length(cuts.nb.exons) & n>1) controls = (which(yy<cuts.nb.exons[n] & yy>=cuts.nb.exons[(n-1)]));
        if(n>length(cuts)) controls = (which(yy>=cuts.nb.exons[(n-1)]))
        
        enrich = c();
        for(n.g in unique(groups))
        {
          # n.g = 1
          finds = (which(groups==n.g))
          positives = intersect(finds, controls)
          enrich = c(enrich, phyper((length(positives)-1), length(controls),  length(groups)-length(controls), length(finds), lower.tail = FALSE));
        }
        test = rbind(test, enrich)
      }
      
      
      #plot((aa$half.life[kk]), aa$translation.effs[kk], log='y')
      pdfname = paste(folder, 'degradation_rate_vs_translation_efficiency.pdf', sep='')
      pdf(pdfname, width=2.4, height=1.6)
      par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
      
      cor(log(1/(aa$half.life[kk])), log(aa$translation.effs[kk]), use = "na.or.complete")
      #lims = c(-4, 5);
      xx = as.numeric(log2(1/(aa$half.life[kk])));yy = log2(aa$translation.effs[kk]);
      jj = which(!is.na(xx)==TRUE & !is.na(yy)==TRUE);
      xx = xx[jj]; yy= yy[jj];
      jj = index.outliers(yy);xx = xx[-jj];yy = yy[-jj];
      #quantiles = quantile(xx, probs = seq(0, 1, by=0.25), na.rm = TRUE);quantiles = quantiles[2:5];
      #groups = rep(NA, length(xx));
      #groups[which(xx<=quantiles[1])] = 1;
      #groups[which(xx<=quantiles[2] & xx>quantiles[1])] = 2;
      #groups[which(xx<=quantiles[3] & xx>quantiles[2])] = 3;
      #groups[which(xx<=quantiles[4] & xx>quantiles[3])] = 4;
      
      #boxplot(yy ~ groups, ylim=c(-3, 4))
      plot(xx,yy, cex=0.1, axes=FALSE, xlim=c(-3, 4))
      fit = lm(log2(aa$translation.effs[kk]) ~ log2(1/(aa$half.life[kk])))
      require('MASS');
      fit2 = rlm(xx[which(!is.na(xx)==TRUE & !is.na(yy)==TRUE)], yy[which(!is.na(xx)==TRUE & !is.na(yy)==TRUE)])
      abline(fit$coefficients[1], fit$coefficients[2], lwd=1.0, col='red')
      axis(1, at=seq(-4, 2, by=1), las=1,cex.axis = 1.0)
      axis(2,at=seq(-4, 5, by=1), cex.axis=1.0)
      box()
      
      dev.off()
      
      #cor(1/(aa$half.life[kk]), log(aa$length.3UTR[kk]), use = "na.or.complete")
      
      sel2 = which(keep$best.model==2 & (keep$cv.splicing.time<cutoff)) #sel2 = which(keep$best.model==2 & keep$cv.half.life<cutoff & is.na(keep$non.identifiability.gamma))
      sel4 = which(keep$best.model==4 & (keep$cv.splicing.time<cutoff))
      kk = c(sel2, sel4);
      cor((aa$splicing.time[kk]), (aa$nb.intron[kk]), use = "na.or.complete")
      plot((aa$nb.intron[kk]), (aa$splicing.time[kk]), log='x')
      
      cor((aa$splicing.time[kk]), log(aa$median.length.intron[kk]), use = "na.or.complete")
      plot(log(aa$median.length.intron[kk]), (aa$splicing.time[kk]), ylim=c(0, 20))
      
      cor((aa$splicing.time[kk]), log(aa$total.length.intron[kk]), use = "na.or.complete")
      plot(log(aa$total.length.intron[kk]), (aa$splicing.time[kk]), ylim=c(0, 20))
      
      cor((aa$splicing.time[kk]), log(aa$length.last.intron[kk]), use = "na.or.complete")
      plot(log(aa$length.last.intron[kk]), (aa$splicing.time[kk]), cex=0.5, ylim = c(0, 20))
      
      cor((aa$splicing.time[kk]), log(aa$length.longest.intron[kk]), use = "na.or.complete")
      plot(log(aa$length.longest.intron[kk]), (aa$splicing.time[kk]), cex=0.6, ylim=c(0, 20))
      
      cor((aa$splicing.time[kk]), log(aa$median.length.transcript[kk]), use = "na.or.complete")
      
    }
  }
  
  
  Half.life_function.analysis.Cedric = FALSE
  if(Half.life_function.analysis.Cedric)
  {
    ### convert gene symbols to ensemble:
    R = read.table('/Users/jiwang/Degradation_Liver/DATA_RNA_Seq/Jingkui_count_IE_Total_040216.txt',sep='\t', header=TRUE)
    #mapping = R[,1]
    xx = unlist(strsplit(as.character(R[,1]),"[|]"))
    mapping = cbind(xx[c(1:nrow(R)*2)], xx[(c(1:nrow(R)*2)-1)])
    #R[,1] = xx
    bb = aa;
    cutoff = 0.4;
    sel2 = which(aa$best.model==2 & (aa$cv.half.life<cutoff|is.na(aa$cv.half.life))) #sel2 = which(keep$best.model==2 & keep$cv.half.life<cutoff & is.na(keep$non.identifiability.gamma))
    sel4 = which(aa$best.model==4 & (aa$cv.half.life<cutoff | is.na(aa$cv.half.life))) #sel4 = which(keep$best.model==4 & keep$cv.hal
    bb = bb[c(sel2, sel4), ];
    
    groups = rep(NA, length(kk)); xx = aa$half.life[c(sel2, sel4)];
    cuts = c(2, 4, 6, 12);
    for(n in 1:(length(cuts)+1)){
      if(n==1) groups[which(xx<=cuts[n])] = 1;
      if(n<=length(cuts) & n>1) groups[which(xx<=cuts[n] & xx>cuts[(n-1)])] = n;
      if(n>length(cuts)) groups[which(xx>cuts[n-1])] = n;
    }
    
    bb$groups = groups
    mm = match(bb$gene, mapping[,1])
    bb$Ensmus = mapping[mm, 2]
    
    ### I'm working with ensembl genes names as rownames of the data frame 
    library("RDAVIDWebService")
    A=list()
    
    david<-DAVIDWebService(email="cedric.gobet@epfl.ch", 
                           url="http://david-d.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
    setTimeOut(david,1000000) ### Increase time out, you can load more genes like that
    #getAllAnnotationCategoryNames(david)
    setAnnotationCategories(david, c("GOTERM_BP_ALL","KEGG_PATHWAY")) # I defined GO BP ALL and KEGG pathway but you can define more
    # You define your background, list of expressed genes, here RPKM.M.model
    bcknd1<-addList(david, bb$Ensmus,idType="ENSEMBL_GENE_ID", listName="Bcknd", listType="Background") 
    
    # GO Enrichment in each model and stock in list A
    for(i in unique(bb$groups)){
      nam = bb$Ensmus[which(bb$groups==i)]
      if(length(nam)!=0){
        addList(david,nam,idType="ENSEMBL_GENE_ID", listName= paste("model",i,sep="_"), listType="Gene")
        setCurrentBackgroundPosition(david,1)
        res=getFunctionalAnnotationChart(david)
        summar=data.frame(pval=as.numeric(res$PValue),qval=as.numeric(res$Benjamini),term=res$Term,Enr=as.numeric(res$Fold.Enrichment),
                          count=as.numeric(res$Count),Pop=as.numeric(res$Pop.Hits),Genes=res$Genes,stringsAsFactors=FALSE)
        A[[i]]=summar
      }
    }
    ########## Parsing of the results , threshold on different value and write in .xls file with different sheets
    library(xlsx)
    pv=0.05; enr=1.5; count=2;  pop=1000;  v=1;
    for(i in unique(bb$groups)){
      summar=A[[i]]
      signi = which(summar$pval < pv & summar$Enr > enr & summar$count > count & summar$Pop < pop)
      summar_signi = summar[signi,]
      summar_signi$term[grep('\\~',summar_signi$term)] = as.character(unlist(sapply(strsplit(as.character(summar_signi$term[grep('\\~',summar_signi$term)]),"~"),"[[",2)))
      summar_signi$term[grep('\\:',summar_signi$term)] = as.character(unlist(sapply(strsplit(as.character(summar_signi$term[grep('\\:',summar_signi$term)]),":"),"[[",2)))
      
      ### write down the results in an Excel table
      if(length(signi)!=0)
      { 
        if(v==1){
          write.xlsx(summar_signi, file="/Users/jiwang/Degradation_Liver/Main_Code_Degradation/Tables/function_analysis_half_life.xlsx",sheetName=paste0("de_model_",i))
        }else{
          write.xlsx(summar_signi, file="/Users/jiwang/Degradation_Liver/Main_Code_Degradation/Tables/function_analysis_half_life.xlsx", sheetName=paste0("de_model_",i), append=TRUE)
        }
        v=v+1
      }
      
      ### make a plot
      Plot.enrichment = FALSE
      if(Plot.enrichment)
      {
        pdfname = paste(folder, 'Function_Enrichment_half_life_group_', i, '.pdf', sep='')
        pdf(pdfname, width=2.5, height=2.5)
        par(cex = 0.7, las = 0, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        palette(c('gray','green3','tomato','black'))
        summar_signi = summar_signi[order(summar_signi$Enr), ]
        counts <- summar_signi$Enr;
        barplot(counts, main=NA, horiz=TRUE, col='gray', space = 1, xlim = c(0, 10))
        for(kk in c(1:nrow(summar_signi)))
        {
          text((counts[length(counts)-kk+1]+2), ((nrow(summar_signi)-kk)*2+1.5), summar_signi$term[kk], cex=0.7, col=i)
        }
        
        dev.off()
      }
    }
    
  }
  
}

FIGURE_4 = FALSE
if(FIGURE_4)
{
  res.version.old = '_total_counts_all_genes_norm_params_v6_OLD';
  load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', 
                  res.version.old, '.Rdata',sep = ''))
  
  #load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_RNA_seq_total_counts_table_sx', res.version, '.Rdata',sep = ''))
  res.version = '_total_counts_all_genes_norm_params_final_v7_All_summary';
  load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
  
  T$BIC.best.model = T$best.model;
  source('functions.R')
  attribute.global.variable.colors();
  folder = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/'

  #########
  #### Phases and Fc of Rhythmic degradation (M3 and M4)
  #########
  #cutoff = 1
  Phases.FCs.RD.M3.M4 = FALSE
  if(Phases.FCs.RD.M3.M4)
  {
    #### select genes in M3 and M4
    #mm = match(examples, keep$gene)
    #cbind(examples, keep$sem.phase.gamma[mm], keep$cv.eps.gamma[mm])
    sel3 = which(keep$best.model==3);sel4 = which(keep$best.model==4)
    cutoff = 0.4;
    sel3 = which(keep$best.model==3 & (keep$sem.phase.gamma<1 | is.na(keep$sem.phase.gamma)) & (keep$cv.eps.gamma<cutoff | is.na(keep$cv.eps.gamma)==TRUE)) #sel3 = which(keep$best.model==3 & keep$cv.half.life<cutoff & is.na(keep$non.identifiability.gamma))
    sel4 = which(keep$best.model==4 & (keep$sem.phase.gamma<1 | is.na(keep$sem.phase.gamma)) & (keep$cv.eps.gamma<cutoff | is.na(keep$cv.eps.gamma)==TRUE)) #se
    
    #### Phase distribution of RD in M3
    phases3 = keep$phase.gamma[sel3];
    phases4 = keep$phase.gamma[sel4];
    #phases = T$phase.gamma.m3[which(T$BIC.best.model==3)]
    pdfname = paste(folder, 'Rhythmic_degradation_phases_M3.pdf', sep='')
    pdf(pdfname, width=2.0, height=1.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.4,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    breaks= seq(0, 24, by=2)
    hist(phases3, breaks = breaks, col=3+4, axes=FALSE, main=NA, xlab=NA, ylab=NA, ylim=c(0, 20));
    #hist(phases3, breaks = breaks, col=3+4, add=TRUE, axes = FALSE);
    axis(2,  at=seq(0, 20, by=5), las = 3, cex.axis=1.0)
    axis(1, at=seq(0, 24, by=4), las=1,cex.axis = 1.0)
    dev.off()
    
    pdfname = paste(folder, 'Rhythmic_degradation_phases_Model_3_subgroups.pdf', sep='')
    pdf(pdfname, width=2., height=2.)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    circular_phase24H_histogram(phases3, col=3+4, breaks=seq(0, 24, by=2), cex.axis=0.5, cex.lab=0.01, lwd=0.5)
    dev.off()
    
    #### Phase distribution of RD in M4
    pdfname = paste(folder, 'Rhythmic_degradation_phases_M4.pdf', sep='')
    pdf(pdfname, width=2.0, height=1.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.4,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    breaks= seq(0, 24, by=2)
    hist(phases4, breaks = breaks, col=4+4, axes=FALSE, main=NA, xlab=NA, ylab=NA, ylim=c(0, 150));
    #hist(phases3, breaks = breaks, col=3+4, add=TRUE, axes = FALSE);
    axis(2,  at=seq(0, 150, by=50), las = 3, cex.axis=1.0)
    axis(1, at=seq(0, 24, by=4), las=1,cex.axis = 1.0)
    dev.off()
    
    pdfname = paste(folder, 'Phases_rhythmic_degradation_Model_4.pdf', sep='')
    pdf(pdfname, width=2., height=2.)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    circular_phase24H_histogram(phases4, col=4+4, breaks=seq(0, 24, by=2),cex.axis=0.5, cex.lab=0.01, lwd=0.5)
    dev.off()
    
    #### FCs of RD in M3 and M4
    #jj = which(keep$best.model==3 & (keep$cv.eps.gamma<cutoff | is.na(keep$cv.eps.gamma)==TRUE))
    fc3 = log2((1+keep$eps.gamma[sel3])/(1-keep$eps.gamma[sel3]));
    #jj = which(keep$best.model==4 & (keep$cv.eps.gamma<cutoff | is.na(keep$cv.eps.gamma)==TRUE))
    fc4 = log2((1+keep$eps.gamma[sel4])/(1-keep$eps.gamma[sel4]));
    
    pdfname = paste(folder, 'FC_rhythmic_degradation_M3.pdf', sep='')
    pdf(pdfname, width=2.0, height=1.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.4,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    breaks = seq(0, 4.5, by=0.5)
    hist(fc3, breaks = breaks, col=3+4, axes=FALSE, main=NA, xlab=NA, ylab=NA, cex.lab=1.0, ylim=c(0, 80))
    #hist(fc3, breaks = breaks, col=3+4, axes=FALSE, add = TRUE)
    axis(2, at=seq(0, 80, by=20), las=3, cex.axis=1.0)
    axis(1, at=seq(0, 4, by=1), las=1,cex.axis = 1.0)
    #box()
    dev.off()
    pdfname = paste(folder, 'FC_rhythmic_degradation_M4.pdf', sep='')
    pdf(pdfname, width=2.0, height=1.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.4,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    breaks = seq(0, 4.5, by=0.5)
    hist(fc4, breaks = breaks, col=4+4, axes=FALSE, main=NA, xlab=NA, ylab=NA, cex.lab=1.0)
    #hist(fc3, breaks = breaks, col=3+4, axes=FALSE, add = TRUE)
    axis(2, at=seq(0, 300, by=100), las=3, cex.axis=1.0)
    axis(1, at=seq(0, 4, by=1), las=1,cex.axis = 1.0)
    #box()
    dev.off()
    
    
    Check.large.FC.M3.M4 = FALSE
    if(Check.large.FC.M3.M4)
    {
     
      #T$BIC.best.model = T$best.model;
      #load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', 
      #                res.version.old, '.Rdata',sep = ''))
      
      test.sel3 = sel3[which(fc3>3)]
      test.sel4 = sel4[which(fc4>3)]
      
      source('functions.R')
      folder = "/Users/jiwang/Dropbox/mRNA_decay/Analysis/M3_M4_largeFC/"
      index = c(match(keep$gene[c(test.sel3, test.sel4)], T.old$gene))
      index = index[which(!is.na(index)==TRUE)]
      index = index[which(T.old$BIC.best.model[index]>2)]
      
      plot.genes.examples(index=index, T=T.old, Figure=TRUE, folder=folder, summary.plot=FALSE);
     
      folder = "/Users/jiwang/Dropbox/mRNA_decay/Analysis/M3_M4_largeFC/"
      pdfname = paste(folder, 'Check_highFC_M3_M4.pdf', sep='')
      pdf(pdfname, width=10., height=8.)
      
      test.index = c(test.sel3, test.sel4)
      plot(keep$eps.gamma[test.index], keep$cv.eps.gamma[test.index])
      plot(keep$eps.gamma[test.index], keep$half.life[test.index],  ylim=c(0, 20))
      plot(keep$eps.gamma[test.index], keep$cv.half.life[test.index],  ylim=c(0, 2)) 
      plot(-log10(keep$pval.log2rpkm.mRNA[test.index]), -log10(keep$pval.log2rpkm.premRNA[test.index]),  ylim=c(0, 20)) 
      plot(keep$amp.log2rpkm.mRNA[test.index], keep$amp.log2rpkm.premRNA[test.index])
      
      dev.off()
      
    }
    
  }
  
  ###################################################f
  ###### Show examples of predicited mRBP acitivities compared with protein or phophoproteins
  ###################################################
  mRBP.motif.activity.vs.protein.abundance = FALSE
  if(mRBP.motif.activity.vs.protein.abundance)
  {
    total = read.table('/Users/jiwang/Dropbox/mRNA_decay/Analysis/Analysis_mRBP_motifs/Table_total_proteins_PNAS.txt', sep='\t', header = TRUE)
    phospho = read.delim('/Users/jiwang/Dropbox/mRNA_decay/Analysis/Analysis_mRBP_motifs/phospho_Robles_2016.txt', header = TRUE)
    robles = read.table('/Users/jiwang/Proteomics_analysis/Nuclear_proteins/Tables_DATA/Table_total_proteins_Robles.txt',sep='\t', header=TRUE)
    
    infer = read.table('/Users/jiwang/Dropbox/mRNA_decay/fromJake/Code/good_results_Jingkui_December_2016.runbyjake.txt')
    load('/Users/jiwang/Dropbox/mRNA_decay/Analysis_mRBP_motifs/RBPs_elstaic_net/mRBPs_positive_controls_rhythmic_proteins.Rdata')
    
    standardization.nona = function(x)
    { xx = (x-mean(x[which(!is.na(x)==TRUE)]))/sd(x[which(!is.na(x)==TRUE)])
    return(xx);}
    
    examples = c('Matr3', 'Syncrip', 'Fxr1', 'Ncl')
    attribute.global.variable.colors();
    
    for(motif in examples)
    {
      motif = examples[3] 
      
      gene = motif
      print(motif)
      xx = infer[grep(motif, rownames(infer)), ]
      pp = as.numeric(xx[c(1,2)])
      amp = as.numeric(xx[c(3:4)])
      kk = which(!is.na(pp)==TRUE)
      if(kk==2) col=4
      if(kk==1) col=3;
      pp = pp[kk];amp=amp[kk];
      
      print(total[grep(gene, total$Gene.names.total), c(grep('phase', colnames(total)), grep('pval', colnames(total)))])
      print(robles[grep(gene, robles$Gene.name.robles), c(grep('phase', colnames(robles)), grep('pval', colnames(robles)))])
      print(phospho[grep(gene, phospho$Gene.names), c(grep('phase', colnames(phospho)), grep('pval', colnames(phospho)))])
      
      
      pdf.name = paste(folder, "Inferred_mRBPs_vs_total_protein_", motif,".pdf", sep = "")
      
      pdf(pdf.name, width=1.8, height=1.5)
      par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
      
      #mm = which(total$Gene.names.total==motif)
      mm = which(robles$Gene.name.robles==motif)
      y0 = standardization.nona(as.numeric(robles[mm, c(2:17)]))
      time = seq(0, 48, by=0.5)
      y1 = amp*cos(2*pi/24*(time-pp))
      y1 = standardization.nona(y1)
      lims = range(c(y0, y1), na.rm=TRUE)
      
      plot(c(0:15)*3, y0, type='l', col="darkblue",lwd=1., main=toupper(gene), cex.main=0.6, ylim=lims, ylab=NA, xlab=NA, xlim=c(0,48), axes=FALSE)
      points(c(0:15)*3, y0, type='p', pch=16, col='darkblue', cex=0.8)
      points(time, y1, type='l', lwd=1.2, col=col)
      
      if(gene=='Fxr1')  axis(1,at=seq(0, 48, by=12),cex.axis =1.0)
      #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
      lims = signif(lims, d=1)
      if(gene=='Hdac3') {
        axis(2,at = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2),las=1,cex.axis = 1.0)
      }else{
        by = signif((lims[2]-lims[1])/5,d=0)
        axis(2,at = seq(-2, 2, by=1),las=1,cex.axis = 1.0)
      }
      box(); abline(h=0,lty=2,lwd=1., col="darkgray");
      dev.off()
    }
  }
  
  
  ######################
  #### Compare estimated phases and amplitudes of RD with PA-test from Sara
  ######################
  Compare.phase.amplitude.with.PA.test = FALSE
  if(Compare.phase.amplitude.with.PA.test)
  {
    cutoff = 0.4;
    sel3 = which(keep$best.model==3 & (keep$sem.phase.gamma<1 | is.na(keep$sem.phase.gamma)) & (keep$cv.eps.gamma<cutoff | is.na(keep$cv.eps.gamma)==TRUE)) #sel3 = which(keep$best.model==3 & keep$cv.half.life<cutoff & is.na(keep$non.identifiability.gamma))
    sel4 = which(keep$best.model==4 & (keep$sem.phase.gamma<1 | is.na(keep$sem.phase.gamma)) & (keep$cv.eps.gamma<cutoff | is.na(keep$cv.eps.gamma)==TRUE)) #se
    sels = c(sel3, sel4)
    
    com2 = read.table('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/Tables/PAtest_result_withEstimated_HL.txt', header = TRUE, sep='\t', as.is = c(1))
    mm = match(keep$gene[sels], com2$gene); 
    jj = which(!is.na(mm)==TRUE & com2$PA.test.p.val[mm]<0.1)
    kk = sels[jj]; mm = mm[jj];
    
    pdfname = paste(folder, 'Sarah_Comparison_phases_v5.pdf', sep='')
    pdf(pdfname, width=1.8, height=1.8)
    par(cex = 0.7, las = 0, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3, pty='s')
    pp1 = keep$phase.gamma[kk]; pp2 = com2$phi.deg[mm];
    lims = range(c(pp1, pp2))
    plot(pp1, pp2, type='p', cex=0.7, pch=1, col='darkblue', xlab=NA, ylab=NA, main=NA, axes=FALSE);
    abline(0, 1, col='darkgray', lwd=1.5);
    axis(1,at=seq(0, 24, by=6),cex.axis =1.0)
    axis(2, at=seq(0, 24, by=6), las=1,cex.axis = 1.0)
    box()
    
    dev.off()
    
    pdfname = paste(folder, 'Sarah_Comparison_amplitudes_v5.pdf', sep='')
    pdf(pdfname, width=1.8, height=1.8)
    par(cex = 0.7, las = 0, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3, pty='s')
    eps1 = keep$eps.gamma[kk]; eps2 = com2$A.deg[mm];
    lims=range(c(eps1, eps2))
    plot(eps1, eps2, type='p', cex=0.7, pch=1, col='darkblue', xlab=NA, ylab=NA, main=NA, axes=FALSE, xlim = lims, ylim = lims);
    abline(0, 1, col='darkgray', lwd=1.0);
    #abline(24, 1, col='blue', lwd=1.0);
    #abline(-24, 1, col='blue', lwd=1.0);
    axis(1,at=seq(0, 1.2, by=0.5),cex.axis =1.0)
    axis(2, at=seq(0, 1.2, by=0.5), las=1,cex.axis = 1.0)
    box()
    
    dev.off()
   
    #com2 = read.table(file='/Users/jiwang/Degradation_Liver/Database_2compare/polyA_length_dynamics.txt', header = TRUE, as.is = c(1, 2))
    #cutoff = 1;jj = which(keep$sem.phase.gamma<cutoff & (is.na(keep$non.identifiability.gamma)==TRUE | keep$non.identifiability.gamma==1))
    #mm = match(keep$gene[jj], com2$gene)
    #kk = which(com2$pval[mm]<0.01)
    #pp1 = keep$phase.gamma[jj]
    #pp2 = com2$phase[mm]
    #plot(pp1, pp2);abline(0, 1, col='blue', lwd=1.0);abline(6, 1, col='blue', lwd=1.0);abline(-18, 1, col='blue', lwd=1.0);
  }
  
  ######################
  #### Phase delays between mRNAs and RD for M3 and M4 (NOT used now)
  ######################
  Phase.delays.mRNA.RD = FALSE
  if(Phase.delays.mRNA.RD)
  {
    sel3 = which(keep$best.model==3 & (keep$sem.phase.gamma<1 | is.na(keep$sem.phase.gamma)) & (keep$cv.eps.gamma<cutoff | is.na(keep$cv.eps.gamma)==TRUE)) #sel3 = which(keep$best.model==3 & keep$cv.half.life<cutoff & is.na(keep$non.identifiability.gamma))
    sel4 = which(keep$best.model==4 & (keep$sem.phase.gamma<1 | is.na(keep$sem.phase.gamma)) & (keep$cv.eps.gamma<cutoff | is.na(keep$cv.eps.gamma)==TRUE)) #se
    #sel3 = which(keep$best.model==3)
    #sel4 = which(keep$best.model==4)
  
    #### Phase distribution of RD in M3
    #phase.delay.m3 = (T$phase.rpkm.mRNA[keep$index[sel3]] - keep$phase.gamma[sel3])%%24;  
    #phase.delay.m4 = (T$phase.rpkm.mRNA[keep$index[sel4]]  - keep$phase.gamma[sel4])%%24;  
    phase.delay.m3 = (keep$phase.mRNA[sel3] - keep$phase.gamma[sel3])%%24;  
    phase.delay.m4 = (keep$phase.mRNA[sel4] - keep$phase.gamma[sel4])%%24;  
    breaks= seq(0, 24, by=1)
    
    #phases = T$phase.gamma.m3[which(T$BIC.best.model==3)]
    pdfname = paste(folder, 'Phase_Delay_mRNA_vs_RD_M3.pdf', sep='')
    pdf(pdfname, width=2.0, height=1.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.4,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    hist(phase.delay.m3, breaks = breaks, col=3+4, axes=FALSE, main=NA, xlab=NA, ylab=NA, ylim=c(0, 20));
    #hist(phases3, breaks = breaks, col=3+4, add=TRUE, axes = FALSE);
    axis(2,  at=seq(0, 40, by=10), las = 3, cex.axis=1.0)
    axis(1, at=seq(0, 24, by=4), las=1,cex.axis = 1.0)
    box()
    dev.off()
  
    #### Phase distribution of RD in M4
    pdfname = paste(folder, 'Phase_Delay_mRNA_vs_RD_M4.pdf', sep='')
    pdf(pdfname, width=2.0, height=1.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.4,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    hist(phase.delay.m4, breaks = breaks, col=4+4, axes=FALSE, main=NA, xlab=NA, ylab=NA, ylim=c(0, 150));
    #hist(phases3, breaks = breaks, col=3+4, add=TRUE, axes = FALSE);
    axis(2,  at=seq(0, 100, by=50), las = 3, cex.axis=1.0)
    axis(1, at=seq(0, 24, by=4), las=1,cex.axis = 1.0)
    box()
    dev.off()
    
  }
  
}

Figure_5 = FALSE
if(Figure_5)
{
  #######################################
  #######################################
  #### Explore advantages of using rhythmic mRNA degradations
  #######################################
  #######################################
  
  ####################
  ##### Prepare the table for this figure
  ####################
  Table.preparation = FALSE
  if(Table.preparation)
  {
    res.version.v6 = '_total_counts_all_genes_norm_params_v6';
    load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_RNA_seq_total_counts_table_sx', res.version.v6, '.Rdata',sep = ''))
    attribute.global.variable.colors();
    folder = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/'
    
    #### Prepare tables and M4 in order to explore advantages of rhythmic degradation 
    sel4 = which(keep$best.model==4 & (keep$sem.phase.gamma<1 | is.na(keep$sem.phase.gamma)) & (keep$cv.eps.gamma<0.4 | is.na(keep$cv.eps.gamma)==TRUE)) #se
    #sel4 = which(keep$best.model==4 & (keep$sem.phase.gamma<1 | is.na(keep$sem.phase.gamma)) & (keep$cv.eps.gamma<0.4 | is.na(keep$cv.eps.gamma)==TRUE) 
    #            & T$qv.rpkm.mRNA[keep$index]<0.05 & T$qv.rpkm.premRNA[keep$inde]<0.05) 
    #sel4 = which(keep$best.model==4 & (keep$sem.phase.gamma<1) & (keep$cv.eps.gamma<0.4)) #se
    #sel4 = which(keep$best.model==4 & (keep$sem.phase.gamma<1.2 | is.na(keep$sem.phase.gamma)) & (keep$cv.eps.gamma<0.5 | is.na(keep$cv.eps.gamma)==TRUE)) #se
    aa = matrix(NA, nrow = length(sel4), ncol=18);  aa = data.frame(aa);
    colnames(aa) = c('gene', 'index', 'half.life', 'cv.half.life', 'eps.gamma', 'cv.eps.gamma', 'phase.gamma', 'sem.phase.gamma',
                     'relamp.premRNA', 'phase.premRNA', 'qv.premRNA', 'relamp.mRNA', 'phase.mRNA',  'qv.mRNA', 'relamp.m22', 'phase.m22',
                     'delay.RD.RS', 'delay.mRNA.RD');
    aa[, c(1:16)] = keep[sel4, match(colnames(aa)[1:16], colnames(keep))];
    #index = aa$index;
    #aa$qv.premRNA = T$qv.rpkm.premRNA[index]; aa$qv.mRNA = T$qv.rpkm.mRNA[index];
    #aa$relamp.premRNA = T$rel.amp.rpkm.premRNA[index];aa$relamp.mRNA = T$rel.amp.rpkm.mRNA[index];
    #aa$phase.premRNA = T$phase.rpkm.premRNA[index];aa$phase.mRNA = T$phase.rpkm.mRNA[index];
    aa$delay.mRNA.RD = (aa$phase.mRNA - aa$phase.gamma)%%24;
    aa$eps.ratio.mRNA.RD = aa$relamp.mRNA/aa$eps.gamma;
    aa$delay.RD.RS = (aa$phase.gamma - aa$phase.premRNA)%%24;
    aa$eps.ratio.RD.RS = aa$eps.gamma/aa$relamp.premRNA;
    aa$delay.m22.RS = (aa$phase.m22 - aa$phase.premRNA)%%24 
    dd = (aa$phase.mRNA - aa$phase.m22)%%24; kk = which(dd>=12); dd[kk] = dd[kk]-24;
    aa$delay.mRNA.m22 = dd;
    aa$delay.mRNA.RS = (aa$phase.mRNA - aa$phase.premRNA)%%24;
    aa$eps.ratio.m22.RS = aa$relamp.m22/aa$relamp.premRNA;
    aa$eps.ratio.mRNA.m22 = aa$relamp.mRNA/aa$relamp.m22;
    aa$eps.ratio.mRNA.RS = aa$relamp.mRNA/aa$relamp.premRNA;
    #gammas = log(2)/aa$half.life;
    #rr2 = aa$relamp.premRNA*gammas/sqrt(w^2+gammas^2);
    #kk = which(aa$cv.half.life>0.4); rr2[kk] = NA; 
    #aa$relamp.mRNA.m2 = rr2
    #qq2 = (atan2(2*pi/24, gammas)/w + aa$phase.premRNA)%%24 
    #kk = which(aa$cv.half.life>0.4); qq2[kk] = NA; 
    #aa$phase.mRNA.m2 = qq2
    ### filter again the genes in M4
    kk = which(aa$qv.mRNA<0.05 & aa$qv.premRNA<0.05 & aa$eps.gamma>0. & aa$relamp.mRNA>0.1 & aa$relamp.premRNA>0.1)
    aa = aa[kk, ]
    
    ### select genes with half-life > 10min
    jj = which(aa$half.life>(15/60))
    aa = aa[jj, ]
  }
   
  Amplifier_Advancer_Delayer.RD = FALSE
  if(Amplifier_Advancer_Delayer.RD)
  {
    kk_1 = which(aa$half.life>=5); #plot(aa$half.life[kk1], aa$delay.RD.RS[kk1], col='darkgreen', cex=cex)
    kk_2 = which(aa$half.life<5 & aa$delay.RD.RS<12); #points(aa$half.life[kk], aa$delay.RD.RS[kk], col='darkorange', cex=cex)
    kk_3 = which(aa$half.life<5 & aa$delay.RD.RS>=12); #points(aa$half.life[kk], aa$delay.RD.RS[kk], col='darkred', cex=cex)
    cols =c('darkgreen','darkorange', 'darkred')
    
    #######
    ### delay and amplitude ratio between RD and RS
    #######
    pdfname = paste(folder, 'Delays_RD_premRNA_Model_4.pdf', sep='')
    pdf(pdfname, width=2.0, height=1.8)
    #par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(0.2,3,0.2,0.4)+0.1, tcl = -0.3)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    breaks = seq(0, 24, by=2)
    dd = aa$delay.RD.RS
    circular_phase24H_histogram(dd, col=4+4, breaks=seq(0, 24, by=2), cex.axis=0.5, cex.lab=0.01, lwd=0.5)
    
    #hist(dd, breaks = breaks, col=4+4, axes=FALSE, main=NA, xlim=c(0, 24), cex.main=1.0, xlab=NA, ylab=NA, cex.lab=1.0, ylim = c(0, 30))
    #axis(2, at=seq(0, 30, by=10), cex.axis=0.8)
    #axis(1, at=seq(0, 24, by=6), cex.axis = 0.8)
    #abline(v=c(median(dd[kk_1]), mean(dd[kk_2]), mean(dd[kk_3])), col=cols, lwd=1.5)
    #box()
    dev.off()
  
    pdfname = paste(folder, 'Relamp_raito_RD_premRNA_Model_4.pdf', sep='')
    pdf(pdfname, width=2.0, height=1.8)
    #par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(0.2,3,0.2,0.4)+0.1, tcl = -0.3)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    #breaks = seq(0, 24, by=2)
    ratios = (aa$eps.gamma/aa$relamp.premRNA)
    ratios = aa$eps.ratio.RD.RS 
    hist(ratios, breaks = 12, col=4+4, axes=FALSE, main=NA, cex.main=1.0, xlab=NA, ylab=NA, cex.lab=1.0)
    #lims = c(0, 1.2)#lims = range(aa$eps.gamma, aa$relamp.premRNA)
    #plot(aa$relamp.premRNA, aa$eps.gamma, col=4+4, xlim=lims, ylim= lims)
    axis(2, at=seq(0, 60, by=20), cex.axis=0.8)
    axis(1, at=seq(0, 6, by=1), cex.axis = 0.8)
    #abline(v=1, col='blue', lty=2, lwd=1.5)
    #box()
    dev.off()
    
    
    #######
    ### delay and amplitude ratios in function of half-life
    #######
    cex = 0.5
    pdfname = paste(folder, 'Phase_delay_RD_RS_M4_half_life.pdf', sep='')
    pdf(pdfname, width=2., height=1.8)
    #par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(0.2,3,0.2,0.4)+0.1, tcl = -0.3)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    plot(aa$half.life, aa$delay.RD.RS, log='x', type='n', xlab=NA, ylab=NA, main=NA, axes = FALSE);
    #attach(mtcars)
    #par(mfrow=c(1,2)) 
    #plot(aa$half.life, aa$delay.RD.RS, log='x');abline(v=5, col='red');abline(h=12, col='green')
    abline(v=5,col='darkgray'); 
    abline(h=12, col='darkgray')
    kk = which(aa$half.life>=5); points(aa$half.life[kk], aa$delay.RD.RS[kk], col='darkgreen', cex=cex)
    kk = which(aa$half.life<5 & aa$delay.RD.RS<12); points(aa$half.life[kk], aa$delay.RD.RS[kk], col='darkorange', cex=cex)
    kk = which(aa$half.life<5 & aa$delay.RD.RS>=12); points(aa$half.life[kk], aa$delay.RD.RS[kk], col='darkred', cex=cex)
    axis(2, at=seq(0, 24, by=6), cex.axis=0.8)
    axis(1, at=c(0.5, 1, 2, 5, 20), cex.axis = 0.8)
    #abline(v=1, col='blue', lty=2, lwd=1.5)
    box()
    
    dev.off()
    
    pdfname = paste(folder, 'Eps_raito_RD_RS_M4_half_life.pdf', sep='')
    pdf(pdfname, width=2.0, height=1.8)
    #par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(0.2,3,0.2,0.4)+0.1, tcl = -0.3)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    plot(aa$half.life, aa$eps.ratio.RD.RS, log='x', type='n', xlab=NA, ylab=NA, main=NA, axes = FALSE);
    abline(v=5,col='darkgray'); 
    abline(h=1, col='darkgray')
    #cex = 0.6
    kk = which(aa$half.life>=5); points(aa$half.life[kk], aa$eps.ratio.RD.RS[kk], col='darkgreen', cex=cex)
    kk = which(aa$half.life<5 & aa$delay.RD.RS<12); points(aa$half.life[kk], aa$eps.ratio.RD.RS[kk], col='darkorange', cex=cex)
    kk = which(aa$half.life<5 & aa$delay.RD.RS>=12); points(aa$half.life[kk], aa$eps.ratio.RD.RS[kk], col='darkred', cex=cex)
    axis(2, at=seq(0, 6, by=1), cex.axis=0.8)
    axis(1, at=c(0.5, 2, 5, 10, 20), cex.axis = 0.8)
    #abline(v=1, col='blue', lty=2, lwd=1.5)
    box()
    
    dev.off()
    
   
    ########
    ### Here we show the RD as amplifier and phase tuner with 3 groups
    ########
    circular.error=function(p,x,y){sum(1-cos(2*pi/24*(y-x-p)));}
    kk_1 = which(aa$half.life>=5); #plot(aa$half.life[kk1], aa$delay.RD.RS[kk1], col='darkgreen', cex=cex)
    kk_2 = which(aa$half.life<5 & aa$delay.RD.RS<12); #points(aa$half.life[kk], aa$delay.RD.RS[kk], col='darkorange', cex=cex)
    kk_3 = which(aa$half.life<5 & aa$delay.RD.RS>=12); #points(aa$half.life[kk], aa$delay.RD.RS[kk], col='darkred', cex=cex)
    cex = 0.4; cols =c('darkgreen','darkorange', 'darkred'); col.bg = 'gray80';lwd.bg=1.2; cex.axis = 0.7;
    roles = c('Amplifier', 'Advancer', 'Delayer')
    for(n in c(1:3))
    {
      #n = 2;
      kk = eval(parse(text=paste('kk_', n, sep='')))
      phases1=aa$phase.m22[kk]; phases2=aa$phase.mRNA[kk]; fit=nlm(f=circular.error, p=2, x=phases1, y=phases2);
      
      lims.phase = c(0, 24);
      lims.relamp = c(0.1, 2)
      mgp = c(1.2, 0.5, 0); mar=c(2,2,1.6,0.8);
      wid=1.3;height=1.3;
      ## relative amplitude mrna vs half-life dampened one
      pdfname = paste(folder, roles[n], '_rel_amp_mRNA_vs_m22.pdf', sep='')
      pdf(pdfname, width=wid, height=height)
      
      par(cex = 0.7, las = 1, mgp = mgp, mar = mar+0.1, tcl = -0.3, pty='s')
      if(n==1) {lims = c(0.02, 2);}else{lims=lims.relamp;}
      plot(aa$relamp.m22[kk], aa$relamp.mRNA[kk], log='xy', type='p', xlab=NA, ylab=NA, main=NA, axes = FALSE, cex=cex, xlim=lims, ylim=lims, col=cols[n]);
      abline(0, 1,col=col.bg, lwd=lwd.bg);box();
      #axis(1, cex.axis = cex.axis);
      axis(1,  las = 1, at = c(0.02, 0.1, 0.5, 2), cex.axis=cex.axis, tck=-0.03, labels = FALSE); 
      axis(1, las=1,labels = c('0.02', '0.1', '0.5', '2'), at = c(0.02, 0.1, 0.5, 2),cex.axis=cex.axis, lwd=0, line = -0.5)
      axis(2,  las = 3, at = c(0.02, 0.1, 0.5, 2), cex.axis=cex.axis, tck=-0.03, labels = FALSE); 
      axis(2, las=3, labels = c('0.02', '0.1', '0.5', '2'), at = c(0.02, 0.1, 0.5, 2),cex.axis=cex.axis, lwd=0, line = -0.3)
      dev.off()
      
      ## phase : mrna vs half-life dampened one
      pdfname = paste(folder, roles[n], '_phase_mRNA_vs_m22.pdf', sep='')
      pdf(pdfname, width=wid, height=height)
      par(cex = 0.7, las = 1, mgp = mgp, mar = mar+0.1, tcl = -0.3, pty='s')
      lims = lims.phase; plot(aa$phase.m22[kk], aa$phase.mRNA[kk], log='', type='p', xlab=NA, ylab=NA, main=NA, axes = FALSE, cex=cex, xlim=lims, ylim=lims,
                            col=cols[n]);
      abline(0, 1,col=col.bg, lwd=lwd.bg);abline(-24, 1,col=col.bg, lwd=lwd.bg);abline(24, 1,col=col.bg, lwd=lwd.bg);
      ats=c(0, 6, 12, 18, 24)
      axis(1,  las = 1, at = ats, cex.axis=cex.axis, tck=-0.03, labels = FALSE); 
      axis(1, las=1, at = ats,cex.axis=cex.axis, lwd=0, line = -0.5)
      axis(2,  las = 3, at = ats, cex.axis=cex.axis, tck=-0.03, labels = FALSE); 
      axis(2, las=3, at = ats,cex.axis=cex.axis, lwd=0, line = -0.3)
      box()
      dev.off()
      
      ## relative ampl mrna vs premRNA
      pdfname = paste(folder,roles[n], '_rel_amp_mRNA_vs_RS.pdf', sep='')
      pdf(pdfname, width=wid, height=height)
      par(cex = 0.7, las = 1, mgp = mgp, mar = mar+0.1, tcl = -0.3, pty='s')
      if(n==1) {lims = c(0.02, 2);}else{lims=lims.relamp;}
      plot(aa$relamp.premRNA[kk], aa$relamp.mRNA[kk], log='xy', type='p', xlab=NA, ylab=NA, main=NA, axes = FALSE, cex=cex, xlim=lims, ylim=lims, col=cols[n]);
      abline(0, 1,col=col.bg, lwd=lwd.bg);
      axis(1,  las = 1, at = c(0.02, 0.1, 0.5, 2), cex.axis=cex.axis, tck=-0.03, labels = FALSE); 
      axis(1, las=1,labels = c('0.02', '0.1', '0.5', '2'), at = c(0.02, 0.1, 0.5, 2),cex.axis=cex.axis, lwd=0, line = -0.5)
      #axis(2,  las = 3, at = c(0.02, 0.1, 0.5, 2), cex.axis=cex.axis, tck=-0.03, labels = FALSE); 
      #axis(2, las=3, at = c(0.02, 0.1, 0.5, 2),cex.axis=cex.axis, lwd=0, line = -0.3)
      box()
      dev.off()
   
      ## phase : mrna vs premRNA
      pdfname = paste(folder, roles[n], '_phase_mRNA_vs_RS.pdf', sep='')
      pdf(pdfname, width=wid, height=height)
      par(cex = 0.7, las = 1, mgp = mgp, mar = mar+0.1, tcl = -0.3, pty='s')
      lims = lims.phase; plot(aa$phase.premRNA[kk], aa$phase.mRNA[kk], log='', type='p', xlab=NA, ylab=NA, main=NA, axes = FALSE, 
                            cex=cex, xlim=lims, ylim=lims, col=cols[n]);
      abline(0, 1,col=col.bg, lwd=lwd.bg);abline(-24, 1,col=col.bg, lwd=lwd.bg);abline(24, 1,col=col.bg, lwd=lwd.bg);
      ats=c(0, 6, 12, 18, 24)
      axis(1,  las = 1, at = ats, cex.axis=cex.axis, tck=-0.03, labels = FALSE); 
      axis(1, las=1, at = ats,cex.axis=cex.axis, lwd=0, line = -0.5)
      #axis(2,  las = 3, at = ats, cex.axis=cex.axis, tck=-0.03, labels = FALSE); 
      #axis(2, las=3, at = ats,cex.axis=cex.axis, lwd=0, line = -0.3)
      box()
      dev.off()
    }
    
    
  }
  
  #######
  ### Choose representive examples
  #######
  Illustrative.examples.for.amplitifer.phase.tuner = FALSE
  if(Illustrative.examples.for.amplitifer.phase.tuner)
  {
    res.version = '_total_counts_all_genes_norm_params_v6';
    #load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
    #res.version = '_total_counts_all_genes_norm_params_v6';
    load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
    source('functions.R')
    BIC = my.model.selection.loglike(T=T, method='BIC', outliers = TRUE)
    T = data.frame(T, BIC, stringsAsFactors = FALSE);
    T = T[which(!is.na(T$BIC.best.model)==TRUE), ]
    Tt = T
    
    dd = (aa$phase.mRNA - aa$phase.premRNA)%%24;jj = which(dd>12);dd[jj] = dd[jj]-24;
    aa$delay.mRNA.RS = dd
    kk = which(aa$delay.RD.RS>10 & aa$delay.RD.RS<14 & aa$relamp.premRNA>0.1 & aa$relamp.mRNA>0.2)
    index = aa$index[kk]
    
    examples = c('Smagp', 'Wee1', 'Slc4a4')
    index = aa$index[match(examples, aa$gene)]
    source('functions.R')
    plot.genes.examples(index=index, T=T, Figure=TRUE, folder='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/Examples', 
                        summary.plot=TRUE, RD.advantage = TRUE);
    
    #source('functions.R')
    #plot.genes.examples(index=index, T=T, Figure=TRUE, folder=folder, summary.plot=TRUE);
    
    index = aa$index[kk_2]
    source('functions.R')
    plot.genes.examples(index=index, T=T, Figure=TRUE, folder='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/Examples/model_4_RD_advantage/phase_tuners_left', 
                        summary.plot=TRUE, RD.advantage = TRUE);
    
    index = aa$index[kk_3]
    source('functions.R')
    plot.genes.examples(index=index, T=T, Figure=TRUE, folder='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Figures/Examples/model_4_RD_advantage/phase_tuners_right', 
                        summary.plot=TRUE, RD.advantage = TRUE);
    
  }
  
}

############
#### Supplementary Tables
############
Make.Tables.Supp = FALSE
if(Make.Tables.Supp)
{
  #### Table S2 
  source('functions.R')
  #library(xlsx)
  library(openxlsx)
  
  Data.read.counts = FALSE
  if(Data.read.counts)
  {
    data.version = '_total_counts_v2';
    load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_sel_alphas', data.version, '.Rdata', sep=''))
    #source("/Users/jiwang/Degradation_Liver/Main_Code_Degradation/f24_modified_1.0.r");
    #attribute.global.variable.colors();
    load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_sel_alphas_log2_RPKM', 
                    data.version, '.Rdata', sep=''))
    #T = R.norm;
    #write.xlsx(summar_signi, file="/Users/jiwang/Degradation_Liver/Main_Code_Degradation/Tables/function_analysis_half_life.xlsx",sheetName=paste0("de_model_",i))
    kk = grep('.count.', colnames(R))
    xx = data.frame(R[, c(1:3, kk)], stringsAsFactors = FALSE)
    jj = grep('log2rpkm', colnames(R.norm))
    xx = data.frame(xx, R.norm[, jj], stringsAsFactors = FALSE)
    
    write.table(xx, file=paste0(supDir, "Table_S1_1.txt"), sep='\t', col.names = TRUE, row.names = FALSE, quote=FALSE)
    
  }
  
  #### Main results
  Make.Table.Main.Results = FALSE
  if(Make.Table.Main.Results)
  {
    source('functions.R')
    
    Make.Summary.Results = FALSE
    if(Make.Summary.Results)
    {
      #res.version = '_total_counts_all_genes_norm_params_v6';
      #load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
      #BIC = my.model.selection.loglike(T=T, method='BIC', outliers = TRUE)
      #T = data.frame(T, BIC, stringsAsFactors = FALSE);
      #T = T[which(!is.na(T$BIC.best.model)==TRUE), ]
      #T.old = T;
      #save(T.old, file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '_OLD.Rdata',sep = ''))
      res.version.old = '_total_counts_all_genes_norm_params_v6';
      res.version = '_total_counts_all_genes_norm_params_final_v7';
      #save(T, file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '_OLD.Rdata',sep = ''))
      #Tt = T
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
      
      #kk = grep('.m', colnames(xx))
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
    
    #keepAll = xx;
    #save(keepAll, 
    #     file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_RNA_seq_total_counts_table_sx_All', res.version, '.Rdata',sep = ''))
    
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
  
  #########
  #### function analysis of different models
  ########
  Functional.Analysis.M2.M3.M4
  if(Functional.Analysis.M2.M3.M4)
  {
    Test.function.analysis.Cedric = FALSE
    if(Test.function.analysis.Cedric)
    {
      ### convert gene symbols to ensemble:
      R = read.table('/Users/jiwang/Degradation_Liver/DATA_RNA_Seq/Jingkui_count_IE_Total_040216.txt',sep='\t', header=TRUE)
      #mapping = R[,1]
      xx = unlist(strsplit(as.character(R[,1]),"[|]"))
      mapping = cbind(xx[c(1:nrow(R)*2)], xx[(c(1:nrow(R)*2)-1)])
      #R[,1] = xx
      mm = match(T$gene, mapping[,1])
      T$Ensmus = mapping[mm, 2]
      
      ### I'm working with ensembl genes names as rownames of the data frame 
      library("RDAVIDWebService")
      A=list()
      
      david<-DAVIDWebService(email="cedric.gobet@epfl.ch", 
                             url="http://david-d.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
      setTimeOut(david,1000000) ### Increase time out, you can load more genes like that
      #getAllAnnotationCategoryNames(david)
      setAnnotationCategories(david, c("GOTERM_BP_ALL","KEGG_PATHWAY")) # I defined GO BP ALL and KEGG pathway but you can define more
      # You define your background, list of expressed genes, here RPKM.M.model
      bcknd1<-addList(david, T$Ensmus,idType="ENSEMBL_GENE_ID", listName="Bcknd", listType="Background") 
      
      # GO Enrichment in each model and stock in list A
      for(i in 2:4){
        nam = T$Ensmus[which(T$BIC.best.model==i & T$qv.log2rpkm.mRNA <0.05)]
        if(length(nam)!=0){
          addList(david,nam,idType="ENSEMBL_GENE_ID", listName= paste("model",i,sep="_"), listType="Gene")
          setCurrentBackgroundPosition(david,1)
          res=getFunctionalAnnotationChart(david)
          summar=data.frame(pval=as.numeric(res$PValue),qval=as.numeric(res$Benjamini),term=res$Term,Enr=as.numeric(res$Fold.Enrichment),
                            count=as.numeric(res$Count),Pop=as.numeric(res$Pop.Hits),Genes=res$Genes,stringsAsFactors=FALSE)
          A[[i]]=summar
        }
      }
      
      ########## Parsing of the results , threshold on different value and write in .xls file with different sheets
      library(xlsx)
      pv=0.01; enr=1.5; count=2;  pop=1000;
      v=1
      for(i in 2:4)
      {
        summar=A[[i]]
        signi = which(summar$pval < pv & summar$Enr > enr & summar$count > count & summar$Pop < pop)
        summar_signi = summar[signi,]
        summar_signi$term[grep('\\~',summar_signi$term)] = as.character(unlist(sapply(strsplit(as.character(summar_signi$term[grep('\\~',summar_signi$term)]),"~"),"[[",2)))
        summar_signi$term[grep('\\:',summar_signi$term)] = as.character(unlist(sapply(strsplit(as.character(summar_signi$term[grep('\\:',summar_signi$term)]),":"),"[[",2)))
        
        ### write down the results in an Excel table
        if(length(signi)!=0)
        { 
          if(v==1){
            write.xlsx(summar_signi, file="/Users/jiwang/Dropbox/mRNA_decay/Manuscript/Tables_supp/Table_S3.xlsx", sheetName=paste0("M",i))
          }else{
            write.xlsx(summar_signi, file="/Users/jiwang/Dropbox/mRNA_decay/Manuscript/Tables_supp/Table_S3.xlsx", sheetName=paste0("M",i), append=TRUE)
          }
          v=v+1
        }
      }
      
    }
  }  
  
  ### summary of predicted mRBPs
  Table.SX.mRBPs = FALSE
  if(Table.SX.mRBPs)
  {
    #load(file='/Users/jiwang/epfl/Motif_analysis/mRNA_decay/mRBPs_positive_controls_rhythmic_proteins.Rdata')
    #load(file='/Users/jiwang/epfl/Motif_analysis/mRNA_decay/mRBPs_positive_controls_rhythmic_phospho_Robles_2016.Rdata')
    total = read.table('/Users/jiwang/Dropbox/mRNA_decay/Analysis/Analysis_mRBP_motifs/Table_total_proteins_PNAS.txt', sep='\t', header = TRUE)
    phospho = read.delim('/Users/jiwang/Dropbox/mRNA_decay/Analysis/Analysis_mRBP_motifs/phospho_Robles_2016.txt', header = TRUE)
    robles = read.table('/Users/jiwang/Proteomics_analysis/Nuclear_proteins/Tables_DATA/Table_total_proteins_Robles.txt',sep='\t', header=TRUE)
    
    res.version = '_total_counts_all_genes_norm_params_final_v7_All_summary';
    load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
    T$BIC.best.model = T$best.model;
    
    infer = read.table('/Users/jiwang/Dropbox/mRNA_decay/Analysis/fromJake/Code/good_results_Jingkui_December_2016.runbyjake.txt', sep='\t', header = TRUE, 
                       row.names = 1)
    
    regulators = matrix(NA, nrow=nrow(infer), ncol=7);
    regulators = data.frame(regulators);
    colnames(regulators) = c('mapping.RBP', 'phase.mRNA', 'pval.mRNA',
                             'phase.protein.Mauvoisin', 'pval.protein.Mauvoision', 
                             'phase.phospho.Robles', 'pval.phospho.Robles')
    
    cyclic = c()
    for(n in 1:nrow(infer))
    {
      #n = 2
      test = unlist(strsplit(as.character(rownames(infer)[n]), "_"))[-c(1:2)]
      regulators[n, 1] =  paste0(test, collapse = ';');
      
      kk = match(test, T$gene);kk = kk[which(!is.na(kk)==TRUE)]
      if(length(kk)>0)
      {
        kk = kk[which(T$pval.log2rpkm.mRNA[kk]==min(T$pval.log2rpkm.mRNA[kk]))]
        regulators[n, 2] = T$phase.log2rpkm.mRNA[kk]
        regulators[n, 3] = T$pval.log2rpkm.mRNA[kk]
      }
      
      kk = match(test, total$Gene.names.total);
      kk = kk[which(!is.na(kk)==TRUE)];
      kk = kk[which(!is.na(total$pval.total[kk])==TRUE)];
      if(length(kk)>0){
        kk = kk[which(total$pval.total[kk]==min(total$pval.total[kk]))]
        regulators[n, 4] = total$phase.total[kk]
        regulators[n, 5] = total$pval.total[kk]
      }
      
      kk = match(test, phospho$Gene.names); kk = kk[which(!is.na(kk)==TRUE)];
      kk = kk[which(!is.na(phospho$pval[kk])==TRUE)];
      if(length(kk)>0){
        kk = kk[which(phospho$pval[kk]==min(phospho$pval[kk]))]
        regulators[n, 6] = phospho$phase[kk]
        regulators[n, 7] = phospho$pval[kk]
      }
      
      #nb.cyclic = length(which(control.prot$pval.total[mm]<0.05)) + length(which(control.prot$pval.robles[mm]<0.05)) + length(which(control.phospho.robles$pval[jj]<0.05))
      #cyclic = c(cyclic, nb.cyclic)
    }
    
    yy = data.frame(rownames(infer), infer, regulators, stringsAsFactors = FALSE)
    colnames(yy)[1] = 'RBP.motif'
    write.table(yy, file="/Users/jiwang/Dropbox/mRNA_decay/Manuscript/Tables_supp/Table_S5.txt", sep='\t', col.names = TRUE, row.names = FALSE, quote=FALSE)
    
    length(which(yy$pval.protein.Mauvoision<0.1 | yy$pval.phospho.Robles<0.1))
    
  }
  
  
  #ptm <- proc.time()
  #source('functions.R')
  #res = t(apply(as.matrix(index), 1, compute.phase.amp.4premrna.mrna.with.fitting.res, T=T));
  #proc.time() - ptm
  
  #xx = data.frame(keep, res, stringsAsFactors = FALSE)
  #yy = matrix(NA, nrow=(nrow(T)-nrow(keep)), ncol=ncol(keep))
  #colnames(yy) = colnames(keep)
  #yy = data.frame(yy, stringsAsFactors = FALSE)
  #yy = rbind(yy, keep)
  
  #yy[1:nrow(keep), ] = keep 
  #yy$index = index
  
  
  #write.xlsx(summar_signi, file="/Users/jiwang/Degradation_Liver/Main_Code_Degradation/Tables/function_analysis_half_life.xlsx",sheetName=paste0("de_model_",i))
}






