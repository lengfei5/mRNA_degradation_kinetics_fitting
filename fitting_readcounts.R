##########################################################################
##########################################################################
## Project: fitting temporal profiles of pre-mRNA and mRNA with a kinetic model to infer the rhythmic transcriptional and post-transcriptional regulation
## Script purpose:
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Fri Jun  1 12:01:00 2018
##########################################################################
##########################################################################

######################################
######################################
## Section: prepare the read count table for fitting
## 1) Import read counts of total RNA-seq data
## 2) filter non-expressed mRNAs
## 3) estimate the dispersion parameters using DESeq2
######################################
######################################
RNA.Seq = FALSE
if(RNA.Seq)
{
    ### Processing the Read Counts Table
    R = read.table('/Users/jiwang/Degradation_Liver/DATA_RNA_Seq/Jingkui_count_IE_Total_040216.txt',sep='\t', header=TRUE)
    gene = R[,1]
    xx = unlist(strsplit(as.character(gene),"[|]"))
    yy =xx[c(1:length(gene)*2)-1]
    xx = xx[c(1:length(gene)*2)]
    R[,1] = xx
    
    R[which(R[,1]=='Npas2'), ]
    R[which(R[,1]=='Bhlhe41'), ]
    
    #### filter expressed genes with RPKM
    filter.expressed.genes.RPKM = FALSE
    if(filter.expressed.genes.RPKM)
    {
      ### filter genes with length of intron or exon smaller than 250 bps
      kk = which(R[,2]>250 & R[,3]>250)
      R = R[kk, ]
      
      ii = grep('Exon', colnames(R))
      jj = grep('Intron', colnames(R))
      exon = as.matrix(R[,ii])
      intron = as.matrix(R[,jj])
      #exon = as.matrix(T[, ii])
      #intron = as.matrix(T[, jj])
      #### Normalized read counts table into RPKM for the selection of expressed genes and also for the downstream analysis
      #xx = diag(ss.norm, nrow=length(ss.norm), ncol=length(ss.norm))
      #ss.norm = xx;
      data = rbind(exon, intron)
      require('edgeR')
      y <- DGEList(counts=data)
      y <- calcNormFactors(y, method=c("TMM"))
      RPKM = rpkm(y, gene.length=c(R$length_exon, R$length_intron), normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25)
      #RPKM = cpm(y, normalized.lib.sizes=TRUE, log=FALSE)
      exon = RPKM[c(1:nrow(exon)), ]
      intron = RPKM[-c(1:nrow(exon)), ]
      source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')
      stat1 = t(apply(exon, 1, f24_R2_alt2, t=c(0:47)*2))
      stat2 = t(apply(intron,1, f24_R2_alt2, t=c(0:47)*2))
      stat1 = cbind(stat1, qvals(stat1[,6]))
      stat2 = cbind(stat2, qvals(stat2[,6]))
    }
    
    #### Filtering expressed genes using read counts
    Selection.read.counts = TRUE
    if(Selection.read.counts)
    {
      ii = grep('Exon', colnames(R))
      jj = grep('Intron', colnames(R))
      exon = as.matrix(R[,ii])
      intron = as.matrix(R[,jj])
      
      fitering.expressed.gene = function(x, threshold=10)
      {
        return(length(which(x>=threshold))>=24)
      }
      ss1 = apply(exon, 2, sum)
      ss2 = apply(intron, 2, sum)
      ss = ss1 + ss2
      ss.norm = ss/min(ss)
      data1 = exon
      data2 = intron
      for(n in 1:length(ss.norm))
      {
        data1[,n] = exon[,n]*ss.norm[n]
        data2[,n] = intron[,n]*ss.norm[n]
      }
      sel.exon = apply(data1, 1, fitering.expressed.gene)
      sel.intron = apply(data2, 1, fitering.expressed.gene)
      
      kk = which(sel.exon==TRUE & sel.intron==TRUE & R[,2]>250 & R[,3]>250)
      
      R = R[kk, ]
    }
    
    examples = c('Npas2', 'Arntl', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry1', 'Cry2', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Nfil3', 'Bhlhe40', 'Bhlhe41', 'Nampt', 'Parp1', 'Prkaa2', 'Prkag1', 'Gsk3a', 'Gsk3b', 'Csnk1d')
    mm = match(examples, R[,1])
    
    data.version = '_total_counts_v2'
    save(R, file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_prepared_expressed_genes', data.version, '.Rdata', sep=''))
    
    ######
    ### estimate dispersion parameter with Deseq or edgeR
    ######
    data.version = '_total_counts_v2'
    load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_prepared_expressed_genes', data.version, '.Rdata', sep=''))
    
    ii = grep('Exon', colnames(R))
    jj = grep('Intron', colnames(R))
    exon = as.matrix(R[,ii])
    intron = as.matrix(R[,jj])
    countData =  rbind(exon, intron)
    
    #### Estimate dispersion parameter with Deseq2
    shared.dispersion.premRNA.mRNA = FALSE
    if(!shared.dispersion.premRNA.mRNA)
    {
      require('DESeq2')
      size.factors = estimateSizeFactorsForMatrix(countData)
      
      alphas = matrix(NA, nrow=nrow(countData), ncol=12)
      
      for(n in 1:ncol(alphas))
      {
        cat('ZT ', (n-1)*2, '\n');
        index = (c(0:3)*12+n);
        countData2 = countData[, index]
        condition <- factor(rep(c(1), 4))
        dds <- DESeqDataSetFromMatrix(countData2, DataFrame(condition), ~1)
        sizeFactors(dds) = size.factors[index]
        dds <- estimateDispersions(dds, maxit = 500)
        alphas[,n] = dispersions(dds);
      }
      
      xx = data.frame(alphas[c(1:nrow(R)),], alphas[-c(1:nrow(R)), ])
      rownames(alphas) = R[,1]
        
      examples = c('Npas2', 'Arntl', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry1', 'Cry2', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Nfil3', 'Bhlhe40', 'Bhlhe41', 'Nampt', 'Parp1', 'Prkaa2', 'Prkag1', 'Gsk3a', 'Gsk3b', 'Csnk1d')
      mm = match(examples, R[,1])
      R[mm, c(1:5)]
      
      n = 9
      plot(c(0:11)*2, alphas[mm[n], ], type='b', ylim=c(0.001, 2), log='y', main=R$gene_name[mm[n]])
      points(c(0:11)*2, countData[mm[n], c(1:12)]/mean(countData[mm[n], c(1:12)]), type='b', col='blue')
      #plot(countData[mm[1], c(1:12)], alphas[mm[1],], log='xy')
      
      condition <- factor(rep(c(1:12), 4))
      dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
      #condition <- factor(rep(1, 48))
      #dds <- DESeqDataSetFromMatrix(countData)
      
      dds <- estimateSizeFactors(dds)
      sizeFactors(dds)
      dds <- estimateDispersions(dds,  fitType = c("local"), maxit = 500)
      
      head(dispersions(dds))
      xx = dispersions(dds)
      #dds <- estimateDispersionsGeneEst(dds)
      yy = mcols(dds)$dispGeneEst
      
      examples = c('Npas2', 'Arntl', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry1', 'Cry2', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Nfil3', 'Bhlhe40', 'Bhlhe41', 'Nampt', 'Parp1', 'Prkaa2', 'Prkag1', 'Gsk3a', 'Gsk3b', 'Csnk1d')
      mm = match(examples, R[,1])
      plot(xx, yy, cex=0.2, log='xy', ylim=c(0.001, 1));
      points(xx[mm], yy[mm], col='blue', cex=1.5)
      points(xx[c(mm+nrow(R))], yy[c(mm+nrow(R))], col='green', cex=1.5)
      text(xx[mm], yy[mm], examples, col='blue', cex=1., offset=4.0)
      abline(0, 1, col='red', lwd=2.0)
      
      
      xx = yy
      
      RR = cbind(R[, c(1:3)], xx[c(1:nrow(R))], xx[-c(1:nrow(R))], exon, intron)
      colnames(RR)[c(4:5)] = c('alpha.mRNA', 'alpha.premRNA')
      R = RR
      
      examples = c('Npas2', 'Arntl', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry1', 'Cry2', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Nfil3', 'Bhlhe40', 'Bhlhe41', 'Nampt', 'Parp1', 'Prkaa2', 'Prkag1', 'Gsk3a', 'Gsk3b', 'Csnk1d')
      mm = match(examples, R[,1])
      R[mm, c(1:5)]
      
      par(mfrow=c(1,1),cex=1.0)
      plot(R$alpha.premRNA, R$alpha.mRNA, cex=0.05, log='xy')
      abline(0, 1, col='darkgray', lwd=2.0)
      mm = match(examples, R[,1])
      abline(h=median(R$alpha.mRNA), col='darkgray', lwd=2.0)
      abline(v=median(R$alpha.premRNA), col='darkgray', lwd=2.0)
      points(R$alpha.premRNA[mm], R$alpha.mRNA[mm], cex=1., col='red')
      text(R$alpha.premRNA[mm], R$alpha.mRNA[mm], examples, cex=0.8, col='blue', offset=0.7, pos=3)
      
    }else{
      require('DESeq2')
      size.factors = estimateSizeFactorsForMatrix(countData)
      
      countData =  cbind(exon, intron)
      condition <- factor(c(rep(c(1:12), 4), rep(c(13:24), 4)))
      dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
      sizeFactors(dds) = rep(size.factors, 2)
      
      dds <- estimateDispersions(dds)
      
      head(dispersions(dds))
      xx = dispersions(dds)
      #dds <- estimateDispersionsGeneEst(dds)
      yy = mcols(dds)$dispGeneEst
      
      RR = cbind(R[, c(1:3)], xx, xx, exon, intron)
      colnames(RR)[c(4:5)] = c('alpha.mRNA', 'alpha.premRNA')
      R = RR
      
      examples = c('Npas2', 'Arntl', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry1', 'Cry2', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Nfil3', 'Bhlhe40', 'Bhlhe41', 'Nampt', 'Parp1', 'Prkaa2', 'Prkag1', 'Gsk3a', 'Gsk3b', 'Csnk1d')
      mm = match(examples, R[,1])
      R[mm, c(1:5)]
      
      par(mfrow=c(1,1),cex=1.0)
      plot(R$alpha.premRNA, R$alpha.mRNA, cex=0.05, log='xy')
      abline(0, 1, col='darkgray', lwd=2.0)
      mm = match(examples, R[,1])
      abline(h=median(R$alpha.mRNA), col='darkgray', lwd=2.0)
      abline(v=median(R$alpha.premRNA), col='darkgray', lwd=2.0)
      points(R$alpha.premRNA[mm], R$alpha.mRNA[mm], cex=1., col='red')
      text(R$alpha.premRNA[mm], R$alpha.mRNA[mm], examples, cex=0.8, col='blue', offset=0.7, pos=3)      
    }
    
    ############
    ##### statistics calculated with RPKM from Deseq2
    ############
    ii = grep('Exon', colnames(R))
    jj = grep('Intron', colnames(R))
    exon = as.matrix(R[,ii])
    intron = as.matrix(R[,jj])
    
    ### rpkm calculated from DESeq2
    require('DESeq2')
    countData =  rbind(exon, intron)
    size.factors = estimateSizeFactorsForMatrix(countData)
    #counts.norm = countData/size.factors
    #exon = counts.norm[c(1:nrow(R)), ]
    #intron = counts.norm[-c(1:nrow(R)), ]
    condition <- rep(factor(colnames(countData)[1:12]), 4)
    dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
    
    mcols(dds)$basepairs <- c(R$length_exon, R$length_intron)
    rpkm = fpkm(dds, robust=TRUE); 
    
    #### rpkm calculated from myself
    ss = apply(countData, 2, sum)
    scale.factors = size.factors*mean(ss/size.factors)
    
    data1 = matrix(NA, ncol=ncol(exon), nrow=nrow(exon))
    data2 = matrix(NA, ncol=ncol(intron), nrow=nrow(intron))
    
    for(n in 1:nrow(exon))
    {
      data1[n, ] = exon[n, ]/scale.factors/R$length_exon[n]*10^9;
      data2[n, ] = intron[n, ]/scale.factors/R$length_intron[n]*10^9;
    }
    
    #ii = grep('Exon', colnames(R))
    #jj = grep('Intron', colnames(R))
    #data1 = R[,ii]
    #data2 = R[,jj]
    source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')
    stat1 = t(apply(data1,1, f24_R2_alt2, t=c(0:47)*2))
    stat2 = t(apply(data2,1, f24_R2_alt2, t=c(0:47)*2))
    
    stat1 = cbind(stat1, qvals(stat1[,6]))
    stat2 = cbind(stat2, qvals(stat2[,6]))
    
    colnames(stat1) = paste(c('nb.timepoints', 'mean', 'amp', 'rel.amp', 'phase','pval', 'qv'), '.rpkm.mRNA', sep='')
    colnames(stat2) = paste(c('nb.timepoints', 'mean', 'amp', 'rel.amp', 'phase','pval', 'qv'), '.rpkm.premRNA', sep='')
    
    R = data.frame(R, stat1, stat2, stringsAsFactors=FALSE)
    R = R[, c(1, 3, 2, 4:115)]
    
    colnames(R)[c(1,2, 3)] = c('gene', 'length.mRNA', 'length.premRNA')
    colnames(R)[6:53] = paste('ZT', c(0:47)*2, '.count.mRNA', sep='')
    colnames(R)[54:101] = paste('ZT', c(0:47)*2, '.count.premRNA', sep='')
    #stat = colnames(R)[50:56]
    #stat = unlist(strsplit(as.character(stat), '.'))
    
    #kk = which(R$mean.rpkm.mRNA>R$mean.rpkm.premRNA)
    #R = R[kk,]
    
    #### save the table for the model selection and parameter estimation
    data.version = '_total_counts_v2'
    save(R, file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis', data.version, '.Rdata', sep=''))
    save(R, file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_sel', data.version, '.Rdata', sep=''))
    
    scaling.factors = scale.factors;
    save(scaling.factors, file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/model_RNA_seq_total/Scaling_factors_48_samples.Rdata')
    save(ss, file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/model_RNA_seq_total/Libary_size_48_samples.Rdata')
    save(size.factors, file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/model_RNA_seq_total/size_factors_48_samples.Rdata')
    #load(file=paste('myRdata/my_genes_RNA_seq_analysis', data.version, '.Rdata', sep=''))
    
    
    test.condition.dependent.alpha = TRUE
    if(test.condition.dependent.alpha)
    {
      Compute.alphas.each.condition = FALSE
      if(Compute.alphas.each.condition){
        ii = grep('.count.mRNA', colnames(T))
        jj = grep('.count.premRNA', colnames(T))
        exon = as.matrix(T[,ii])
        intron = as.matrix(T[,jj])
        countData =  rbind(exon, intron)
        require('DESeq2')
        load(file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/model_RNA_seq_total/size_factors_48_samples.Rdata')
        
        ### estimate dispersion parameters with Deseq2 for each gene and each condition
        alphas = matrix(NA, nrow=nrow(countData), ncol=12)
        alphas.genes = matrix(NA, nrow=nrow(countData), ncol=12)
        
        for(n in 1:ncol(alphas))
        {
          cat('ZT ', (n-1)*2, '\n');
          index = (c(0:3)*12+n);
          countData2 = countData[, index];
          condition <- factor(rep(c(1), 4))
          dds <- DESeqDataSetFromMatrix(countData2, DataFrame(condition), ~1);
          sizeFactors(dds) = size.factors[index]
          dds <- estimateDispersions(dds, maxit = 500)
          alphas[,n] = dispersions(dds);
          alphas.genes[,n] = mcols(dds)$dispGeneEst
        }
        
        xx = data.frame(alphas[c(1:nrow(T)),], alphas[-c(1:nrow(T)), ])
        colnames(xx) = c(paste('alpha.mRNA.ZT', c(0:11)*2, sep=''),  paste('alpha.premRNA.ZT', c(0:11)*2, sep=''))
        alphas = xx;
        
        T = data.frame(T[,c(1:3)], alphas, T[, -c(1:5)], stringsAsFactors = FALSE)
        
        T$mRNA.outlier = '';
        T$premRNA.outlier = '';
        
        
        #save(alphas, alphas.genes, file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/Deseq2_estimation_alpha_conditions', 
        #                       data.version, '.Rdata', sep=''))
        
        load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/Deseq2_estimation_alpha_conditions', 
                        data.version, '.Rdata', sep=''))
        xx = data.frame(alphas[c(1:nrow(T)),], alphas[-c(1:nrow(T)), ])
        colnames(xx) = c(paste('alpha.mRNA.ZT', c(0:11)*2, sep=''),  paste('alpha.premRNA.ZT', c(0:11)*2, sep=''))
        alphas = xx;
        
        T$mRNA.outlier = '';
        T$premRNA.outlier = '';
        
        kk = grep('alpha', colnames(T))
        T = data.frame(T[,c(1:3)], alphas, T[, -kk], stringsAsFactors = FALSE)
        T = T[, -c(28:30)]
        R = T;
        data.version = '_total_counts_v2'
        #save(R, file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_sel_alphas', data.version, '.Rdata', sep=''))
        
      }
    }  
    
}

######################################
######################################
## Section: fitting the read count from RNA-seq data
######################################
######################################
Real.Data.Fitting = FALSE
if(Real.Data.Fitting)
{
  cat('TEST THE FIT on THE REAL DATA\n')
  data.version = '_total_counts_v2'
  load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_sel_alphas', data.version, '.Rdata', sep=''))
  source("/Users/jiwang/Degradation_Liver/Main_Code_Degradation/f24_modified_1.0.r")
  T = R;
  
  Save.table = FALSE
  if(Save.table)
  {
    ZT.int = grep('.count.premRNA', colnames(T))
    ZT.ex = grep('.count.mRNA', colnames(T))
    exon = as.matrix(T[, ZT.ex])
    intron = as.matrix(T[, ZT.int])
    #xx = data.frame(T[, c(1:3)], T[, ZT.ex], T[, ZT.int], stringsAsFactors = FALSE)
    #write.table(xx,  file = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/Tables/Raw_data_total_RNA_seq.txt', 
    #            quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
    load(file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/model_RNA_seq_total/Scaling_factors_48_samples.Rdata')
    load(file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/model_RNA_seq_total/Libary_size_48_samples.Rdata')
    #save(size.factors, file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/model_RNA_seq_total/size_factors_48_samples.Rdata')
    #ss = apply(countData, 2, sum)
    #scale.factors = size.factors*mean(ss/size.factors)
    
    data1 = matrix(NA, ncol=ncol(exon), nrow=nrow(exon))
    data2 = matrix(NA, ncol=ncol(intron), nrow=nrow(intron))
    
    for(n in 1:nrow(exon))
    {
      data1[n, ] = exon[n, ]/scaling.factors/R$length.mRNA[n]*10^9;
      data2[n, ] = intron[n, ]/scaling.factors/R$length.premRNA[n]*10^9;
    }
    
    colnames(data1) = paste('ZT', seq(0, 94, by=2), '.rpkm.mRNA', sep='')
    colnames(data2) = paste('ZT', seq(0, 94, by=2), '.rpkm.premRNA', sep='')
    xx = data.frame(T[, c(1:3)], data1, data2, stringsAsFactors = FALSE)
    #write.table(xx,  file = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/Tables/Raw_data_total_RNA_seq_RPKM_normalized.txt', 
    #                      quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  }
  
  examples = c('Per2', 'Cry1', 'Per3','Per1', 'Npas2', 'Arntl', 'Clock', 'Cry2', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Rorc','Rora',    
               'Nfil3', 'Bhlhe40', 'Bhlhe41', 'Nampt', 'Parp1', 'Gsk3a','Csnk1d', 
               'Por', 'Tymp', 'Ppara', 'Hnrnpc', 'Elmo3', 'Camta1', 'Abcb11', 'Tfrc', 'Loxl4', 'Rcan1', 'Cdkn1a', 'Tubb2a', 'Mfsd2', 'Ppard',
               'Nedd4l', 'Cbs', 'Gsta3', 'Eps8l2', 'Insig2', 'Rbm3')
  mm = match(unique(examples), T[,1])
  #mm = mm[order(-T[mm, 5])]
  mm = mm[which(!is.na(mm)==TRUE)]
  T[mm, c(1:3)]
  
  ZT.int = grep('.count.premRNA', colnames(T))
  ZT.ex = grep('.count.mRNA', colnames(T))
  zt = seq(0,94,by = 2)
  
  source('functions.R')
  
  keep5 = c()
  #for(j in match(controls, T[,1]))
  for(j in mm)
  {
    cat(T[j, 1], '\n');
    
    gg = 'Npas2'
    j = which(T$gene==gg)
    source('functions.R')
    #T$mRNA.outlier[j] = paste(outlier.m, sep='', collapse = ';')
    #T$premRNA.outlier[j] = paste(outlier.s, sep='', collapse = ';')
    ptm <- proc.time()
    param.fits.results = make.fits.with.all.models.for.one.gene.remove.outliers(T = T, gene.index = j, debug = TRUE,
                                                                                zt = zt, i.ex = ZT.ex, i.int = ZT.int, outliers = TRUE);
    proc.time() - ptm
    
    index = match(c('outlier.m', 'outlier.s'), names(param.fits.results))
    res.fit = as.numeric(param.fits.results[-index])
    names(res.fit) = names(param.fits.results)[-index]
    
    missing.data = length(unlist(strsplit(param.fits.results[index[1]], ';'))) + length(unlist(strsplit(param.fits.results[index[2]], ';')))
    test1 = c(res.fit, my.model.selection.one.gene.loglike(res.fit, nb.data = (96-missing.data), method = 'BIC'))
    print(test1);
    m = test1[which(names(test1)=='BIC.best.model')]
    name1 = paste('gamma.m', m, sep='');name2 = paste('gamma.stderr.m', m, sep='')
    eval(parse(text=paste('print(c(half.lfie=log(2)/test1[which(names(test1)==name1)], test1[which(names(test1)==name2)]/test1[which(names(test1)==name1)]))', sep='')));
    
  }
  
  keep5 = data.frame(T[mm[c(1:14)], 1], keep5, stringsAsFactors = FALSE)
  keep5[, c(1, grep('best.model', colnames(keep5)))]
  #keep4 = keep4[c(1:4),]
  
  
  
  
  #save(keep4, file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_core_clock_fitting_results_alpha_conditions_MLE_outliers_v4', data.version, '.Rdata', sep=''))
  #rownames(keep2) = T[mm, 1]
  
  keep1 = c()
  for(j in mm)
  {
    ptm <- proc.time()
    param.fits.results = make.fits.with.all.models.for.one.gene(T = T, gene.index = j, debug = TRUE, zt = zt, i.ex = ZT.ex, i.int = ZT.int);
    proc.time() - ptm
    #c(param.fits.results, my.BIC.this.gene(param.fits.results))
    c(param.fits.results, my.AIC.this.gene.loglike(param.fits.results),
      my.AIC.this.gene.loglike(param.fits.results, correction=TRUE), 
      my.BIC.this.gene.loglike(param.fits.results))
    keep1 = rbind(keep1, c(T[, c(1:5)], param.fits.results, my.AIC.this.gene.loglike(param.fits.results),
                           my.AIC.this.gene.loglike(param.fits.results, correction=TRUE), 
                           my.BIC.this.gene.loglike(param.fits.results)))
  }
  
  xx = keep1
  xx = data.frame(T[mm, c(1:3)], keep1[, -c(1:5)], stringsAsFactors = FALSE)
  keep1 = xx
  
  keep1[, c(1, match(c('AIC.best.model', 'AICc.best.model', 'BIC.best.model'), colnames(keep1)))]
  
  #save(keep1, file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_core_clock_fitting_results_alpha_conditions', data.version, '.Rdata', sep=''))
  

  #### check fitting quality
  i = j
  zt = seq(0, 94, by=2)
  source('functions.R')
  zt.p = seq(0, 94, by=2)
  model = 4;
  gamma = param.fits.results[27];eps.gamma = param.fits.results[28];phase.gamma = param.fits.results[29];splicing.k =  param.fits.results[30];
  Min = param.fits.results[31]; Amp=param.fits.results[32]; phase=param.fits.results[33]; beta=param.fits.results[34];
  m1 = compute.m.beta(zt.p, gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
  s1 = compute.s.beta(zt.p, Min, Amp, phase, beta);
  #read.m4 = convert.nb.reads(m1, R$length.mRNA[i])
  #read.s4 = convert.nb.reads(s1, R$length.premRNA[i])
  read.s4 = s1;
  read.m4 = m1;
  
  gamma = param.fits.results[3];eps.gamma = 0;phase.gamma = 12;splicing.k =  param.fits.results[4];
  Min = param.fits.results[5]; Amp=param.fits.results[6]; phase=param.fits.results[7]; beta=param.fits.results[8];
  m2 = compute.m.beta(zt.p, gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
  s2 = compute.s.beta(zt.p, Min, Amp, phase, beta);
  #read.m2 = convert.nb.reads(m2, R$length.mRNA[i])
  #read.s2 = convert.nb.reads(s2, R$length.premRNA[i])
  read.s2 = s2;
  read.m2 = m2;
  
  L.m = T$length.mRNA[i];
  L.s = T$length.premRNA[i];
  M = norm.RPKM(as.numeric(T[i,ZT.ex]), L.m)
  S = norm.RPKM(as.numeric(T[i,ZT.int]), L.s)
  
  par(mfcol=c(1,3),cex=1.0)
  lims = range(S, read.s4, read.s2)
  plot(zt, as.numeric(S), col='darkblue',type='b',ylim=lims, lty=1, main=paste(T$gene[i], ', premRNA', sep=''), log='y')
  points(zt.p, read.s4, col='violet', type='l',lty=1, lwd=2.0)
  points(zt.p, read.s2, col='green', type='l',lty=2, lwd=2.0)
  #abline(h=1.0,col='gray',lwd=2.0)
  lims = range(M, read.m4, read.m2)
  plot(zt, as.numeric(M), col='darkblue',type='b',lty=1, ylim=lims, main=paste(T$gene[i], ', mRNA', sep=''), log='y')
  points(zt.p, read.m4, col='violet', type='l',lty=1, lwd=2.0)
  points(zt.p, read.m2, col='green', type='l',lty=2, lwd=2.0)        
  
  M = M/mean(M)
  S = S/mean(S)
  
  lims = range(S, M)
  plot(zt, S, col='green', type='b', ylim=lims, main=paste(T$gene[i], ', premRNA', sep=''), log='')
  points(zt, M, col='blue', type='b',lty=1, lwd=2.0)
  abline(h=1, col='darkgray', lwd=2.0)
  #points(zt.p, s2, col='violet', type='l',lty=1, lwd=2.0)
  #points(zt.p, s1, col='green', type='l',lty=2, lwd=2.0)
  
  #lims = range(M, m2, m1)
  #plot(zt, M, col='darkblue', type='b', ylim=lims, main=paste(T$gene[i], ', mRNA', sep=''), log='y')
  #points(zt.p, m2, col='violet', type='l',lty=1, lwd=2.0)
  #points(zt.p, m1, col='green', type='l',lty=2, lwd=2.0)
  #eval(parse(text = paste('gamma = param.fits.results[which(names(param.fits.results)=='')];
  #T$gamma[i];eps.gamma = T$eps.gamma[i];phase.gamma = T$phase.gamma[i];splicing.k = T$splicing.k[i]; 
  #mean=T$mean.int[i];fold.change = T$fold.change.int[i]; phase.int = T$phase.int[i]; beta = T$beta.int[i]', sep = '')))  
  
  
  ### check differences between RPKM and Deseq2 normalization
  Check.normalization = FALSE
  if(Check.normalization)
  {
    gg = 'Camta1'
    #my_genes_RNA_seq_analysis_total_counts_v1.Rdata
    data.version = '_total_v1'
    load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_sel', data.version, '.Rdata', sep=''))
    T0 = R
    ZT.int = grep('.abs.premRNA', colnames(T0))
    ZT.ex = grep('.abs.mRNA', colnames(T0))
    zt = seq(0,94,by = 2)
    j0 = which(T0$gene==gg)
    #print(T0[j0,])
    print(T0[j0, match(c('rel.amp.mRNA', 'rel.amp.premRNA'), colnames(T0))])
    #ZT.int = grep('.count.premRNA', colnames(T0))
    #ZT.ex = grep('.count.mRNA', colnames(T0))
    #zt = seq(0,94,by = 2)
    M = (as.numeric(T0[j0,ZT.ex]))
	  S = (as.numeric(T0[j0,ZT.int]))
	  M = M/mean(M)
	  S = S/mean(S)
    lims = range(S, M)
    
    par(mfcol=c(1,3),cex=1.0)
	  plot(zt, S, col='green', type='b', ylim=lims, main=paste(T0$gene[j0], ', RPKM', sep=''), log='')
	  points(zt, M, col='blue', type='b',lty=1, lwd=2.0)
	  abline(h=1, col='darkgray', lwd=2.0)
	  #points(zt.p, s2, col='violet', type='l',lty=1, lwd=2.0)
	  
    data.version = '_total_counts_v2'
    load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_sel', data.version, '.Rdata', sep=''))
    source("/Users/jiwang/Degradation_Liver/Main_Code_Degradation/f24_modified_1.0.r")
    T = R;
    source('functions.R')
    #best.model = model.check
    ZT.int = grep('.count.premRNA', colnames(T))
    ZT.ex = grep('.count.mRNA', colnames(T))
    zt = seq(0,94,by = 2)
    j = which(T$gene==gg)
    print(T[j, match(c('rel.amp.rpkm.mRNA', 'rel.amp.rpkm.premRNA'), colnames(T))])
    
    L.m = T$length.mRNA[j];
	  L.s = T$length.premRNA[j];
    
    M = norm.RPKM.libary.size(as.numeric(T[j,ZT.ex]), L.m)
	  S = norm.RPKM.libary.size(as.numeric(T[j,ZT.int]), L.s)
	  M = M/mean(M)
	  S = S/mean(S)
    #lims = range(S, M)
	  plot(zt, S, col='green', type='b', ylim=lims, main=paste(T$gene[j], ', RPKM #2', sep=''), log='')
	  points(zt, M, col='blue', type='b',lty=1, lwd=2.0)
	  abline(h=1, col='darkgray', lwd=2.0)
    
	  M = norm.RPKM(as.numeric(T[j,ZT.ex]), L.m)
	  S = norm.RPKM(as.numeric(T[j,ZT.int]), L.s)
	  M = M/mean(M)
	  S = S/mean(S)
    #lims = range(S, M)
	  plot(zt, S, col='green', type='b', ylim=lims, main=paste(T$gene[j], ', Scaling factors', sep=''), log='')
	  points(zt, M, col='blue', type='b',lty=1, lwd=2.0)
	  abline(h=1, col='darkgray', lwd=2.0)
    
	  par(mfcol=c(1,1),cex=1.0)
	  source('functions.R')
	  load(file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/model_RNA_seq_total/Libary_size_48_samples.Rdata')
    set.scaling.factors()
	  plot(c(0:47)*2, 1/(scaling.factors/mean(scaling.factors)), type='b', col='violet')
	  points(c(0:47)*2, 1/(ss/mean(ss)), type='b', col='black')
	  abline(h=1, col='red', lwd=2.0)
  }
   
}
