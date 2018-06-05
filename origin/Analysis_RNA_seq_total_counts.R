####################### Front:: HERE we started again this project and we really hope that we can finish it this time !!!
### install some packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("fdrtool", "circular", "preprocessCore", "gtools", "biomaRt", "numDeriv", "Matrix")
ipak(packages)

#############
####### Import read counts of total RNA-seq data
#############
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
      
      #countData2 =  cbind(exon, intron)
      #counts.norm = countData/size.factors
      #exon = counts.norm[c(1:nrow(R)), ]
      #intron = counts.norm[-c(1:nrow(R)), ]
      #cnts <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10)
      #cond <- factor(rep(1:2, each=5))
      # object construction
      #dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)
      #condition <- factor(c(rep(c(1:12), 4), rep(c(13:24), 4)))
      
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
   
    Compare.Deseq.edgeR = FALSE
    if(Compare.Deseq.edgeR)
    {
      require('edgeR')
      #xx1 = xx
      #yy1 = yy
      #y <- matrix(rnbinom(1000, mu=10, size=5), ncol=4)
      group <- factor(rep(c(1:12), 4))
      design <- model.matrix(~group-1)
      norm.factors = calcNormFactors(countData)
      
      #ss = apply(countData, 2, sum)
      
      #par(mfrow=c(1,3),cex=1.0)
      #plot(ss/mean(ss), norm.factors*ss/mean(ss), cex=0.7);abline(0, 1, col='red', lwd=2.0)
      #plot(ss/mean(ss), size.factors, cex=0.7);abline(0, 1, col='red', lwd=2.0)
      #plot(norm.factors*ss/mean(ss), size.factors, cex=0.7);abline(0, 1, col='red', lwd=2.0)
      #d <- DGEList(counts=countData, norm.factors = norm.factors, group=group)
      d <- DGEList(counts=intron, norm.factors = norm.factors, group=group)
      #d1 <- estimateDisp(d)
      d1 <- estimateDisp(d, design, trend.method="locfit")
      
      d <- DGEList(counts=exon, norm.factors = norm.factors, group=group)
      #d1 <- estimateDisp(d)
      d2 <- estimateDisp(d, design, trend.method="locfit")
      
      
      xx1 = d1$tagwise.dispersion
      d2 = estimateDisp(d, design, trend.method="locfit", robust=TRUE)
      xx2 = d2$tagwise.dispersion
      
      d3 <- estimateGLMTrendedDisp(d, design, min.n=10)
      d3 <- estimateGLMTagwiseDisp(d3, design)
      xx3 = d3$tagwise.dispersion
      
      d4 = estimateGLMRobustDisp(d, design)
      xx4 = d4$tagwise.dispersion
      
      xx2 = d2$tagwise.dispersion
      yy2 = d2$common.dispersion
      zz2 = d2$trended.dispersion
      
      par(mfrow=c(1,3),cex=1.0)
      plot(xx1, xx2, cex=0.2, log='xy', main='edgeR');abline(0, 1, col='red')
      plot(xx1, xx3, cex=0.2, log='xy', main='edgeR');abline(0, 1, col='red')
      plot(xx1, xx4, cex=0.2, log='xy', main='edgeR');abline(0, 1, col='red')
      
      xx = zz2
      
      par(mfrow=c(2,2),cex=1.0)
      plot(xx1[c(1:nrow(R))], xx1[-c(1:nrow(R))], cex=0.4, log='xy', main='Deseq2');abline(0, 1, col='red')
      plot(xx2[c(1:nrow(R))], xx2[-c(1:nrow(R))], cex=0.4, log='xy', main='edgeR');abline(0, 1, col='red')
      plot(xx1[c(1:nrow(R))], xx2[c(1:nrow(R))], cex=0.4, log='xy', main='Deseq2 vs edgeR (mRNA)');abline(0, 1, col='red')
      plot(xx1[-c(1:nrow(R))], xx2[-c(1:nrow(R))], cex=0.4, log='xy', main='Deqseq2 vs edgeR (premRNA)');abline(0, 1, col='red')
      
      RR = cbind(R[, c(1:3)], xx4[c(1:nrow(R))], xx4[-c(1:nrow(R))], exon, intron)
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
    Prepare.NF.KO = FALSE
    if(Prepare.NF.KO)
    {
        prepare.NF.KO.tables();
    }
    
}

################### (TO DO)
##### Data Control and Validation for the total RNA-seq data by comparing with qPCR (David Gatefield), ExonArray, nascent-RNA-seq
###################
Data.Control = FALSE
if(Data.Control)
{
    data.version = '_total_v2'
    load(file=paste('myRdata/my_genes_RNA_seq_analysis_sel', data.version, '.Rdata', sep=''))
    
    source("f24_modified_1.0.r")
    source('model_RNA_seq_total/functions.R')
    
    ###### Summary of data: expression of premRNAs and mRNAs, phases of rhythmic premRNAs and mRNAs
    Tt = R;
   
    ii = grep('abs.mRNA', colnames(Tt))
    jj = grep('abs.premRNA', colnames(Tt))
    data1 = log2(as.matrix(Tt[,ii]))
    data2 = log2(as.matrix(Tt[,jj]))
    
    stat1 = t(apply(data1,1, f24_R2_alt2, t=c(0:47)*2))
    stat2 = t(apply(data2,1, f24_R2_alt2, t=c(0:47)*2))
    stat1 = cbind(stat1, qvals(stat1[,6]))
    stat2 = cbind(stat2, qvals(stat2[,6]))
    
    
    sel = which(stat1[,2]>-1)
    Tt = Tt[sel, ]
    
    data.version = '_total_v1'
    save(R, Tt, file=paste('myRdata/my_genes_RNA_seq_analysis_sel', data.version, '.Rdata', sep=''))
    
    
    pdf('myplots/premRNAs_expression_log2.pdf', width=2.2, height=2.2)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    x = stat2[sel,2]
    hist(x,breaks=25,xlab=NA,ylab=NA, main=NA,col='gray',cex.lab=1.0,xlim=c(-5, 15), axes=FALSE)
    abline(v=median(x), lwd=1.5, col='navajowhite4')
    axis(1,at=seq(-5,15,by=3),cex.axis =0.7)
    axis(2,at = seq(0,1000, by=200),las=1,cex.axis = 0.7)
    
    dev.off()
    
    pdf('myplots/mRNAs_expression_log2.pdf', width=2.2, height=2.2)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    x = stat1[sel,2]
    hist(x,breaks=50,xlab=NA,ylab=NA, main=NA,col='gray',cex.lab=1.0,xlim=c(-5, 15), axes=FALSE)
    abline(v=median(x), lwd=1.5, col='navajowhite4')
    axis(1,at=seq(-5,15,by=3),cex.axis =0.7)
    axis(2,at = seq(0,800, by=200),las=1,cex.axis = 0.7)
    
    dev.off()
    
    library(plotrix)
    library("circular")
    make_circ_coord = function(t,x,ttot=24)
    {
        dt=(t[2]-t[1])*.45
        a=(rep(t,rep(4,length(t)))+rep(c(-dt,-dt,dt,dt),length(t)))*2*pi/ttot
        h=rep(x,rep(4,length(x)))*rep(c(0,1,1,0),length(t))
        list(angles=a,heights=h)
    }
    circular_phase24H_histogram = function(x,color_hist = rgb(0.6,0,0.2), cex.axis=0.5, cex.lab=0.5, lwd=0.5)
    {
        #color_DHS = rgb(0.6,0,0.2)
        par(lwd=lwd,cex.axis=cex.axis, cex.main=0.1,cex.lab=cex.lab)
        #par(mfrow=c(1,1),mar=c(4.5,4.5,1,.5)+.1,las=1)
        br=0:24
        h=hist(x, br=br,plot=FALSE)
        co=make_circ_coord(br[-1],h$counts)
        radial.plot(co$heights,co$angles,br[-1]-br[2], clockwise=TRUE,start=pi/2,main=NA, rp.type='p',poly.col=color_hist)
    }
    
    kk = which(stat2[,7]<0.05)
    phases = stat2[kk, 5]
    
    pdf('myplots/Phases_rhythmic_premRNAs.pdf', width=2.2, height=2.2)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    circular_phase24H_histogram(phases, col=rgb(0.0,0.5,0.2), cex.axis=0.7, cex.lab=0.01, lwd=0.5)
    
    dev.off()
    
    kk = which(stat1[,7]<0.05)
    phases = stat1[kk, 5]
    
    pdf('myplots/Phases_rhythmic_RNAs.pdf', width=2.2, height=2.2)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    circular_phase24H_histogram(phases, col=rgb(0.7,0.,0.2), cex.axis=0.7, cex.lab=0.01, lwd=0.5)
    
    dev.off()
    
    ### compare with half-lives
    half.life.comparison = FALSE
    if(half.life.comparison)
    {
        shar = read.table('/Users/jiwang/Degradation_Liver/_x_MicroArray/complementary_data/half-lives/Sharova.txt', stringsAsFactors = FALSE, sep = '\t', header = TRUE)
        fried = read.table('/Users/jiwang/Degradation_Liver/_x_MicroArray/complementary_data/half-lives/Friedel.txt', stringsAsFactors = FALSE, sep = '\t', header = TRUE)
        schwa = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Annotations/Schwarnhausser_2011.csv', header=TRUE)
        ii = match(c("Gene.Names", "mRNA.half.life.replicate..h."), colnames(schwa))
        schwa = schwa[,ii]
        schwa = schwa[which(!is.na(schwa[,2])==TRUE), ]
        gene = c()
        
        for(n in 1:nrow(schwa))
        {
            test = schwa[n, 1]
            test = unlist(strsplit(as.character(test), ';'))
            mm = match(test, Tt$gene)
            if(length(which(!is.na(mm)==TRUE))>0) {gene = c(gene, as.character(test[which(!is.na(mm)==TRUE)][1]))
            }else{
                gene = c(gene, as.character(schwa[n, 1]))
            }
        }
        schwa[, 1] = gene
        
        colnames(shar) = c('gene', 'half.life')
        colnames(fried) = c('gene', 'probeset.ID','hlmin','half.life')
        shar = shar[which(shar[,2]>0),]
        fried = fried[which(fried[,4]>0),]
        fried = fried[, c(1,4)]
        colnames(schwa) = c('gene', 'half.life')
        #save(shar, fried, schwa, file='myRdata/mRNAs_half_lives_databases.Rdata')
        
        overlap2 = intersect(Tt$gene, shar[,1])
        overlap2 = intersect(overlap2, fried[,1])
        overlap2 = intersect(overlap2, schwa[,1] )
        
        ii = match(overlap2, Tt$gene)
        
        pdf('myplots/Compare_half-lives_polyA_Total_Exon_Array_Nascent.pdf', width=16.0, height=6.0)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        par(mfrow=c(1,3))
        
        xx1 = as.numeric(Tt$mean.mRNA[ii]/Tt$mean.premRNA[ii])
        mm = match(overlap2, fried[,1])
        yy1 = as.numeric(fried[mm, 2])
        
        cex = 0.7
        plot(xx1, yy1, cex=cex, xlab='total', ylab='fried', main=paste('R = ', signif(cor(xx1, yy1), d=2), sep=''), log='xy')
        cor(xx1, yy1)
        
        mm = match(overlap2, shar[,1])
        yy1 = as.numeric(shar[mm, 2])
        plot(xx1, yy1, cex=cex, xlab='total', ylab='shar', main=paste('R = ', signif(cor(xx1, yy1), d=2), sep=''), log='xy')
        cor(xx1, yy1)
        
        mm = match(overlap2, schwa[,1])
        yy1 = as.numeric(schwa[mm, 2])
        plot(xx1, yy1, cex=cex, xlab='total', ylab='schwa', main=paste('R = ', signif(cor(xx1, yy1), d=2), sep=''), log='xy')
        cor(xx1, yy1)
        
        dev.off()
    }
}

###########################
######## FITS WITH REAL DATA (test part)
###########################
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
    
    #ptm <- proc.time()
    #param.fits.results = make.fits.with.all.models.for.one.gene(T = T, gene.index = j, debug = TRUE, zt = zt, i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE);
    #proc.time() - ptm
    #c(param.fits.results, my.BIC.this.gene(param.fits.results))
    #test1 = c(param.fits.results, my.model.selection.one.gene.loglike(param.fits.results, method = 'BIC'),
    #  my.model.selection.one.gene.loglike(param.fits.results, method = 'AICc'), my.model.selection.one.gene.loglike(param.fits.results, method = 'AIC'))
    #keep5 = rbind(keep5, test1)
    
  }
  
  keep5 = data.frame(T[mm[c(1:14)], 1], keep5, stringsAsFactors = FALSE)
  keep5[, c(1, grep('best.model', colnames(keep5)))]
  #keep4 = keep4[c(1:4),]
  
  
  test.condition.dependent.alpha = TRUE
  if(test.condition.dependent.alpha)
  {
    Compute.alphas.each.condition = FALSE
    if(Compute.alphas.each.condition)
    {
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
   
   
    #### check functions for mRNA premRNA profiles
    Check.functions = FALSE
    if(Check.functions)
    {
      source('functions.R')
      zt = seq(0, 94, by=2)
      Min = 5
      Amp = 1.0
      beta = 5.0
      phase = 6;
      s1 = compute.s.beta(zt, Min, Amp, phase, beta);
      fold.change = Amp/Min+1;
      #mean = Min*(1+(fold.change-1)/4^beta*gamma(1+2*beta)/gamma(1+beta)^2)
      s2 = compute.s.beta.v1(zt, mean, fold.change, phase, beta)
    
      gamma = log(2)/5;
      eps.gamma = 0.2;
      phase.gamma = 15;
      splicing.k = log(2)/(5/60);
      m1 = compute.m.beta(zt, gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta, simulation.only=FALSE)
      mm1 = compute.m.beta(zt, gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta, simulation.only=TRUE)
    
      source('functions_RPKM_v1.R')
      m2 = compute.m.beta(zt, gamma*(1-eps.gamma), 2*gamma*eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
    }
    
}

##########################################################################################
############################################################## Method validation
##########################################################################################

##################
######## Generate simulated data with real data
##################
my.do.step = TRUE
data.version = '_total_counts_v2'
fake.version = '_total_counts_s5_all'

if(my.do.step)
{
  cat('GENERATE THE FAKE DATA\n')
  
  #data.version = '_total_counts_v2'
  load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_sel_alphas', data.version, '.Rdata', sep=''))
  source("/Users/jiwang/Degradation_Liver/Main_Code_Degradation/f24_modified_1.0.r")
  T = R; Tt = R;
  
  ZT.int = grep('.count.premRNA', colnames(T))
  ZT.ex = grep('.count.mRNA', colnames(T))
  zt = seq(0,94,by = 2)  
  
  #load(file=paste('myRdata/my_genes_RNA_seq_analysis_sel', data.version, '.Rdata', sep=''))
  #source("f24_modified_1.0.r")
  # T = Tt;
  
  #### 500 fake genes for each model
  X = 200
  F = T[1:(4*X),];
  source('functions.R')
  
  for(model in c(1:4))
  {
    cat('model : ',model, '\n')
    # model = 4;
    Fm = generate.fake.read.counts(T = F[((model-1)*X+1):(model*X),], Tt=Tt, X = X, model = model, zt=zt, test.noise.identifiability=FALSE)
    
    if(model == 1)
    {
      FF = Fm
    }else
    {
      FF = rbind(FF,Fm)
    }    
  }
  F = FF;
  
  #### add statistics calculated with RPKM from Deseq2
  #R = T;
  ii = grep('count.mRNA', colnames(F)); jj = grep('count.premRNA', colnames(F))
  exon = as.matrix(F[,ii]);intron = as.matrix(F[,jj]);
  
  data1 = matrix(NA, ncol=ncol(exon), nrow=nrow(exon))
  data2 = matrix(NA, ncol=ncol(intron), nrow=nrow(intron))
  for(n in 1:nrow(exon))
  {
    data1[n, ] = norm.RPKM(exon[n, ], F$length.mRNA[n]);
    data2[n, ] = norm.RPKM(intron[n, ], F$length.premRNA[n])
  }
  source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')
  stat1 = t(apply(data1,1, f24_R2_alt2, t=c(0:47)*2))
  stat2 = t(apply(data2,1, f24_R2_alt2, t=c(0:47)*2))
  stat1 = cbind(stat1, qvals(stat1[,6]))
  stat2 = cbind(stat2, qvals(stat2[,6]))
  
  F[, grep('rpkm.mRNA', colnames(F))] = stat1
  F[, grep('rpkm.premRNA', colnames(F))] = stat2
  
  F$true.model = ceiling(c(1:nrow(F))/200)
  
  ### estimate dispersion parameter with Deseq2
  Compute.alphas.each.condition = TRUE
  if(Compute.alphas.each.condition)
  {
    ss1 = apply(exon, 1, sum);  ss2 = apply(intron, 1, sum)
    kk = which(ss1<10^10 & ss2<10^10);
    F = F[kk, ];
    countData =  rbind(exon[kk,], intron[kk, ])
    require('DESeq2')
    #size.factors = estimateSizeFactorsForMatrix(countData)
    load(file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/model_RNA_seq_total/size_factors_48_samples.Rdata')
    
    ### estimate dispersion parameters with Deseq2 for each gene and each condition
    alphas = matrix(NA, nrow=nrow(countData), ncol=12)
    #alphas.genes = matrix(NA, nrow=nrow(countData), ncol=12)
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
      #alphas.genes[,n] = mcols(dds)$dispGeneEst
    }
    xx = data.frame(alphas[c(1:nrow(F)),], alphas[-c(1:nrow(F)), ])
    yy = F[, grep('alpha', colnames(F))];
    F[, grep('alpha', colnames(F))] = xx;
    colnames(yy) = c(paste('alpha.mRNA.ZT', c(0:11)*2, '.true', sep=''),  paste('alpha.premRNA.ZT', c(0:11)*2, '.true', sep=''))
    #alphas = xx;
    #colnames(alphas) = c(paste('a.m.ZT', c(0:11)*2, '.estimated', sep=''),  paste('a.s.ZT', c(0:11)*2, '.estimated', sep=''))
    
    F = data.frame(F, yy, stringsAsFactors = FALSE)
    
    Check.alpha.estimation = FALSE
    if(Check.alpha.estimation){
      par(mfrow=c(1,2))
      lims = range(c(F$alpha.mRNA.ZT0, F$alpha.mRNA.ZT0.true));
      plot(F$alpha.mRNA.ZT0, F$alpha.mRNA.ZT0.true, log='xy', xlim=lims, ylim=lims, cex=0.7, main='Alpha for mRNA at ZT0', xlab='true values', 
           ylab='estimated values');
      abline(0, 1, lwd=2.0, col='red'); abline(log(2), 1, col='red'); abline(log(0.5), 1, col='red')
      
      plot(T$alpha.premRNA.ZT0, T$a.s.ZT0estimated, log='xy', xlim=c(0.001, 2), ylim=c(0.001, 2), cex=0.7, main = 'Alpha for premRNA at ZT0', xlab='true values', 
           ylab='estimated values');
      abline(0, 1, lwd=2.0, col='red'); abline(log(2), 1, col='red'); abline(log(0.5), 1, col='red');
    }
  }
  
  #load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_simulated_data', fake.version, '_without_noise.Rdata', sep=''))
  #F0$true.model = ceiling(c(1:nrow(F0))/200)
  #mm = match(F[,1], F0[,1])
  #F0 = F0[mm, ]
  save(F, file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_simulated_data', fake.version, '.Rdata', sep=''))
  
  ### test the model selection and parameter estimation
  Test.Fitting = FALSE
  if(Test.Fitting)
  {
    cat('TEST THE FIT on THE FAKE DATA\n')
    #fake.version = '_total_counts_s1'
    fake.version = '_total_counts_s5_all'
    load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_simulated_data', fake.version, '.Rdata', sep=''))
    
    #load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_simulated_data', fake.version, '.Rdata', sep=''))
    #source("f24_modified_1.0.r")
    T = F;
    
    #source('functions.R')
    #source('functions_read_counts_gamma_amp_gamma_v2.R')
    #best.model = model.check
    ZT.ex = grep('.count.mRNA', colnames(T))
    ZT.int = grep('.count.premRNA', colnames(T))
    zt = seq(0,94,by = 2)
    
    test.model = 4;
    w = 2*pi/24
    jj = which(T$true.model==test.model);
    j = jj[5]
    source('functions.R')
    #T$mRNA.outlier[j] = paste(outlier.m, sep='', collapse = ';')
    #T$premRNA.outlier[j] = paste(outlier.s, sep='', collapse = ';')
    ptm <- proc.time()
    param.fits.results = make.fits.with.all.models.for.one.gene.remove.outliers(T = T, gene.index = j, debug = TRUE,
                                                                                zt = zt, i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE);
    proc.time() - ptm
    
    index = match(c('outlier.m', 'outlier.s'), names(param.fits.results))
    res.fit = as.numeric(param.fits.results[-index])
    names(res.fit) = names(param.fits.results)[-index]
    #missing.data = length(unlist(strsplit(param.fits.results[index[1]], ';'))) + length(unlist(strsplit(param.fits.results[index[2]], ';')))
    test1 = c(param.fits.results, my.model.selection.one.gene.loglike(res.fit, nb.data = (96), method = 'BIC'))
    print(test1)
    print(T[j, c(4:27, 142, 146, 147, 141, 140, 143:145)])
    
    ### fit simulated data with weak noise
    T = F0;
    ptm <- proc.time()
    param.fits.results = make.fits.with.all.models.for.one.gene.remove.outliers(T = T, gene.index = j, debug = TRUE,
                                                                                zt = zt, i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE);
    proc.time() - ptm
    
    index = match(c('outlier.m', 'outlier.s'), names(param.fits.results))
    res.fit = as.numeric(param.fits.results[-index])
    names(res.fit) = names(param.fits.results)[-index]
    #missing.data = length(unlist(strsplit(param.fits.results[index[1]], ';'))) + length(unlist(strsplit(param.fits.results[index[2]], ';')))
    test0 = c(param.fits.results, my.model.selection.one.gene.loglike(res.fit, nb.data = (96), method = 'BIC'))
    
    print(test1)
    print(test0)
    print(F[j, c(142, 146, 147, 141, 140, 143:145)])
    print(T[j, c(142, 146, 147, 141, 140, 143:145)])
    
    #### test noise effects
    test.noise.identifiability = FALSE
    if(test.noise.identifiability)
    {
      load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_simulated_data', fake.version, '_without_noise.Rdata', sep=''))
      T = F;
      
      T[, grep('alpha.mRNA.ZT', colnames(T))] = 0.0001;
      T[, grep('alpha.premRNA.ZT', colnames(T))] = 0.0001;
      
      ZT.int = grep('.count.premRNA', colnames(T))
      ZT.ex = grep('.count.mRNA', colnames(T))
      zt = seq(0,94,by = 2)
      
      kk = c(401:600)
      model = 3
      #gg = 'Per2'
      j = kk[50]
      source('functions.R')
      #T$mRNA.outlier[j] = paste(outlier.m, sep='', collapse = ';')
      #T$premRNA.outlier[j] = paste(outlier.s, sep='', collapse = ';')
      ptm <- proc.time()
      param.fits.results = make.fits.with.all.models.for.one.gene.remove.outliers(T = T, gene.index = j, debug = TRUE,
                                                                                  zt = zt, i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE);
      proc.time() - ptm
      
    }
    #c(param.fits.results, my.BIC.this.gene(param.fits.results))
    source('functions.R')
    c(param.fits.results, my.AIC.this.gene.loglike(param.fits.results),
      my.AIC.this.gene.loglike(param.fits.results, correction=TRUE),
      my.BIC.this.gene.loglike(param.fits.results))
    
    i = j;
    zt = seq(0, 94, by=2)
    m0 = compute.m.beta(zt, T$gamma[i], T$eps.gamma[i], T$phase.gamma[i], T$splicing.k[i], T$Min.int[i], T$Amp.int[i], T$phase.int[i], T$beta.int[i])
    s0 = compute.s.beta(zt, T$Min.int[i], T$Amp.int[i], T$phase.int[i], T$beta.int[i]);
    read.m0 = convert.nb.reads(m0, T[j, 2]);
    read.s0 = convert.nb.reads(s0, T[j, 3]);
    
    R.m=T[j, ZT.ex]; R.s=T[j, ZT.int]; alpha.m=T[j, 4]; alpha.s=T[j, 5]; mu.m=read.m0; mu.s=read.s0;
    cat(NB.error(R.m=T[j, ZT.ex], R.s=T[j, ZT.int], alpha.m=T[j, 4], alpha.s=T[j, 5], read.m0, read.s0), '\n');
    cat(-2*sum(R.m*log(mu.m*alpha.m/(1.0+mu.m*alpha.m)) - 1.0/alpha.m*log(1.0+alpha.m*mu.m)
               + lgamma(R.m+1.0/alpha.m) -lgamma(1.0/alpha.m) - lfactorial(R.m)), '\n');    
    cat(-2*sum(R.s*log(mu.s*alpha.s/(1.0+mu.s*alpha.s)) - 1.0/alpha.s*log(1.0+alpha.s*mu.s)
               + lgamma(R.s+1.0/alpha.s) -lgamma(1.0/alpha.s) - lfactorial(R.s)), '\n')
    
    gamma = param.fits.results[27];eps.gamma = param.fits.results[28];phase.gamma = param.fits.results[29];splicing.k =  param.fits.results[30];
    Min = param.fits.results[31]; Amp=param.fits.results[32]; phase=param.fits.results[33]; beta=param.fits.results[34];
    m4 = compute.m.beta(zt, gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
    s4 = compute.s.beta(zt, Min, Amp, phase, beta);
    read.m4 = convert.nb.reads(m4, as.numeric(T[j,2]))
    read.s4 = convert.nb.reads(s4, T$length.premRNA[i])
    R.m=T[j, ZT.ex]; R.s=T[j, ZT.int]; alpha.m=T[j, 4]; alpha.s=T[j, 5]; mu.m=read.m4; mu.s=read.s4;
    #cat(NB.error(R.m=T[j, ZT.ex], R.s=T[j, ZT.int], alpha.m=T[j, 4], alpha.s=T[j, 5], read.m0, read.s0), '\n');
    cat(-2*sum(R.m*log(mu.m*alpha.m/(1.0+mu.m*alpha.m)) - 1.0/alpha.m*log(1.0+alpha.m*mu.m)
               + lgamma(R.m+1.0/alpha.m) -lgamma(1.0/alpha.m) - lfactorial(R.m)), '\n');    
    cat(-2*sum(R.s*log(mu.s*alpha.s/(1.0+mu.s*alpha.s)) - 1.0/alpha.s*log(1.0+alpha.s*mu.s)
               + lgamma(R.s+1.0/alpha.s) -lgamma(1.0/alpha.s) - lfactorial(R.s)), '\n')
    
    ### check fitted temporal profiles
    i = j
    zt = seq(0, 94, by=2)
    source('functions.R')
    zt.p = seq(0, 94, by=0.2)
    
    model = 0;
    #gamma = T$gamma[i];eps.gamma = 0;phase.gamma = 12;splicing.k =  param.fits.results[4];
    #Min = param.fits.results[5]; Amp=param.fits.results[6]; phase=param.fits.results[7]; beta=param.fits.results[8];
    m0 = compute.m.beta(zt.p, T$gamma[i], T$eps.gamma[i], T$phase.gamma[i], T$splicing.k[i], T$Min.int[i], T$Amp.int[i], T$phase.int[i], T$beta.int[i])
    s0 = compute.s.beta(zt.p, T$Min.int[i], T$Amp.int[i], T$phase.int[i], T$beta.int[i]);
    read.s0 = s0;
    read.m0 = m0;
    
    model = 4;
    gamma = param.fits.results[27];eps.gamma = param.fits.results[28];phase.gamma = param.fits.results[29];splicing.k =  param.fits.results[30];
    Min = param.fits.results[31]; Amp=param.fits.results[32]; phase=param.fits.results[33]; beta=param.fits.results[34];
    m4 = compute.m.beta(zt.p, gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
    s4 = compute.s.beta(zt.p, Min, Amp, phase, beta);
    #read.m4 = convert.nb.reads(m1, R$length.mRNA[i])
    #read.s4 = convert.nb.reads(s1, R$length.premRNA[i])
    read.s4 = s4;
    read.m4 = m4;
    
    model = 3;
    gamma = param.fits.results[16];eps.gamma = param.fits.results[17];phase.gamma = param.fits.results[18];splicing.k =  param.fits.results[19];
    Min = param.fits.results[20]; Amp=0; phase=12; beta=1;
    m3 = compute.m.beta(zt.p, gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
    s3 = compute.s.beta(zt.p, Min, Amp, phase, beta);
    read.s3 = s3;
    read.m3 = m3;
    
    model = 2;
    gamma = param.fits.results[3];eps.gamma = 0;phase.gamma = 12;splicing.k =  param.fits.results[4];
    Min = param.fits.results[5]; Amp=param.fits.results[6]; phase=param.fits.results[7]; beta=param.fits.results[8];
    m2 = compute.m.beta(zt.p, gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
    s2 = compute.s.beta(zt.p, Min, Amp, phase, beta);
    #read.m2 = convert.nb.reads(m2, R$length.mRNA[i])
    #read.s2 = convert.nb.reads(s2, R$length.premRNA[i])
    read.s2 = s2;
    read.m2 = m2;
    
    model = 0;
    #gamma = T$gamma[i];eps.gamma = 0;phase.gamma = 12;splicing.k =  param.fits.results[4];
    #Min = param.fits.results[5]; Amp=param.fits.results[6]; phase=param.fits.results[7]; beta=param.fits.results[8];
    m0 = compute.m.beta(zt.p, T$gamma[i], T$eps.gamma[i], T$phase.gamma[i], T$splicing.k[i], T$Min.int[i], T$Amp.int[i], T$phase.int[i], T$beta.int[i])
    s0 = compute.s.beta(zt.p, T$Min.int[i], T$Amp.int[i], T$phase.int[i], T$beta.int[i]);
    read.s0 = s0;
    read.m0 = m0;
    
    L.m = T$length.mRNA[i];
    L.s = T$length.premRNA[i];
    M = norm.RPKM(as.numeric(T[i,ZT.ex]), L.m)
    S = norm.RPKM(as.numeric(T[i,ZT.int]), L.s)
    
    par(mfcol=c(1,3),cex=1.0)
    
    lims = range(S, read.s4, read.s2)
    plot(zt, as.numeric(S), col='darkblue',type='b',ylim=lims, lty=1, main=paste(T$gene[i], ', premRNA', sep=''), log='y', ylab='rpkm')
    points(zt.p, read.s4, col='magenta', type='l',lty=2, lwd=2.0)
    points(zt.p, read.s3, col='red', type='l',lty=2, lwd=2.0)
    points(zt.p, read.s2, col='green', type='l',lty=2, lwd=2.0)
    points(zt.p, read.s0, col='blue', type='l',lty=1, lwd=1.5)
    #abline(h=1.0,col='gray',lwd=2.0)
    
    lims = range(M, read.m4, read.m2)
    plot(zt, as.numeric(M), col='darkblue',type='b',lty=1, ylim=lims, main=paste(T$gene[i], ', mRNA', sep=''), log='y', ylab='rpkm')
    points(zt.p, read.m4, col='magenta', type='l',lty=2, lwd=2.0)
    points(zt.p, read.m3, col='red', type='l',lty=2, lwd=2.0)
    points(zt.p, read.m2, col='green', type='l',lty=2, lwd=2.0)        
    points(zt.p, read.m0, col='blue', type='l',lty=1, lwd=1.5)
    
    M = M/mean(M)
    S = S/mean(S)
    lims = range(S, M)
    plot(zt, S, col='green', type='b', ylim=lims, main=paste(T$gene[i], ', premRNA & mRNA', sep=''), log='', ylab='relative rpkm')
    points(zt, M, col='blue', type='b',lty=1)
    abline(h=1, col='darkgray', lwd=2.0)
    
  }
  #####################
  ### control the temporal profiles
  ####################
  Control.profiles = FALSE
  if(Control.profiles)
  {
    #pdf('myplots/Total/simulated_data_control.pdf', width=6, height=6)
    #par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,3,0.8)+0.1, tcl = -0.3)
    
    F$true.model = ceiling(c(1:nrow(F))/500)
    
    par(mfrow=c(2,2))
    lims = range(c(F$alpha.mRNA, F$alpha.mRNA.estimated))
    plot(F$alpha.mRNA, F$alpha.mRNA.estimated, ylim=lims, xlim=lims, type='p', cex=0.7, log='xy',xlab='true value', ylab='estimated value', 
         main='dispersion of mRNAs')
    abline(0, 1, col='red', lwd=2.0)
    
    lims = range(c(F$alpha.premRNA, F$alpha.premRNA.estimated))
    plot(F$alpha.premRNA, F$alpha.premRNA.estimated, ylim=lims, xlim=lims, type='p', log='xy', cex=0.7, main='dispersion of premRNAs')
    abline(0, 1, col='red', lwd=2.0)
    
    #boxplot(log10(F$alpha.mRNA.estimated/F$alpha.mRNA) ~ F$true.model)
    #boxplot(log10(F$alpha.premRNA.estimated/F$alpha.premRNA) ~ F$true.model)
    
    mm = which(F$pval.rpkm.mRNA<0.001 & F$rel.amp.rpkm.mRNA>0.05)
    plot(F$rel.amp.rpkm.mRNA[mm], F$alpha.mRNA.estimated[mm]/F$alpha.mRNA[mm], log='xy', cex=0.7)
    abline(h=1, col='red', lwd=2.0)
    mm = which(F$pval.rpkm.premRNA<0.001 & F$rel.amp.rpkm.premRNA>0.05)
    plot(F$rel.amp.rpkm.premRNA[mm], F$alpha.premRNA.estimated[mm]/F$alpha.premRNA[mm],  log='xy', cex=0.7)
    abline(h=1, col='red', lwd=2.0)
    
    kk = c(101:400)
    model = c(rep(2, 100), rep(3, 100), rep(4, 100))
    boxplot(-log10(F$pval.mRNA[kk]) ~ model, col=c('green', 'red', 'gray'), xlab='Model', ylab='-log10(pval)')
    boxplot((F$amp.mRNA[kk]) ~ model, col=c('green', 'red', 'gray'), xlab='Model', ylab='log2 fold change')
    
    dev.off()
    
  }
  
  dispersion.Deseq.vs.edgeR = FALSE
  if(dispersion.Deseq.vs.edgeR)
  {
    fake.version = '_total_counts_s2'
    load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_simulated_data', fake.version, '.Rdata', sep=''))
    
    #### check the estimation of dispersion for individual examples
    kk = which(F$qv.rpkm.mRNA<0.01 & F$qv.rpkm.premRNA<0.01 & F$rel.amp.rpkm.mRNA>0.5 & F$rel.amp.rpkm.premRNA>0.5)
    #examples = c('Npas2', 'Arntl', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry1', 'Cry2', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2', 
    #'Dbp', 'Hlf', 'Tef', 'Nfil3', 'Bhlhe40', 'Bhlhe41', 'Nampt', 'Parp1', 'Prkaa2', 'Prkag1', 'Gsk3a', 'Gsk3b', 'Csnk1d')
    #mm = match(examples, R[,1])
    mm = kk[order(-F[kk, 5])]
    mm = kk[which(!is.na(mm)==TRUE)]
    #T[mm, c(1:5)]
    source('functions.R')
    plot.dispersion.examples(T=F, index=mm)
    
    
    ii = grep('count.mRNA', colnames(F))
    jj = grep('count.premRNA', colnames(F))
    exon = as.matrix(F[,ii])
    intron = as.matrix(F[,jj])
    countData =  rbind(exon, intron)
    
    ### estimate dispersion from Deseq2
    require('DESeq2')
    #countData =  rbind(exon, intron)
    #size.factors = estimateSizeFactorsForMatrix(countData)
    #condition <- factor(rep(c(1:12), 4))
    kk = c(0:3)*12 +1
    
    countData = exon[, kk]
    condition <- factor(rep(c(1, 2), 2));
    dds <- DESeqDataSetFromMatrix(countData, colData=DataFrame(condition), design = ~ 1)
    
    #load(file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/model_RNA_seq_total/size_factors_48_samples.Rdata')
    dds <- estimateSizeFactors(dds)
    #sizeFactors(dds) = size.factors;
    #sizeFactors(dds)
    dds <- estimateDispersions(dds)
    xx1 = mcols(dds)$dispersion
    xx = xx1
    #FF = cbind(F, xx[c(1:nrow(F))], xx[-c(1:nrow(F))])
    #colnames(FF)[c(124:125)] = c('alpha.mRNA.estimated', 'alpha.premRNA.estimated')
    F$alpha.mRNA.estimated = xx[c(1:nrow(F))] 
    F$alpha.premRNA.estimated = xx[-c(1:nrow(F))] 
    
    
    ### estimate dispersion with edgeR
    require('edgeR')
    kk = c(0:3)*12 +1
    countData = exon[,kk]
    
    group <- factor(rep(c(1), 4))
    #design <- model.matrix(~group)
    norm.factors = calcNormFactors(countData)
    
    d <- DGEList(counts=countData, norm.factors = norm.factors, group=group)
    #d <- DGEList(counts=intron, norm.factors = norm.factors, group=group)
    #d1 <- estimateDisp(d)
    d1 <- estimateDisp(d, design=NULL, trend.method="locfit")
    
    #d <- DGEList(counts=exon, norm.factors = norm.factors, group=group)
    #d1 <- estimateDisp(d)
    #d2 <- estimateDisp(d, design, trend.method="locfit")
    xx1 = d1$tagwise.dispersion
    #d2 = estimateDisp(d, design, trend.method="locfit", robust=TRUE)
    #xx2 = d2$tagwise.dispersion
    #d3 <- estimateGLMTrendedDisp(d, design, min.n=10)
    #d3 <- estimateGLMTagwiseDisp(d3, design)
    #xx3 = d3$tagwise.dispersion
    #d4 = estimateGLMRobustDisp(d, design)
    #xx4 = d4$tagwise.dispersion
    xx = xx1
    #FF = cbind(F, xx[c(1:nrow(F))], xx[-c(1:nrow(F))])
    #colnames(FF)[c(124:125)] = c('alpha.mRNA.estimated', 'alpha.premRNA.estimated')
    F$alpha.mRNA.estimated = xx[c(1:nrow(F))] 
    F$alpha.premRNA.estimated = xx[-c(1:nrow(F))] 
  }
    
}

##################
######## Results of fitting simulated data
##################
my.do.step = FALSE
if(my.do.step)
{
  Fitting.results  = FALSE
  if(Fitting.results)
  {
    #### Model selection for simulated data
    norm.parameters.v2 = FALSE
    if(norm.parameters.v2)
    {
      xx = read.table('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/Vital-IT_fits_optim/all_genes_fits_results_with_header_RNA_seq_beta_total_counts_s1_fitting_norm_parameters_v2.txt',
                      sep='\t', header=FALSE)
      fake.version = '_total_counts_s1'
      load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_simulated_data', fake.version, '.Rdata', sep=''));
      gene_list = F$gene;
      
      m = match(gene_list, xx[,1])
      xx = xx[m,]
      
      names = colnames(T)[-c(1:172)]
      names = c('gene', names[1:2], 'rank.hess.m2', names[3:15], 'rank.hess.m3', names[16:26], 'rank.hess.m4', names[27:44])
      colnames(xx) = names
      #kk = match(gene_list, T$gene)
      #Tt = Tt[kk,]
      res.version = '_total_counts_s1_fitting_norm_parameters_v2';
      T = data.frame(F, xx[,-1], stringsAsFactors = FALSE)	
      save(T, file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos',
                           res.version, '.Rdata',sep = ''))
      
    }
    #res.version = '_total_counts_s1_fitting_v1';
    #res.version = '_total_counts_s1_fitting_norm_parameters_v2';
    res.version = '_total_counts_s1_fitting_norm_parameters_constrained_gamma_v4';
    load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
    source('functions.R')
    
    BIC = my.model.selection.loglike(T=T, method='BIC', outliers = TRUE)
    AIC = my.model.selection.loglike(T=T, method='AIC', outliers = TRUE)
    AICc = my.model.selection.loglike(T=T, method='AICc', outliers = TRUE)
    T = data.frame(T, AIC, AICc, BIC, stringsAsFactors = FALSE);
    
    Fitting.Test = FALSE
    if(Fitting.Test)
    {
      missing = T$gene[which(is.na(T$BIC.best.model))] 
      #save(missing, file = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos_norm_params_v3.Rdata')
      j = which(T$gene=='fake_m4_157')
      
      source('functions.R')
      ZT.ex = grep('.count.mRNA', colnames(T))
      ZT.int = grep('.count.premRNA', colnames(T))
      zt = seq(0,94,by = 2)
      ptm <- proc.time()
      param.fits.results = make.fits.with.all.models.for.one.gene.remove.outliers(T = T, gene.index = j, debug = TRUE,
                                                                                  zt = zt, i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE);
      proc.time() - ptm
      index = match(c('outlier.m', 'outlier.s'), names(param.fits.results))
      res.fit = as.numeric(param.fits.results[-index])
      names(res.fit) = names(param.fits.results)[-index]
      #missing.data = length(unlist(strsplit(param.fits.results[index[1]], ';'))) + length(unlist(strsplit(param.fits.results[index[2]], ';')))
      test0 = c(res.fit, my.model.selection.one.gene.loglike(res.fit, nb.data = (96), method = 'BIC'))
      #print(test1)
      print(test0)
      print(T[j, c(142, 146, 147, 141, 140, 143:145)])
    }
    
    T = T[which(!is.na(T$BIC.best.model)==TRUE), ]
    
    Tt = T
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
    }
    
    #### Plot selected mdoels and estimated parameters
    Plot.Results = TRUE
    if(Plot.Results)
    {
      source('functions.R')
      library(stringr)
      #selecting.approach = c('BIC.best.model','AICc.best.model','AIC.best.model')
      selecting.approach = c('BIC.best.model')
      #sd.noise = c(0.1, 0.25, 0.5)
      #set.global.sigma();
      #filter = TRUE
      source('functions.R')
      pdf.name = paste(file = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Total/Total_counts/my_simulated_data_with_fits_results_BIC_AIC_LRT', res.version, '.pdf', sep='')
      pdf(pdf.name, width = 10, height = 7.5)
      par(mfrow = c(3,4), cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
      
      for(select.method in selecting.approach)
      {
        #select.method = 'BIC.best.model'
        #T$true.model = ceiling(c(1:400)/100)
        eval(parse(text = paste('T$best.model=T$', select.method, sep='')))
        kk = which(!is.na(T$best.model)==TRUE)
        Tt = T[kk,]
        
        compare.parameters.RNA.seq.total.fake.by.model(Tt, select.method)
      }
      
      dev.off()
    }
    
  }
  
  ##################
  ################## Learn identifiablity from sythetic data fitting
  ##################
  Search.robust.parameter.combination = FALSE
  if(Search.robust.parameter.combination)
  {
    ###### start with Model 3
    kk = which(T$true.model==3 & T$BIC.best.model==3)
    gamma.true = c(T$gamma[kk])
    gamma.estimate = c(T$gamma.m3[kk])
    eps.true = c(T$eps.gamma[kk])
    eps.estimate = c(T$eps.gamma.m3[kk])
    #eps.true = T$amp.gamma[kk]/2;
    #eps.estimate = T$amp.gamma.m3[kk]/2
    ratio1 = gamma.estimate/gamma.true;
    ratio2 = eps.estimate/eps.true;
    
    yy = cbind(gamma.true, gamma.estimate, ratio1, eps.true, eps.estimate, ratio2)
    #yy = yy[order(-yy[,2]), ]
    jj1 = which(ratio1>5 | ratio1<0.2)
    #jj2 = which(ratio2>2. | ratio2<0.5)
    #jj3 = intersect(jj1, jj2)
    #jj1 = setdiff(jj1, jj3)
    #jj2 = setdiff(jj2, jj3)
    
    plot(yy[, 4:5], log='');abline(0, 1, col='red');abline(0, 2, col='red');abline(0, 0.5, col='red');
    
    plot(log(2)/yy[, 1:2], log='xy');abline(0, 1, col='red');abline(log(2), 1, col='red');abline(log(1/2), 1, col='red');
    plot(w/gamma.true/(2*pi/log(2)), w/gamma.estimate/(2*pi/log(2)), log='xy');
    abline(0, 1, col='red');abline(log(2), 1, col='red');abline(log(1/2), 1, col='red');
    plot(atan2(w, gamma.true), atan2(w,gamma.estimate), log='xy');abline(0, 1, col='red');abline(log(2), 1, col='red');abline(log(1/2), 1, col='red');
    #plot(1/sqrt(1+w^2/gamma.true^2), 1/sqrt(1+w^2/gamma.estimate^2),  log='xy');abline(0, 1, col='red');abline(log(2), 1, col='red');abline(log(1/2), 1, col='red');
    
    plot(eps.true/sqrt(1+w^2/gamma.true^2), eps.estimate/sqrt(1+w^2/gamma.estimate^2),  log='xy');abline(0, 1, col='red');abline(log(2), 1, col='red');abline(log(1/2), 1, col='red');
    
    test = c(jj3, jj1, jj2)
    plot(gamma.true[test]*eps.true[test], gamma.estimate[test]*eps.estimate[test],  log='xy');abline(0, 1, col='red');abline(log(2), 1, col='red');abline(log(1/2), 1, col='red');
    
    ######################################################
    ######## good combination to further consider
    ######################################################
    #plot(gamma.true*eps.true, gamma.estimate*eps.estimate,  log='xy');abline(0, 1, col='red');abline(log(2), 1, col='red');abline(log(1/2), 1, col='red');
    plot((eps.true/sqrt(1+w^2/gamma.true^2)), (eps.estimate/sqrt(1+w^2/gamma.estimate^2)),  cex=0.7,log='');abline(0, 1, col='red');abline(log(2), 1, col='red');abline(log(1/2), 1, col='red');
    plot((phase.true+atan2(w, gamma.true)/w), (phase.estimate+atan2(w,gamma.estimate)/w), log='', cex=0.7);
    abline(0, 1, col='red'); abline(24, 1, col='red');abline(-24, 1, col='red');
    
    ####### test 
    w = 2*pi/24
    kk0 = which(T$true.model==2)
    kk1 = which(T$true.model==3)
    #kk1 = c()
    kk2 = which(T$true.model==4)
    #kk2 = c()
    #kk = which(T$true.model==3 & T$BIC.best.model==3)
    gamma.true = c(T$gamma[kk1], T$gamma[kk2])
    gamma.estimate = c(T$gamma.m3[kk1], T$gamma.m4[kk2])
    eps.true = c(T$eps.gamma[kk1], T$eps.gamma[kk2])
    eps.estimate = c(T$eps.gamma.m3[kk1], T$eps.gamma.m4[kk2])
    phase.true = c(T$phase.gamma[kk1], T$phase.gamma[kk2])
    phase.estimate = c(T$phase.gamma.m3[kk1], T$phase.gamma.m4[kk2])
    ratio.true = c(T$splicing.k[kk0]/T$gamma[kk0], T$splicing.k[kk1]/T$gamma[kk1], T$splicing.k[kk2]/T$gamma[kk2])
    ratio.estimate = c(T$splicing.k.m2[kk0]/T$gamma.m2[kk0], T$splicing.k.m3[kk1]/T$gamma.m3[kk1], T$splicing.k.m4[kk2]/T$gamma.m4[kk2])
    
    #pdf.name = paste("/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Total/Total_counts/new_parameter_combination_Model_4.pdf", sep='')
    #pdf(pdf.name, width=10, height=4.0)
    #par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3, pty='s')
    
    par(mfrow=c(1,3),cex=1.0)
    plot(gamma.true, gamma.estimate, log='xy', cex=0.7);abline(0, 1, col='red', cex=0.7);
    abline(log(2), 1, col='red');abline(log(1/2), 1, col='red');
    abline(0, 1, col='red');abline(log(2), 1, col='red');abline(log(1/2), 1, col='red');
    plot((eps.true/sqrt(1+w^2/gamma.true^2)), (eps.estimate/sqrt(1+w^2/gamma.estimate^2)), cex=0.7, log='xy');
    abline(0, 1, col='red');abline(v=0.1, col='red', lwd=2.0)
    abline(log(2), 1, col='red');abline(log(1/2), 1, col='red');
    #plot(phase.true, phase.estimate, log='', cex=0.7);abline(0, 1, col='red'); abline(6, 1, col='red');abline(-6, 1, col='red');
    
    plot((phase.true+atan2(w, gamma.true)/w), (phase.estimate+atan2(w,gamma.estimate)/w), log='', cex=0.7);
    abline(0, 1, col='red'); abline(24, 1, col='red');abline(-24, 1, col='red');
    
    #dev.off()
    #plot(atan2(w, gamma.true), atan2(w,gamma.estimate), log='xy');abline(0, 1, col='red');abline(log(2), 1, col='red');abline(log(1/2), 1, col='red');
    #plot((phase.true+atan2(gamma.true, w))[-c(1:length(kk1))], (phase.estimate+atan2(gamma.estimate, w))[-c(1:length(kk1))], log='xy', cex=0.7);abline(0, 1, col='red');abline(log(2), 1, col='red');abline(log(1/2), 1, col='red');
    
    plot(ratio.true, ratio.estimate, log='xy');abline(0, 1, col='red', lwd=2.0)
    abline(log(2), 1, col='red');abline(log(1/2), 1, col='red');
  }
  
  ##################
  ################## Identifiablity analysis for Model 3
  ##################
  Check.identifiability = FALSE
  if(Check.identifiability)
  {
    ###### Check Model 3 and try to find conditions in which estimation is vulnerable to non-identifiability
    kk = which(T$true.model==3 & T$BIC.best.model==3)
    gamma.true = c(T$gamma[kk])
    gamma.estimate = c(T$gamma.m3[kk])
    eps.true = c(T$eps.gamma[kk])
    eps.estimate = c(T$eps.gamma.m3[kk])
    #eps.true = T$amp.gamma[kk]/2;
    #eps.estimate = T$amp.gamma.m3[kk]/2
    ratio1 = gamma.estimate/gamma.true;
    ratio2 = eps.estimate/eps.true;
    
    yy = cbind(gamma.true, gamma.estimate, ratio1, eps.true, eps.estimate, ratio2)
    #yy = yy[order(-yy[,2]), ]
    jj1 = which(ratio1>5); jj2 = which(ratio1<0.2);jj3 = which(ratio1<=5 & ratio1>=0.2);
    #jj2 = which(ratio2>2. | ratio2<0.5)
    pdf.name = paste("/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Total/Total_counts/Parameter_non_identifiable_M3_norm_parameters.pdf", 
                     sep='')
    pdf(pdf.name, width=6, height=4.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    par(mfrow=c(1,1),cex=1.0)
    plot(gamma.true[jj3], eps.true[jj3], log='x', cex=0.8, xlim = range(c(gamma.true, gamma.estimate)), xlab='true gamma', ylab='true eps.norm')
    points(gamma.true[jj1], eps.true[jj1], col='red', pch=16, cex=1.)
    points(gamma.estimate[jj1], eps.estimate[jj1], col='red', pch=1, cex=0.8)
    arrows(gamma.true[jj1], eps.true[jj1], gamma.estimate[jj1], eps.estimate[jj1], col='gray', length = 0.1, lwd=0.5)
    points(gamma.true[jj2], eps.true[jj2], col='cadetblue4', pch=16, cex=1.)
    points(gamma.estimate[jj2], eps.estimate[jj2], col='cadetblue4', pch=1, cex=0.8)
    arrows(gamma.true[jj2], eps.true[jj2], gamma.estimate[jj2], eps.estimate[jj2], col='gray', length = 0.1, lwd=0.5)
    
    w = 2*pi/24; 
    #abline(v=w); abline(v=2*w)
    rr = lseq(log(2)/24, log(2)/(10/60), length=100)
    points(rr, rr/sqrt(w^2 + rr^2), col='orange', type='l', lwd=1.5, lty=2);
    abline(h=0.4, col='darkorange', lwd=1.5, lty=2)
    
    dev.off()
    
    aa = matrix(NA, nrow=max(c(length(jj1), length(jj2), length(jj3)))*12, ncol=3)
    aa[1:(length(jj1)*12), 1] = -log10(as.numeric(unlist(T[kk[jj1], c(4:15)])))
    aa[1:(length(jj2)*12), 2] = -log10(as.numeric(unlist(T[kk[jj2], c(4:15)])))
    aa[1:(length(jj3)*12), 3] = -log10(as.numeric(unlist(T[kk[jj3], c(4:15)])))
   
    boxplot(aa)
    
    mean(as.matrix(T[kk[jj1], c(4:15)]))
    mean(as.matrix(T[kk[jj2], c(4:15)]))
    mean(as.matrix(T[kk[jj3], c(4:15)]))
    
    #####################
    ###### Identifiability analysis
    #####################
    jj1 = which(ratio1>5); jj2 = which(ratio1<0.2);jj3 = which(ratio1<=5 & ratio1>=0.2);
    
    j = kk[jj1[1]]
    Identifiablity.Analysis.by.Profile.Likelihood = TRUE
    if(Identifiablity.Analysis.by.Profile.Likelihood)
    {
      source('functions.R')
      ZT.ex = grep('.count.mRNA', colnames(T))
      ZT.int = grep('.count.premRNA', colnames(T))
      zt = seq(0,94,by = 2)
      
      for(j in c(kk[jj1], kk[jj2]))
      {
        cat(j, '\n');
        #### Detail analysis with individual examples with plot
        ptm <- proc.time()
        Identifiablity.analysis.M3.with.plot.result(T = T, gene.index = j, debug = TRUE, zt = zt, i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE); 
        proc.time() - ptm
      }
      
    }
    
    ####################
    ### Check Refitting
    ####################
    j = kk[jj1[1]]
    #j = 701
    j = which(T$gene=='fake_m4_120')
    source('functions.R')
    ZT.ex = grep('.count.mRNA', colnames(T))
    ZT.int = grep('.count.premRNA', colnames(T))
    zt = seq(0,94,by = 2)
    ptm <- proc.time()
    param.fits.results = make.fits.with.all.models.for.one.gene.remove.outliers(T = T, gene.index = j, debug = TRUE,
                                                                                zt = zt, i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE);
    proc.time() - ptm
    index = match(c('outlier.m', 'outlier.s'), names(param.fits.results))
    res.fit = as.numeric(param.fits.results[-index])
    names(res.fit) = names(param.fits.results)[-index]
    #missing.data = length(unlist(strsplit(param.fits.results[index[1]], ';'))) + length(unlist(strsplit(param.fits.results[index[2]], ';')))
    test0 = c(res.fit, my.model.selection.one.gene.loglike(res.fit, nb.data = (96), method = 'BIC'))
    #print(test1)
    print(test0)
    print(T[j, c(142, 146, 147, 141, 140, 143:145)])
    
    #### contour plot of objective function
    #source('functions.R')
    #contour.plot.objective.funciton.parameters(T=T, j=j, zt=zt, i.ex=ZT.ex, i.int=ZT.int)
    
    #####################
    #### Find regions of non-identifiability in parameter space by simulaiton
    #####################
    Searching.for.non.identifiability.parameter.regions = FALSE
    if(Searching.for.non.identifiability.parameter.regions)
    {
      #### Generate simulated data for all parameter space (epsilon, gamma) 
      ####at different noise controlled by dispersion parameter alpha
      res.version = '_total_counts_all_genes_norm_params_v4';
      load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
      source('functions.R')
      BIC = my.model.selection.loglike(T=T, method='BIC', outliers = TRUE)
      T = data.frame(T, BIC, stringsAsFactors = FALSE);
      T = T[which(!is.na(T$BIC.best.model)==TRUE), ]
      
      ZT.int = grep('.count.premRNA', colnames(T))
      ZT.ex = grep('.count.mRNA', colnames(T))
      zt = seq(0,94,by = 2)  
      aa = as.matrix(T[, grep('alpha', colnames(T))])
      hist(log10(aa), breaks=200);
      a.dispersion = c(10^-4, 10^-3, 10^-2);
      nb.realization = 50;
      #gene.template = 'Per1'
      model = 3;
      
      #F0 = c()
      source('functions.R')
      for(n in 1:length(a.dispersion))
      {
        xx = generate.fake.read.counts.4identifiability.test(T = T, alpha = a.dispersion[n], nb.realization = nb.realization, model = model, zt=zt)
        xx$gene = paste('alpha_', n, '_', xx$gene, sep='')
        if(n==1)
        {
          F0 = xx;
        }else{
          F0 = rbind(F0, xx);
        }
      }
      
      fake.version = '_total_counts_s4_nonident_model3'
      save(F0, file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_simulated_data', fake.version, '.Rdata', sep=''))
      
      Test.fitting = TRUE
      if(Test.fitting)
      {
        fake.version = '_total_counts_s4_nonident_model3'
        load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_simulated_data', fake.version, '.Rdata', sep='')) 
        T = F0;
        
        j = which(T$gene=='alpha_1_fake_m3_100_r1')
        
        source('functions.R')
        ZT.ex = grep('.count.mRNA', colnames(T))
        ZT.int = grep('.count.premRNA', colnames(T))
        zt = seq(0,94,by = 2)
        ptm <- proc.time()
        param.fits.results = make.fits.with.all.models.for.one.gene.remove.outliers.4specific.model(T = T, gene.index = j, model=3, debug = FALSE,
                                                                                    zt = zt, i.ex = ZT.ex, i.int = ZT.int);
        proc.time() - ptm
        print(param.fits.results)
        print(T[j, c(124:131)])
        
      }
      
      Search.for.regions = FALSE
      if(Search.for.regions)
      {
        #res.version = '_total_counts_s3_nonident_model3_v1';
        res.version = '_total_counts_s4_nonident_model3_v1';
        load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', 
                          res.version, '.Rdata',sep = ''))
        #T = F0;
        normalize.parameteres = FALSE
        if(normalize.parameteres)
        {
          w = 2*pi/24;
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
        Plot.Non.Identifiability.Results = TRUE
        if(Plot.Non.Identifiability.Results)
        {
          alphas = a.dispersion;
          source('functions.R')
          #par(mfrow = c(1,length(alphas)))
          #par(mfrow=c(2,2))
          for(m0 in 1:length(alphas))
          {
            pdf.name = paste(file = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Total/Total_counts/non_ident_simulated_data_', 
                             res.version, '_alpha_', alphas[m0], '.pdf', sep='')
            pdf(pdf.name, width = 5, height = 3)
            par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
            
            m0 = 1;
            cat(m0, '\n');
            index0 = grep(paste('alpha_', m0, sep=''), T$gene);
            index.par = unique(unlist(strsplit(as.character(T$gene[index0]), '_'))[5+6*c(0:(length(index0)-1))])
            grids = matrix(NA, nrow = length(index.par), ncol = 3);
            for(m1 in 1:length(index.par))
            {
              ii = index.par[m1];
              kk = grep(paste('alpha_', m0, '_fake_m', model,'_', ii, '_r', sep=''), T$gene)
              eps.true = T$eps.gamma[kk[1]];
              gamma.true = T$gamma[kk[1]];
              ## error meric
              errs = T$gamma.m3[kk]/gamma.true
              errs = errs[which(!is.na(errs)==TRUE)];
              #err.estimate = mean(abs(T$gamma.m3[kk]-gamma.true)/gamma.true)
              err.estimate = mean(errs);
              grids[m1, ] = c(gamma.true, eps.true, err.estimate)
            }
            grids = data.frame(grids)
            colnames(grids) = c('gamma.true', 'eps.true', 'err.gamma.estimate')
            
            grids[,3] = -log2(grids[,3]); 
            
            grids[,1] = log10(log(2)/(grids[,1]))
            #grids[,1] = log10(grids[,1])
            
            w = 2*pi/24; A = 0.2;
            rr = seq(log(2)/24, log(2)/(10/60), length.out=1000)
            #zz = A*sqrt(w^2/rr^2+1)
            zz = A/rr
            library(lattice)
            levelplot(err.gamma.estimate~gamma.true*eps.true, grids, xlab="true half.life",  ylab="true eps.gamma", main=paste("Err of gamma estimation : alpha =  ", alphas[m0], sep=''),
                      sub=NA, region = TRUE, contour=FALSE, panel = function(...){
                        panel.levelplot(...)
                        panel.abline(v = log10(2))
                        #panel.abline(v = 0.4)
                        #panel.xyplot(rr, zz, col='blue', type='l', lwd=1.5)
                        },
                      #colorkey = list(at=seq(6, -4, by=1))
                      )
            #plot(c(1:10), c(1:10), add=TRUE)
            #abline(v=2*pi/24, col='red')
            dev.off()
            
          }
        }
      }
    
    }
  }
  
}

########################
###### Plot for Daniel of animal request
########################
Plot.4.Animal.Request = FALSE
if(Plot.4.Animal.Request){
  aa = read.table(file = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/Tables/Raw_data_total_RNA_seq_RPKM_normalized.txt', sep='\t', header = TRUE)
  
  data1 = log2(as.matrix(aa[, c(4:51)])); 
  data2 = data1[, c(1:36)]
  
  source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')
  
  stat1 = t(apply(data1,1, f24_R2_alt2, t=c(0:47)*2))
  stat1[,4] = t(apply(2^data1,1, f24_R2_alt2, t=c(0:47)*2))[,4]
  stat1 = cbind(stat1, qvals(stat1[,6]))
  
  stat2 = t(apply(data2,1, f24_R2_alt2, t=c(0:35)*2))
  stat2[,4] = t(apply(2^data2,1, f24_R2_alt2, t=c(0:35)*2))[,4]
  stat2 = cbind(stat2, qvals(stat2[,6]))
  
  cutoff = seq(-6, 0, by=0.2)
  nbs = c()
  for(n in 1:length(cutoff)){
    nbs = rbind(nbs, c(length(which(stat1[,6]<10^(cutoff[n]))), length(which(stat2[,6]<10^(cutoff[n])))))
  }
  
  pdfname = paste('/Users/jiwang/Dropbox/GachonProt/Cumulant_genes.pdf', sep='')
  pdf(pdfname, width=1.8, height=1.8)
  par(cex = 0.7, las = 1, mgp = c(2.0,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
  
  plot(cutoff, nbs[, 2], log='', type='l', col='red', lwd=1.5, xlab='log10(p-val)', ylab='nb of rhythmic genes', cex.lab=0.7, axes=FALSE)
  points(cutoff, nbs[,1], type='l', col='green', lwd=1.5)
  abline(v=-2, col='orange', lwd=1., lty=2)
  abline(h=nbs[which(cutoff==-2),1], lwd=0.7, lty=1, col='gray')
  abline(h=nbs[which(cutoff==-2),2], lwd=0.7, lty=1, col='gray')
  box()
  axis(1, at=seq(-8, 0, by=1), cex.axis=0.7)
  axis(2, las=1, at=seq(0, 12000, by=2000), cex.axis=0.7, tick = FALSE)
  axis(2, las=1, at=seq(0, 12000, by=1000), cex.axis=0.7, labels = FALSE)
  dev.off()
  #colnames(stat1) = paste(c('nb.timepoints', 'mean', 'amp', 'rel.amp', 'phase','pval', 'qv'), '.log2rpkm.mRNA', sep='')
  #colnames(stat2) = paste(c('nb.timepoints', 'mean', 'amp', 'rel.amp', 'phase','pval', 'qv'), '.log2rpkm.premRNA', sep='')
  
  jj = which(stat1[, 6]<0.01 & stat2[, 6]>0.01)
  xx = data.frame(jj, aa$gene[jj], stat1[jj,6], stat2[jj, 6], stringsAsFactors = FALSE)
  plot(xx[, c(3, 4)], log='xy');abline(0, 1, col='red')
  text(xx[, (3)], xx[,4], xx[,2], offset=0.2, pos=2)
  
  kk = which(aa$gene=='Bag2')
  
  zt = seq(0, 22, by=2)
  err.S = mean.err(data1[kk,], interval = 2)[2,];mean.S = mean.err(data1[kk,], interval = 2)[1,];
  err.M = mean.err(data2[kk,], interval = 2)[2,];mean.M = mean.err(data2[kk,], interval = 2)[1,];
  xlim = c(0, 24); 
  ylim = range(c((mean.S+err.S), (mean.S-err.S), s, (mean.M+err.M),(mean.M-err.M)), na.rm = TRUE);
  #ylim = range(c((mean.S+err.S), (mean.S-err.S), s, (mean.M+err.M),(mean.M-err.M), m, dd, m22))
  #if(ylim[1]<10^-6){ylim = range(c(s, m, dd));}

  pdfname = paste('/Users/jiwang/Dropbox/GachonProt/gene_examples.pdf', sep='')
  pdf(pdfname, width=1.8, height=1.6)
  par(cex = 0.7, las = 1, mgp = c(2.0,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
  #### relative RPKM with fitting for one day
  plot(zt, mean.S, col = 'green', type = 'n', lwd = 1.0, ylim=ylim, xlim=xlim, 
       main = aa$gene[kk], cex.main=1.0, xlab=NA, log='', ylab=NA, axes=TRUE);
  points(zt, mean.S, col='green', type='p', cex=0.6, pch=16);
  points(zt, mean.S, col='green', type='l', cex=0.6, pch=16);
  arrows(zt, mean.S-err.S, zt, mean.S+err.S, length=0.03, angle=90, code=3, col='green', lwd=1, cex=1.0)
  #points(zt, M, type='l', lwd=0.8, col='steelblue');
  #points(zt, mean.M, type='p', cex=0.6, col='steelblue', pch=16);
  #points(zt, mean.M, type='l', cex=0.6, col='steelblue', pch=16);
  #arrows(zt, mean.M-err.M, zt, mean.M+err.M, length=0.03, angle=90, code=3, col=col.ex, lwd=1.0, cex=1.0)
  text(5, 0.4, paste('pval = ', signif(stat1[kk, 6], d=1), sep=''), cex=0.7)
  dev.off()
  
  pdfname = paste('/Users/jiwang/Dropbox/GachonProt/gene_examples_2.pdf', sep='')
  pdf(pdfname, width=1.8, height=1.6)
  par(cex = 0.7, las = 1, mgp = c(2.0,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
  #### relative RPKM with fitting for one day
  plot(zt, mean.S, col = 'green3', type = 'n', lwd = 1.0, ylim=ylim, xlim=xlim, 
       main = aa$gene[kk], cex.main=1.0, xlab=NA, log='', ylab=NA, axes=TRUE);
  #points(zt, mean.S, col='green3', type='p', cex=0.6, pch=16);
  #points(zt, mean.S, col='green3', type='l', cex=0.6, pch=16);
  #arrows(zt, mean.S-err.S, zt, mean.S+err.S, length=0.03, angle=90, code=3, col='green3', lwd=1, cex=1.0)
  #points(zt, M, type='l', lwd=0.8, col='steelblue');
  points(zt, mean.M, type='p', cex=0.6, col='red', pch=16);
  points(zt, mean.M, type='l', cex=0.6, col='red', pch=16);
  arrows(zt, mean.M-err.M, zt, mean.M+err.M, length=0.03, angle=90, code=3, col='red', lwd=1.0, cex=1.0)
  text(5, 0.4, paste('pval = ', signif(stat2[kk, 6], d=1), sep=''), cex=0.7)
  dev.off()
  
  #bb = data.frame(aa, stat1, stat2, stringsAsFactors=FALSE)
}

#############
### Search for candidates for Clemence
##############
Candidates4Clemence = FALSE
if(Candidates4Clemence)
{
  data.version = '_total_counts_v2';
  load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_
                  analysis_sel_alphas', data.version, '.Rdata', sep=''))
  #source("/Users/jiwang/Degradation_Liver/Main_Code_Degradation/f24_modified_1.0.r");
  #attribute.global.variable.colors();
  load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_sel_alphas_log2_RPKM', 
                  data.version, '.Rdata', sep=''))
  T = R.norm;
  #T = 2^T;
  source('functions.R')
  
  ZT.ex = grep('.rpkm.mRNA', colnames(T))[1:48]
  ZT.int = grep('.rpkm.premRNA', colnames(T))[1:48]
  #zt = seq(0,94,by = 2)
  
  xx = T[, c(1:3, ZT.ex, ZT.int)]
  
  write.table(xx, file='total_RNA_seq_mRNA_premRNA_expressed_linear.txt', col.names = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
  
  
  min(as.numeric(T[which(T$gene=='Arntl'),ZT.ex]))
  
  data = as.matrix(T[, ZT.ex])
  max = apply(data, 1, max); min = apply(data, 1, min);fc = T$amp.log2rpkm.mRNA;
  kk = which(T$qv.log2rpkm.mRNA<10^-3 & min<20 & fc>1  & max>10)
  
  plot(T$amp.log2rpkm.mRNA[kk], log10(max[kk]), log='', cex=0.5)
  abline(h=c(2, 2.5), col='green')
  text(T$amp.log2rpkm.mRNA[kk], log10(max[kk]), T$gene[kk], cex=0.7, pos=3)
  
}



