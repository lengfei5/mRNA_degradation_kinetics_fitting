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
    
    gg = 'Per2'
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
  
   
}
