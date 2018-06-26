##########################################################################
##########################################################################
## Project:
## Script purpose: 
## preprocess the data table and estimate parameters required afterwards but needs the whole dataset to obtain
## 0) import libraries
## 1) processing the expresson table 
## 2) calculate the size factor and dispersion parameters using DESeq2
## 3) calculate variances using voom
## 4) Important Note !!!!!! 
##   
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Tue Jun  5 11:31:45 2018
##########################################################################
##########################################################################
library(emdbook)
library(deSolve)
library(fdrtool)
library(circular)
library(preprocessCore)
library(gtools)
library(biomaRt)
library(numDeriv)
library(Matrix)
library(DESeq2)

## creat a S3 class MDfitDataSet to store data matrix, (P and M for pre-mRNA and mRNA), time points (zt), length.pre-mRNA and length.mRNA
## scaling factors for rpkm calculation, dispersion parameters
## some simple functions associated to extract these parameters, which will be used as global parameters in the model fitting
MDfitDataSet = function(P, M, length.P, length.M, zt=zt)
{
  cat("creat a S3 class MDfitDataSet to store the tables and also necessary parameters after processing ...")
  
  mds = list();
  class(mds) = "MDfitDataSet"
  
  # store the input tables
  mds$P = P;
  mds$M = M;
  
  # lengths of pre-mRNA and mRNA
  mds$length.P = length.P
  mds$length.M = length.M
  
  # associated time points
  mds$zt = zt;
  
  # scaling factor for each time points calculated by DESeq2 
  mds$scaling.factors = calculate.scaling.factors.DESeq2(P, M, zt)
  
  # dispersion parameter for each time point and for pre-mRNA and mRNA respectivley, calculated by DESeq2
  # becasue 4 replicates for each time point ennables us to calculate the dispersion in such way; and also because
  # large differences in gene expression mean between time points requires this. 
  estimateDispersions = calculate.dispersions.for.each.time.point.DESeq2(P, M, zt)
  mds$dispersions.P = estimateDispersions$alphas.P 
  mds$dispersions.M = estimateDispersions$alphas.M
  
  return(mds)
}

calculate.SizeFactors.DESeq2 = function(countData)
{
  require('DESeq2')
  size.factors = estimateSizeFactorsForMatrix(countData) 
  
  return(size.factors)
}

calculate.scaling.factors.DESeq2 = function(P, M, zt)
{
  # P = T[, ZT.int]; M = T[, ZT.ex]
  countData =  rbind(as.matrix(M), as.matrix(P))
  size.factors = calculate.SizeFactors.DESeq2(countData)
  
  # here similar but different from DESeq2 for the scaling factor, because a constant mean(ss/size.factors) was multipled with size.factors
  # However in DESeq2, when fpkm is calculated, a slight different constant is used, which is genometric mean of library sizes.
  ss = apply(countData, 2, sum)
  scaling.factors = size.factors*mean(ss/size.factors)
  names(scaling.factors) = paste0('ZT', zt)
  
  ## here is the old scaling factor used in the paper
  #scaling.factors = c(39920608, 42250245, 38121270, 45609244, 41511752, 45781196, 43722568, 39638552, 30496638, 30573333, 54950572, 47158379,
  #                      31722765, 39931646, 36317783, 35382708, 47293167, 42408985, 39842283, 40230336, 43691685, 39237518, 51051196, 44778546,
  #                      43858841, 42791401, 42357301, 49782402, 44628140, 44561463, 43485553, 47853067, 43318817, 45055723, 30180984, 46825671,
  #                      43270558, 37496344, 40971385, 45828360, 37065376, 35776330, 45025514, 43026714, 43116633, 35173387, 28538212, 36707156);
  
  return(scaling.factors)
  
}

calculate.dispersions.for.each.time.point.DESeq2 = function(P, M, zt, fitType = 'local')
{
  require('DESeq2')
  # P = T[, ZT.int]; M = T[, ZT.ex]
  countData = rbind(as.matrix(M), as.matrix(P))
  size.factors = calculate.SizeFactors.DESeq2(countData)
  
  # estimate dispersion parameters with Deseq2 for each gene and each condition
  zt.24 = zt%%24
  zt.uniq = unique(zt.24)
  alphas = matrix(NA, nrow=nrow(countData), ncol=length(zt.uniq))
    
  condition <- factor(rep(c(1), 4))
  rownames(countData) = c(1:nrow(countData))
  
  for(n in 1:ncol(alphas))
  {
    index.reps = which(zt.24==zt.uniq[n])
    cat('ZT ',  zt.uniq[n], "--", length(index.reps), "replicates", '\n');
    
    #countData2 = countData[, index];
    dds = DESeqDataSetFromMatrix(countData[, index.reps], DataFrame(condition), ~1);
    sizeFactors(dds) = size.factors[index.reps]
    
    dds <- estimateDispersions(dds, maxit = 500, fitType=fitType)
    alphas[,n] = dispersions(dds);
    
    plotDispEsts(dds)
    
    plot(alphas[c(1:nrow(T)), 1], T$alpha.mRNA.ZT0, log='xy', cex=0.2);abline(0, 1, lwd=2.0, col='red')
    
    #alphas.genes[,n] = mcols(dds)$dispGeneEst
  }
  
  
  xx = data.frame(alphas[c(1:nrow(T)),], alphas[-c(1:nrow(T)), ])
  colnames(xx) = c(paste('alpha.mRNA.ZT', c(0:11)*2, sep=''),  paste('alpha.premRNA.ZT', c(0:11)*2, sep=''))
  alphas = xx;
  
  # T = data.frame(T[,c(1:3)], alphas, T[, -c(1:5)], stringsAsFactors = FALSE)
  #Disps = 1;
  return(list(alphas.P = alphas.P, alphas.M = alphas.M))
  
}

test.make.S3.class = function()
{
  j = list (name = "Joe", salary = 2000, union = T)
  class(j) = "employee"
  print.employee = function(wrkr)
  {
    cat(wrkr$name, "\n")
    cat("salary", wrkr$salary, "\n")
    cat("union member", wrkr$union, "\n")
  }
  
  
}