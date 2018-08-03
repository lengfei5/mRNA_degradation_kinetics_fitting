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
library(limma)

## creat a S3 class MDfitDataSet to store data matrix, (P and M for pre-mRNA and mRNA), time points (zt), length.pre-mRNA and length.mRNA
## scaling factors for rpkm calculation, dispersion parameters
## some simple functions associated to extract these parameters, which will be used as global parameters in the model fitting
MDfitDataSet = function(P, M, length.P=c(10000), length.M=c(1000), zt=zt, mode = "NB", 
                        fitType.dispersion = "local", fitType.var = "separate")
{
  cat("creat a S3 class MDfitDataSet to store the tables and also necessary parameters after processing ...\n")
  
  mds = list();
  class(mds) = "MDfitDataSet"
  
  # store the input tables
  mds$P = data.frame(P);
  mds$M = data.frame(M);
  
  # associated time points
  mds$zt = zt;
  
  mds$mode = mode;
  
  if(mode == "NB")
  {
    # lengths of pre-mRNA and mRNA
    mds$length.P = length.P
    mds$length.M = length.M
    
    # scaling factor for each time points calculated by DESeq2 
    mds$scaling.factors = calculate.scaling.factors.DESeq2(P, M, zt)
    
    # dispersion parameter for each time point and for pre-mRNA and mRNA respectivley, calculated by DESeq2
    # becasue 4 replicates for each time point ennables us to calculate the dispersion in such way; and also because
    # large differences in gene expression mean between time points requires this. 
    estimatedDispersions = calculate.dispersions.for.each.time.point.DESeq2(P, M, zt,  fitType.dispersion = fitType.dispersion)
    
    mds$dispersions.P = estimatedDispersions$alphas.P
    mds$dispersions.M = estimatedDispersions$alphas.M
    
  }else{
    if(mode == "logNormal"){
      ####################
      # logNormal mode
      # employ limma package to estimate variances for eahc gene and for each time points (empirical Bayes methods implemented in limma
      # and the shrinked variance estimation is more robust in the context of small samples)
      # in addtion, here only the sample for the same time points considered as biological replicates 
      # to have the gene-wise and time-wise variance estiamtion, because for circadian genes, the expression can change dramatically and 
      # variance depends on the expression. Thus commonly used only gene-wise variance (e.g. DESeq2 and limma) is not appropriate any more.
      # After checking the data, it turns out that the pre-mRNAs have different variance dependance on the log expression
      # the variance was estiamted separatetly for pre-mRNAs and mRNAs
      ####################
      #estimateDispersions = estimateVariances.for.each.time.point.limma(P, M, zt)
      estimatedVars = estimateVariances.for.each.time.point.limma(P = P, M = M, zt, fitType.var = fitType.var);
      
      mds$var.P = data.frame(estimatedVars$var.P)
      mds$var.M = data.frame(estimatedVars$var.M)
      
    }else{
      stop("current function supports only 'NB' and 'logNormal' two modes")
    }
    
  }
  
  return(mds)
}

print.MDfitDataSet = function(mds)
{
  cat("MDfitDataSet object (S3 class) \n")
  
  for(n in 1:length(mds))
  {
    cat("$", names(mds)[n], "\n")
    if(is.vector(mds[[n]])) {
      max.nb = min(length(mds[[n]]), length(mds$zt))
      cat("\t", mds[[n]][1:max.nb], "\n")
    }
    if(is.data.frame(mds[[n]])) print(mds[[n]][c(1:2), ])
  }
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

calculate.dispersions.for.each.time.point.DESeq2 = function(P, M, zt,  fitType.dispersion = "local") 
{
  # require('DESeq2')
  # P = T[, ZT.int]; M = T[, ZT.ex]; fitType.dispersion = "local"
  countData = rbind(as.matrix(M), as.matrix(P))
  size.factors = calculate.SizeFactors.DESeq2(countData)
  
  # estimate dispersion parameters with Deseq2 for each gene and each condition
  zt.24 = zt%%24
  zt.uniq = unique(zt.24)
  alphas = matrix(NA, nrow=nrow(countData), ncol=length(zt.uniq))
  
  rownames(countData) = c(1:nrow(countData))
  
  for(n in 1:ncol(alphas))
  {
    # n = 1
    index.reps = which(zt.24==zt.uniq[n])
    cat('ZT',  zt.uniq[n], "--", length(index.reps), "replicates", '\n');
    
    #countData2 = countData[, index];
    condition <- factor(rep(c(1), length(index.reps)))
    dds = DESeqDataSetFromMatrix(countData[, index.reps], DataFrame(condition), ~1);
    sizeFactors(dds) = size.factors[index.reps]
    
    dds <- estimateDispersions(dds, maxit = 500, fitType=fitType.dispersion)
    alphas[,n] = dispersions(dds);
    
    #plotDispEsts(dds)
    #plot(alphas[c(1:nrow(T)), 1], T$alpha.mRNA.ZT0, log='xy', cex=0.2);abline(0, 1, lwd=2.0, col='red')
  }
  
  alphas.all = alphas[, match(zt.24, zt.uniq)]
  alphas.M = data.frame(alphas.all[c(1:nrow(M)), ])
  alphas.P = data.frame(alphas.all[-c(1:nrow(M)), ])  
  
  #colnames(xx) = c(paste('alpha.mRNA.ZT', c(0:11)*2, sep=''),  paste('alpha.premRNA.ZT', c(0:11)*2, sep=''))
  colnames(alphas.P) = paste0('alpha.premRNA.ZT', zt)
  colnames(alphas.M) = paste0('alpha.mRNA.ZT', zt)
  
  return(list(alphas.P = alphas.P, alphas.M = alphas.M))
  
}

estimateVariances.for.each.time.point.limma = function(P, M, zt, robust = TRUE, fitType.var = "separate")
{
  # P = T[, ZT.int]; M = T[, ZT.ex]; robust = TRUE; 
  if(fitType.var == "pool"){
    cat("estimate variance for pre-mRNA and mRNA---\n")
    
    normData = rbind(as.matrix(P), as.matrix(M));
    var.PM = data.frame(squeezeVar.fromMatrix(normData, zt))
    var.P = var.PM[c(1:nrow(P)), ]
    var.M = var.PM[(c(nrow(P)+1):nrow(var.PM)), ]
    
  }else{
    cat("estimate variance for pre-mRNA ---\n")
    var.P = data.frame(squeezeVar.fromMatrix(P, zt))
    
    cat("estimate variance for mRNA ---\n")
    var.M = data.frame(squeezeVar.fromMatrix(M, zt))
  }
  
  return(list(var.P = var.P, var.M = var.M))
  
}

squeezeVar.fromMatrix = function(normData, zt, robust = TRUE)
{
  normData = as.matrix(normData)
  normData[which(normData==0 | is.na(normData))] = 10^-10; 
  normData = log(normData)
  # estimate dispersion parameters with Deseq2 for each gene and each condition
  zt.24 = zt%%24
  zt.uniq = unique(zt.24)
  
  #condition <- factor(rep(c(1), 4))
  rownames(normData) = c(1:nrow(normData))
  out.vars = matrix(NA, nrow=nrow(normData), ncol=length(zt.uniq))
  
  ## pool vars and mean for all time points
  vars = c();
  covariate = c();
  df = c()
  for(n in 1:length(zt.uniq))
  {
    # n = 1
    index.reps = which(zt.24==zt.uniq[n])
    cat('\tZT',  zt.uniq[n], "--", length(index.reps), "replicates", '\n');
    norm.sel = normData[, index.reps];
    
    vars = c(vars, apply(norm.sel, 1, var));
    df = c(df, rep((length(index.reps)-1), nrow(norm.sel)));
    covariate = c(covariate, apply(norm.sel, 1, mean));
    out <- squeezeVar(vars, df, covariate=covariate, robust=robust);
    #covariate.P = covariate;
    #out.P = out
    #plot(covariate, out$var.prior^(1/4),col = 'red', cex=1., ylim = c(0, 1.5))
    #jj = c(1:(length(covariate)/2))
    #points(covariate[jj], vars[jj]^(1/4), cex=0.05, col = 'black' )
    #points(covariate[-jj], vars[-jj]^(1/4), cex=0.05, col = 'blue' )
    #points(covariate.P, out.P$var.prior^(1/4), col = 'blue', cex=1.)
    #points(covariate, vars^(1/4), cex=0.05, col = 'black' )
    #points(covariate, out$var.post^(1/4), col='blue', cex=0.4)
  }
  
  out <- squeezeVar(vars, df, covariate=covariate, robust=robust);
  
  out.vars = out$var.post
  dim(out.vars) = c(nrow(normData), length(zt.uniq))
  out.vars.all = out.vars[, match(zt.24, zt.uniq)]
  
  plot.EB.Trend.shunkage = FALSE
  if(plot.EB.Trend.shunkage){
    plot(covariate, out$var.prior^(1/4),col = 'red', cex=1., ylim = c(0, 4))
    points(covariate, vars^(1/4), cex=0.05, col = 'black' )
    points(covariate, out$var.post^(1/4), col='blue', cex=0.4)
  }
  
  colnames(out.vars.all) = paste0("variance.ZT", zt)
  
  #colnames(alphas.P) = paste0('alpha.premRNA.ZT', zt)
  #colnames(alphas.M) = paste0('alpha.mRNA.ZT', zt)
  return(out.vars.all)
  
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

Make.data.example.RPKM = function(mds)
{
  xx = matrix(NA, nrow = nrow(mds$P), ncol = ncol(mds$P))
  yy = xx;
  
  for(n in 1:nrow(xx))
  {
    xx[n, ] = norm.RPKM(unlist(mds$P[n,]), mds$length.P[n]) 
    yy[n, ] = norm.RPKM(unlist(mds$M[n,]), mds$length.M[n]) 
  }
  colnames(xx) = paste0('ZT', mds$zt, ".rpkm.premRNA")
  colnames(yy) = paste0('ZT', mds$zt, ".rpkm.mRNA")
  
  R = cbind(xx, yy)
  return(R);
  
}

Test.limma.eBayes = function()
{
  #  Simulate gene expression data,
  #  6 microarrays and 100 genes with one gene differentially expressed
  set.seed(2016)
  sigma2 <- 0.05 / rchisq(100, df=10) * 10
  y <- matrix(rnorm(100*6,sd=sqrt(sigma2)),100,6)
  design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))
  y[1,4:6] <- y[1,4:6] + 1
  fit <- lmFit(y,design)
  
  ## part of source code from eBayes function in limma (https://rdrr.io/bioc/limma/src/R/ebayes.R)
  trend = TRUE;
  robust = TRUE;
  coefficients <- fit$coefficients
  stdev.unscaled <- fit$stdev.unscaled
  sigma <- fit$sigma
  df.residual <- fit$df.residual
  if(is.null(coefficients) || is.null(stdev.unscaled) || is.null(sigma) || is.null(df.residual)) stop("No data, or argument is not a valid lmFit object")
  if(all(df.residual==0)) stop("No residual degrees of freedom in linear model fits")
  if(all(!is.finite(sigma))) stop("No finite residual standard deviations")
  if(trend) {
    covariate <- fit$Amean
    if(is.null(covariate)) stop("Need Amean component in fit to estimate trend")
  } else {
    covariate <- NULL
  }
  
  #	Moderated t-statistic
  out <- squeezeVar(sigma^2, df.residual, covariate=covariate, robust=robust)
  out$s2.prior <- out$var.prior
  out$s2.post <- out$var.post
  out$var.prior <- out$var.post <- NULL
  out$t <- coefficients / stdev.unscaled / sqrt(out$s2.post)
  df.total <- df.residual + out$df.prior
  df.pooled <- sum(df.residual,na.rm=TRUE)
  df.total <- pmin(df.total,df.pooled)
  out$df.total <- df.total
  
}

