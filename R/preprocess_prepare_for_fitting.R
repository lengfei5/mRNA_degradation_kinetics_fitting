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
MDfitDataSet = function(T=T, zt=zt, i.ex. = ZT.ex, i.int = ZT.int)
{
  
  
  return(mds)
}
set.scaling.factors = function()
{
    scaling.factors <<- c(39920608, 42250245, 38121270, 45609244, 41511752, 45781196, 43722568, 39638552, 30496638, 30573333, 54950572, 47158379,
                        31722765, 39931646, 36317783, 35382708, 47293167, 42408985, 39842283, 40230336, 43691685, 39237518, 51051196, 44778546,
                        43858841, 42791401, 42357301, 49782402, 44628140, 44561463, 43485553, 47853067, 43318817, 45055723, 30180984, 46825671,
                        43270558, 37496344, 40971385, 45828360, 37065376, 35776330, 45025514, 43026714, 43116633, 35173387, 28538212, 36707156);
    
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