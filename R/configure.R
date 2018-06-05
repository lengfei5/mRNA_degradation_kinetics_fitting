##########################################################################
##########################################################################
## Project: 
## Script purpose: install R pakcages and make configuration 
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Tue Jun  5 11:28:20 2018
##########################################################################
##########################################################################
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

library(emdbook)
library(deSolve)
library(fdrtool)
library(circular)
library(preprocessCore)
library(gtools)
library(biomaRt)
#library(plotrix)
library(numDeriv)
library('Matrix')
#library('matrixcalc')
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
